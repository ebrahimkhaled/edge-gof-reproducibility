## grid_asymmetry.R — EXPERIMENT 2.2: the LINK-ASYMMETRY SWEEP (SPECKIT §5 Exp 2.2).
## DGP: y ~ inv_stukel(0.6*x + 0.5*d, a1=-1, a2=a), a in seq(-1,1,0.25); fit y~x+d.
##   a=-1  => inv_stukel(.,-1,-1) (symmetric heavy tails, a "parity" corner)
##   a=+1  => inv_stukel(.,-1, 1) (the asymmetric stukel_asym departure)
## n in {500,1000}, B=5000, curated rival set (A3 -- NO slow rivals), G=10, seed 20260709.
##
## Primary output = Delta = reject(EF) - reject(HL) with +/-2*mcse_diff bands, where
##   mcse_diff = sqrt((p10 + p01)/B)  (paired McNemar cells from the stored per-rep p-values,
##   at alpha=0.05); PLUS Stukel's own reject curve.
##
## Contract: source _harness.R then source(DGP). Uses ek_cluster/ek_run_cell/ek_append.
## Store per-rep p.<test> AND s.<test> so size-adjusted power / alpha quantiles recompute.
## STOP: Delta must be ~0 at a=-1 and rise monotone to > +0.03 at a=+1; flat/decreasing =>
##       append to _STOP_TRIGGERS.log + cat() (never alter numbers to avoid a trigger).
## ---- archive paths (inserted by make_archive.py; the convention is in README.md, "How to run") -------------
EDGE_ARCHIVE_ROOT <- local({
  r <- Sys.getenv("EDGE_ARCHIVE_ROOT")
  if (!nzchar(r)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    r <- if (length(f)) file.path(dirname(normalizePath(f[1], winslash = "/")), "..", "..") else getwd()
  }
  r <- normalizePath(r, winslash = "/", mustWork = FALSE)
  Sys.setenv(EDGE_ARCHIVE_ROOT = r)        # so worker processes started from here resolve the same root
  r
})
edge_path <- function(...) file.path(EDGE_ARCHIVE_ROOT, ...)
edge_battery <- function(...) {           # the author's battery/ folder: analysis/ -> results/analysis, rest -> results/blocks
  if (...length() == 0L) return(edge_path("results", "blocks"))
  p <- file.path(...)
  ifelse(p == "analysis" | startsWith(p, "analysis/"), edge_path("results", p), edge_path("results", "blocks", p))
}
edge_out <- function(...) { d <- edge_path("output", ...); dir.create(d, showWarnings = FALSE, recursive = TRUE); d }
## ---------------------------------------------------------------------------------------------------------------

suppressMessages(library(parallel))

SIMDIR <- edge_path("code/simulations")
source(file.path(SIMDIR, "_harness.R"))   # SIMDIR, DGP, EK_NCORES, ek_reps, ek_cluster, ek_run_cell, ek_append, ek_mcse
source(DGP)                                # inv_stukel, run_curated_full, EK_CURATED (also loaded in workers by ek_cluster)

## ---------------------------------------------------------------------------
## Paths & constants
## ---------------------------------------------------------------------------
PV_PATH   <- file.path(SIMDIR, "sim_asym_sweep_pvalues.csv")  # per-rep p AND statistic (wide)
SUM_PATH  <- file.path(SIMDIR, "sim_asym_sweep.csv")          # a,n,test,reject_rate,mcse (+alpha)
DELTA_PATH<- file.path(SIMDIR, "sim_asym_delta.csv")          # the Delta = reject(EF)-reject(HL) series
STOP_LOG  <- file.path(SIMDIR, "_STOP_TRIGGERS.log")

GRID_SEED  <- 20260709L
B_SWEEP    <- 5000L
ALPHAS     <- c(0.01, 0.05, 0.10)
ALPHA_STAR <- 0.05           # alpha at which Delta / STOP trigger are evaluated
G_FIX      <- 10L

## ---------------------------------------------------------------------------
## Cell grid: a in seq(-1,1,0.25) (9 values) x n in {500,1000}. Distinct seed_base
## per cell = GRID_SEED + cell_id*1e5. cell_id fixed on the FULL grid so a subset
## keeps identical seeds.
## ---------------------------------------------------------------------------
A_SEQ <- seq(-1, 1, 0.25)
N_SEQ <- c(500L, 1000L)

cells <- expand.grid(a = A_SEQ, n = N_SEQ, G = G_FIX,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
cells$cell_id <- seq_len(nrow(cells))

## Optional self-test subset: EK_ASYM_CELLS="1,9,10,18" runs only those (stable) cell_ids.
.subset <- Sys.getenv("EK_ASYM_CELLS", "")
if (nzchar(.subset)) {
  keep  <- as.integer(strsplit(.subset, ",")[[1]])
  cells <- cells[cells$cell_id %in% keep, , drop = FALSE]
  cat("SUBSET run: cell_ids", paste(keep, collapse = ","), "\n")
}

## ---------------------------------------------------------------------------
## one_rep: returns a NAMED numeric vector (identical names every call).
##   Custom Stukel-asymmetry link DGP: y ~ inv_stukel(0.6*x + 0.5*d, -1, a); fit y~x+d.
##   (inv_stukel is sourced from the DGP library in every worker.)
## ---------------------------------------------------------------------------
one_rep <- function(rep, cell) {
  x <- runif(cell$n, -3, 3)
  d <- rbinom(cell$n, 1, 0.5)
  y <- rbinom(cell$n, 1, inv_stukel(0.6 * x + 0.5 * d, a1 = -1, a2 = cell$a))
  fit <- suppressWarnings(glm(y ~ x + d, family = binomial()))
  run_curated_full(fit, G = cell$G)   # named: p.<test> then s.<test>, curated order = EK_CURATED
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------
B <- ek_reps(B_SWEEP)
cat(sprintf("grid_asymmetry.R: %d cells (a x n), B=%d (REPS_SCALE=%s), EK_NCORES=%d\n",
            nrow(cells), B, Sys.getenv("REPS_SCALE", "1"), EK_NCORES))

## Fresh outputs each run (this script fully regenerates its CSVs).
for (p in c(PV_PATH, SUM_PATH, DELTA_PATH)) if (file.exists(p)) file.remove(p)

## NOTE: do NOT use on.exit() here -- at the top level of a source()'d script it is not
## bound to a function frame and can fire immediately (killing the cluster before the loop).
## Register cleanup via reg.finalizer / an explicit stopCluster at the end instead.
cl <- ek_cluster(GRID_SEED)

pcols <- paste0("p.", EK_CURATED)
scols <- paste0("s.", EK_CURATED)

## Accumulate the (a,n) reject rates for EF/HL/Stukel + the paired McNemar cells for Delta.
delta_rows <- list()

for (i in seq_len(nrow(cells))) {
  cell      <- as.list(cells[i, ])
  seed_base <- GRID_SEED + cell$cell_id * 1e5
  t0 <- proc.time()[["elapsed"]]

  df <- ek_run_cell(cl, B, seed_base, one_rep, cell)   # B rows x (p.* , s.*)

  ## ---- per-rep p-value/statistic CSV (wide): family,a,n,G,rep,seed, p.*, s.* ----
  reps <- seq_len(nrow(df))
  pv <- data.frame(family = "asym", a = cell$a, n = cell$n, G = cell$G,
                   rep = reps, seed = seed_base + reps,
                   check.names = FALSE, stringsAsFactors = FALSE)
  pv <- cbind(pv, df[, pcols, drop = FALSE], df[, scols, drop = FALSE])
  ek_append(pv, PV_PATH)

  ## ---- summary rollup per (a,n,test,alpha): reject_rate + mcse ----
  sum_rows <- list()
  for (tst in EK_CURATED) {
    pval <- df[[paste0("p.", tst)]]
    nval <- sum(!is.na(pval))
    for (a in ALPHAS) {
      reject_rate <- if (nval > 0) mean(pval < a, na.rm = TRUE) else NA_real_
      mcse        <- ek_mcse(reject_rate, B)
      sum_rows[[length(sum_rows) + 1L]] <- data.frame(
        family = "asym", a = cell$a, n = cell$n, G = cell$G, test = tst, alpha = a,
        reject_rate = reject_rate, mcse = mcse, B = B,
        stringsAsFactors = FALSE)
    }
  }
  ek_append(do.call(rbind, sum_rows), SUM_PATH)

  ## ---- Delta = reject(EF) - reject(HL) with paired McNemar mcse_diff at ALPHA_STAR ----
  ef_p <- df[["p.EF"]]; hl_p <- df[["p.HL"]]; st_p <- df[["p.Stukel"]]
  paired <- !is.na(ef_p) & !is.na(hl_p)                 # reps where BOTH computed
  rej_ef <- ef_p[paired] < ALPHA_STAR
  rej_hl <- hl_p[paired] < ALPHA_STAR
  Bp     <- sum(paired)
  reject_ef <- if (Bp > 0) mean(rej_ef) else NA_real_
  reject_hl <- if (Bp > 0) mean(rej_hl) else NA_real_
  p10 <- if (Bp > 0) mean(rej_ef & !rej_hl) else NA_real_   # EF rejects, HL not
  p01 <- if (Bp > 0) mean(!rej_ef & rej_hl) else NA_real_   # HL rejects, EF not
  delta      <- reject_ef - reject_hl
  mcse_diff  <- if (Bp > 0) sqrt((p10 + p01) / Bp) else NA_real_
  reject_stk <- if (any(!is.na(st_p))) mean(st_p < ALPHA_STAR, na.rm = TRUE) else NA_real_

  delta_rows[[length(delta_rows) + 1L]] <- data.frame(
    a = cell$a, n = cell$n, G = cell$G, alpha = ALPHA_STAR,
    reject_EF = reject_ef, reject_HL = reject_hl, reject_Stukel = reject_stk,
    delta = delta, p10 = p10, p01 = p01, mcse_diff = mcse_diff,
    delta_lo = delta - 2 * mcse_diff, delta_hi = delta + 2 * mcse_diff,
    B_paired = Bp,
    stringsAsFactors = FALSE)

  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("  cell %d/%d done: a=%+.2f n=%d  Delta=%+.4f (+/-%.4f)  Stukel=%.3f  (%.1fs)\n",
              i, nrow(cells), cell$a, cell$n,
              ifelse(is.na(delta), NA, delta), ifelse(is.na(mcse_diff), NA, 2 * mcse_diff),
              ifelse(is.na(reject_stk), NA, reject_stk), el))
}

delta_df <- do.call(rbind, delta_rows)
ek_append(delta_df, DELTA_PATH)

## ---------------------------------------------------------------------------
## STOP trigger (SPECKIT §5 Exp 2.2): Delta must be ~0 at a=-1 and rise MONOTONE
## to > +0.03 at a=+1. Flat/decreasing => halt-and-report (parity claim unsupported).
## Evaluated per n on the Delta series at ALPHA_STAR. Never alter numbers to avoid it.
## "~0 at a=-1": |Delta(a=-1)| <= 2*mcse_diff (a tie with zero).
## "monotone": non-decreasing in a within 2*mcse_diff tolerance (no real drop).
## ---------------------------------------------------------------------------
stop_hits <- list()
for (nn in sort(unique(delta_df$n))) {
  s <- delta_df[delta_df$n == nn, ]
  s <- s[order(s$a), ]
  d_m1 <- s$delta[which.min(abs(s$a - (-1)))]
  d_p1 <- s$delta[which.min(abs(s$a - ( 1)))]
  md_m1 <- s$mcse_diff[which.min(abs(s$a - (-1)))]

  ## (1) ~0 at a=-1
  near_zero_m1 <- is.finite(d_m1) && abs(d_m1) <= 2 * md_m1
  ## (2) rises to > +0.03 at a=+1
  rises_at_p1  <- is.finite(d_p1) && d_p1 > 0.03
  ## (3) monotone non-decreasing in a (allow noise slack of 2*mcse_diff on each step)
  dd   <- diff(s$delta)
  slack<- 2 * s$mcse_diff[-1]
  monotone <- all(is.finite(dd)) && all(dd >= -slack)

  ok <- near_zero_m1 && rises_at_p1 && monotone
  if (!ok) {
    detail <- sprintf(paste0("n=%d alpha=%.2f Delta(a=-1)=%.4f (|.|<=%.4f? %s) ",
                             "Delta(a=+1)=%.4f (>0.03? %s) monotone? %s"),
                      nn, ALPHA_STAR, d_m1, 2 * md_m1, near_zero_m1,
                      d_p1, rises_at_p1, monotone)
    stop_hits[[length(stop_hits) + 1L]] <- detail
  }
}

if (length(stop_hits) > 0) {
  for (d in stop_hits) {
    line <- paste("grid_asymmetry.R", "ASYM_DELTA_NOT_MONOTONE_RISE", d, sep = "\t")
    cat(line, "\n", file = STOP_LOG, append = TRUE, sep = "")
    cat("STOP TRIGGER:", line, "\n")
  }
} else {
  cat("No STOP triggers fired (Delta ~0 at a=-1 and rises monotone to >+0.03 at a=+1, per n).\n")
}

stopCluster(cl)

cat("grid_asymmetry.R DONE. Wrote:\n  ", PV_PATH, "\n  ", SUM_PATH, "\n  ", DELTA_PATH, "\n")
