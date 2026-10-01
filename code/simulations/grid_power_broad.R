## grid_power_broad.R — EXPERIMENT 2.1 (severity x n power surfaces) + 2.0 size gate
## + 2.7 G-sensitivity of power (SPECKIT §5 Exp 2.1 / 2.0 / 2.7).
##
## THE HEADLINE POWER STUDY (A1). Curated ~10-test rival set only (A3). NO combined-basis
## column (A5). Pre-specified default basis = DEF.poly3 (headline EDGE).
##
## Contract: source _harness.R then source(DGP). Uses ek_cluster/ek_run_cell/ek_append.
## one_rep(rep,cell): d<-dgp_alt(family,param,n); fit<-glm(d$f,data=d$d,binomial); run_curated_full(fit,G).
## Store per-rep p.<test> AND s.<test> (wide) -> reshaped to long sim_power_broad_pvalues.csv.
##
## Size gate (2.0): every power row carries
##   size_at_cell    = this test's realized null_size at matching (family,n,G,alpha) from sim_null.csv
##   power_size_adj  = reject rate at the per-test empirical (1-alpha) quantile of s.<test> from
##                     sim_null_pvalues.csv for the matching (family,n,G) (upper-tail convention).
## DEPENDS ON sim_null.csv + sim_null_pvalues.csv (grid_null.R). If absent -> STOP (spec).
##
## Outputs:
##   sim_power_broad_pvalues.csv  (per rep,test):  family,param,n,G,rep,seed,test,statistic,p_value
##   sim_power_broad.csv          (summary):       family,param,n,G,alpha,test,reject_rate,
##                                                 power_size_adj,size_at_cell,B,mcse
##   sim_g_sensitivity.csv        (2.7 sweep):     family,param,n,G,test,reject_rate,mcse
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
source(DGP)                                # dgp_alt, run_curated_full, EK_CURATED

## ---------------------------------------------------------------------------
## Paths & constants
## ---------------------------------------------------------------------------
PV_PATH    <- file.path(SIMDIR, "sim_power_broad_pvalues.csv")
SUM_PATH   <- file.path(SIMDIR, "sim_power_broad.csv")
GSENS_PATH <- file.path(SIMDIR, "sim_g_sensitivity.csv")
NULL_SUM   <- file.path(SIMDIR, "sim_null.csv")
NULL_PV    <- file.path(SIMDIR, "sim_null_pvalues.csv")
STOP_LOG   <- file.path(SIMDIR, "_STOP_TRIGGERS.log")

GRID_SEED   <- 20260708L
G_DEFAULT   <- 10L
ALPHAS      <- c(0.01, 0.05, 0.10)
B_BASE      <- 5000L    # B at every cell ...
B_SMALLN    <- 10000L   # ... except n=200 rows of EVERY family (spec 2.1 B rule)
N_SMALL     <- 200L

## G-sensitivity (2.7): sweep G on {link=cloglog, quad J=0.05} x n=1000.
GSENS_SEED  <- GRID_SEED           # shares the power stream but distinct cell_ids -> distinct seed_base
GSENS_GVALS <- c(6L, 8L, 10L, 12L, 14L, 20L)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

ek_stop <- function(trigger, detail) {
  line <- paste("grid_power_broad.R", trigger, detail, sep = "\t")
  cat(line, "\n", file = STOP_LOG, append = TRUE, sep = "")
  cat("STOP TRIGGER:", line, "\n")
}

## ---------------------------------------------------------------------------
## Cell grid — 2.1 severity x n surfaces (SPECKIT §5 Exp 2.1 table)
##   Each cell = (family, param, n). Distinct seed_base = GRID_SEED + cell_id*1e5.
## ---------------------------------------------------------------------------
N_LADDER <- c(200L, 500L, 1000L, 2000L, 5000L)

build_family_cells <- function(family, params, n_grid) {
  expand.grid(family = family, param = as.character(params), n = n_grid, G = G_DEFAULT,
              stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
}

cells_quad    <- build_family_cells("quad",    c(0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.40), N_LADDER)
cells_binint  <- build_family_cells("binint",  c(0.1, 0.2, 0.3, 0.5, 0.7),                   N_LADDER)
cells_contint <- build_family_cells("contint", c(0.1, 0.3, 0.5, 0.7),                        N_LADDER)
cells_link    <- build_family_cells("link",    c("probit", "cloglog", "stukel_heavy", "stukel_light", "stukel_asym"), N_LADDER)
cells_rough   <- build_family_cells("rough",   c("osc2", "osc4", "sawtooth", "bump"),         c(500L, 1000L, 2000L))

cells <- rbind(cells_quad, cells_binint, cells_contint, cells_link, cells_rough)

## 2.7 G-sensitivity cells appended AFTER the 2.1 cells so their cell_ids are stable & disjoint.
gsens_cells <- rbind(
  build_family_cells("link", "cloglog", 1000L)[rep(1, length(GSENS_GVALS)), , drop = FALSE],
  build_family_cells("quad", "0.05",    1000L)[rep(1, length(GSENS_GVALS)), , drop = FALSE]
)
gsens_cells$G <- rep(GSENS_GVALS, times = 2)
gsens_cells$is_gsens <- TRUE

cells$is_gsens <- FALSE
all_cells <- rbind(cells, gsens_cells)
all_cells$cell_id <- seq_len(nrow(all_cells))   # fixed on the FULL grid => seeds stable under subsetting

## per-cell B: 10000 at n=200 (2.1 rows only), else 5000. G-sensitivity uses B_BASE.
cell_B <- function(cell) {
  if (isTRUE(cell$is_gsens)) return(B_BASE)
  if (cell$n == N_SMALL) B_SMALLN else B_BASE
}

## Optional self-test subset: EK_POWER_CELLS="1,50" runs only those stable cell_ids.
.subset <- Sys.getenv("EK_POWER_CELLS", "")
if (nzchar(.subset)) {
  keep <- as.integer(strsplit(.subset, ",")[[1]])
  all_cells <- all_cells[all_cells$cell_id %in% keep, , drop = FALSE]
  cat("SUBSET run: cell_ids", paste(keep, collapse = ","), "\n")
}

## ---------------------------------------------------------------------------
## one_rep: returns a NAMED numeric vector (identical names every call).
##   p.<10 tests>, s.<10 tests> from run_curated_full (curated order = EK_CURATED).
## ---------------------------------------------------------------------------
one_rep <- function(rep, cell) {
  fam <- cell$family
  ## quad/binint/contint take numeric severities; link/rough take scenario names.
  param <- if (fam %in% c("quad", "binint", "contint")) as.numeric(cell$param) else cell$param
  d   <- dgp_alt(fam, param, cell$n)
  fit <- suppressWarnings(glm(d$f, data = d$d, family = binomial()))
  run_curated_full(fit, G = cell$G)   # named: p.<test> then s.<test>
}

## ---------------------------------------------------------------------------
## Size-gate helper tables from grid_null.R outputs.
##   sim_null.csv        : family,n,G,test,alpha,null_size,mcse,size_ok
##   sim_null_pvalues.csv: family,n,G,rep,seed, p.<test>..., s.<test>...  (WIDE)
## We precompute, per (family,n,G,test,alpha), the empirical (1-alpha) upper-tail critical
## value of the NULL statistic s.<test>, so size-adjusted power = mean(s_power > crit).
## ---------------------------------------------------------------------------
if (!file.exists(NULL_SUM) || !file.exists(NULL_PV)) {
  ek_stop("MISSING_NULL_INPUT",
          sprintf("sim_null.csv exists=%s sim_null_pvalues.csv exists=%s (run grid_null.R first)",
                  file.exists(NULL_SUM), file.exists(NULL_PV)))
  stop("grid_power_broad.R: sim_null.csv / sim_null_pvalues.csv missing — run grid_null.R first.")
}

null_sum <- utils::read.csv(NULL_SUM, stringsAsFactors = FALSE)
null_pv  <- utils::read.csv(NULL_PV, stringsAsFactors = FALSE, check.names = FALSE)

## size_at_cell lookup: named by "family|n|G|test|alpha" -> null_size
.size_key <- with(null_sum, paste(family, n, G, test, alpha, sep = "|"))
SIZE_AT <- setNames(null_sum$null_size, .size_key)

## Empirical upper-tail critical value for size-adjusted power.
## crit(family,n,G,test,alpha) = (1-alpha) quantile of NULL s.<test> for that (family,n,G).
null_crit <- function(family, n, G, test, alpha) {
  scol <- paste0("s.", test)
  if (!scol %in% names(null_pv)) return(NA_real_)
  sub <- null_pv[[scol]][null_pv$family == family & null_pv$n == n & null_pv$G == G]
  sub <- sub[is.finite(sub)]
  if (length(sub) < 10) return(NA_real_)
  as.numeric(stats::quantile(sub, probs = 1 - alpha, type = 7, names = FALSE))
}

## Cache critical values (unique family,n,G,test,alpha) so we don't recompute per row.
CRIT_CACHE <- new.env(parent = emptyenv())
get_crit <- function(family, n, G, test, alpha) {
  key <- paste(family, n, G, test, alpha, sep = "|")
  if (!is.null(CRIT_CACHE[[key]])) return(CRIT_CACHE[[key]])
  v <- null_crit(family, n, G, test, alpha)
  CRIT_CACHE[[key]] <- v
  v
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------
cat(sprintf("grid_power_broad.R: %d cells (%d power + %d g-sens), REPS_SCALE=%s, EK_NCORES=%d\n",
            nrow(all_cells), sum(!all_cells$is_gsens), sum(all_cells$is_gsens),
            Sys.getenv("REPS_SCALE", "1"), EK_NCORES))

## Fresh outputs each run (this script fully regenerates its CSVs).
for (p in c(PV_PATH, SUM_PATH, GSENS_PATH)) if (file.exists(p)) file.remove(p)

cl <- ek_cluster(GRID_SEED)
on.exit(stopCluster(cl), add = TRUE)

pcols <- paste0("p.", EK_CURATED)
scols <- paste0("s.", EK_CURATED)

for (i in seq_len(nrow(all_cells))) {
  cell      <- as.list(all_cells[i, ])
  seed_base <- GRID_SEED + cell$cell_id * 1e5
  B         <- ek_reps(cell_B(cell))
  t0        <- proc.time()[["elapsed"]]

  df <- ek_run_cell(cl, B, seed_base, one_rep, cell)   # B rows x (p.* , s.*)
  reps <- seq_len(nrow(df))

  if (isTRUE(cell$is_gsens)) {
    ## ---- 2.7 G-sensitivity: reject_rate + mcse per test at nominal alpha=0.05 ----
    gs_rows <- list()
    for (tst in EK_CURATED) {
      pval <- df[[paste0("p.", tst)]]
      rr   <- if (sum(!is.na(pval)) > 0) mean(pval < 0.05, na.rm = TRUE) else NA_real_
      gs_rows[[length(gs_rows) + 1L]] <- data.frame(
        family = cell$family, param = cell$param, n = cell$n, G = cell$G, test = tst,
        reject_rate = rr, mcse = ek_mcse(rr, B), stringsAsFactors = FALSE)
    }
    ek_append(do.call(rbind, gs_rows), GSENS_PATH)

    el <- proc.time()[["elapsed"]] - t0
    cat(sprintf("  cell %d/%d [g-sens] family=%s param=%s n=%d G=%d  (%.1fs)\n",
                i, nrow(all_cells), cell$family, cell$param, cell$n, cell$G, el))
    next
  }

  ## ---- per-rep p-value/statistic CSV (LONG: one row per rep x test) ----
  pv_list <- vector("list", length(EK_CURATED))
  for (j in seq_along(EK_CURATED)) {
    tst <- EK_CURATED[j]
    pv_list[[j]] <- data.frame(
      family = cell$family, param = cell$param, n = cell$n, G = cell$G,
      rep = reps, seed = seed_base + reps, test = tst,
      statistic = df[[paste0("s.", tst)]], p_value = df[[paste0("p.", tst)]],
      stringsAsFactors = FALSE)
  }
  ek_append(do.call(rbind, pv_list), PV_PATH)

  ## ---- summary rollup per (family,param,n,G,alpha,test) ----
  sum_rows <- list()
  for (tst in EK_CURATED) {
    pval <- df[[paste0("p.", tst)]]
    sval <- df[[paste0("s.", tst)]]
    nval <- sum(!is.na(pval))
    for (a in ALPHAS) {
      reject_rate <- if (nval > 0) mean(pval < a, na.rm = TRUE) else NA_real_

      ## size_at_cell: matching null_size from sim_null.csv (family,n,G,test,alpha)
      skey <- paste(cell$family, cell$n, cell$G, tst, a, sep = "|")
      size_at_cell <- unname(SIZE_AT[skey]) %||% NA_real_

      ## power_size_adj: reject at the empirical (1-alpha) upper-tail NULL critical value
      crit <- get_crit(cell$family, cell$n, cell$G, tst, a)
      psa  <- if (is.finite(crit)) {
        sf <- sval[is.finite(sval)]
        if (length(sf) > 0) mean(sf > crit) else NA_real_
      } else NA_real_

      sum_rows[[length(sum_rows) + 1L]] <- data.frame(
        family = cell$family, param = cell$param, n = cell$n, G = cell$G, alpha = a,
        test = tst, reject_rate = reject_rate, power_size_adj = psa,
        size_at_cell = size_at_cell, B = B, mcse = ek_mcse(reject_rate, B),
        stringsAsFactors = FALSE)
    }
  }
  ek_append(do.call(rbind, sum_rows), SUM_PATH)

  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("  cell %d/%d done: family=%s param=%s n=%d G=%d B=%d  (%.1fs)\n",
              i, nrow(all_cells), cell$family, cell$param, cell$n, cell$G, B, el))
}

cat("grid_power_broad.R DONE. Wrote:\n  ", PV_PATH, "\n  ", SUM_PATH, "\n  ", GSENS_PATH, "\n")
