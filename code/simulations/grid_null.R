## grid_null.R — EXPERIMENT 0: the NULL GRID / size gate (SPECKIT §5 Exp 0).
## Correctly-specified null of each power family => realized rejection == realized SIZE.
## Produces sim_null_pvalues.csv (per-rep p AND statistic) + sim_null.csv (size summary
## per family,n,G,test,alpha) that §2.0's size-before-power gate consumes.
##
## Contract: source _harness.R then source(DGP). Uses ek_cluster/ek_run_cell/ek_append.
## one_rep(rep,cell): d<-dgp_null(family,n); fit<-glm(d$f,data=d$d,binomial); run_curated_full(fit,G).
## Store per-rep p.<test> and s.<test> so size-adjusted power / empirical-alpha quantiles recompute.
## STOP: any DEF.poly2/DEF.poly3/DEF.stukel/EF with null_size outside the symmetric +/-3*mcse_nom
##       band at G=10 & n>=1000 -> append to _STOP_TRIGGERS.log + cat(). PH conservative => flag not STOP.
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
source(file.path(SIMDIR, "_harness.R"))   # defines SIMDIR, DGP, EK_NCORES, ek_reps, ek_cluster, ek_run_cell, ek_append, ek_mcse
source(DGP)                                # dgp_null, run_curated_full, EK_CURATED (also loaded in workers by ek_cluster)

## ---------------------------------------------------------------------------
## Paths & constants
## ---------------------------------------------------------------------------
PV_PATH    <- file.path(SIMDIR, "sim_null_pvalues.csv")
SUM_PATH   <- file.path(SIMDIR, "sim_null.csv")
STOP_LOG   <- file.path(SIMDIR, "_STOP_TRIGGERS.log")
STALE_TEST <- file.path(SIMDIR, "_test_null_pvalues.csv")

GRID_SEED <- 20260705L
B_NULL    <- 10000L
ALPHAS    <- c(0.01, 0.05, 0.10)

## STOP-trigger set (mis-calibrated null here invalidates every downstream power claim)
STOP_TESTS <- c("DEF.poly2", "DEF.poly3", "DEF.stukel", "EF")

## Clean up any stale test file from a previous smoke run (spec instruction).
if (file.exists(STALE_TEST)) { file.remove(STALE_TEST); cat("Removed stale", STALE_TEST, "\n") }

## ---------------------------------------------------------------------------
## Cell grid
##   Main gate: families x n in {200,500,1000,2000,5000} at G=10.
##   G-sensitivity size block: families {link,quad} x n=1000 x G in {6,8,12,14,20}.
## Distinct seed_base per cell = GRID_SEED + cell_id*1e5 (spec).
## ---------------------------------------------------------------------------
FAMILIES <- c("quad", "binint", "contint", "link", "rough")
N_MAIN   <- c(200L, 500L, 1000L, 2000L, 5000L)

main_grid <- expand.grid(family = FAMILIES, n = N_MAIN, G = 10L,
                         stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

gsens_grid <- expand.grid(family = c("link", "quad"), n = 1000L, G = c(6L, 8L, 12L, 14L, 20L),
                          stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

cells <- rbind(main_grid, gsens_grid)
cells$cell_id <- seq_len(nrow(cells))   # cell_id fixed on the FULL grid => seeds stable if subset

## Optional self-test subset: EK_NULL_CELLS="1,6" runs only those (stable) cell_ids.
.subset <- Sys.getenv("EK_NULL_CELLS", "")
if (nzchar(.subset)) {
  keep <- as.integer(strsplit(.subset, ",")[[1]])
  cells <- cells[cells$cell_id %in% keep, , drop = FALSE]
  cat("SUBSET run: cell_ids", paste(keep, collapse = ","), "\n")
}

## ---------------------------------------------------------------------------
## one_rep: returns a NAMED numeric vector (identical names every call).
##   p.<10 tests>, s.<10 tests> from run_curated_full (curated order = EK_CURATED).
## ---------------------------------------------------------------------------
one_rep <- function(rep, cell) {
  d   <- dgp_null(cell$family, cell$n)
  fit <- suppressWarnings(glm(d$f, data = d$d, family = binomial()))
  run_curated_full(fit, G = cell$G)   # named: p.<test> then s.<test>
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------
B <- ek_reps(B_NULL)
cat(sprintf("grid_null.R: %d cells, B=%d (REPS_SCALE=%s), EK_NCORES=%d\n",
            nrow(cells), B, Sys.getenv("REPS_SCALE", "1"), EK_NCORES))

## Fresh outputs each run (this script fully regenerates the two CSVs).
if (file.exists(PV_PATH))  file.remove(PV_PATH)
if (file.exists(SUM_PATH)) file.remove(SUM_PATH)

cl <- ek_cluster(GRID_SEED)
on.exit(stopCluster(cl), add = TRUE)

pcols <- paste0("p.", EK_CURATED)
scols <- paste0("s.", EK_CURATED)

stop_hits <- list()

for (i in seq_len(nrow(cells))) {
  cell      <- as.list(cells[i, ])
  seed_base <- GRID_SEED + cell$cell_id * 1e5
  t0 <- proc.time()[["elapsed"]]

  df <- ek_run_cell(cl, B, seed_base, one_rep, cell)   # B rows x (p.* , s.*)

  ## ---- per-rep p-value/statistic CSV (wide) ----
  reps <- seq_len(nrow(df))
  pv <- data.frame(family = cell$family, n = cell$n, G = cell$G,
                   rep = reps, seed = seed_base + reps,
                   check.names = FALSE, stringsAsFactors = FALSE)
  pv <- cbind(pv, df[, pcols, drop = FALSE], df[, scols, drop = FALSE])
  ek_append(pv, PV_PATH)

  ## ---- summary rollup per (family,n,G,test,alpha) ----
  sum_rows <- list()
  for (tst in EK_CURATED) {
    pcol <- paste0("p.", tst)
    pval <- df[[pcol]]
    nval <- sum(!is.na(pval))
    for (a in ALPHAS) {
      null_size <- if (nval > 0) mean(pval < a, na.rm = TRUE) else NA_real_
      mcse      <- ek_mcse(null_size, B)                 # SE of the realized rate
      mcse_nom  <- ek_mcse(a, B)                          # nominal-alpha SE => fixed band
      size_ok   <- !is.na(null_size) &&
                   null_size >= (a - 3 * mcse_nom) &&
                   null_size <= (a + 3 * mcse_nom)        # symmetric +/-3 band
      sum_rows[[length(sum_rows) + 1L]] <- data.frame(
        family = cell$family, n = cell$n, G = cell$G, test = tst, alpha = a,
        null_size = null_size, mcse = mcse, size_ok = size_ok,
        stringsAsFactors = FALSE)

      ## ---- STOP trigger check: DEF.*/EF outside band at G=10 & n>=1000 ----
      if (tst %in% STOP_TESTS && cell$G == 10L && cell$n >= 1000L && !size_ok && !is.na(null_size)) {
        detail <- sprintf("test=%s family=%s n=%d G=%d alpha=%.2f null_size=%.4f band=[%.4f,%.4f]",
                          tst, cell$family, cell$n, cell$G, a, null_size,
                          a - 3 * mcse_nom, a + 3 * mcse_nom)
        stop_hits[[length(stop_hits) + 1L]] <- detail
      }
    }
  }
  ek_append(do.call(rbind, sum_rows), SUM_PATH)

  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("  cell %d/%d done: family=%s n=%d G=%d  (%.1fs)\n",
              i, nrow(cells), cell$family, cell$n, cell$G, el))
}

## ---------------------------------------------------------------------------
## Honor STOP triggers: log + cat (never alter numbers to avoid a trigger).
## ---------------------------------------------------------------------------
if (length(stop_hits) > 0) {
  for (d in stop_hits) {
    line <- paste("grid_null.R", "NULL_SIZE_OUT_OF_BAND", d, sep = "\t")
    cat(line, "\n", file = STOP_LOG, append = TRUE, sep = "")
    cat("STOP TRIGGER:", line, "\n")
  }
} else {
  cat("No STOP triggers fired (DEF.*/EF null sizes within +/-3*mcse_nom band at G=10, n>=1000).\n")
}

## PH conservative-flag (not a STOP): report PH's realized null size at G=10.
ph_flag <- tryCatch({
  s <- utils::read.csv(SUM_PATH, stringsAsFactors = FALSE)
  ph <- s[s$test == "Pigeon-Heyse" & s$G == 10L & s$alpha == 0.05, ]
  if (nrow(ph) > 0) cat(sprintf("PH caveat (flag, not STOP): mean PH null_size @alpha=.05,G=10 = %.4f\n",
                                mean(ph$null_size, na.rm = TRUE)))
  invisible(TRUE)
}, error = function(e) invisible(FALSE))

cat("grid_null.R DONE. Wrote:\n  ", PV_PATH, "\n  ", SUM_PATH, "\n")
