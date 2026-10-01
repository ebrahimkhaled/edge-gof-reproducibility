## bench_time_T.R -- block T: what each test costs. Design note: paper_EDGE/theory/PREDECLARATION_blockT_computation.md
## (sha256 85441ae0...), frozen before any timing. One core, nothing else running; one data set per n.
##
##   Rscript bench_time_T.R            the declared grid
##   Rscript bench_time_T.R --smoke    n = 300 only, fast tests and le Cessie (a code check; writes to battery/T_smoke)
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
SIMDIR <- edge_path("code/simulations")
PKG    <- Sys.getenv("EBRAHIM_GOF_SRC")
SMOKE  <- "--smoke" %in% commandArgs(TRUE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
suppressPackageStartupMessages({
  library(ebrahim.gof)                              # the development 2.8.0: edge.gof, cubic.calib.gof, ...  # archive: was load_all() of the author's development tree; install ebrahim.gof (>= 2.9.0) from CRAN instead
  library(data.table)
})
source(file.path(SIMDIR, "_proj_test.R"))
OUT <- edge_battery(if (SMOKE) "T_smoke" else "T")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "timing_log.txt")
lg <- function(...) { s <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", sprintf(...)); cat(s, "\n"); cat(s, "\n", file = LOG, append = TRUE) }

N_FAST <- if (SMOKE) 300L else c(500L, 1000L, 2000L, 5000L, 20000L, 100000L)
N_LC   <- if (SMOKE) 300L else c(500L, 1000L, 2000L, 5000L, 10000L)
N_PROJ <- if (SMOKE) integer(0) else c(500L, 1000L)
N_BAG  <- if (SMOKE) integer(0) else c(500L, 1000L, 2000L)

make_data <- function(n) {
  set.seed(20260920 + n)
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(0.6 * x + 0.5 * d))
  dat <- data.frame(y = y, x = x, d = d)
  list(dat = dat, fit = glm(y ~ x + d, family = binomial(), data = dat))
}

## median seconds per call over `reps` repetitions; a fast call is repeated K times inside each repetition so that
## a repetition lasts about 0.2 s, well above the timer's resolution
time_it <- function(f, reps) {
  t1 <- system.time(f())[["elapsed"]]
  K <- if (t1 >= 0.2) 1L else as.integer(min(1000, ceiling(0.2 / max(t1, 1e-4))))
  s <- vapply(seq_len(reps), function(i) {
    t0 <- proc.time()[["elapsed"]]; for (k in seq_len(K)) f(); (proc.time()[["elapsed"]] - t0) / K
  }, numeric(1))
  c(median = stats::median(s), min = min(s), reps = reps, inner = K)
}
peak_mb <- function(f) {
  invisible(gc(reset = TRUE)); f(); g <- gc()
  sum(g[, ncol(g)])                                                  # "max used (Mb)" over Ncells and Vcells
}

rows <- list()
add <- function(test, n, tt, note = "", mem = NA_real_) {
  rows[[length(rows) + 1]] <<- data.frame(test = test, n = n, median_sec = tt[["median"]], min_sec = tt[["min"]],
                                          reps = tt[["reps"]], inner_loop = tt[["inner"]], peak_mb = mem, note = note)
  lg("%-20s n=%-6d median %.4g s  (reps %d, inner loop %d)%s", test, n, tt[["median"]], tt[["reps"]], tt[["inner"]],
     if (is.finite(mem)) sprintf("  peak memory %.0f MB", mem) else "")
}
lg("block T: one core; grid fast %s | le Cessie %s | projection %s | BAGofT %s%s", paste(N_FAST, collapse = ","),
   paste(N_LC, collapse = ","), paste(N_PROJ, collapse = ","), paste(N_BAG, collapse = ","), if (SMOKE) " [SMOKE]" else "")

## ---- the fast tests, every n ------------------------------------------------------------------------------------------
for (n in N_FAST) {
  D <- make_data(n); fit <- D$fit
  add("EDGE",             n, time_it(function() edge.gof(fit, G = "auto"), 11))
  add("Hosmer-Lemeshow",  n, time_it(function() run.all.gof(fit, tests = "HL", G = 10), 11))
  add("Stukel joint",     n, time_it(function() run.all.gof(fit, tests = "Stukel"), 11))
  add("cubic calib. LR",  n, time_it(function() cubic.calib.gof(fit), 11))
  add("GiViTI",           n, time_it(function() suppressWarnings(givitiR::givitiCalibrationTest(
                                        o = fit$y, e = fitted(fit), devel = "internal")), 11))
}
## ---- le Cessie: n x n kernel ------------------------------------------------------------------------------------------
for (n in N_LC) {
  D <- make_data(n); fit <- D$fit
  f <- function() suppressMessages(run.all.gof(fit, tests = "le-Cessie"))
  mem <- peak_mb(f)
  add("le Cessie", n, time_it(f, if (n <= 2000) 3L else 1L), mem = mem)
}
if (!SMOKE) add("le Cessie", 20000L, c(median = NA, min = NA, reps = 0, inner = 0),
                note = "not run: an n x n kernel at n = 20000 needs an estimated 16 GB or more")
## ---- the Liu projection test and BAGofT --------------------------------------------------------------------------------
for (n in N_PROJ) {
  D <- make_data(n); fit <- D$fit
  add("Liu projection", n, time_it(function() proj_pvalue(fit$y, model.matrix(fit), B = 250), if (n <= 500) 3L else 1L))
}
for (n in N_BAG) {
  D <- make_data(n)
  add("BAGofT", n, time_it(function() suppressWarnings(suppressMessages(BAGofT::BAGofT(
    testModel = BAGofT::testGlmBi(formula = y ~ x + d, link = "logit"), data = D$dat))), 1L))
}
R <- rbindlist(rows)
fwrite(R, file.path(OUT, "timing.csv"))
lg("block T finished: %d rows written to %s", nrow(R), file.path(OUT, "timing.csv"))
