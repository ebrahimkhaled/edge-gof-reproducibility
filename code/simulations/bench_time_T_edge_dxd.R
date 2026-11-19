## bench_time_T_edge_dxd.R -- block T again, EDGE only: the released code against the code that never forms
## the G x G matrix Omega. Same protocol as bench_time_T.R (design note paper_EDGE/theory/PREDECLARATION_blockT_computation.md):
## its data generator, its n grid for the fast tests, 11 repetitions, the median of the repetitions with the inner loop
## that makes a repetition last about 0.2 s, edge.gof(fit, G = "auto"), one core, nothing else running.
## The two versions share one namespace, so each runs in its own R process, one after the other.
##   EDGE_old = ebrahim.gof 3.1.0 (the CRAN release)          EDGE_dxd = branch edge-fast-omega (d x d and G x p products only)
##
##   Rscript bench_time_T_edge_dxd.R            both versions, writes battery/T/timing_edge_dxd.csv
##   Rscript bench_time_T_edge_dxd.R --child <label> <pkg> <out.csv>     (internal: one version)
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
## The two package sources compared: OLD is ebrahim.gof 3.1.0 from CRAN, NEW is the source in
## code/package/ (the same package with the d x d computation). Unpack both and point these at them.
OLD    <- Sys.getenv("EBGOF_OLD_SRC", edge_path("ebrahim.gof-3.1.0"))
NEW    <- Sys.getenv("EBGOF_NEW_SRC", edge_path("ebrahim.gof-dxd"))
OUT    <- edge_battery("T")
args   <- commandArgs(TRUE)

N_FAST <- c(500L, 1000L, 2000L, 5000L, 20000L, 100000L)
REPS   <- 11L

make_data <- function(n) {                                           # as in bench_time_T.R
  set.seed(20260920 + n)
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(0.6 * x + 0.5 * d))
  dat <- data.frame(y = y, x = x, d = d)
  list(dat = dat, fit = glm(y ~ x + d, family = binomial(), data = dat))
}
time_it <- function(f, reps) {                                       # as in bench_time_T.R
  t1 <- system.time(f())[["elapsed"]]
  K <- if (t1 >= 0.2) 1L else as.integer(min(1000, ceiling(0.2 / max(t1, 1e-4))))
  s <- vapply(seq_len(reps), function(i) {
    t0 <- proc.time()[["elapsed"]]; for (k in seq_len(K)) f(); (proc.time()[["elapsed"]] - t0) / K
  }, numeric(1))
  c(median = stats::median(s), min = min(s), reps = reps, inner = K)
}

if (length(args) >= 4 && args[1] == "--child") {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
  suppressPackageStartupMessages(pkgload::load_all(args[3], quiet = TRUE))
  rows <- list()
  for (n in N_FAST) {
    fit <- make_data(n)$fit
    tt  <- time_it(function() edge.gof(fit, G = "auto"), REPS)
    rows[[length(rows) + 1]] <- data.frame(test = args[2], n = n, median_sec = tt[["median"]], reps = tt[["reps"]])
    cat(sprintf("%-9s n=%-6d median %.4g s (inner loop %d)\n", args[2], n, tt[["median"]], tt[["inner"]]))
  }
  write.csv(do.call(rbind, rows), args[4], row.names = FALSE)
  quit(save = "no")
}

rs  <- file.path(R.home("bin"), "Rscript")
me  <- file.path(SIMDIR, "bench_time_T_edge_dxd.R")
tmp <- tempfile(c("old", "new"), fileext = ".csv")
system2(rs, c(shQuote(me), "--child", "EDGE_old", shQuote(OLD), shQuote(tmp[1])))
system2(rs, c(shQuote(me), "--child", "EDGE_dxd", shQuote(NEW), shQuote(tmp[2])))
R <- rbind(read.csv(tmp[1]), read.csv(tmp[2]))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write.csv(R, file.path(OUT, "timing_edge_dxd.csv"), row.names = FALSE)
print(R)
