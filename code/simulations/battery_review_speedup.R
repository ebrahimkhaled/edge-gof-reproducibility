## battery_review_speedup.R -- parallel throughput of battery_one() on warm workers (review item f), within the review
## caps: at most 8 workers and at most 50 replicates per cell. Four cells x 50 replicates at n = 5000 and at n = 1000,
## run serially in this process and on 4 and 8 socket workers with one chunk per worker, BLAS on one thread.
## Writes battery/_review/review_speedup.log.
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
OUT <- edge_battery("_review")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
suppressPackageStartupMessages({ library(parallel); library(data.table) })
source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
sink(file.path(OUT, "review_speedup.log"), split = TRUE)
cat("battery_review_speedup.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
RNGkind("L'Ecuyer-CMRG")
Cells <- battery_cells()
SETS <- list(n5000 = c("null_quad_n5000", "null_binint_n5000", "null_contint_n5000", "null_link_n5000"),
             n1000 = c("null_binint_n1000", "null_contint_n1000", "null_rough_n1000", "null_quad_n500"))
jobs_of <- function(set) do.call(c, lapply(SETS[[set]], function(cn) {
  ce <- as.list(Cells[Cells$block == "1a" & Cells$cell == cn, ]); lapply(1:50, function(r) list(r = r, ce = ce)) }))
run_jobs <- function(js) for (j in js) battery_one(j$r, j$ce)

for (set in names(SETS)) {
  js <- jobs_of(set)
  invisible(battery_one(99L, js[[1]]$ce))
  ts <- system.time(run_jobs(js))[["elapsed"]]
  cat(sprintf("\n%s: %d replicates serially in %.1f s (%.3f s/rep)\n", set, length(js), ts, ts / length(js)))
  for (w in c(4L, 8L)) {
    cl <- makeCluster(w)
    parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
    clusterCall(cl, function(f) { Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
      source(f); RNGkind("L'Ecuyer-CMRG"); NULL }, file.path(SIMDIR, "_battery_tests.R"))
    invisible(clusterCall(cl, function(ce) { battery_one(98L, ce); NULL }, js[[1]]$ce))          # warm every worker
    chunks <- splitIndices(length(js), w)
    tp <- system.time(parLapply(cl, chunks, function(ix, js) { for (i in ix) battery_one(js[[i]]$r, js[[i]]$ce); NULL }, js = js))[["elapsed"]]
    ## the driver's own chunking: W * 4 chunks, load-balanced, results returned
    chunks2 <- splitIndices(length(js), w * 4L)
    tp2 <- system.time(parLapplyLB(cl, chunks2, function(ix, js) do.call(rbind, lapply(ix, function(i) battery_one(js[[i]]$r, js[[i]]$ce))), js = js))[["elapsed"]]
    stopCluster(cl)
    cat(sprintf("  %d workers: %.2f s, speed-up %.2f (efficiency %.2f); with the driver's %d load-balanced chunks: %.2f s, speed-up %.2f\n",
                w, tp, ts / tp, ts / tp / w, w * 4L, tp2, ts / tp2))
  }
}
sink()
