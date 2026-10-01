## run_M_block9R3.R -- block 9R3: does the tolerance rule of EDGE-FR hold at 2% and 4% corruption?
## Contract: paper_EDGE/theory/PREDECLARATION_block9R3_tolerance.md (sha256 421217c0...), frozen before any scenario ran.
##
##   Rscript run_M_block9R3.R [--workers 20] [--reps N] [--out battery/9R3]
##
## Fresh seeds, disjoint from every other block (390000000 + j * 10000), because these corruption rates are outside
## block 9's grid and there is nothing stored to regenerate. The variants are block 9R2's, unchanged.
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
suppressPackageStartupMessages({ library(parallel); library(data.table) })
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))
source(file.path(SIMDIR, "_block9_contam.R"))
source(file.path(SIMDIR, "_block9R2_hybrid.R"))
BAT <- edge_battery()

r93_cells <- function() {
  P <- data.frame(
    cell       = c("fr_clean_n1000", "fr_C1_k020_n1000", "fr_C1_k040_n1000",
                   "fr_clean_n5000", "fr_C1_k100_n5000", "fr_C1_k200_n5000"),
    n          = c(1000L, 1000L, 1000L, 5000L, 5000L, 5000L),
    corruption = c("clean", "C1", "C1", "clean", "C1", "C1"),
    rate       = c(0, 0.02, 0.04, 0, 0.02, 0.04),
    stringsAsFactors = FALSE)
  P$seed_base <- 390000000 + seq_len(nrow(P)) * 10000
  P$reps <- 1000L
  P
}

r93_one <- function(rep, ce) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  set.seed(ce$seed_base + rep)
  g <- b9_gen("logit", ce$n, ce$corruption, ce$rate)
  d <- g$d; ev <- sum(d$y)
  fq <- if (min(ev, nrow(d) - ev) == 0) NULL else bt_fit(list(d = d, f = g$f))
  c(rep = rep, n = nrow(d), events = ev, k = length(g$corrupt), fr2_pvalues(fq))
}
r93_task <- function(rep, ce) r93_one(rep, ce)

OPT <- list(); a <- commandArgs(TRUE); i <- 1L
while (i <= length(a)) { OPT[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
W <- if (is.null(OPT$workers)) 20L else as.integer(OPT$workers)
REPS <- if (is.null(OPT$reps)) NA_integer_ else as.integer(OPT$reps)
OUT <- normalizePath(if (!is.null(OPT$out)) OPT$out else edge_battery(if (is.na(REPS)) "9R3" else "9R3_test"),
                     winslash = "/", mustWork = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lg <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
                      cat(s, "\n"); cat(s, "\n", file = file.path(OUT, "_progress.log"), append = TRUE) }

P <- r93_cells()
lg("block 9R3: %d scenarios, %d workers, out %s", nrow(P), W, OUT)
cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl))
invisible(clusterCall(cl, function(simdir) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  SIMDIR <<- simdir
  source(file.path(simdir, "_battery_tests.R")); source(file.path(simdir, "_battery_cells.R"))
  source(file.path(simdir, "_block9_contam.R")); source(file.path(simdir, "_block9R2_hybrid.R"))
  RNGkind("L'Ecuyer-CMRG"); NULL
}, SIMDIR))
clusterExport(cl, c("r93_one", "r93_task"))

for (i in seq_len(nrow(P))) {
  ce <- as.list(P[i, ]); R <- if (is.na(REPS)) ce$reps else min(REPS, ce$reps)
  f <- file.path(OUT, sprintf("%s_variants.csv.gz", ce$cell))
  if (file.exists(f)) { lg("%s: present, skipped", ce$cell); next }
  t0 <- Sys.time()
  M <- as.data.table(do.call(rbind, clusterApplyLB(cl, seq_len(R), r93_task, ce = ce)))
  fwrite(M, f, compress = "gzip")
  lg("[%d/%d] %-18s n=%-5d k=%-4d R=%d  V0 %.3f HYB05 %.3f HYB10 %.3f HYB20 %.3f G10 %.3f  %.1f min",
     i, nrow(P), ce$cell, ce$n, as.integer(M$k[1]), R,
     mean(M$V0 <= 0.05, na.rm = TRUE), mean(M$HYB05 <= 0.05, na.rm = TRUE), mean(M$HYB10 <= 0.05, na.rm = TRUE),
     mean(M$HYB20 <= 0.05, na.rm = TRUE), mean(M$G10 <= 0.05, na.rm = TRUE),
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
lg("block 9R3: run finished")
