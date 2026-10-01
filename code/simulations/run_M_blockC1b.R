## run_M_blockC1b.R -- block Z Part A: a sign error in one covariate.
## Contract: paper_EDGE/theory/PREDECLARATION_blockZ_signerror_and_robust5000.md (sha256 f5e1a7f1...),
## frozen before any scenario ran.
##
##   Rscript run_M_blockC1b.R [--workers 20] [--reps 1000] [--out battery/C1b]
##
## Ten cells, the battery's own tests, fresh seeds. The generator is checked against block 9's before
## anything is computed: with the corruption set to C1 the two must agree exactly.
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
source(file.path(SIMDIR, "_blockC1b_contam.R"))
BAT <- edge_battery()

OPT <- list(); a <- commandArgs(TRUE); i <- 1L
while (i <= length(a)) { OPT[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
W <- if (is.null(OPT$workers)) 20L else as.integer(OPT$workers)
REPS <- if (is.null(OPT$reps)) 1000L else as.integer(OPT$reps)
OUT <- normalizePath(if (!is.null(OPT$out)) OPT$out else edge_battery("C1b"),
                     winslash = "/", mustWork = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lg <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
                      cat(s, "\n"); cat(s, "\n", file = file.path(OUT, "_progress.log"), append = TRUE) }

c1b_identity_check()
lg("block C1b: generator identity gate passed (C1 reproduces b9_gen exactly)")

C <- c1b_cells()
lg("block C1b: %d cells, %d workers, %d replicates, out %s", nrow(C), W, REPS, OUT)
cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl))
invisible(clusterCall(cl, function(simdir) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  SIMDIR <<- simdir
  source(file.path(simdir, "_battery_tests.R")); source(file.path(simdir, "_battery_cells.R"))
  source(file.path(simdir, "_blockC1b_contam.R"))
  RNGkind("L'Ecuyer-CMRG"); NULL
}, SIMDIR))
clusterExport(cl, "c1b_one")

for (i in seq_len(nrow(C))) {
  ce <- as.list(C[i, ])
  f <- file.path(OUT, sprintf("%s_pvalues.csv.gz", ce$cell))
  if (file.exists(f)) { lg("%s: present, skipped", ce$cell); next }
  t0 <- Sys.time()
  M <- as.data.table(do.call(rbind, clusterApplyLB(cl, seq_len(REPS), function(r, cell) c1b_one(r, cell),
                                                   cell = ce)))
  fwrite(M, f, compress = "gzip")
  ed <- mean(M[["EDGE.poly3.u.Grule"]] <= 0.05, na.rm = TRUE)
  sk <- mean(M[["Stk.joint"]] <= 0.05, na.rm = TRUE)
  cu <- mean(M[["Cubic.LR"]] <= 0.05, na.rm = TRUE)
  gv <- mean(M[["GiViTI"]] <= 0.05, na.rm = TRUE)
  lg("[%2d/%d] %-18s n=%-5d k=%-3d  EDGE %.3f  Stukel %.3f  cubic %.3f  GiViTI %.3f  %.1f min",
     i, nrow(C), ce$cell, ce$n, ce$k, ed, sk, cu, gv,
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
lg("block C1b: run finished")
