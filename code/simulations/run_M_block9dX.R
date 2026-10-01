## run_M_block9dX.R -- block Z Part B: the robust fit at the paper's own boundary, n = 5000.
## Contract: paper_EDGE/theory/PREDECLARATION_blockZ_signerror_and_robust5000.md (sha256 f5e1a7f1...),
## frozen before any scenario ran.
##
##   Rscript run_M_block9dX.R [--workers 20] [--reps 1000] [--out battery/9dX]
##
## Block 9d's machinery unchanged -- the same two estimators on the same data sets, the same replacement
## of bt_fit, the same battery -- with the grid moved to n = 5000, the logistic truth and k = 25, 50, 100,
## which is where Section 6.5 says the directed test fails. Only b9d_configs() is overridden.
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
options(block9d.source_only = TRUE)
suppressPackageStartupMessages({ library(parallel); library(data.table) })
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))
source(file.path(SIMDIR, "_block9_contam.R"))
source(file.path(SIMDIR, "_block9d_robust.R"))
BAT <- edge_battery()

B9DX_SEED0 <- 650000000
b9d_configs <- function() {
  C <- data.frame(truth = "logit", k = c(0L, 25L, 50L, 100L), stringsAsFactors = FALSE)
  C$config_id <- seq_len(nrow(C))
  C$seed_base <- B9DX_SEED0 + C$config_id * 10000
  C$n <- 5000L
  C
}

OPT <- list(); a <- commandArgs(TRUE); i <- 1L
while (i <= length(a)) { OPT[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
W <- if (is.null(OPT$workers)) 20L else as.integer(OPT$workers)
REPS <- if (is.null(OPT$reps)) 1000L else as.integer(OPT$reps)
OUT <- normalizePath(if (!is.null(OPT$out)) OPT$out else edge_battery("9dX"),
                     winslash = "/", mustWork = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lg <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
                      cat(s, "\n"); cat(s, "\n", file = file.path(OUT, "_progress.log"), append = TRUE) }

K <- b9d_configs()
lg("block 9dX: %d configurations (%d cells), %d workers, %d replicates, out %s",
   nrow(K), 2L * nrow(K), W, REPS, OUT)
cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl))
invisible(clusterCall(cl, function(simdir, seed0) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  SIMDIR <<- simdir
  options(block9d.source_only = TRUE)
  suppressPackageStartupMessages(library(robustbase))
  source(file.path(simdir, "_battery_tests.R")); source(file.path(simdir, "_battery_cells.R"))
  source(file.path(simdir, "_block9_contam.R")); source(file.path(simdir, "_block9d_robust.R"))
  B9DX_SEED0 <<- seed0
  b9d_configs <<- function() {
    C <- data.frame(truth = "logit", k = c(0L, 25L, 50L, 100L), stringsAsFactors = FALSE)
    C$config_id <- seq_len(nrow(C)); C$seed_base <- B9DX_SEED0 + C$config_id * 10000; C$n <- 5000L
    C
  }
  RNGkind("L'Ecuyer-CMRG"); NULL
}, SIMDIR, B9DX_SEED0))

## One replicate draws one data set and fits it twice, so the loop runs over CONFIGURATIONS and writes
## one file per estimator, exactly as block 9d does.
for (i in seq_len(nrow(K))) {
  cfg  <- as.list(K[i, ])
  outs <- file.path(OUT, sprintf("%s_k%03d_%s_pvalues.csv.gz", cfg$truth, cfg$k, c("ML", "robust")))
  if (all(file.exists(outs))) { lg("%s k=%d: present, skipped", cfg$truth, cfg$k); next }
  t0 <- Sys.time()
  rows <- clusterApplyLB(cl, seq_len(REPS), function(r, config) b9d_one(r, config), config = cfg)
  for (j in 1:2) {
    est <- c("ML", "robust")[j]
    M <- as.data.frame(do.call(rbind, lapply(rows, `[[`, est)))
    M <- M[order(M$rep), , drop = FALSE]
    fwrite(M, outs[j], compress = "gzip")
    lg("[%d/%d] %-6s k=%-3d %-6s  EDGE %.3f  HL %.3f  Stukel %.3f  cubic %.3f  %.1f min",
       i, nrow(K), cfg$truth, cfg$k, est,
       mean(M[["EDGE.poly3.u.Grule"]] <= 0.05, na.rm = TRUE),
       mean(M[["HL.Grule"]] <= 0.05, na.rm = TRUE),
       mean(M[["Stk.joint"]] <= 0.05, na.rm = TRUE),
       mean(M[["Cubic.LR"]] <= 0.05, na.rm = TRUE),
       as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  Mr <- as.data.frame(do.call(rbind, lapply(rows, `[[`, "robust")))
  lg("        robust converged %d of %d", sum(Mr$robust_converged %in% 1), REPS)
}
lg("block 9dX: run finished")
