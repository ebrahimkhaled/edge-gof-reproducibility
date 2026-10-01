## run_M_block9d.R -- block 9d: does robust estimation solve the problem instead?
## Contract: paper_EDGE/theory/PREDECLARATION_block9d_robust.md (sha256 799a85aa...).
##
##   Rscript run_M_block9d.R [--workers 14] [--reps N] [--out battery/9d]
##
## 8 data configurations x 1000 replicates. Each replicate is one data set fitted twice, by maximum
## likelihood and by robustbase::glmrob (Mqle), so both estimators see identical data; the rows are
## written to 16 per-cell files named <truth>_k<k>_<estimator>_pvalues.csv.gz. Gated by
## block9d_selftest.R (42 checks). Refuses to write into any finished block's folder.
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
suppressMessages({ library(parallel); library(data.table) })
SIMDIR <- edge_path("code/simulations")
ROOT   <- edge_battery()
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_block9_contam.R"))
source(file.path(SIMDIR, "_block9d_robust.R"))

opt <- list(); a <- commandArgs(TRUE); i <- 1L
while (i <= length(a)) { opt[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
W    <- if (is.null(opt$workers)) 14L else as.integer(opt$workers)
REPS <- if (is.null(opt$reps)) 1000L else as.integer(opt$reps)
OUT  <- normalizePath(if (is.null(opt$out)) edge_battery("9d") else opt$out, winslash = "/", mustWork = FALSE)
busy <- normalizePath(edge_battery(c(as.character(0:9), "9c", "dryrun", "analysis")), winslash = "/", mustWork = FALSE)
if (any(startsWith(paste0(OUT, "/"), paste0(busy, "/")))) stop("--out cannot be a finished block's folder: ", OUT)
if (!is.null(opt$reps) && identical(OUT, normalizePath(edge_battery("9d"), winslash = "/", mustWork = FALSE)))
  stop("a test run (--reps) cannot write to the real block 9d folder; pass --out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "_progress.log")
say <- function(m) { l <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "", m)
                     cat(l, "\n", sep = ""); cat(l, "\n", sep = "", file = LOG, append = TRUE) }

K <- b9d_configs()
say(sprintf("block 9d: %d configurations x %d replicates, 2 estimators each, %d workers, out %s",
            nrow(K), REPS, W, OUT))

cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, RNGkind("L'Ecuyer-CMRG")))
clusterExport(cl, "SIMDIR", envir = environment())
invisible(clusterEvalQ(cl, { suppressMessages(requireNamespace("robustbase"))
                   source(file.path(SIMDIR, "_battery_tests.R"))
                   source(file.path(SIMDIR, "_block9_contam.R"))
                   source(file.path(SIMDIR, "_block9d_robust.R")); TRUE }))

t0 <- Sys.time()
for (i in seq_len(nrow(K))) {
  cfg  <- as.list(K[i, ])
  outs <- file.path(OUT, sprintf("%s_k%02d_%s_pvalues.csv.gz", cfg$truth, cfg$k, c("ML", "robust")))
  if (all(file.exists(outs))) { say(sprintf("%s k=%d: finished files present, skipped", cfg$truth, cfg$k)); next }
  ## the extra argument must not be a prefix of clusterApplyLB's own formals (cl, x, fun): R's
  ## partial matching read an argument named 'c' as 'cl' and replaced the cluster itself
  rows <- clusterApplyLB(cl, seq_len(REPS), function(r, config) b9d_one(r, config), config = cfg)
  for (j in 1:2) {
    est <- c("ML", "robust")[j]
    M <- as.data.frame(do.call(rbind, lapply(rows, `[[`, est)))
    M <- M[order(M$rep), , drop = FALSE]
    tmp <- paste0(outs[j], ".tmp"); fwrite(M, tmp, compress = "gzip")
    if (file.exists(outs[j])) file.remove(outs[j])
    if (!file.rename(tmp, outs[j])) stop("could not rename ", tmp)
  }
  Mr <- as.data.frame(do.call(rbind, lapply(rows, `[[`, "robust")))
  say(sprintf("[%d/%d] %-6s k=%-2d written both estimators; robust converged %d of %d, errors %d  elapsed %.1f min",
              i, nrow(K), cfg$truth, cfg$k, sum(Mr$robust_converged %in% 1), REPS,
              sum(Mr$flag.b9d_error %in% 1), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
say("block 9d finished")
