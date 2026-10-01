## run_M_block9b.R -- block 9b: are the slow rivals robust to a corrupted record?
## Contract: paper_EDGE/theory/PREDECLARATION_block9b_slow_rivals_robustness.md (sha256 72d703a9...).
##
##   Rscript run_M_block9b.R [--workers 20] [--reps N] [--out battery/9b]
##
## 3 cells x 100 replicates; BAGofT at the package defaults dominates the cost (~8 minutes a data set
## under load), so replicates are written to <cell>_pvalues.csv.gz.part after every batch and a
## restarted run resumes from it. Refuses to write into a finished block's folder.
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
source(file.path(SIMDIR, "_proj_test.R"))
source(file.path(SIMDIR, "_block9b_rivals.R"))

opt <- list(); a <- commandArgs(TRUE); i <- 1L
while (i <= length(a)) { opt[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
W    <- if (is.null(opt$workers)) 20L else as.integer(opt$workers)
REPS <- if (is.null(opt$reps)) NA_integer_ else as.integer(opt$reps)
OUT  <- normalizePath(if (is.null(opt$out)) edge_battery("9b") else opt$out, winslash = "/", mustWork = FALSE)
busy <- normalizePath(edge_battery(c(as.character(0:9), "9c", "9d", "dryrun", "analysis")),
                      winslash = "/", mustWork = FALSE)
if (any(startsWith(paste0(OUT, "/"), paste0(busy, "/")))) stop("--out cannot be a finished block's folder: ", OUT)
if (!is.na(REPS) && identical(OUT, normalizePath(edge_battery("9b"), winslash = "/", mustWork = FALSE)))
  stop("a test run (--reps) cannot write to the real block 9b folder; pass --out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "_progress.log")
say <- function(m) { l <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "", m)
                     cat(l, "\n", sep = ""); cat(l, "\n", sep = "", file = LOG, append = TRUE) }

C <- b9b_cells()
say(sprintf("block 9b: %d cells, B = %s, %d workers, out %s", nrow(C),
            if (is.na(REPS)) "100" else paste0(REPS, " (TEST RUN)"), W, OUT))

cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, RNGkind("L'Ecuyer-CMRG")))
clusterExport(cl, "SIMDIR", envir = environment())
invisible(clusterEvalQ(cl, { suppressMessages({ requireNamespace("BAGofT"); requireNamespace("randomForest") })
                             source(file.path(SIMDIR, "_battery_tests.R"))
                             source(file.path(SIMDIR, "_proj_test.R"))
                             source(file.path(SIMDIR, "_block9b_rivals.R")); TRUE }))

read_part <- function(p) {
  for (f in c(p, paste0(p, ".tmp"))) if (file.exists(f)) {
    x <- tryCatch(as.data.frame(fread(f)), error = function(e) NULL); if (!is.null(x)) return(x) }
  NULL
}

t0 <- Sys.time()
for (i in seq_len(nrow(C))) {
  ce  <- as.list(C[i, ]); R <- if (is.na(REPS)) as.integer(ce$B) else REPS
  out <- file.path(OUT, paste0(ce$cell, "_pvalues.csv.gz")); part <- paste0(out, ".part")
  if (file.exists(out)) { say(sprintf("%s: finished file present, skipped", ce$cell)); next }
  M <- read_part(part)
  todo <- setdiff(seq_len(R), M$rep)
  say(sprintf("%s: %d of %d replicates to run%s", ce$cell, length(todo), R,
              if (is.null(M)) "" else sprintf(" (%d resumed)", nrow(M))))
  for (b in split(todo, ceiling(seq_along(todo) / W))) {
    ## the extra argument must not be a prefix of clusterApplyLB's formals (cl, x, fun)
    rows <- clusterApplyLB(cl, b, function(r, cellrow) b9b_one(r, cellrow), cellrow = ce)
    M <- rbind(M, as.data.frame(do.call(rbind, rows)))
    M <- M[order(M$rep), , drop = FALSE]; rownames(M) <- NULL
    tmp <- paste0(part, ".tmp"); fwrite(M, tmp, compress = "gzip")
    if (file.exists(part)) file.remove(part)
    if (!file.rename(tmp, part)) stop("could not rename ", tmp)
    say(sprintf("  %-5s %d/%d replicates; median BAGofT %.0f s, projection %.0f s; elapsed %.1f min",
                ce$cell, nrow(M), R, median(M$sec_bagoft, na.rm = TRUE), median(M$sec_proj, na.rm = TRUE),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  if (nrow(M) != R || any(M$rep != seq_len(R))) stop("internal: ", part, " is incomplete")
  if (!file.rename(part, out)) stop("could not rename ", part)
  say(sprintf("[%d/%d] %s written: errors %d, BAGofT missing %d, projection missing %d",
              i, nrow(C), ce$cell, sum(M$flag.b9b_error %in% 1), sum(!is.finite(M$BAGofT)),
              sum(!is.finite(M$proj))))
}
say("block 9b finished")
