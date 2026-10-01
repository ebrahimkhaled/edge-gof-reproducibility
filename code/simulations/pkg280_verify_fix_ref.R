## pkg280_verify_fix_ref.R -- runs pkg280_reference.R in "after2 devtree" mode against a clean clone of the fixed
## branch, writing to a scratch folder instead of the theory folder, so the fixer's pkg280_reference_after2.rds is
## not overwritten. The script's own after2 block compares against pkg280_reference_after.rds (2.8.0 before the
## review fixes). This runner then also compares the untouched sections with pkg280_reference_before.rds (2.7.0).
## Run: Rscript pkg280_verify_fix_ref.R after2 devtree > ../paper_EDGE/theory/pkg280_verify_fix_ref.log 2>&1
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

VF   <- file.path(tempdir(), "vfix")
SIM  <- edge_path("code/simulations")
THEO <- edge_path("declarations")
CLONE <- file.path(VF, "head")
RDIR  <- file.path(VF, "ref")
dir.create(RDIR, showWarnings = FALSE)
stopifnot(file.copy(file.path(THEO, "pkg280_reference_after.rds"), file.path(RDIR, "pkg280_reference_after.rds"),
                    overwrite = TRUE))
src <- readLines(file.path(SIM, "pkg280_reference.R"))
i1 <- grep("^PKG <- ", src); i2 <- grep("^OUT <- ", src)
stopifnot(length(i1) == 1, length(i2) == 1)
src[i1] <- sprintf('PKG <- "%s"', CLONE)
src[i2] <- sprintf('OUT <- "%s"', RDIR)
eval(parse(text = src, keep.source = FALSE), envir = globalenv())

## ---- extra: the sections that 2.8.0 must not move, against 2.7.0 ----
before <- readRDS(file.path(THEO, "pkg280_reference_before.rds"))
cat("\n==== sections against pkg280_reference_before.rds (", before$meta$version, "HEAD", before$meta$head, ") ====\n")
flat <- function(x, path = "") {
  if (is.list(x)) {
    out <- list()
    for (i in seq_along(x)) out <- c(out, flat(x[[i]], paste0(path, "$", if (is.null(names(x))) i else names(x)[i])))
    out
  } else setNames(list(x), path)
}
for (sec in c("def", "edge", "ensemble", "features", "features_right", "cdef", "calm", "shrink", "legoft",
              "legoft_localize")) {
  a <- flat(before[[sec]]); b <- flat(ref[[sec]])
  same_names <- identical(names(a), names(b))
  maxd <- 0; nbad <- 0; nother <- 0
  if (same_names) for (k in names(a)) {
    x <- a[[k]]; y <- b[[k]]
    if (is.numeric(x) && is.numeric(y) && length(x) == length(y) && identical(is.na(x), is.na(y))) {
      ok <- !is.na(x); dd <- abs(x[ok] - y[ok])
      if (length(dd)) { maxd <- max(maxd, dd); nbad <- nbad + sum(dd > 1e-10 * pmax(1, abs(x[ok]))) }
    } else if (!identical(x, y)) nother <- nother + 1
  }
  cat(sprintf("  %-16s leaves %4d | same structure %s | numbers beyond 1e-10: %d | non-numeric differences: %d | max abs diff %.3g\n",
              sec, length(a), same_names, nbad, nother, maxd))
}
cat("  battery rows common to 2.7.0, apart from Stukel:\n")
for (nm in names(before$battery)) {
  A <- before$battery[[nm]]; B <- ref$battery[[nm]]
  common <- setdiff(intersect(A$Test, B$Test), "Stukel")
  A <- A[match(common, A$Test), ]; B <- B[match(common, B$Test), ]
  d <- abs(c(A$Statistic, A$df, A$p_value) - c(B$Statistic, B$df, B$p_value))
  cat(sprintf("    %-10s rows %d -> %d | common %d | max abs diff %.3g | NA pattern same %s | Notes same %s | new rows: %s\n",
              nm, nrow(before$battery[[nm]]), nrow(ref$battery[[nm]]), length(common), max(d, na.rm = TRUE),
              identical(is.na(c(A$Statistic, A$df, A$p_value)), is.na(c(B$Statistic, B$df, B$p_value))),
              identical(A$Note, B$Note), paste(setdiff(ref$battery[[nm]]$Test, before$battery[[nm]]$Test), collapse = ",")))
}
cat("done\n")
