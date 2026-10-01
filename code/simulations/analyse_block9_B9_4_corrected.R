## analyse_block9_B9_4_corrected.R -- the reformulation of declared criterion B9.4, deposited.
##
## WHY THIS FILE EXISTS. B9.4 was declared as a RELATIVE change in each fitted coefficient, with the
## criterion max(|d_intercept|, |d_x|, |d_d|) < 0.10. The design's true intercept is zero, so a
## relative change in it divides by a quantity estimated near zero: in battery/9/analysis/B9_4_fit.csv
## that column reaches 61.6, and the criterion cannot be read as declared. The quantity the criterion
## was about -- whether the corruption moves the fitted model -- is well defined for the intercept on
## the linear-predictor scale, and that is what B9_4_fit_corrected.csv reports.
##
## The reformulation was made after the numbers were computed. It is recorded here, and in the paper's
## declaration-conformance note, rather than left as an undocumented second file beside the first: a
## study whose argument is that thresholds were fixed in advance cannot leave a criterion silently
## restated. Both files are kept; the paper quotes the relative change for the slopes and reports the
## intercept on its own scale.
##
##   Rscript analyse_block9_B9_4_corrected.R
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
suppressPackageStartupMessages(library(data.table))
setwd(SIMDIR)
FITS <- edge_battery("9", "analysis", "fits.csv")
OUT  <- edge_battery("9", "analysis", "B9_4_fit_corrected.csv")

F <- fread(FITS)
need <- c("cell", "truth", "n", "corruption", "rate", "k", "b_intercept", "b_x", "b_d")
if (!all(need %in% names(F))) stop("fits.csv does not carry ", paste(setdiff(need, names(F)), collapse = ", "))

## the clean fit of the same truth at the same sample size is the reference for every contaminated cell
clean <- F[corruption == "clean", .(truth, n, c_intercept = b_intercept, c_x = b_x, c_d = b_d)]
D <- merge(F[corruption != "clean"], clean, by = c("truth", "n"), all.x = TRUE)
D[, abs_intercept := b_intercept - c_intercept]          # the intercept on the linear-predictor scale
D[, rel_x := (b_x - c_x) / c_x]                          # the slopes keep their declared relative reading
D[, rel_d := (b_d - c_d) / c_d]
D[, worst_slope := pmax(abs(rel_x), abs(rel_d))]
D[, holds := worst_slope < 0.10]
setorder(D, truth, n, corruption, rate)
fwrite(D[, .(cell, truth, n, corruption, rate, k, b_intercept, b_x, b_d,
             abs_intercept, rel_x, rel_d, worst_slope, holds)], OUT)

C1 <- D[corruption == "C1"]
cat("contaminated scenarios:", nrow(D), " of which C1:", nrow(C1), "\n")
cat("max |relative change| in the slope on the corrupted covariate (C1): ",
    sprintf("%.4f", max(abs(C1$rel_x))), "\n")
cat("max |relative change| in the slope on the second covariate (C1):    ",
    sprintf("%.4f", max(abs(C1$rel_d))), "\n")
cat("max |change| in the intercept, linear-predictor scale (C1):         ",
    sprintf("%.4f", max(abs(C1$abs_intercept))), "\n")
cat("C1 scenarios meeting the slope criterion: ", sum(C1$holds), " of ", nrow(C1), "\n", sep = "")
cat("\nwrote ", OUT, "\n", sep = "")
