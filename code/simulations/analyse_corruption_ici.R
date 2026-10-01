## analyse_corruption_ici.R -- how miscalibrated are the corrupted data sets? (prompted by the author, 2026-10-01)
## For the corruptions of Section 6 (c1b_gen: x ~ U(-3,3), d ~ Bern(.5), logistic truth eta = 0.6x + 0.5d; C1 x4 and
## C1b x(-4) on k random records), fit the model to the recorded data and compare each record's fitted risk with its
## true risk:
##   ICI_all   = mean |p_true - p_fitted| over all n records (the calibration of the data set as recorded);
##   ICI_clean = the same over the uncorrupted records only (what the fit did to everyone else);
##   ICI_bad   = the same over the k corrupted records.
## 200 replicates a setting; the same generator as blocks 9 and C1b.
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
suppressMessages(library(data.table))
source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
source(file.path(SIMDIR, "_block9_contam.R")); source(file.path(SIMDIR, "_blockC1b_contam.R"))
S <- CJ(corruption = c("clean", "C1", "C1b"), n = c(1000L, 5000L), k = c(0L, 1L, 10L, 50L))
S <- S[(corruption == "clean") == (k == 0)][!(n == 1000 & k == 50)]
R <- rbindlist(lapply(seq_len(nrow(S)), function(j) {
  s <- S[j]
  v <- t(vapply(1:200, function(r) {
    set.seed(5500000 + 1000 * j + r)
    g <- c1b_gen("logit", s$n, s$corruption, s$k / s$n)
    f <- suppressWarnings(glm(g$f, data = g$d, family = binomial()))
    ad <- abs(g$p_true - fitted(f)); bad <- seq_len(s$n) %in% g$corrupt
    c(all = mean(ad), clean = mean(ad[!bad]), bad = if (any(bad)) mean(ad[bad]) else NA_real_)
  }, numeric(3)))
  data.table(s, ICI_all = mean(v[, 1]), ICI_clean = mean(v[, 2]), ICI_bad = mean(v[, 3], na.rm = TRUE))
}))
print(R, digits = 3)
fwrite(R, edge_battery("analysis", "corruption_ici.csv"))
