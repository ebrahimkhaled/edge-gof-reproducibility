## analyse_stukel_marginal_size.R -- size of the Stukel form distributed in LogisticDx (sum of the two marginal score
## components referred to chi-square(2), column Stk.marg) over the null scenarios of the main battery, beside the joint
## score the paper uses (referee minor comment 11, 2026-10-01). Reads the battery's block summaries.
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
fs <- edge_battery(c("0", "1a", "1b", "2", "3", "4", "5", "6", "7"), "_summary.csv")
S <- rbindlist(lapply(fs[file.exists(fs)], fread), fill = TRUE)
for (t in c("Stk.marg", "Stk.joint")) {
  x <- S[test == t & role == "null" & alpha == 0.05 & is.finite(rejection)]
  cat(sprintf("%-9s %3d null cells: median %.4f, mean %.4f, quartiles %.4f to %.4f\n", t, nrow(x),
              median(x$rejection), mean(x$rejection), quantile(x$rejection, 0.25), quantile(x$rejection, 0.75)))
}
