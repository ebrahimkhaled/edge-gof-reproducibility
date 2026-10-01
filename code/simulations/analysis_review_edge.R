## analysis_review_edge.R -- edge cases of analyse_M_battery.R on a copy of the synthetic set1 root, treated as a "real root"
## (bat = the copy, launch.log records "launch finished"): what happens when a family cell HAS a summary but no
## size-adjusted power (matched null missing), for one cell and for a whole family. No value of the real battery is read.
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
suppressPackageStartupMessages(library(data.table))
SIMDIR <- edge_path("code/simulations")
RV <- edge_battery("_review", "analysis", "review")
src <- file.path(RV, "set1")
sink(file.path(RV, "review_edge.log"), split = TRUE)
BE <- new.env()
sys.source(file.path(SIMDIR, "analyse_M_battery.R"), envir = BE)

make_copy <- function(dir, blank) {
  unlink(dir, recursive = TRUE); dir.create(dir, recursive = TRUE)
  file.copy(file.path(src, "cells.csv"), file.path(dir, "cells.csv"))
  for (b in c("1a", "1b", "2", "3", "4", "5", "6", "7")) {
    S <- fread(file.path(src, b, "_summary.csv"), colClasses = list(character = c("block", "cell", "role", "null_cell", "subset", "test", "status")))
    hit <- paste(S$block, S$cell) %in% blank & startsWith(S$test, "EDGE.")
    S[hit, `:=`(size_adj_power = NA_real_, status = "matched null not run yet")]
    dir.create(file.path(dir, b), showWarnings = FALSE)
    fwrite(S, file.path(dir, b, "_summary.csv"))
  }
  writeLines(c("2026-09-14 08:28:48  launch from step 1 of 19", "2026-09-15 02:00:00  launch finished"), file.path(dir, "launch.log"))
}
try_run <- function(dir) {
  res <- NULL
  msg <- tryCatch({ txt <- capture.output(res <- BE$an_main(c("--root", dir, "--out", file.path(dir, "out")), bat = dir)); "ran" },
                  error = function(e) conditionMessage(e))
  list(msg = msg, res = res)
}

cat("case 1: one family 1 cell (3 link_probit_n1000) has a summary but no EDGE size-adjusted power\n")
d1 <- file.path(RV, "edge_one_cell"); make_copy(d1, "3 link_probit_n1000")
r1 <- try_run(d1)
cat("  outcome:", r1$msg, "\n")
if (!is.null(r1$res)) {
  print(r1$res$families[family %in% c("1", "macro") & basis == "poly3", .(form, family, cells, with_value, mean_power)])
  cat("  decision reason poly3:", r1$res$decision[basis == "poly3"]$reason, "\n")
  cat("  H1 rows:\n"); print(r1$res$hypotheses[hypothesis == "H1", .(form, cells, cells_with_value, statistic, verdict)])
}

cat("\ncase 2: both family 6 cells (crossover_n1000, crossover_n2000) have summaries but no EDGE size-adjusted power\n")
d2 <- file.path(RV, "edge_family6"); make_copy(d2, c("4 crossover_n1000", "4 crossover_n2000"))
r2 <- try_run(d2)
cat("  outcome:", r2$msg, "\n")
if (!is.null(r2$res)) {
  print(r2$res$families[family %in% c("6", "macro") & basis == "poly3", .(form, family, cells, with_value, mean_power)])
  cat("  decision reason poly3:", r2$res$decision[basis == "poly3"]$reason, "\n")
  cat("  macro-average check: mean of the five remaining family means =",
      all(abs(r2$res$families[basis == "poly3" & family == "macro"]$mean_power -
              r2$res$families[basis == "poly3" & family %in% as.character(1:5), .(m = mean(mean_power)), by = form]$m) < 1e-12), "\n")
}
unlink(c(file.path(d1, "out"), file.path(d2, "out")), recursive = TRUE)
sink()
