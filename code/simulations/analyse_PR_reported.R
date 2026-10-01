## analyse_PR_reported.R -- the Pulkstenis-Robinson test, run in every battery scenario but not reported in paper 3
## (referee, 2026-10-01). Reads the per-replicate p-values already deposited; nothing is re-run.
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
rate <- function(f, col) { d <- fread(f); c(rej = mean(!is.na(d[[col]]) & d[[col]] <= 0.05), na = mean(is.na(d[[col]]))) }
B9 <- edge_battery("9"); C1B <- edge_battery("C1b")
rows <- list()
for (n in c(1000, 5000)) for (r in c("001", "002", "005", "010")) {
  f <- file.path(B9, sprintf("logit_C1_r%s_n%d_pvalues.csv.gz", r, n))
  if (file.exists(f)) rows[[length(rows) + 1]] <- data.table(error = "x4", n = n, k = as.numeric(r) / 1000 * n,
                                                             t(rate(f, "PR")), EDGE = rate(f, "EDGE.poly3.u.Grule")[["rej"]])
}
for (f in list.files(C1B, "pvalues.csv.gz$", full.names = TRUE)) {
  d <- fread(f); if (!"PR" %in% names(d)) next
  rows[[length(rows) + 1]] <- data.table(error = "x-4", n = NA, k = NA, file = basename(f), t(rate(f, "PR")))
}
R <- rbindlist(rows, fill = TRUE)
print(R)
CEN <- fread(edge_battery("analysis", "census_extended.csv"))
pw <- CEN[test %in% c("PR", "EDGE.poly3.u.Grule", "HL.G10"), .(mean_size_adjusted_power = mean(power, na.rm = TRUE), cells = .N), by = test]
print(pw)
fwrite(R, edge_battery("analysis", "PR_corrupted.csv")); fwrite(pw, edge_battery("analysis", "PR_power.csv"))
