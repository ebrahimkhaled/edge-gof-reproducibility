## analyse_blockRQ.R -- the declared analysis of block RQ (declaration sha256 f7354147): false alarms, size-adjusted
## power at the clean logistic null, and exact McNemar tests of the directed test against the robust Stukel tests.
##   Rscript analyse_blockRQ.R   -> battery/RQ/_summary.csv
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
suppressMessages(library(data.table))
B <- edge_battery("RQ")
T <- c("RQD.Huber", "RQD.Mallows", "EDGE.default", "EDGE.G10")
rd <- function(cell) fread(file.path(B, sprintf("%s_pvalues.csv.gz", cell)))
rej <- function(p, crit = 0.05) !is.na(p) & p <= crit
null <- rd("RQ_logit_clean_n1000")
crit <- vapply(T, function(t) as.numeric(quantile(null[[t]], 0.05, na.rm = TRUE, type = 1)), 0)
cells <- sub("_pvalues.csv.gz", "", list.files(B, pattern = "_pvalues.csv.gz$"))
S <- rbindlist(lapply(cells, function(cl) {
  M <- rd(cl)
  rbindlist(lapply(T, function(t) data.table(cell = cl, test = t, rejection = mean(rej(M[[t]])),
    size_adj = if (grepl("cloglog", cl)) mean(!is.na(M[[t]]) & M[[t]] <= crit[[t]]) else NA_real_,
    missing = sum(is.na(M[[t]])))))
}))
## paired: EDGE at ten groups against the Mallows robust test, in every cell
mc <- rbindlist(lapply(cells, function(cl) {
  M <- rd(cl); a <- rej(M$EDGE.G10); b <- rej(M$RQD.Mallows)
  n10 <- sum(a & !b); n01 <- sum(!a & b)
  data.table(cell = cl, edge_only = n10, mallows_only = n01,
             mcnemar_p = if (n10 + n01 > 0) binom.test(n10, n10 + n01)$p.value else 1)
}))
fwrite(S, file.path(B, "_summary.csv")); fwrite(mc, file.path(B, "_paired.csv"))
print(dcast(S, cell ~ test, value.var = "rejection"), digits = 3)
cat("\nsize-adjusted power, cloglog n = 1000:\n"); print(S[grepl("cloglog", cell), .(test, size_adj)], digits = 3)
cat("\npaired, EDGE (ten groups) vs RQD-Mallows:\n"); print(mc, digits = 3)
