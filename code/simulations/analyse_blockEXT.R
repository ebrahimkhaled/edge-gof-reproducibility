## analyse_blockEXT.R -- the analysis fixed in PREDECLARATION_blockEXT_external_validation.md (sha256 754cf309...).
##   level: rejection at 5% in cells 1-3; holds within 0.05 +/- 3 SE (0.029 to 0.071 at 1000 replicates)
##   power: size-adjusted at the 5% quantile of each test's p-values in cell 2 (null, n = 500)
##   protection: false-alarm rate in cells 14-21; a test holds where it is at most 0.10
## A replicate with no p-value counts as a non-rejection.
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
IN <- edge_battery("EXT")
fs <- sort(list.files(IN, pattern = "_pvalues\\.csv\\.gz$", full.names = TRUE))
stopifnot(length(fs) == 21L)
P <- rbindlist(lapply(fs, function(f) {
  d <- fread(f); m <- regmatches(basename(f), regexec("c(\\d+)_([a-z]+)_C([-0-9.]+)_n(\\d+)_k(\\d+)", basename(f)))[[1]]
  d[, `:=`(fam = m[3], C = as.numeric(m[4]), n = as.integer(m[5]), k = as.integer(m[6]))]; d }))
TESTS <- c("EDGE_auto", "EDGE_G10", "cox", "spiegelhalter", "giviti", "HL", "stukel")
L <- melt(P, id.vars = c("id", "rep", "fam", "C", "n", "k"), measure.vars = TESTS, variable.name = "test", value.name = "p")
L[, test := as.character(test)]
cat(sprintf("missing p-values by test: %s\n", paste(L[, .(m = sum(is.na(p))), by = test][, sprintf("%s %d", test, m)], collapse = ", ")))
crit <- L[id == 2, .(crit = as.numeric(quantile(p, 0.05, type = 1, na.rm = TRUE))), by = test]
L <- merge(L, crit, by = "test")
S <- L[, .(raw = mean(!is.na(p) & p < 0.05), adj = mean(!is.na(p) & p <= crit), reps = .N), by = .(id, fam, C, n, k, test)]
S[, role := ifelse(fam == "null", "level", ifelse(k > 0, "corrupted", "alternative"))]
S[, holds := ifelse(role == "level", raw >= 0.029 & raw <= 0.071, ifelse(role == "corrupted", raw <= 0.10, NA))]
fwrite(S, file.path(IN, "_summary.csv"))

wide <- function(D, v) dcast(D, id + fam + C + n + k ~ test, value.var = v)[, c("id", "fam", "C", "n", "k", TESTS), with = FALSE]
cat("\n--- level (raw rejection at 5%) ---\n");            print(wide(S[role == "level"], "raw"), digits = 3)
cat("\n--- power, size-adjusted at cell 2 ---\n");        print(wide(S[role == "alternative"], "adj"), digits = 3)
cat("\n--- mean size-adjusted power over the ten alternatives ---\n")
print(S[role == "alternative", .(mean_adj = round(mean(adj), 3)), by = test][order(-mean_adj)])
cat("\n--- false alarms, correct model, k corrupted records (x2 recorded as C * x2), n = 1000 ---\n")
print(wide(S[role == "corrupted"], "raw"), digits = 3)
cat("\n--- how many corrupted cells each test holds (<= 0.10), of 8 ---\n")
print(S[role == "corrupted", .(holds = sum(holds)), by = test][order(-holds)])
cat("\n--- GiViTI and Cox: share of replicates with identical p-values ---\n")
print(P[, .(identical = mean(abs(giviti - cox) < 1e-12, na.rm = TRUE)), by = id])
