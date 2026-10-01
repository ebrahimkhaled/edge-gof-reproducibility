## analyse_block9R2.R -- block 9R2 read against its declaration (sha256 9943e551...) and beside block 9R.
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
DIR <- edge_battery("9R2")
V <- c("V0", "HYB05", "HYB10", "HYB20", "G10")
BATT <- c("cauchit_n1000", "t4_n1000", "loglog_n1000", "stk_short_n1000", "stk_long_n1000", "stk_asym_n1000")

rd <- function(cell) fread(file.path(DIR, sprintf("%s_variants.csv.gz", cell)))
cells <- sub("_variants.csv.gz$", "", basename(Sys.glob(file.path(DIR, "*_variants.csv.gz"))))
R <- rbindlist(lapply(cells, function(cl) {
  M <- rd(cl)
  data.table(cell = cl, variant = V, rate = vapply(V, function(v) mean(M[[v]] <= 0.05, na.rm = TRUE), 0),
             na = vapply(V, function(v) sum(!is.finite(M[[v]])), 0L))
}))
W <- dcast(R, cell ~ variant, value.var = "rate")[, c("cell", V), with = FALSE]

cat("\n== every scenario, rejection at 0.05 (1000 replicates) ==\n")
print(W, digits = 3)

cat("\n== the three declared criteria ==\n")
lev <- W[cell %in% c("logit_clean_n1000", "logit_clean_n5000")]
cat("  1. level, clean logistic (band 0.029-0.071):\n"); print(lev, digits = 3)
prot <- W[cell %in% c("logit_C1_r005_n5000", "logit_C1_r010_n5000")]
cat("  2. protection, n = 5000 at k = 25 and k = 50 (bar 0.10):\n"); print(prot, digits = 3)
pw <- W[cell %in% BATT]
mp <- colMeans(pw[, V, with = FALSE])
cat("  3. mean power over the six battery scenarios:\n")
print(round(mp, 4))
cat("     loss against the default:\n"); print(round(mp - mp[["V0"]], 4))

## paired comparisons inside each battery scenario: HYB05 against the default and against the coarse partition
cat("\n== paired differences on the battery scenarios (exact McNemar, Holm within each family) ==\n")
pair <- function(a, b) {
  D <- rbindlist(lapply(BATT, function(cl) {
    M <- rd(cl); ra <- M[[a]] <= 0.05; rb <- M[[b]] <= 0.05
    ra[is.na(ra)] <- FALSE; rb[is.na(rb)] <- FALSE
    n10 <- sum(ra & !rb); n01 <- sum(!ra & rb)
    data.table(cell = cl, a = mean(ra), b = mean(rb), diff = mean(ra) - mean(rb), win = n10, loss = n01,
               p = if (n10 + n01 == 0) 1 else stats::binom.test(n10, n10 + n01, 0.5)$p.value)
  }))
  D[, p_holm := p.adjust(p, "holm")][]
}
cat(sprintf("\n  HYB05 (a) against the default V0 (b):\n")); print(pair("HYB05", "V0"), digits = 3)
cat(sprintf("\n  HYB05 (a) against the coarse partition G10 (b):\n")); print(pair("HYB05", "G10"), digits = 3)
cat(sprintf("\n  HYB10 (a) against the coarse partition G10 (b):\n")); print(pair("HYB10", "G10"), digits = 3)

cat("\n== the other power scenarios ==\n")
print(W[cell %in% c("probit_clean_n1000", "probit_clean_n5000", "cloglog_clean_n1000", "cloglog_clean_n5000",
                    "probit_C1_r005_n1000", "cloglog_C1_r005_n1000")], digits = 3)
cat("\n== false alarms, logistic truth with a corrupted covariate ==\n")
print(W[grepl("^logit_C1", cell)], digits = 3)
cat(sprintf("\nreplicates with no p-value: %d of %d\n", sum(R$na), 22L * 1000L * length(V)))
