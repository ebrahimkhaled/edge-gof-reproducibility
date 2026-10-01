## check_slope_signerror.R -- how far the maximum-likelihood fit moves under the sign error (block C1b's generator,
## logistic truth). Block C1b stored p-values but not the fitted coefficients, and the paper-3 referee round asked
## whether the model under test is still the clean one when the error reverses predictions. The true slope on x is
## 0.6. 400 replicates a cell, seeds disjoint from the block's.
##
##   Rscript check_slope_signerror.R   -> battery/analysis/slope_signerror.csv
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
source(file.path(SIMDIR, "_blockC1b_contam.R"))
RNGkind("L'Ecuyer-CMRG")
R <- list()
for (n in c(1000L, 5000L)) for (cor in c("C1", "C1b")) for (rate in c(0, 0.001, 0.002, 0.005, 0.01)) {
  b <- vapply(1:400, function(rep) {
    set.seed(990000000 + n + 1e6 * rate + rep)
    g <- c1b_gen("logit", n, if (rate == 0) "clean" else cor, rate)
    coef(glm(g$f, data = g$d, family = binomial()))[["x"]]
  }, 0)
  R[[length(R) + 1]] <- data.table(n = n, error = c(C1 = "exaggeration", C1b = "sign error")[[cor]],
                                   k = round(rate * n), mean_slope = mean(b), change_pct = 100 * (mean(b) / 0.6 - 1))
}
R <- rbindlist(R)
fwrite(R, edge_battery("analysis", "slope_signerror.csv"))
print(R, digits = 3)
