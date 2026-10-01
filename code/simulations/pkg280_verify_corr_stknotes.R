## pkg280_verify_corr_stknotes.R -- the Note of every Stukel NA with at least one event in pkg280_verify_corr_size.R,
## regenerated from its seed (same generator and seed rule).
## Run: Rscript pkg280_verify_corr_stknotes.R > ../paper_EDGE/theory/pkg280_verify_corr_stknotes.log 2>&1
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
SIM <- edge_path("code/simulations")
OUT <- edge_path("declarations")
options(width = 200)
suppressPackageStartupMessages(library(ebrahim.gof))
source(file.path(SIM, "_dgp_library.R"))
P <- read.csv(file.path(OUT, "pkg280_verify_corr_size_pvalues.csv"))
K <- P[is.na(P$Stukel.joint) & P$events > 0, c("n", "r", "events")]
for (k in seq_len(nrow(K))) {
  set.seed(78000000 + 1000 * K$n[k] + K$r[k])
  g <- gen_sparse_link(K$n[k], intercept = -4.9)
  fit <- suppressWarnings(glm(g$f, family = binomial(), data = g$d))
  stopifnot(sum(fit$y) == K$events[k])
  st <- run.all.gof(fit, tests = "Stukel", install = "no")
  cat(sprintf("n %d r %4d events %d slope %9.5f | %s\n", K$n[k], K$r[k], K$events[k], coef(fit)[2], st$Note))
}
cat("done\n")
