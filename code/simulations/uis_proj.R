## Run OUR (corrected) proj implementation on the UIS linear-NDRGTX model,
## to (a) validate it reproduces Liu 2024's proj rejecting the linear model,
## (b) give the direct proj-vs-EDGE contrast for the real-data vignette.
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
suppressMessages(library(stats))
SIM <- edge_path("code/simulations")
source(file.path(SIM, "_proj_test.R"))
set.seed(20260707)

## ---- load UIS (persist it this time) ----
rds <- file.path(SIM, "uis_data.rds")
if (file.exists(rds)) {
  uis <- readRDS(rds)
} else if (requireNamespace("LogisticDx", quietly=TRUE)) {
  data(uis, package="LogisticDx"); saveRDS(uis, rds)
} else {
  tmp <- tempfile(fileext=".tar.gz")
  download.file("https://cran.r-project.org/src/contrib/Archive/LogisticDx/LogisticDx_0.3.tar.gz",
                tmp, mode="wb", quiet=TRUE)
  ex <- tempfile(); dir.create(ex)
  untar(tmp, files="LogisticDx/data/uis.rda", exdir=ex)
  load(file.path(ex, "LogisticDx/data/uis.rda"))   # -> uis
  saveRDS(uis, rds)
}
cat("uis dim:", dim(uis), "\n")

d <- uis
d$y <- as.integer(d$DFREE == "no")     # event = returned to drug use (147), rate ~0.256
cat("event rate:", round(mean(d$y),4), "\n")

fit_lin  <- glm(y ~ AGE + NDRGTX + IVHX + RACE + TREAT + SITE, data=d, family=binomial)
suppressMessages(library(splines))
fit_rich <- glm(y ~ AGE + ns(NDRGTX,3) + IVHX + RACE + TREAT + SITE, data=d, family=binomial)

Xlin  <- model.matrix(fit_lin)
Xrich <- model.matrix(fit_rich)

cat("\n--- proj on LINEAR-NDRGTX model (Liu model 17), B=1000 ---\n")
t0 <- Sys.time()
pl <- proj_pvalue(fit_lin$y, Xlin, B=1000, seed=101)
cat("proj  T =", signif(pl$stat,5), "  p =", round(pl$p_value,4),
    "  (", round(as.numeric(difftime(Sys.time(),t0,units="secs")),1), "s )\n")

cat("\n--- proj on RICHER spline-NDRGTX model (control) ---\n")
pr <- proj_pvalue(fit_rich$y, Xrich, B=1000, seed=102)
cat("proj  T =", signif(pr$stat,5), "  p =", round(pr$p_value,4), "\n")

cat("\n=== SUMMARY (proj vs the index-based tests already computed) ===\n")
cat(sprintf("LINEAR model:  proj p=%.4f   [EDGE-poly3 p=0.478, HL p=0.674, Stukel p=0.386]\n", pl$p_value))
cat(sprintf("SPLINE control: proj p=%.4f\n", pr$p_value))
cat("PROJ_DONE\n")
