## check_theta.R -- theta, the probability that a corrupted record's outcome contradicts its recorded prediction,
## under the two covariate errors of the simulation design (x ~ U(-3, 3), d ~ Bernoulli(0.5), eta = 0.6 x + 0.5 d,
## logistic truth). "Contradicts" means the outcome differs from the class the recorded prediction points to,
## y != 1{recorded eta > 0}. Paper 3 quotes about 0.3 for x -> 4x and close to 0.8 for x -> -4x; the referee round
## found no deposited computation behind either number, so this is it.
##
##   Rscript check_theta.R   -> battery/analysis/theta.csv
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
set.seed(20260930)
n <- 2e6
x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
y <- rbinom(n, 1, plogis(0.6 * x + 0.5 * d))
th <- function(xr) { er <- 0.6 * xr + 0.5 * d; mean(y != as.integer(er > 0)) }
ex <- function(xr) { er <- 0.6 * xr + 0.5 * d; top <- er > quantile(er, 0.96); mean(y[top] != 1) }
out <- data.frame(error = c("exaggeration", "sign error"),
                  theta_all = c(th(4 * x), th(-4 * x)),
                  theta_top_group = c(ex(4 * x), ex(-4 * x)))
print(out, digits = 3)
write.csv(out, edge_path("results/analysis/theta.csv"), row.names = FALSE)
