## run_I_p3_ici.R -- effect sizes for the miscalibration of Section 7 (referee M7, 2026-10-01): the integrated calibration
## index (ICI), E50, E90 and Emax of Austin and Steyerberg (Statistics in Medicine 2019), from a loess calibration curve
## of the observed outcome on the predicted risk, on the validation half of the cohort.
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
setwd(SIMDIR)
source("_cohort_p3.R")
CO <- cohort_p3(SIMDIR)
fit <- stats::glm(CO$f, data = CO$Ddev, family = stats::binomial())
p <- as.numeric(stats::predict(fit, newdata = CO$Dval, type = "response"))
y <- CO$Dval$y
## Austin and Steyerberg's recipe: loess of y on p with the default span and degree 2, evaluated at each prediction
lo <- stats::loess(y ~ p, degree = 2)
pc <- stats::predict(lo, newdata = data.frame(p = p))
ad <- abs(pc - p)
out <- data.frame(n = length(y), ICI = mean(ad), E50 = stats::median(ad), E90 = as.numeric(stats::quantile(ad, 0.9)),
                  Emax = max(ad), mean_pred = mean(p), mean_obs = mean(y))
print(out, digits = 4)
utils::write.csv(out, edge_path("results/cohort/runI_p3_ici.csv"), row.names = FALSE)
