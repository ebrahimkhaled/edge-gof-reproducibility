## run_I_p3_calslope.R -- the calibration intercept and slope on the paper-3 validation half, which the Cox/Miller
## line of run_I_bigdata_p3.R prints as NA (c() renames coef(f1)[1] to "a.(Intercept)", so sl["a"] finds nothing).
## Also the calibration-in-the-large gap and the model's discrimination, for the text.
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
suppressPackageStartupMessages(library(stats))
source("_cohort_p3.R")
CO <- cohort_p3(getwd())
fit <- glm(CO$f, data = CO$Ddev, family = binomial())
p <- as.numeric(predict(fit, newdata = CO$Dval, type = "response")); y <- CO$Dval$y
lp <- qlogis(p)
cal <- coef(glm(y ~ lp, family = binomial()))                                   # slope, with free intercept
citl <- coef(glm(y ~ offset(lp), family = binomial()))[1]                       # calibration in the large
auc <- { r <- rank(p); (sum(r[y == 1]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(1 - y)) }
cat(sprintf("validation n=%d events=%d | mean predicted %.4f observed %.4f\n", length(y), sum(y), mean(p), mean(y)))
cat(sprintf("calibration slope %.3f (intercept %.3f) | calibration-in-the-large %.3f | AUC %.3f\n",
            cal[2], cal[1], citl, auc))
cat(sprintf("coefficients: num_lab_procedures %.5f, number_inpatient %.4f; sd of linear predictor %.3f\n",
            coef(fit)[["num_lab_procedures"]], coef(fit)[["number_inpatient"]], sd(lp)))

## the numbers the figures and Section 7 print, written once so that nothing is typed into a plot
ll0 <- sum(y * log(p) + (1 - y) * log(1 - p))
cox <- 2 * (as.numeric(logLik(glm(y ~ lp, family = binomial()))) - ll0)
V <- read.csv(edge_path("results/cohort/runI_p3_val.csv")); Dv <- read.csv(edge_path("results/cohort/runI_p3_dev.csv"))
out <- data.frame(cox_stat = cox, cox_p = pchisq(cox, 2, lower.tail = FALSE), cal_slope = cal[2], cal_intercept = cal[1],
                  citl = citl, auc = auc, n_val = length(y), n_dev = nrow(CO$Ddev),
                  edge_val_p_min = min(V$poly_p), edge_val_p_max = max(V$poly_p),
                  edge_val_S_min = min(V$poly_stat), edge_val_S_max = max(V$poly_stat),
                  edge_dev_p_min = min(Dv$EDGE.poly3), edge_dev_p_max = max(Dv$EDGE.poly3))
write.csv(out, edge_path("results/cohort/runI_p3_summary.csv"), row.names = FALSE)
print(t(out))

## the GiViTI belt's external-mode test on the same frozen predictions, deposited rather than printed only
suppressPackageStartupMessages(library(givitiR))
gv <- tryCatch(givitiCalibrationTest(y, p, devel = "external")$p.value, error = function(e) NA_real_)
out$giviti_external_p <- gv
write.csv(out, edge_path("results/cohort/runI_p3_summary.csv"), row.names = FALSE)
cat(sprintf("GiViTI external p = %.3g\n", gv))
