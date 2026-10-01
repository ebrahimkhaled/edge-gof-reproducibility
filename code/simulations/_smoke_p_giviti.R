## Smoke test before the three real runs: confirm the test API, the returned test names,
## and that givitiR can be driven from a fitted glm the same way the package tests are.
suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)})
source("_dgp_library.R")

set.seed(20260910)

cat("--- run.all.gof on a well-specified fit, p*=4 ---\n")
g <- gen_bench(500, p = 4)
fit <- glm(g$f, data = g$d, family = binomial())
r <- run.all.gof(fit, G = 10, include_slow = FALSE)
print(r)
cat("\nclass:", class(r), "\n")
cat("names:", paste(names(r), collapse = " | "), "\n")

cat("\n--- p*=10 (six nuisance covariates) ---\n")
g10 <- gen_bench(500, p = 10)
fit10 <- glm(g10$f, data = g10$d, family = binomial())
cat("coefficients estimated:", length(coef(fit10)), "\n")
r10 <- run.all.gof(fit10, G = 10, include_slow = FALSE)
print(r10)

cat("\n--- p*=20 with G=10: G < p*, expect trouble ---\n")
g20 <- gen_bench(1000, p = 20)
fit20 <- glm(g20$f, data = g20$d, family = binomial())
cat("coefficients estimated:", length(coef(fit20)), "\n")
r20 <- tryCatch(run.all.gof(fit20, G = 10, include_slow = FALSE), error = function(e) paste("ERROR:", conditionMessage(e)))
print(r20)

cat("\n--- givitiR on the same fit ---\n")
o <- fit$y
e <- fitted(fit)
gt <- tryCatch(givitiCalibrationTest(o, e, devel = "internal"),
               error = function(er) paste("ERROR:", conditionMessage(er)))
print(gt)
cat("\nstr of giviti result:\n"); str(gt, max.level = 1)
