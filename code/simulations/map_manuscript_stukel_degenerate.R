## map_manuscript_stukel_degenerate.R -- READ-ONLY check made while mapping the manuscript (2026-09-13).
## When every fitted risk lies below 0.5, one Stukel carrier is identically zero. What do
## (1) ebrahim.gof 2.6.0 gof_stukel (statmod path), (2) the joint score, (3) the harness LR refit return?
## Prints only; writes nothing.
suppressPackageStartupMessages({ library(ebrahim.gof); library(statmod) })
set.seed(20260913)
n <- 1500; x <- runif(n, -3, 3)
y <- rbinom(n, 1, plogis(-2.5 + 0.4 * x))
fit <- glm(y ~ x, family = binomial())
ph <- fitted(fit); eta <- predict(fit)
cat(sprintf("max fitted risk %.3f (all below 0.5: %s); events %d\n", max(ph), all(ph < 0.5), sum(y)))
za <- 0.5 * eta^2 * (ph >= 0.5); zb <- -0.5 * eta^2 * (ph < 0.5)
cat(sprintf("za identically zero: %s\n", all(za == 0)))

## (1) the package function, through its internal context
g <- tryCatch(run.all.gof(fit, tests = "Stukel", include_slow = FALSE, install = "no"), error = function(e) conditionMessage(e))
cat("\n(1) run.all.gof(tests = 'Stukel'):\n"); print(g)
z <- tryCatch(glm.scoretest(fit, cbind(za, zb)), error = function(e) conditionMessage(e))
cat("    statmod::glm.scoretest z's:", format(z), "\n")

## (2) joint score with a one-df fallback
zz <- glm.scoretest(fit, zb)
cat(sprintf("\n(2) one-df score on the non-zero carrier: z = %.3f, p = %.4f\n", zz, 2 * pnorm(-abs(zz))))

## (3) harness LR refit (stukel_p logic, without the separation guard)
d <- data.frame(y = y, x = x, za = za, zb = zb)
fa <- suppressWarnings(glm(y ~ x + za + zb, family = binomial(), data = d))
lr <- deviance(fit) - deviance(fa)
cat(sprintf("\n(3) LR refit: coefficients %s\n    LR = %.3f; p on 2 df = %.4f; p on 1 df = %.4f; rank of augmented fit = %d (base %d)\n",
    paste(names(coef(fa)), signif(coef(fa), 3), collapse = ", "), lr,
    pchisq(lr, 2, lower.tail = FALSE), pchisq(lr, 1, lower.tail = FALSE), fa$rank, fit$rank))
