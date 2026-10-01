## pkg280_verify_num_repro.R -- minimal reproducible examples for the findings of the 2.8.0 numerical verification,
## using the INSTALLED ebrahim.gof. Each block runs on its own.
## Run: Rscript pkg280_verify_num_repro.R > ../paper_EDGE/theory/pkg280_verify_num_repro.log 2>&1

suppressPackageStartupMessages(library(ebrahim.gof))
stopifnot(as.character(packageVersion("ebrahim.gof")) == "2.8.0")
options(width = 160, warn = 1)
cat("ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "|", R.version.string, "\n")

## ---- R1. joint Stukel is NA on a rank-deficient glm (anova Rao and form = "lr" are fine) ---------------------------
cat("\n==== R1 ====\n")
set.seed(3)
gd <- data.frame(xa = runif(300, -3, 3))
gd$out <- rbinom(300, 1, plogis(0.5 * gd$xa))
gd$xa2 <- 2 * gd$xa
fit_rd <- glm(out ~ xa + xa2, family = binomial(), data = gd)       # coefficient of xa2 is NA
fit_fr <- glm(out ~ xa, family = binomial(), data = gd)             # the same model, full rank
print(as.data.frame(run.all.gof(fit_rd, tests = "Stukel")))
print(as.data.frame(run.all.gof(fit_rd, tests = "Stukel", control = list(Stukel = list(form = "lr")))))
print(as.data.frame(run.all.gof(fit_fr, tests = "Stukel")))
eta <- predict(fit_rd); ph <- fitted(fit_rd)
za <- 0.5 * eta^2 * (ph >= 0.5); zb <- -0.5 * eta^2 * (ph < 0.5)
print(anova(fit_rd, glm(out ~ xa + xa2 + za + zb, family = binomial(), data = gd), test = "Rao"))

## ---- R2. unit form stops with a singular solve when one group mean risk is a hair above 0.5 (also 2.7.0) ------------
cat("\n==== R2 ====\n")
set.seed(4)
n <- 500
ph <- c(sort(runif(475, 0.02, 0.45)), rep(0.5001, 25))              # G = 20: the top group has mean risk 0.5001
X <- cbind(1, qlogis(ph))
y <- rbinom(n, 1, ph)
r_unit  <- tryCatch(def.gof(y, ph, X = X, G = 20, basis = "stukel"), error = function(e) conditionMessage(e))
r_score <- tryCatch(def.gof(y, ph, X = X, G = 20, basis = "stukel", weights = "score"), error = function(e) conditionMessage(e))
print(r_unit); print(r_score)
cat("qlogis(0.5001)^2 =", qlogis(0.5001)^2, "(kept by the colSums(abs(Z)) > 1e-8 rule)\n")

## ---- R3. form = "lr" is NA under quasi-separation of the augmented fit (joint is fine) -----------------------------
cat("\n==== R3 ====\n")
set.seed(142)
dq <- data.frame(xa = runif(400, -3, 3), db = rbinom(400, 1, 0.5))
dq$out <- rbinom(400, 1, plogis(-2.2 + 0.6 * dq$xa))
fq <- glm(out ~ xa + db, family = binomial(), data = dq)
cat("fitted >= 0.5:", sum(fitted(fq) >= 0.5), "with y =", fq$y[fitted(fq) >= 0.5], "\n")
print(as.data.frame(run.all.gof(fq, tests = "Stukel", control = list(Stukel = list(form = "lr")))))
print(as.data.frame(run.all.gof(fq, tests = "Stukel")))

## ---- R4. weights = "score" on a probit fit: not the probit Rao test, but still a calibrated test ---------------------
cat("\n==== R4 ====\n")
set.seed(3)
dp <- data.frame(xa = runif(300, -3, 3)); dp$out <- rbinom(300, 1, pnorm(0.5 * dp$xa))
fp <- glm(out ~ xa, family = binomial("probit"), data = dp, control = glm.control(epsilon = 1e-13, maxit = 100))
nG <- 10
php <- pmin(pmax(fitted(fp), 1e-6), 1 - 1e-6)
grp <- pmin(ceiling(rank(php, ties.method = "first") / (300 / nG)), nG)
pbar <- tapply(php, grp, mean); e <- qlogis(pbar)
S <- (e * abs(e))[grp]
fp1 <- glm(out ~ xa + S, family = binomial("probit"), data = dp, start = c(coef(fp), 0),
           control = glm.control(epsilon = 1e-13, maxit = 100))
cat("def.gof sym score on the probit fit:", def.gof(fp, G = nG, basis = "sym", weights = "score")$Test_Statistic,
    "| probit anova Rao for the same step covariate:", anova(fp, fp1, test = "Rao")$Rao[2], "\n")
set.seed(20260914)
B <- 1000; pv <- matrix(NA_real_, B, 2, dimnames = list(NULL, c("sym.score", "poly3.score")))
for (b in 1:B) {
  xa <- runif(1000, -3, 3); db <- rbinom(1000, 1, 0.5); yb <- rbinom(1000, 1, pnorm(0.6 * xa + 0.5 * db))
  fb <- suppressWarnings(glm(yb ~ xa + db, family = binomial("probit")))
  pv[b, 1] <- def.gof(fb, G = "auto", basis = "sym", weights = "score")$p_value
  pv[b, 2] <- def.gof(fb, G = "auto", basis = "poly3", weights = "score")$p_value
}
cat(sprintf("probit truth, probit fit, n 1000, G auto, B %d: size05 sym score %.4f, poly3 score %.4f (MCSE %.4f)\n",
            B, mean(pv[, 1] <= 0.05), mean(pv[, 2] <= 0.05), sqrt(0.05 * 0.95 / B)))
cat("\ndone\n")
