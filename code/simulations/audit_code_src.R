## audit_code_src.R -- code-lens audit (read-only): print the package implementations the claims rest on,
## and run three small deterministic checks. Writes nothing.
suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR); library(mvtnorm)})

cat("\n==== ebrahim.gof:::gof_stukel ====\n"); print(ebrahim.gof:::gof_stukel)

cat("\n==== LogisticDx gof.glm, lines near Stukel ====\n")
g <- tryCatch(getS3method("gof", "glm", envir = asNamespace("LogisticDx")), error = function(e) NULL)
if (is.null(g)) g <- tryCatch(get("gof.glm", envir = asNamespace("LogisticDx")), error = function(e) NULL)
if (!is.null(g)) { f <- deparse(g); i <- grep("Sst|stukel|Stukel|za|zb", f)
  idx <- sort(unique(unlist(lapply(i, function(k) max(1, k - 10):min(length(f), k + 10)))))
  cat(paste(idx, f[idx]), sep = "\n") } else cat("gof.glm not found\n")

cat("\n==== givitiR internals ====\n")
print(givitiR:::givitiCalibrationTestComp); print(givitiR:::polynomialLogRegrFw); print(givitiR:::givitiStatCdf)

cat("\n==== does pmvnorm (dim 3) consume the RNG? ====\n")
R3 <- matrix(.5, 3, 3); diag(R3) <- 1
set.seed(7); a <- runif(1); set.seed(7); invisible(pmvnorm(lower = rep(-2, 3), upper = rep(2, 3), corr = R3)); b <- runif(1)
cat("runif after pmvnorm identical to without:", identical(a, b), "\n")
set.seed(7); p1 <- 1 - pmvnorm(lower = rep(-2.5, 3), upper = rep(2.5, 3), corr = R3)
set.seed(8); p2 <- 1 - pmvnorm(lower = rep(-2.5, 3), upper = rep(2.5, 3), corr = R3)
cat(sprintf("max3 p at T=2.5 under two seeds: %.7f %.7f (error est %.2e)\n", p1, p2, attr(p1, "error")))
cat("pmvnorm at very large T (p could be < 0?):", 1 - as.numeric(pmvnorm(lower = rep(-9, 3), upper = rep(9, 3), corr = R3)), "\n")

cat("\n==== does givitiCalibrationTest consume the RNG? ====\n")
set.seed(3); x <- runif(2000, -3, 3); y <- rbinom(2000, 1, plogis(0.6 * x)); ph <- fitted(glm(y ~ x, family = binomial()))
set.seed(9); a <- runif(1); set.seed(9); invisible(givitiCalibrationTest(y, ph, devel = "internal")); b <- runif(1)
cat("runif after GiViTI identical to without:", identical(a, b), "\n")

cat("\n==== algebra: score-form grouped probe == ungrouped score test of the step function z_{g(i)} ====\n")
set.seed(11); n <- 4000; G <- 160
x <- runif(n, -3, 3); d <- rbinom(n, 1, .5); X <- cbind(1, x, d); y <- rbinom(n, 1, plogis(2 * (0.6 * x + 0.5 * d)))
fit <- glm.fit(X, y, family = binomial()); ph <- fit$fitted.values; w <- ph * (1 - ph)
g <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
o <- as.numeric(tapply(y, g, sum)); e <- as.numeric(tapply(ph, g, sum)); V <- as.numeric(tapply(w, g, sum))
pb <- as.numeric(tapply(ph, g, mean)); eb <- qlogis(pb); r <- (o - e) / sqrt(V)
U <- rowsum(w * X, g) / sqrt(V); XWX <- crossprod(X, w * X)
ZOZ <- function(Z) crossprod(Z) - crossprod(Z, U) %*% solve(XWX, crossprod(U, Z))
z <- eb * abs(eb); Zt <- matrix(z * sqrt(V)); zi <- z[g]
lhs_u <- drop(crossprod(Zt, r)); rhs_u <- sum(zi * (y - ph))
lhs_I <- drop(ZOZ(Zt)); zW <- crossprod(zi, w * X); rhs_I <- sum(w * zi^2) - drop(zW %*% solve(XWX, t(zW)))
cat(sprintf("Z'r %.10f vs sum z_g(i)(y-p) %.10f | Z'OmegaZ %.8f vs ungrouped post-fit info %.8f\n", lhs_u, rhs_u, lhs_I, rhs_I))

cat("\n==== grouped Omega by brute force: Monte Carlo covariance of r at fixed X, refit every draw ====\n")
set.seed(12); n <- 1500; G <- 60; x <- runif(n, -3, 3); d <- rbinom(n, 1, .5); X <- cbind(1, x, d); p0 <- plogis(2 * (0.6 * x + 0.5 * d))
zs <- function(y) { f <- glm.fit(X, y, family = binomial()); ph <- f$fitted.values; w <- ph * (1 - ph)
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  o <- as.numeric(tapply(y, g, sum)); e <- as.numeric(tapply(ph, g, sum)); V <- as.numeric(tapply(w, g, sum))
  pb <- as.numeric(tapply(ph, g, mean)); eb <- qlogis(pb); r <- (o - e) / sqrt(V)
  U <- rowsum(w * X, g) / sqrt(V); XWX <- crossprod(X, w * X)
  ZOZ <- function(Z) crossprod(Z) - crossprod(Z, U) %*% solve(XWX, crossprod(U, Z))
  Zu <- matrix(eb * abs(eb)); Zs <- Zu * sqrt(V)
  c(unit = drop(crossprod(Zu, r)) / sqrt(drop(ZOZ(Zu))), score = drop(crossprod(Zs, r)) / sqrt(drop(ZOZ(Zs)))) }
M <- t(replicate(3000, zs(rbinom(n, 1, p0))))
cat(sprintf("z mean (unit, score) %.3f %.3f | var %.3f %.3f (SE of var ~ %.3f) | P(|z|>1.96) %.4f %.4f\n",
  mean(M[, 1]), mean(M[, 2]), var(M[, 1]), var(M[, 2]), sqrt(2 / 3000), mean(abs(M[, 1]) > 1.96), mean(abs(M[, 2]) > 1.96)))
