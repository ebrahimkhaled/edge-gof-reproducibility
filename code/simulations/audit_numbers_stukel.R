## audit_numbers_stukel.R -- independent replication (new seed, 5 workers) of the Stukel implementation check:
## size of the marginal-sum statistic (ebrahim.gof / LogisticDx "SstBoth"), the joint score u'I^-1 u and the
## likelihood-ratio refit, post-fit correlation of the two score directions, and refit failure rates.
suppressPackageStartupMessages(library(parallel))
one <- function(b, design, n) {
  if (design == "base") { x <- runif(n, -3, 3); d <- rbinom(n, 1, .5); eta <- 0.6 * x + 0.5 * d; X <- cbind(1, x, d) }
  else { x <- as.numeric(scale(rchisq(n, 4))); eta <- -4 + 0.9 * x; X <- cbind(1, x) }
  y <- rbinom(n, 1, plogis(eta))
  fit <- suppressWarnings(glm.fit(X, y, family = binomial()))
  ph <- fit$fitted.values; eh <- drop(X %*% fit$coefficients); W <- ph * (1 - ph)
  za <- 0.5 * eh^2 * (eh >= 0); zb <- -0.5 * eh^2 * (eh < 0); Z <- cbind(za, zb)
  u <- drop(crossprod(Z, y - ph)); WZX <- crossprod(X, W * Z)
  I <- crossprod(Z, W * Z) - t(WZX) %*% solve(crossprod(X, W * X), WZX)
  ok <- all(is.finite(I)) && det(I) > 1e-12
  marg <- if (ok) sum(u^2 / diag(I)) else NA
  joint <- if (ok) drop(t(u) %*% solve(I, u)) else NA
  rho <- if (ok) I[1, 2] / sqrt(I[1, 1] * I[2, 2]) else NA
  sep <- FALSE
  fa <- tryCatch(withCallingHandlers(glm.fit(cbind(X, Z), y, family = binomial(), control = list(maxit = 50)),
        warning = function(w) { if (grepl("did not converge|numerically 0 or 1", conditionMessage(w))) sep <<- TRUE; invokeRestart("muffleWarning") }),
        error = function(e) NULL)
  lr <- if (is.null(fa) || sep) NA else fit$deviance - fa$deviance
  c(marg = marg, joint = joint, lr = lr, rho = rho, events = sum(y))
}
cl <- makeCluster(5); clusterExport(cl, "one"); clusterSetRNGStream(cl, 20260913)
wil <- function(k, m) { ci <- binom.test(k, m)$conf.int; sprintf("%.4f [%.4f, %.4f]", k / m, ci[1], ci[2]) }
for (cfg in list(c("base", 1000), c("base", 5000), c("sparse", 500), c("sparse", 2000))) {
  M <- t(parSapply(cl, 1:4000, function(b, d, n) one(b, d, n), d = cfg[1], n = as.integer(cfg[2])))
  M <- M[M[, "events"] >= 10, , drop = FALSE]
  cr <- qchisq(.95, 2)
  f <- function(v) { v <- v[is.finite(v)]; wil(sum(v > cr), length(v)) }
  cat(sprintf("%-6s n=%-5s reps=%d | marginal-sum %s | joint %s | LR %s | score NA %.4f  refit NA %.4f | corr(za,zb) mean %.3f sd %.3f\n",
    cfg[1], cfg[2], nrow(M), f(M[, "marg"]), f(M[, "joint"]), f(M[, "lr"]),
    mean(is.na(M[, "joint"])), mean(is.na(M[, "lr"])), mean(M[, "rho"], na.rm = TRUE), sd(M[, "rho"], na.rm = TRUE)))
}
stopCluster(cl)
## population post-fit correlation at the base design (N = 400,000, logistic truth)
set.seed(7); N <- 400000; x <- runif(N, -3, 3); d <- rbinom(N, 1, .5); X <- cbind(1, x, d); eta <- 0.6 * x + 0.5 * d
W <- plogis(eta) * (1 - plogis(eta)); Z <- cbind(0.5 * eta^2 * (eta >= 0), -0.5 * eta^2 * (eta < 0))
WZX <- crossprod(X, W * Z); I <- crossprod(Z, W * Z) - t(WZX) %*% solve(crossprod(X, W * X), WZX)
cat(sprintf("population post-fit corr(za,zb), base design: %.4f\n", I[1, 2] / sqrt(I[1, 1] * I[2, 2])))
## implied size of the marginal sum when corr = r: X = z1^2 + z2^2 with corr r, eigenvalues 1 +- |r|
for (r in c(-0.714, -0.70)) { lam <- c(1 + abs(r), 1 - abs(r)); z <- matrix(rnorm(2e6), ncol = 2)
  cat(sprintf("implied limit size of marginal sum at corr %.3f: %.4f\n", r, mean(lam[1] * z[, 1]^2 + lam[2] * z[, 2]^2 > qchisq(.95, 2)))) }
