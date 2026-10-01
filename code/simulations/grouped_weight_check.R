## grouped_weight_check.R -- is it the WEIGHTING? Population share of the detectable misfit (rho)
## captured by one-column grouped shapes in the paper's unit-weighted form (z'r) and in the score form
## (z*sqrt(V))'r, against the ungrouped score probe, at G = 160 equal groups.
set.seed(11); N <- 400000L; G <- 160L
x <- runif(N, -3, 3); d <- rbinom(N, 1, 0.5); X <- cbind(1, x, d)
aoc <- function(e, p) ifelse(e < -30, 0, 1 - log1p(exp(e)) / p)
for (S in c(1, 2, 3)) for (lk in c("probit", "cauchit", "cloglog")) {
  et <- S * (0.6 * x + 0.5 * d)
  p <- switch(lk, probit = pnorm(et), cauchit = pcauchy(et), cloglog = 1 - exp(-exp(et)))
  fit <- suppressWarnings(glm.fit(X, p, family = quasibinomial()))
  ph <- fit$fitted.values; eh <- drop(X %*% fit$coefficients); w <- ph * (1 - ph); dv <- p - ph
  D <- mean(dv^2 / w); XWX <- crossprod(X, w * X)
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (N / G)), G)
  V <- as.numeric(tapply(w, g, sum)); mu <- as.numeric(tapply(dv, g, sum)) / sqrt(V)
  pb <- as.numeric(tapply(ph, g, mean)); eb <- qlogis(pb); U <- rowsum(w * X, g) / sqrt(V)
  rho_g <- function(z) { z <- matrix(z); zO <- crossprod(z) - crossprod(z, U) %*% solve(XWX, crossprod(U, z)); drop(crossprod(z, mu))^2 / drop(zO) / N / D }
  rho_u <- function(z) { zW <- crossprod(z, w * X); I <- sum(w * z^2) - drop(zW %*% solve(XWX, t(zW))); sum(z * dv)^2 / I / N / D }
  sh <- list(sym = function(e, q) e * abs(e), ao = function(e, q) aoc(e, q), cub = function(e, q) e^3)
  cat(sprintf("%-8s s=%g (fitted logit %5.1f..%4.1f, V_g max/min %6.0f) ", lk, S, min(eh), max(eh), max(V) / min(V)))
  for (nm in names(sh)) cat(sprintf("| %s unit %.3f score %.3f ungrouped %.3f ", nm, rho_g(sh[[nm]](eb, pb)),
                                   rho_g(sh[[nm]](eb, pb) * sqrt(V)), rho_u(sh[[nm]](eh, ph))))
  cat("\n")
}
