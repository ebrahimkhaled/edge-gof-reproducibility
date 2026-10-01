## grouped_weight_poly3.R -- the paper's DEFAULT basis (cubic in the group-mean risk, eigen form): how much of the
## detectable misfit does it keep in the unit-weighted form vs the score form (columns times sqrt(V_g))?
## Two-column version P2, P3 (the tilt P1 is nearly absorbed by the fit), non-centrality u'I^-1 u.
set.seed(11); N <- 400000L; G <- 160L
x <- runif(N, -3, 3); d <- rbinom(N, 1, 0.5); X <- cbind(1, x, d)
for (S in c(1, 2, 3)) for (lk in c("probit", "cauchit", "cloglog", "loglog")) {
  et <- S * (0.6 * x + 0.5 * d)
  p <- switch(lk, probit = pnorm(et), cauchit = pcauchy(et), cloglog = 1 - exp(-exp(et)), loglog = exp(-exp(-et)))
  fit <- suppressWarnings(glm.fit(X, p, family = quasibinomial()))
  ph <- fit$fitted.values; w <- ph * (1 - ph); dv <- p - ph; D <- mean(dv^2 / w); XWX <- crossprod(X, w * X)
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (N / G)), G)
  V <- as.numeric(tapply(w, g, sum)); mu <- as.numeric(tapply(dv, g, sum)) / sqrt(V)
  pb <- as.numeric(tapply(ph, g, mean)); U <- rowsum(w * X, g) / sqrt(V)
  rho <- function(Z) { Z <- as.matrix(Z); I <- crossprod(Z) - crossprod(Z, U) %*% solve(XWX, crossprod(U, Z))
    u <- crossprod(Z, mu); drop(t(u) %*% solve(I, u)) / N / D }
  P <- as.matrix(poly(pb, 3))
  cat(sprintf("%-8s s=%g | poly(P2,P3) unit %.3f score %.3f | poly(P1..P3) unit %.3f score %.3f\n", lk, S,
    rho(P[, 2:3]), rho(P[, 2:3] * sqrt(V)), rho(P), rho(P * sqrt(V))))
}
