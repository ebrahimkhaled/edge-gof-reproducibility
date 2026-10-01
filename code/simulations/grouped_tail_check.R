## grouped_tail_check.R -- why do the grouped logit-scale probes (EDGE-stk, g.sym, g.ao) collapse on a
## probit truth at high discrimination (s=2: EDGE-stk size-adjusted power 0.004, g.sym 0.267 against the
## ungrouped u.sym 0.697 at n=4000)? Population computation: the grouped signal and its post-fit
## variance for one-column shapes evaluated (a) at logit(mean fitted risk), as the paper defines the
## group position, and (b) at the mean fitted logit, for G = n/25 at n = 4000 (160 equal groups).
set.seed(11); N <- 400000L; G <- 160L
x <- runif(N, -3, 3); d <- rbinom(N, 1, 0.5); X <- cbind(1, x, d)
for (S in c(1, 2)) for (lk in c("probit", "cauchit")) {
  et <- S * (0.6 * x + 0.5 * d); p <- if (lk == "probit") pnorm(et) else pcauchy(et)
  fit <- suppressWarnings(glm.fit(X, p, family = quasibinomial()))
  ph <- fit$fitted.values; eh <- drop(X %*% fit$coefficients); w <- ph * (1 - ph)
  D <- mean((p - ph)^2 / w)
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (N / G)), G)
  V <- as.numeric(tapply(w, g, sum)); mu <- as.numeric(tapply(p - ph, g, sum)) / sqrt(V)
  pb <- as.numeric(tapply(ph, g, mean)); eb <- qlogis(pb); em <- as.numeric(tapply(eh, g, mean))
  U <- rowsum(w * X, g) / sqrt(V); XWX <- crossprod(X, w * X)
  ncp1 <- function(z) { z <- matrix(z); zO <- crossprod(z) - crossprod(z, U) %*% solve(XWX, crossprod(U, z))
    drop(crossprod(z, mu))^2 / drop(zO) / N }
  shapes <- list(sym = function(e) e * abs(e), stkpos = function(e) e^2 * (e >= 0), cub = function(e) e^3,
                 ao = function(e) 1 - log1p(exp(e)) / plogis(e))
  cat(sprintf("\n%s s=%g  D=%.2e  |eta| range of fitted logit %.1f to %.1f; tail groups: logit(mean p) vs mean logit  g1 %.2f vs %.2f, gG %.2f vs %.2f\n",
    lk, S, D, min(eh), max(eh), eb[1], em[1], eb[G], em[G]))
  for (nm in names(shapes)) {
    ra <- ncp1(shapes[[nm]](eb)) / D; rb <- ncp1(shapes[[nm]](em)) / D
    cat(sprintf("   %-7s rho at logit(mean p) %.3f | at mean logit %.3f\n", nm, ra, rb))
  }
  cat(sprintf("   grouped matched (mu itself): rho %.3f\n", ncp1(mu / sqrt(1)) / D))
}
