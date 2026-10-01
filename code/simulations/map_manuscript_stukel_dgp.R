## map_manuscript_stukel_dgp.R -- READ-ONLY check made while mapping the manuscript (2026-09-13).
## Question: are the "Stukel link" DGPs of the paper (inv_stukel in _dgp_library.R) members of
## Stukel's (1988) generalised logistic h-family, as sections/05_simulation.tex and 09_appendix.tex say?
## No simulation, no files written; prints a table only.

## verbatim copy of the harness generator (_dgp_library.R lines 20-23)
inv_stukel <- function(ev, a1, a2) { z <- numeric(length(ev)); pos <- ev >= 0
  if (any(pos))  z[pos]  <- if (abs(a1) < 1e-12) ev[pos]  else (-1 + sqrt(pmax(0, 1 + 2 * a1 * ev[pos])))  / a1
  if (any(!pos)) z[!pos] <- if (abs(a2) < 1e-12) ev[!pos] else (-1 + sqrt(pmax(0, 1 + 2 * a2 * ev[!pos]))) / a2
  plogis(z) }

## Stukel (1988) h-family, p = plogis(h(eta))
stukel_h <- function(ev, a1, a2) {
  h <- numeric(length(ev)); pos <- ev >= 0; ae <- abs(ev)
  hp <- function(a, t) if (abs(a) < 1e-12) t else if (a > 0) (exp(a * t) - 1) / a else -log(1 - a * t) / a
  h[pos]  <- hp(a1, ae[pos])
  h[!pos] <- -hp(a2, ae[!pos])
  plogis(h) }

eta <- c(-2.3, -1.5, -1, -0.5, 0, 0.25, 0.5, 1, 1.5, 2.3)
cells <- list(heavy = c(-1, -1), light = c(1, 1), asym = c(-1, 1))
for (nm in names(cells)) {
  a <- cells[[nm]]
  cat(sprintf("\n== stukel_%s (a1, a2) = (%g, %g) ==\n", nm, a[1], a[2]))
  cat(sprintf("%6s %10s %10s %12s %12s\n", "eta", "logit(p)", "harness p", "Stukel88 lp", "Stukel88 p"))
  ph <- inv_stukel(eta, a[1], a[2]); ps <- stukel_h(eta, a[1], a[2])
  for (i in seq_along(eta)) cat(sprintf("%6.2f %10.3f %10.3f %12.3f %12.3f\n", eta[i], qlogis(ph[i]), ph[i], qlogis(ps[i]), ps[i]))
}

## the base design's range of eta = 0.6 x + 0.5 d, x ~ U(-3, 3): [-1.8, 2.3]
set.seed(1); x <- runif(2e5, -3, 3); d <- rbinom(2e5, 1, 0.5); e <- 0.6 * x + 0.5 * d
for (nm in names(cells)) {
  a <- cells[[nm]]; p <- inv_stukel(e, a[1], a[2])
  cat(sprintf("\nharness stukel_%s on base design: p range [%.3f, %.3f]; share of patients with p at its cap = %.3f; event rate %.3f\n",
      nm, min(p), max(p), mean(abs(p - max(p)) < 1e-12 | abs(p - min(p)) < 1e-12), mean(p)))
  ## local second-order carrier coefficients of logit(p) - eta, by least squares on za, zb within |eta| <= 0.4
  k <- abs(e) <= 0.4; za <- 0.5 * e^2 * (e >= 0); zb <- -0.5 * e^2 * (e < 0)
  cf <- coef(lm(I(qlogis(p) - e) ~ 0 + za + zb, subset = k))
  cat(sprintf("   local Stukel carrier coefficients near eta = 0: alpha1 %.2f, alpha2 %.2f\n", cf[1], cf[2]))
}
