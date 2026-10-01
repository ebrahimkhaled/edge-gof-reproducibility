## clinical_visibility.R -- how VISIBLE is each link misfit on the risk scale, beside how many patients a test needs?
## Population (N = 400,000): pseudo-true logistic fit pi*; the calibration curve E[p_F | pi*] in 100 quantile bins;
## ICI = mean |p_F - pi*| (Austin & Steyerberg's integrated calibration index), E90 and Emax of the binned curve,
## and the worst relative error of the risk in the lowest decile. n80 for the symmetric / AO probes and the matched
## direction come from the same local formula as scout_envelope.R (lambda80 / (rho D)); GiViTI is NOT given (its
## scout column used an unconditional critical value -- audit bug B1).
set.seed(20260920); N <- 400000L
lam80 <- uniroot(function(l) pchisq(qchisq(.95, 1), 1, ncp = l, lower.tail = FALSE) - .8, c(.01, 200))$root
lnconv <- function(e, s) { K <- 40L; w <- qnorm((seq_len(K) - 0.5) / K); o <- numeric(length(e)); for (k in 1:K) o <- o + plogis(e + s * w[k]); o / K }
LINKS <- list(probit = pnorm, t4 = function(e) pt(e, 4), cauchit = pcauchy, lnorm2 = function(e) lnconv(e, 2),
              cloglog = function(e) 1 - exp(-exp(e)), loglog = function(e) exp(-exp(-e)))
aoc <- function(e, p) ifelse(e < -30, 0, 1 - log1p(exp(e)) / p)
out <- list()
for (des in list(c(1, 0, 0), c(2, 0, 0), c(1, -2, 0), c(1, 0, 1))) {
  s <- des[1]; c0 <- des[2]; normx <- des[3] == 1
  x <- if (normx) rnorm(N, 0, 1.5) else runif(N, -3, 3); d <- rbinom(N, 1, 0.5); X <- cbind(1, x, d)
  for (lk in names(LINKS)) {
    et <- c0 + s * (0.6 * x + 0.5 * d); p <- LINKS[[lk]](et)
    fit <- suppressWarnings(glm.fit(X, p, family = quasibinomial()))
    ph <- fit$fitted.values; eh <- drop(X %*% fit$coefficients); w <- ph * (1 - ph); dv <- p - ph
    D <- mean(dv^2 / w); XWX <- crossprod(X, w * X)
    rho1 <- function(z) { zW <- crossprod(z, w * X); I <- sum(w * z^2) - drop(zW %*% solve(XWX, t(zW))); sum(z * dv)^2 / I / N / D }
    b <- cut(ph, unique(quantile(ph, seq(0, 1, 0.01))), include.lowest = TRUE)
    cal <- tapply(p, b, mean) - tapply(ph, b, mean)
    dec1 <- ph <= quantile(ph, 0.1)
    out[[length(out) + 1]] <- data.frame(design = sprintf("s=%g c0=%g x=%s", s, c0, if (normx) "normal" else "uniform"),
      link = lk, event = round(mean(p), 3), ICI = mean(abs(dv)), E90 = unname(quantile(abs(dv), .9)), Emax_binned = max(abs(cal)),
      rel_err_low_decile = mean(p[dec1]) / mean(ph[dec1]) - 1,
      n80_matched = lam80 / D, n80_sym = lam80 / (rho1(eh * abs(eh)) * D), n80_ao = lam80 / (rho1(aoc(eh, ph)) * D),
      n80_aoM = lam80 / (rho1(aoc(-eh, 1 - ph)) * D))
  }
}
R <- do.call(rbind, out); write.csv(R, "clinical_visibility.csv", row.names = FALSE)
for (dsg in unique(R$design)) { cat("\n--", dsg, "--\n")
  z <- R[R$design == dsg, ]
  for (i in seq_len(nrow(z))) with(z[i, ], cat(sprintf("  %-8s event %.2f | ICI %.4f  E90 %.4f  Emax %.4f  low-decile risk off by %+5.1f%% | n80 matched %7.0f  sym %8.0f  ao %8.0f  aoM %8.0f\n",
    link, event, ICI, E90, Emax_binned, 100 * rel_err_low_decile, n80_matched, n80_sym, n80_ao, n80_aoM))) }
