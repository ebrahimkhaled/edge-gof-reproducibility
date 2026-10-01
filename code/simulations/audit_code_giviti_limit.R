## audit_code_giviti_limit.R -- code-lens audit (read-only, writes nothing).
## scout_envelope.R models GiViTI in the Gaussian limit as: statistic T = z2^2 + z3^2*1(gate3) + z4^2*1(gate3 & gate4),
## rejected above its UNCONDITIONAL 95% null quantile. givitiR instead reports p = 1 - F_m(T) with F_m the null cdf
## CONDITIONAL on the selected degree m (givitiStatCdf, internal mode), and rejects when p <= 0.05. This script puts both
## rules on the scout's own population u and I and compares them with the finite-sample GiViTI power from run L.
## It also reports the population post-fit correlation of Stukel's two columns under the logistic null and the
## limiting size of the marginal-sum statistic, and the GiViTI nuisance-space issue (projection on X vs on (1, eta)).
suppressPackageStartupMessages(library(givitiR))
set.seed(20260920); N <- 400000L
x <- runif(N, -3, 3); dd <- rbinom(N, 1, 0.5); X <- cbind(1, x, dd); lin <- 0.6 * x + 0.5 * dd
FL <- list(logit = plogis, probit = pnorm, t4 = function(e) pt(e, 4), cauchit = pcauchy,
           cloglog = function(e) 1 - exp(-exp(e)), loglog = function(e) exp(-exp(-e)))
k <- qchisq(.95, 1)
q3 <- uniroot(function(t) 1 - givitiR:::givitiStatCdf(t, 3, "internal", .95) - .05, c(k + 1e-6, 60))$root
q4 <- uniroot(function(t) 1 - givitiR:::givitiStatCdf(t, 4, "internal", .95) - .05, c(2 * k + 1e-3, 80))$root
Ec <- matrix(rnorm(3 * 2e6), ncol = 3)
g3 <- Ec[, 2]^2 > k; g4 <- g3 & Ec[, 3]^2 > k
cat(sprintf("conditional 95%% points: m=2 %.3f | m=3 %.3f (MC check %.3f) | m=4 %.3f (MC check %.3f)\n", k, q3,
  quantile((Ec[, 1]^2 + Ec[, 2]^2)[g3 & !g4], .95), q4, quantile(rowSums(Ec^2)[g4], .95)))
cat(sprintf("conditional rule, overall null size in the limit: %.4f\n",
  mean((!g3 & Ec[, 1]^2 > k) | (g3 & !g4 & Ec[, 1]^2 + Ec[, 2]^2 > q3) | (g4 & rowSums(Ec^2) > q4))))

E0 <- matrix(rnorm(3 * 200000), ncol = 3)
cellfun <- function(s, c0, lk) {
  eta <- c0 + s * lin; p <- FL[[lk]](eta)
  fit <- suppressWarnings(glm.fit(X, p, family = quasibinomial()))
  pis <- fit$fitted.values; eh <- drop(X %*% fit$coefficients); w <- pis * (1 - pis); dv <- p - pis
  infoA <- function(Z, A) { Z <- as.matrix(Z); ZW <- crossprod(Z, w * A); (crossprod(Z, w * Z) - ZW %*% solve(crossprod(A, w * A), t(ZW))) / N }
  list(eh = eh, w = w, dv = dv, infoA = infoA, D = mean(dv^2 / w))
}
rules <- function(u, I) {
  L <- t(chol(I))
  stats <- function(n) { S <- sweep(E0 %*% t(L), 2, sqrt(n) * u, "+"); Zw <- t(solve(L, t(S))); z2 <- Zw^2
    g3 <- z2[, 2] > k; g4 <- g3 & z2[, 3] > k; T <- z2[, 1] + z2[, 2] * g3 + z2[, 3] * g4
    list(T = T, m = 2 + g3 + g4) }
  s0 <- stats(0); crit <- quantile(s0$T, .95)
  pw_scout <- function(n) mean(stats(n)$T > crit)
  pw_cond <- function(n) { z <- stats(n); mean((z$m == 2 & z$T > k) | (z$m == 3 & z$T > q3) | (z$m == 4 & z$T > q4)) }
  n80 <- function(pw) { if (pw(1e8) < .8) return(Inf); exp(uniroot(function(t) pw(exp(t)) - .8, c(log(5), log(1e8)))$root) }
  wu <- drop(solve(L, u)); share <- wu^2 / sum(wu^2)
  list(pw_scout = pw_scout, pw_cond = pw_cond, n80_scout = n80(pw_scout), n80_cond = n80(pw_cond), share = share,
       size_scout = pw_scout(0), size_cond = pw_cond(0))
}

cat("\n=== GiViTI in the Gaussian limit: scout rule (unconditional critical value) vs givitiR rule (conditional on degree) ===\n")
CELLS <- list(list(1, 0, "probit", c(10000, 16000), c(0.243, 0.476)), list(1, 0, "cauchit", c(2500, 4000), c(0.237, 0.386)),
              list(2, 0, "cauchit", c(600, 900), c(0.254, 0.376)), list(2, 0, "t4", c(6000, 9000), c(0.221, 0.378)),
              list(2, 0, "probit", c(4000, 6500), c(0.269, 0.561)), list(1, -2, "probit", c(5000, 8000), c(0.708, 0.880)),
              list(1, 0, "cloglog", c(600, 1000), c(0.661, 0.859)), list(1, 0, "loglog", c(600, 1000), c(0.583, 0.786)))
for (cc in CELLS) {
  P <- cellfun(cc[[1]], cc[[2]], cc[[3]]); eh <- P$eh
  Z <- cbind(eh^2, eh^3, eh^4); u <- drop(crossprod(Z, P$dv)) / N
  IX <- P$infoA(Z, X); I1 <- P$infoA(Z, cbind(1, eh))
  R <- rules(u, IX)
  cat(sprintf("s=%g c0=%+g %-8s | n80 scout-rule %6.0f  givitiR-rule %6.0f | limit size %.3f / %.3f | whitened signal share deg2 %.2f deg3|2 %.2f deg4|3 %.2f | I_X/I_(1,eta) deg2 %.3f\n",
    cc[[1]], cc[[2]], cc[[3]], R$n80_scout, R$n80_cond, R$size_scout, R$size_cond, R$share[1], R$share[2], R$share[3], IX[1, 1] / I1[1, 1]))
  for (j in 1:2) cat(sprintf("      n=%-6d limit power scout-rule %.3f  givitiR-rule %.3f  | run L GiViTI (size-adjusted) %.3f\n",
    cc[[4]][j], R$pw_scout(cc[[4]][j]), R$pw_cond(cc[[4]][j]), cc[[5]][j]))
}

cat("\n=== Stukel columns under the LOGISTIC null (population): post-fit correlation and limiting size of the marginal sum ===\n")
for (dz in list(c(1, 0), c(1, -2), c(2, 0))) {
  P <- cellfun(dz[1], dz[2], "logit"); eh <- P$eh
  Zs <- cbind(eh^2 * (eh >= 0), -eh^2 * (eh < 0)); I <- P$infoA(Zs, X); rr <- cov2cor(I)[1, 2]
  zz <- cbind(E0[, 1], rr * E0[, 1] + sqrt(1 - rr^2) * E0[, 2])
  cat(sprintf("  s=%g c0=%+g  post-fit corr(za, zb) %.3f | limiting size at .05 of z1^2+z2^2 vs chi2_2: %.4f  (share of eta>=0: %.2f)\n",
    dz[1], dz[2], rr, mean(rowSums(zz^2) > qchisq(.95, 2)), mean(eh >= 0)))
}
