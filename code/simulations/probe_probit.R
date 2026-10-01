## probe_probit.R -- IS probit-vs-logit detectable in principle, or is it an identifiability wall?
##
## The paper reports every test at nominal size against a probit truth and calls it an
## identifiability limit. Before building a probit-specific basis it is worth knowing WHICH of two
## worlds we are in:
##
##   (A) the signal is genuinely ~zero  -> no basis can help, the wall is real;
##   (B) the signal is small but lives in a SHAPE the poly3 basis does not span (e.g. 4th order)
##       -> a targeted basis could find it, and the current failure is a basis choice, not a law.
##
## This script answers that by computing the POPULATION calibration deviation, with beta refitted
## at its logit pseudo-true value, then decomposing it onto orthonormal polynomials in the
## calibration position. No sampling noise: this is the signal that would be there at infinite n.

suppressPackageStartupMessages(library(stats))

## --- the true DGP: probit link on eta = 0.6 x + 0.5 d, x ~ U(-3,3), d ~ Bern(.5) ------------
## Everything is done on a fine deterministic grid, so this is exact up to quadrature.
NG <- 40001
x <- seq(-3, 3, length.out = NG)
wx <- rep(1 / NG, NG)

grid <- rbind(cbind(x, 0), cbind(x, 1))
w <- c(wx, wx) / 2
eta_true <- 0.6 * grid[, 1] + 0.5 * grid[, 2]

truep <- function(scn, e) switch(scn,
  probit  = pnorm(e),
  cloglog = 1 - exp(-exp(e)),
  loglog  = exp(-exp(-e)),
  logit   = plogis(e),
  stop(scn))

## --- the pseudo-true logit fit: the beta a logistic model converges to under this truth -------
## Maximises the expected log-likelihood, i.e. solves E[w (p_true - plogis(Xb)) X] = 0.
pseudo_true <- function(p_true) {
  X <- cbind(1, grid[, 1], grid[, 2])
  b <- c(0, 0.6, 0.5)
  for (it in 1:200) {
    mu <- plogis(as.vector(X %*% b))
    W <- w * mu * (1 - mu)
    score <- crossprod(X, w * (p_true - mu))
    H <- crossprod(X * W, X)
    step <- solve(H, score)
    b <- b + as.vector(step)
    if (max(abs(step)) < 1e-12) break
  }
  b
}

## --- orthonormal polynomials in the calibration position, in the w-weighted L2 -----------------
orthopoly <- function(pos, w, K) {
  B <- matrix(0, length(pos), K)
  u <- 2 * (pos - min(pos)) / (max(pos) - min(pos)) - 1     # map to [-1,1]
  for (k in 1:K) {
    v <- u^k
    if (k > 1) for (j in 1:(k - 1)) v <- v - sum(w * v * B[, j]) * B[, j]
    v <- v - sum(w * v)                                      # orthogonal to the constant
    B[, k] <- v / sqrt(sum(w * v^2))
  }
  B
}

cat("=====================================================================\n")
cat(" POPULATION calibration deviation after the logit refit, decomposed\n")
cat(" onto orthonormal polynomial shapes of the fitted risk.\n")
cat(" delta(pi) = p_true - p_logit_pseudotrue.  Units: probability points.\n")
cat("=====================================================================\n")

K <- 6
res <- list()
for (scn in c("probit", "cloglog", "loglog")) {
  pt <- truep(scn, eta_true)
  b <- pseudo_true(pt)
  pf <- plogis(as.vector(cbind(1, grid[, 1], grid[, 2]) %*% b))
  delta <- pt - pf

  ## weight by the Fisher information of the fit, which is what a residual test actually sees
  wv <- w * pf * (1 - pf)
  wv <- wv / sum(wv)
  B <- orthopoly(pf, wv, K)
  coef <- as.vector(crossprod(B * wv, delta / sqrt(pf * (1 - pf))))
  tot <- sqrt(sum(coef^2))

  cat(sprintf("\n--- %s truth, logit fit ---\n", toupper(scn)))
  cat(sprintf("  pseudo-true beta : (%.4f, %.4f, %.4f)   [truth was (0, 0.6, 0.5)]\n", b[1], b[2], b[3]))
  cat(sprintf("  max |delta|      : %.5f  (%.3f percentage points)\n", max(abs(delta)), 100 * max(abs(delta))))
  cat(sprintf("  total signal norm: %.6f\n", tot))
  cat("  share of the signal by polynomial order:\n")
  for (k in 1:K)
    cat(sprintf("     P%d %-9s %8.5f   %5.1f%%  %s\n", k,
                c("(tilt)", "(bow)", "(S)", "(W)", "(5th)", "(6th)")[k],
                coef[k], 100 * coef[k]^2 / sum(coef^2),
                strrep("#", round(60 * coef[k]^2 / sum(coef^2)))))
  cat(sprintf("  captured by a poly3 basis (P1..P3): %5.1f%%\n", 100 * sum(coef[1:3]^2) / sum(coef^2)))
  cat(sprintf("  captured by a poly4 basis (P1..P4): %5.1f%%\n", 100 * sum(coef[1:4]^2) / sum(coef^2)))
  res[[scn]] <- list(coef = coef, tot = tot, maxd = max(abs(delta)))
}

cat("\n=====================================================================\n")
cat(" THE COMPARISON THAT MATTERS: how big is probit's signal relative to\n")
cat(" a departure every test DOES detect?\n")
cat("=====================================================================\n")
cat(sprintf("  cloglog total signal norm : %.6f\n", res$cloglog$tot))
cat(sprintf("  loglog  total signal norm : %.6f\n", res$loglog$tot))
cat(sprintf("  probit  total signal norm : %.6f\n", res$probit$tot))
cat(sprintf("\n  probit is %.1fx weaker than cloglog, %.1fx weaker than loglog.\n",
            res$cloglog$tot / res$probit$tot, res$loglog$tot / res$probit$tot))
cat(sprintf("  Non-centrality scales with the SQUARE of the norm, so at equal n the\n"))
cat(sprintf("  probit non-centrality is about %.0fx smaller than cloglog's.\n",
            (res$cloglog$tot / res$probit$tot)^2))
