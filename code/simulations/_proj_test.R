## _proj_test.R -- Escanciano/Liu projection-based GOF test for logistic regression.
##
## PURPOSE (paper EDGE): a faithful, self-contained R implementation of the projection
## test ("proj") of Liu, Li, Chen, Haerdle & Liang (2024, Stat. Comput. 34:175), which is
## Escanciano's (2006) residual-marked-empirical-process (RMEP) omnibus test carried over
## to logistic regression with a model-based bootstrap (MBB) critical value.
##
## We implement it ONLY to obtain an honest wall-clock TIMING point for the compute-power
## frontier (Fig 5). We do NOT run its power -- Liu et al. (2024) already establish that proj
## is the most powerful test in their study, including the purely off-index quadratic (their
## Setting 2) that every grouped/Stukel test misses. proj is slow (their MBB is ~2.5 h per
## simulation setting) precisely because it RE-ESTIMATES beta inside every bootstrap replicate;
## that cost is what Fig 5 records.
##
## STATISTIC (Liu et al. 2024, eqs. following their Thm 1):
##   T_proj = (1/n^2) * sum_i sum_j eps_i eps_j A_ij,   A_ij = sum_l A_ijl,
##   A_ijl  = C_p * measure{ w in S^p : X_i'w <= X_l'w  and  X_j'w <= X_l'w },
## the sphere measure of directions along which both projected covariates fall below the
## anchor l (Escanciano 2006, eq. for the indicator-product integral). For centered vectors
## u = X_i - X_l and v = X_j - X_l at angle theta, that half-space-overlap measure is
## proportional to (pi - theta), i.e.
##   A_ijl \propto 1 - arccos( u'v / (||u|| ||v||) ) / pi.
## This "(pi - angle)/pi" overlap form is the positive-semidefinite kernel that makes T_proj a
## genuine Cramer-von Mises norm (verified: A has no negative eigenvalues; the raw "angle/pi"
## print of the paper is its complement and is not PSD). eps_i = y_i - fitted_i are the raw
## (response) residuals. The constant C_p depends only on the covariate dimension; because the
## MBB p-value compares T_proj to bootstrap replicates computed with the SAME fixed A, C_p
## cancels exactly, so we set C_p = 1.
##
## The matrix A is a FIXED function of the design X only (not of y), so it is built ONCE and
## reused across all B bootstrap replicates. The refit-in-bootstrap cost lives entirely in the
## glm() re-estimation, not in A -- faithful to Liu's algorithm.
##
## No external packages are required (base R only).

## ---- build the fixed n x n matrix A from the design X (with intercept) -------------------
## Vectorized over l: for each anchor point l we center all rows (C = X - X_l), form the
## Gram matrix G = C C^T, normalize by the outer product of row norms to get cosines, clamp
## to [-1, 1], and accumulate the overlap kernel 1 - acos(cos)/pi into A. Degenerate terms
## (i == l or a zero-norm centered vector, where the cosine is undefined) contribute the
## neutral value 0.5 (= 1 - arccos(0)/pi), the "orthogonal / no information" limit.
proj_build_A <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  A <- matrix(0.0, n, n)
  for (l in seq_len(n)) {
    C  <- sweep(X, 2L, X[l, ], "-")      # C_i = X_i - X_l   (n x p)
    nrm <- sqrt(rowSums(C * C))          # ||C_i||           (length n)
    G  <- tcrossprod(C)                  # G_ij = C_i . C_j  (n x n)
    den <- outer(nrm, nrm)               # ||C_i|| ||C_j||
    cosv <- G / den                      # cosine of the angle at anchor l
    cosv[!is.finite(cosv)] <- 0.0        # zero-norm rows (incl. i==l): undefined -> 0 -> 0.5
    cosv[cosv >  1] <-  1                 # clamp for numerical safety
    cosv[cosv < -1] <- -1
    A <- A + (1 - acos(cosv) / pi)       # (pi - angle)/pi overlap kernel; C_p = 1 (cancels)
  }
  A
}

## ---- the projection test with a model-based bootstrap p-value ----------------------------
## y : 0/1 response (length n)
## X : model matrix WITH intercept (n x p) -- the same design passed to glm
## B : number of MBB replicates (default 1000, per Liu et al.)
## Returns list(stat, p_value, B, T_boot).
proj_pvalue <- function(y, X, B = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  n <- length(y)
  stopifnot(nrow(X) == n)

  ## FIXED marking matrix -- built once, reused for every bootstrap replicate.
  A <- proj_build_A(X)

  ## fit the null (assumed) logistic model on the design as given (X already has intercept).
  fit0   <- suppressWarnings(glm.fit(X, y, family = binomial()))
  mu_hat <- fit0$fitted.values
  e      <- y - mu_hat
  Tstat  <- as.numeric(crossprod(e, A %*% e)) / n^2

  ## model-based bootstrap: resample y* ~ Bernoulli(mu_hat), REFIT, recompute the statistic
  ## with the SAME fixed A (Liu et al. 2024, MBB Steps 1-3).
  Tboot <- numeric(B)
  for (b in seq_len(B)) {
    ystar <- rbinom(n, 1L, mu_hat)
    fb    <- suppressWarnings(glm.fit(X, ystar, family = binomial()))
    eb    <- ystar - fb$fitted.values
    Tboot[b] <- as.numeric(crossprod(eb, A %*% eb)) / n^2
  }

  list(stat = Tstat,
       p_value = mean(Tboot >= Tstat),
       B = B,
       T_boot = Tboot)
}
