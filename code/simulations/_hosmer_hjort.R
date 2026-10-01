## _hosmer_hjort.R -- the weighted grouped ("decile-of-risk") goodness-of-fit tests of
## Hosmer & Hjort (2002), Statistics in Medicine 21:2723-2738, as a comparator arm.
##
## WHAT IS IMPLEMENTED (equation numbers: published paper first, the Oslo technical-report
## preprint of the same paper in brackets).
##
##   Process, grouped over g ordered risk groups of the fitted logit x'b:
##     W_{n,j} = (1/n) sum_i I(r_{j-1} < x_i'b <= r_j) w_i (y_i - pi_i) = (o_j - e_j)/n
##                                                                         eq (A1) [eq 2]
##   Covariance of the g-vector:  Omega = D - B' J^{-1} B,  with
##     d_j = (1/n) sum_{I_j} w_i^2 v_i,   b_j = (1/n) sum_{I_j} w_i v_i x_i,
##     J   = (1/n) X'VX,                   v_i = pi_i (1 - pi_i).
##   HLw  = sum_j (o_j - e_j)^2 / (n d_j)              eq (A4) [eq 5]  -> chi^2(g - 2)
##   X2w  = n W' Omega^- W  (full covariance)          eq (A2)/(A3) [eqs 3/4]
##                                                     -> chi^2(df), df = rank(Omega)
##   Weights  w = (I - H) z,  H = X (X'VX)^{-1} X'V    eq (2) [eq 7], i.e. the residuals of
##   the V-weighted linear regression of z on X (the added-variable-plot co-ordinate).
##
## WHICH WEIGHTS. The paper's general-purpose (no specific omitted covariate) weight is
## z_i = pi_i ln(pi_i), the score direction of the one-parameter generalisation
## pi^(1+gamma) of the logistic model (eq 3) [eq 8]. Its recommendation for practice (Section
## 6, item 3) is the optimally weighted test HLop, which needs a user-chosen omitted
## covariate z (e.g. x^2 or x*d). So hh_test returns
##   HL1, X2_1    w = 1           (Table I rows 4-5; HL1 ~ the ordinary C-hat)
##   HLnp, X2np   z = pi ln pi    (Table I rows 6-7)
##   HLop, X2op   z supplied      (Table I rows 10-11) -- only when `z` is given.
## The partial-sum test Ww (eq A5) is NOT implemented: it is ungrouped and needs a refitting
## parametric bootstrap (M = 80 refits), which is outside the "well under a second at
## n = 5000" budget for an arm and is not a grouped test.
##
## REFERENCE DISTRIBUTIONS AND THE PAPER'S AMBIGUITIES (how each is resolved here).
##  1. HLw uses chi^2(g - 2), as the paper does (Appendix, after eq A4): the paper takes it
##     from the Hosmer-Lemeshow C-hat theory and shows by simulation that it is close.
##  2. X2w: the paper says "each weight function yields G of full rank, thus chi^2(g)". That
##     holds for z = pi ln pi (Omega is then non-singular), but NOT for w = 1: with an
##     intercept in the model the grouped raw residuals sum to exactly 0 (the score
##     equation), so Omega 1 = 0 and rank(Omega) = g - 1. We therefore follow the general
##     statement of eq (A2) -- a generalised inverse with df = rank(Omega) -- computed from an
##     eigen-decomposition with a relative tolerance. This gives chi^2(g) for np/op weights
##     and chi^2(g - 1) for w = 1.
##  3. The printed closed form (A3) has a scaling slip (the second term lacks a factor n in
##     the published version; the preprint divides by sqrt(n d_j)). We avoid it by evaluating
##     the quadratic form (A2) directly on the un-scaled sums: r = o - e and
##     n*Omega = diag(sum_{I_j} w^2 v) - Bt (X'VX)^{-1} Bt', Bt_j = sum_{I_j} w v x.
##     (The Woodbury form with G = J - B D^{-1} B' agrees with this to machine precision when
##     Omega is non-singular; checked in the validation script.)
##  4. Groups: the cutpoint r_j is the logit of the (n j / g)-th ordered fitted value. We
##     order ascending (the grouping is symmetric either way) and use floor(n j / g). Groups
##     are value-defined (r_{j-1} < x'b <= r_j), so tied fitted values never straddle a cut;
##     a cut that repeats leaves an empty group, which is dropped and the df shrink with it.
##  5. The weights are treated as fixed (their dependence on beta-hat is ignored), exactly as
##     in the paper's Omega.
##
## Base R + stats only. `fit` must be a binomial glm with the logit link and a 0/1 response.

## ---- one weighted grouped test pair (HLw, X2w) for a given weight vector -----------------
## r   : grouped weighted residual sums o_j - e_j          (length G)
## dn  : n * d_j = sum_{I_j} w^2 v                          (length G)
## Bt  : n * b_j' = sum_{I_j} w v x'                        (G x p)
## XtVX: n * J = X'VX                                       (p x p)
hh_grouped_pair <- function(w, res, v, X, grp, XtVX_inv, tol = 1e-8) {
  r   <- rowsum(w * res, grp, reorder = TRUE)[, 1]
  dn  <- rowsum(w * w * v, grp, reorder = TRUE)[, 1]
  Bt  <- rowsum((w * v) * X, grp, reorder = TRUE)
  G   <- length(r)

  ## HLw (A4): diagonal-variance version; chi^2(G - 2) as in the paper.
  HL    <- sum(r^2 / dn)
  p_HL  <- pchisq(HL, df = G - 2, lower.tail = FALSE)

  ## X2w (A2): full covariance n*Omega = D - B'J^{-1}B, generalised inverse, df = rank.
  Om  <- diag(dn, G) - Bt %*% XtVX_inv %*% t(Bt)
  Om  <- (Om + t(Om)) / 2                        # symmetrise against round-off
  eg  <- eigen(Om, symmetric = TRUE)
  keep <- eg$values > tol * max(eg$values)       # drops the exact zero of the w = 1 case
  u   <- crossprod(eg$vectors[, keep, drop = FALSE], r)
  X2  <- sum(u^2 / eg$values[keep])
  df  <- sum(keep)
  p_X2 <- pchisq(X2, df = df, lower.tail = FALSE)

  c(p_HL = p_HL, p_X2 = p_X2, HL = HL, X2 = X2, df_X2 = df, G = G)
}

## ---- main entry ----------------------------------------------------------------------------
## fit     : fitted glm(binomial(link = "logit")) with a 0/1 response
## g       : number of risk groups (paper: 10)
## z       : optional "omitted covariate" (length n) for the optimally weighted HLop/X2op,
##           e.g. x^2 or x*d. Must be aligned with the rows used in `fit`.
## details : if TRUE return a list with the statistics and df as well
## Returns a named numeric vector of p-values:
##   HL1, X2_1, HLnp, X2np [, HLop, X2op]
hh_test <- function(fit, g = 10, z = NULL, details = FALSE) {
  if (!inherits(fit, "glm") || fit$family$family != "binomial" ||
      fit$family$link != "logit")
    stop("hh_test: `fit` must be a binomial glm with the logit link")
  y <- fit$y
  if (any(y != 0 & y != 1) || any(fit$prior.weights != 1))
    stop("hh_test: needs a 0/1 response with unit prior weights (strictly binary case)")

  X   <- model.matrix(fit)
  X   <- X[, !is.na(coef(fit)), drop = FALSE]    # drop aliased columns
  pi_ <- fit$fitted.values
  n   <- length(y)
  v   <- pi_ * (1 - pi_)
  res <- y - pi_

  ## Risk groups on the fitted values (equivalent to the logit, a monotone map).
  ps   <- sort(pi_)
  cuts <- ps[floor(n * seq_len(g - 1) / g)]
  grp  <- findInterval(pi_, cuts, left.open = TRUE) + 1L   # r_{j-1} < . <= r_j

  XtVX_inv <- solve(crossprod(X, v * X))

  ## Optimal weights (2): residuals of the v-weighted LS regression of z on X.
  opt_w <- function(zz) lm.wfit(X, zz, v)$residuals

  w1  <- rep(1, n)
  wnp <- opt_w(pi_ * log(pi_))

  a1  <- hh_grouped_pair(w1,  res, v, X, grp, XtVX_inv)
  anp <- hh_grouped_pair(wnp, res, v, X, grp, XtVX_inv)
  out <- c(HL1 = a1[["p_HL"]], X2_1 = a1[["p_X2"]],
           HLnp = anp[["p_HL"]], X2np = anp[["p_X2"]])
  stats <- rbind(w1 = a1, np = anp)

  if (!is.null(z)) {
    if (length(z) != n) stop("hh_test: `z` must have one value per fitted observation")
    aop <- hh_grouped_pair(opt_w(as.numeric(z)), res, v, X, grp, XtVX_inv)
    out <- c(out, HLop = aop[["p_HL"]], X2op = aop[["p_X2"]])
    stats <- rbind(stats, op = aop)
  }

  if (details) return(list(p = out, stats = stats))
  out
}
