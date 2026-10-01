## _pstar_giviti_harness.R -- shared harness for the three 2026-09-10 runs:
##   A. the p* sweep  (Reviewer 1, point 3: the number of covariates is never varied)
##   B. GiViTI added to the comparison (Reviewer 1, point 2: a rival cited 10x and never run)
##   C. the EF gate   (does a second-order n x G story exist at all?)
##
## Tests are called individually, not through run.all.gof: the full battery costs 2.7 s/call
## because of the covariate-space tests, the individual functions cost 0.1-0.3 ms.
## HL and HL_w are implemented here and checked against ResourceSelection::hoslem.test in
## _harness_selfcheck.R before any run is trusted.

suppressPackageStartupMessages({
  library(ebrahim.gof)
  library(givitiR)
})

## ---- Hosmer-Lemeshow on G equal-frequency groups of fitted risk -------------------------
hl_stat <- function(y, p, G = 10, equal_width = FALSE) {
  br <- if (equal_width) seq(0, 1, length.out = G + 1)
        else stats::quantile(p, probs = seq(0, 1, length.out = G + 1), na.rm = TRUE)
  br[1] <- -Inf; br[length(br)] <- Inf
  g <- cut(p, breaks = unique(br), include.lowest = TRUE)
  if (nlevels(droplevels(g)) < 3) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  o <- tapply(y, g, sum); e <- tapply(p, g, sum); n <- tapply(y, g, length)
  o <- o[!is.na(o)]; e <- e[!is.na(e)]; n <- n[!is.na(n)]
  v <- e * (1 - e / n)
  keep <- v > 1e-8
  X2 <- sum((o[keep] - e[keep])^2 / v[keep])
  df <- sum(keep) - 2
  if (df < 1) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  c(stat = X2, df = df, p = stats::pchisq(X2, df, lower.tail = FALSE))
}

## ---- Stukel's score test: refit with the two sign-split eta^2 terms ---------------------
## Returns NA when the auxiliary fit separates or fails -- the failure IS the datum.
stukel_p <- function(fit) {
  e <- stats::predict(fit)
  d <- fit$model
  d$.za <- 0.5 * e^2 * (e >= 0)
  d$.zb <- -0.5 * e^2 * (e < 0)
  sep <- FALSE
  fa <- tryCatch(withCallingHandlers(
      stats::glm(stats::update(stats::formula(fit), . ~ . + .za + .zb), data = d,
                 family = stats::binomial(), control = list(maxit = 50)),
      warning = function(w) {
        if (grepl("did not converge|numerically 0 or 1", conditionMessage(w))) sep <<- TRUE
        invokeRestart("muffleWarning")
      }), error = function(er) NULL)
  if (is.null(fa) || sep) return(NA_real_)
  lr <- stats::deviance(fit) - stats::deviance(fa)
  if (!is.finite(lr) || lr < 0) return(NA_real_)
  stats::pchisq(lr, df = 2, lower.tail = FALSE)
}

## ---- GiViTI calibration test (internal validation = the development data) ---------------
giviti_p <- function(fit) {
  out <- tryCatch(givitiR::givitiCalibrationTest(fit$y, stats::fitted(fit), devel = "internal"),
                  error = function(e) NULL, warning = function(w) NULL)
  if (is.null(out)) NA_real_ else as.numeric(out$p.value)
}

## ---- one replication: every p-value we care about, from ONE fit -------------------------
## Returns a named numeric vector. NA means the test declined to run on this sample; that is
## recorded, never silently dropped.
one_rep <- function(dat, G = 10) {
  fit <- tryCatch(stats::glm(dat$f, data = dat$d, family = stats::binomial(),
                             control = list(maxit = 50)),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(EDGE.poly2 = NA, EDGE.poly3 = NA, EDGE.stk = NA, HL_F = NA,
                             HL = NA, HL_w = NA, Stukel = NA, GiViTI = NA, pstar = NA))
  gp <- function(basis) tryCatch(edge.gof(fit, G = G, basis = basis)$p_value,
                                 error = function(e) NA_real_)
  efp <- tryCatch(ef.gof(fit, G = G)$p_value, error = function(e) NA_real_)
  y <- fit$y; p <- stats::fitted(fit)
  c(EDGE.poly2 = gp("poly2"), EDGE.poly3 = gp("poly3"), EDGE.stk = gp("stukel"),
    HL_F = efp,
    HL   = unname(hl_stat(y, p, G, FALSE)["p"]),
    HL_w = unname(hl_stat(y, p, G, TRUE)["p"]),
    Stukel = stukel_p(fit),
    GiViTI = giviti_p(fit),
    pstar = length(stats::coef(fit)))
}

TESTS <- c("EDGE.poly2", "EDGE.poly3", "EDGE.stk", "HL_F", "HL", "HL_w", "Stukel", "GiViTI")

## ---- size-adjusted rejection: use the null run's empirical alpha-quantile of p ----------
## Reports BOTH raw and size-adjusted power, because a comparison on raw p-values against a
## test whose realised size differs is not a comparison.
size_adjusted <- function(p_alt, p_null, alpha = 0.05) {
  p_null <- p_null[is.finite(p_null)]
  if (length(p_null) < 50) return(NA_real_)
  crit <- stats::quantile(p_null, probs = alpha, na.rm = TRUE, type = 1)
  mean(p_alt[is.finite(p_alt)] <= crit)
}

raw_power <- function(p, alpha = 0.05) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  mean(p <= alpha)
}

na_rate <- function(p) mean(!is.finite(p))

mcse <- function(rate, B) ifelse(is.na(rate), NA_real_, sqrt(rate * (1 - rate) / B))
