## run_H_overconfidence.R -- can EDGE detect an OVERCONFIDENT classifier, and where does it beat
## the standard tools? EXTERNAL-VALIDATION mode: the predictions are frozen, nothing is refitted.
##
## Why this is the right experiment. The probit win showed EDGE-stk owns the S-shaped (cubic /
## sign-split) distortion. An S-shaped calibration curve is not exotic: it is the signature of an
## overconfident classifier, which pushes predictions toward 0 and 1 so that observed rates are
## LESS extreme than predicted at both ends. Modern ML classifiers do this routinely (Guo et al.
## 2017), and the standard fix, temperature scaling, is a symmetric stretch of the logit.
##
## The comparator that MUST be beaten, or the whole idea is dead: the Cox (1958) / Miller
## calibration-slope test -- fit y ~ a + b*logit(p), test (a,b) = (0,1). It is THE standard in
## clinical external validation, it is a 2-df LRT, and it is the exact likelihood-ratio test for
## a pure temperature distortion (logit(p') = logit(p)/T  <=>  slope b = T). So on pure
## temperature scaling the slope test should WIN and EDGE cannot beat an LRT for its own family.
##
## EDGE's opening is therefore where the S is NOT a pure logit rescaling:
##   (a) ASYMMETRIC overconfidence -- overconfident on one tail only. This is what class imbalance
##       does to a classifier: confident on the majority class, not on the minority. The slope
##       test is blind to it because a single slope cannot bend two ways.
##   (b) probit-type compression (a genuine S, not a rescaling).
## The Stukel basis is a sign-split quadratic in the logit, i.e. it is built for (a).
##
## In external mode there is no fitting correction: Omega = I_G, so EDGE's null is a PLAIN
## chi-squared_k. The constant and linear directions are NOT absorbed (nothing was fitted), so
## the basis includes them: k = 4. Note that the slope test IS the k = 2 sub-basis [1, eta].
##
## Every cell is size-adjusted against T = 1 (the null) at the same n, G, design.

suppressPackageStartupMessages({library(parallel); library(givitiR)})

SEED0 <- 20260916L
NS <- c(500L, 1000L, 2000L, 5000L)
B  <- 2000L

## ---- external-validation tests on frozen predictions p and outcomes y ----------------------
grp <- function(p, G) {
  br <- stats::quantile(p, probs = seq(0, 1, length.out = G + 1), type = 1)
  br[1] <- -Inf; br[length(br)] <- Inf
  as.integer(cut(p, breaks = unique(br), include.lowest = TRUE))
}
resid_ext <- function(y, p, G) {
  g <- grp(p, G); K <- max(g)
  o <- tapply(y, g, sum); e <- tapply(p, g, sum); v <- tapply(p * (1 - p), g, sum)
  pb <- tapply(p, g, mean)
  list(r = as.numeric((o - e) / sqrt(v)), eta = as.numeric(stats::qlogis(pb)), K = K)
}
edge_ext <- function(y, p, G, basis) {
  z <- resid_ext(y, p, G); eta <- z$eta; r <- z$r
  Z <- switch(basis,
    slope  = cbind(1, eta),                                          # = Cox/Miller directions
    poly3  = cbind(1, eta, eta^2, eta^3),
    stukel = cbind(1, eta, eta^2 * (eta >= 0), -eta^2 * (eta < 0)))
  Z <- Z[, apply(Z, 2, function(c) stats::sd(c) > 1e-10 | all(c == 1)), drop = FALSE]
  Q <- qr(Z); if (Q$rank < ncol(Z)) Z <- Z[, Q$pivot[seq_len(Q$rank)], drop = FALSE]
  S <- sum(qr.fitted(qr(Z), r)^2)
  stats::pchisq(S, df = ncol(Z), lower.tail = FALSE)
}
hl_ext <- function(y, p, G) { z <- resid_ext(y, p, G)
  stats::pchisq(sum(z$r^2), df = z$K, lower.tail = FALSE) }
slope_lrt <- function(y, p) {                                       # Cox/Miller: H0 a=0,b=1
  lp <- stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  f1 <- tryCatch(stats::glm(y ~ lp, family = stats::binomial()), error = function(e) NULL)
  if (is.null(f1)) return(NA_real_)
  ll0 <- sum(y * log(p) + (1 - y) * log(1 - p))
  stats::pchisq(2 * (as.numeric(stats::logLik(f1)) - ll0), df = 2, lower.tail = FALSE)
}
giviti_ext <- function(y, p) {
  o <- tryCatch(givitiR::givitiCalibrationTest(y, p, devel = "external"),
                error = function(e) NULL, warning = function(w) NULL)
  if (is.null(o)) NA_real_ else as.numeric(o$p.value)
}
one <- function(y, p, G) c(
  EDGE.stk4 = edge_ext(y, p, G, "stukel"), EDGE.poly4 = edge_ext(y, p, G, "poly3"),
  SlopeLRT = slope_lrt(y, p), EDGE.slope2 = edge_ext(y, p, G, "slope"),
  HL.ext = hl_ext(y, p, G), GiViTI.ext = giviti_ext(y, p))
TESTS <- c("EDGE.stk4", "EDGE.poly4", "SlopeLRT", "EDGE.slope2", "HL.ext", "GiViTI.ext")

## ---- the distortions: what the frozen classifier REPORTS vs what is TRUE --------------------
## truth: logit p = 0.6x + 0.5d, balanced. The classifier reports p' = f(p). y ~ Bern(p_true).
distort <- function(p, kind, s) {
  e <- stats::qlogis(p)
  switch(kind,
    null        = p,
    temp        = stats::plogis(e / s),                              # symmetric: s<1 overconfident
    asym_hi     = stats::plogis(ifelse(e >= 0, e / s, e)),           # overconfident on top tail only
    asym_lo     = stats::plogis(ifelse(e <  0, e / s, e)),           # overconfident on bottom only
    probit_comp = stats::pnorm(e * 0.588 / s),                       # probit-type compression
    stop(kind))
}
gen <- function(n) { x <- runif(n, -3, 3); d <- rbinom(n, 1, .5); stats::plogis(0.6 * x + 0.5 * d) }

CELLS <- expand.grid(kind = c("temp", "asym_hi", "asym_lo", "probit_comp"), s = c(0.85, 0.7),
                     n = NS, stringsAsFactors = FALSE)
cat(sprintf("run H: %d alternative cells + %d null cells, B=%d\n", nrow(CELLS), length(NS), B))

ncores <- max(1L, min(20L, parallel::detectCores() - 2L))
cl <- parallel::makeCluster(ncores); on.exit(parallel::stopCluster(cl), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {suppressPackageStartupMessages(library(givitiR)); NULL}))
parallel::clusterExport(cl, c("grp", "resid_ext", "edge_ext", "hl_ext", "slope_lrt", "giviti_ext",
                              "one", "distort", "gen"))

Gof <- function(n) max(6L, as.integer(round(n / 25)))              # the m ~ 25 rule from run F
t0 <- Sys.time(); nul <- list(); alt <- list()
for (n in NS) {
  parallel::clusterSetRNGStream(cl, SEED0 + n)
  nul[[as.character(n)]] <- t(parallel::parSapply(cl, seq_len(B), function(b, n, G) {
    p <- gen(n); y <- rbinom(n, 1, p); one(y, p, G) }, n = n, G = Gof(n)))
  cat(sprintf("  null n=%5d G=%3d (%s)\n", n, Gof(n), format(round(difftime(Sys.time(), t0, units="mins"),1)))); flush.console()
}
for (i in seq_len(nrow(CELLS))) {
  ce <- CELLS[i, ]
  parallel::clusterSetRNGStream(cl, SEED0 + 10000L + i)
  alt[[i]] <- t(parallel::parSapply(cl, seq_len(B), function(b, n, G, kind, s) {
    p <- gen(n); y <- rbinom(n, 1, p); one(y, distort(p, kind, s), G) },
    n = ce$n, G = Gof(ce$n), kind = ce$kind, s = ce$s))
  cat(sprintf("  [%2d/%2d] %-12s s=%.2f n=%5d (%s)\n", i, nrow(CELLS), ce$kind, ce$s, ce$n,
              format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); flush.console()
}

adj <- function(pa, pn, a = .05) { pn <- pn[is.finite(pn)]; cr <- stats::quantile(pn, a, type = 1)
  mean(pa[is.finite(pa)] <= cr) }
out <- list()
for (i in seq_len(nrow(CELLS))) { ce <- CELLS[i, ]; N <- nul[[as.character(ce$n)]]
  for (t in TESTS) out[[length(out) + 1]] <- data.frame(kind = ce$kind, s = ce$s, n = ce$n, G = Gof(ce$n),
    test = t, size = mean(N[, t] <= .05, na.rm = TRUE), power_adj = adj(alt[[i]][, t], N[, t]),
    na = mean(!is.finite(alt[[i]][, t]))) }
R <- do.call(rbind, out)
utils::write.csv(R, "runH_overconfidence_summary.csv", row.names = FALSE)

cat("\n=== SIZE in external mode (nominal .05) ===\n")
for (t in TESTS) cat(sprintf("  %-12s %s\n", t, paste(sprintf("%.3f", sapply(NS, function(n) mean(nul[[as.character(n)]][, t] <= .05, na.rm = TRUE))), collapse = "  ")))
for (k in c("temp", "asym_hi", "asym_lo", "probit_comp")) for (s in c(0.85, 0.7)) {
  cat(sprintf("\n=== %s  s=%.2f  (size-adjusted power; cols n=%s) ===\n", k, s, paste(NS, collapse = ",")))
  for (t in TESTS) cat(sprintf("  %-12s %s\n", t, paste(sprintf("%.3f", sapply(NS, function(n) R$power_adj[R$kind == k & R$s == s & R$n == n & R$test == t])), collapse = "  ")))
}
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
