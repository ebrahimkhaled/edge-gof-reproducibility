## scout_envelope.R -- the population landscape of link misfit (n -> infinity).
##
## For a true link F and a design, pi* is the pseudo-true logistic fit and
##     D = E[ (p_F - pi*)^2 / (pi*(1 - pi*)) ]
## is the largest local non-centrality per patient that ANY test can reach against F: the one-df
## test in the matched direction (p_F - pi*)/w attains it. A test built on the directions Z
## captures the share rho = u'I^{-1}u / D, with u = E[Z(p_F - pi*)] and I the post-fit
## information of Z, and pays k = ncol(Z) degrees of freedom, so its 80%-power sample size is
##     n80 = lambda80(k) / (rho * D).
## GiViTI's forward selection and the max over link probes have no closed-form lambda80, so their
## power is computed by Monte Carlo in the same Gaussian limit. Everything here is ungrouped; the
## growing-partition theorem says the grouped statistics reach these limits as G grows.
suppressPackageStartupMessages(library(stats))
GIVITI_START <- as.integer(Sys.getenv("GIVITI_START", "2"))   # internal-mode starting degree; checked against givitiR source

N <- 400000L
set.seed(20260920)
x <- runif(N, -3, 3); dd <- rbinom(N, 1, 0.5); X <- cbind(1, x, dd); lin <- 0.6 * x + 0.5 * dd
lnconv <- function(e, s) { K <- 40L; w <- qnorm((seq_len(K) - 0.5) / K)
  out <- numeric(length(e)); for (k in seq_len(K)) out <- out + plogis(e + s * w[k]); out / K }
LINKS <- list(
  probit  = function(e) pnorm(e),
  t4      = function(e) pt(e, 4),          # robit: heavier tails than the logistic
  cauchit = function(e) pcauchy(e),
  lnorm2  = function(e) lnconv(e, 2),      # an omitted N(0, 2^2) random intercept, marginalised
  cloglog = function(e) 1 - exp(-exp(e)),
  loglog  = function(e) exp(-exp(-e)))
DESIGNS <- expand.grid(s = c(1, 2, 3), c0 = c(0, -2))

lam80 <- sapply(1:10, function(k) uniroot(function(l)
  pchisq(qchisq(.95, k), k, ncp = l, lower.tail = FALSE) - .8, c(0.01, 200))$root)
E0 <- matrix(rnorm(3 * 100000), ncol = 3)               # common random numbers for every MC power
auc <- function(p, eta) { o <- order(eta); p <- p[o]; q <- 1 - p
  sum(p * (cumsum(q) - q)) / (sum(p) * sum(q)) }

one_cell <- function(s, c0, lk) {
  eta <- c0 + s * lin; p <- LINKS[[lk]](eta)
  fit <- suppressWarnings(glm.fit(X, p, family = quasibinomial()))
  pis <- fit$fitted.values; eh <- drop(X %*% fit$coefficients); w <- pis * (1 - pis); dv <- p - pis
  D <- mean(dv^2 / w)
  XtW <- t(X * w); A <- solve(XtW %*% X)
  info <- function(Z) { ZW <- t(Z * w); (ZW %*% Z - ZW %*% X %*% A %*% (XtW %*% Z)) / N }
  score <- function(Z) { Z <- as.matrix(Z); u <- drop(crossprod(Z, dv)) / N; I <- info(Z); list(u = u, I = I) }
  ## a column that is identically zero on this design (a Stukel half with no fitted risks on its
  ## side) or not finite is dropped; a basis left with no column returns NA rather than a crash
  keep_cols <- function(Z) { Z <- as.matrix(Z); ok <- apply(Z, 2, function(z) all(is.finite(z)) && sd(z) > 1e-10)
    Z[, ok, drop = FALSE] }
  closed <- function(Z) { Z <- keep_cols(Z); k <- ncol(Z)
    if (!k) return(c(rho = NA, n80 = NA, k = 0))
    sc <- score(Z); ncp <- tryCatch(drop(t(sc$u) %*% solve(sc$I, sc$u)), error = function(e) NA_real_)
    c(rho = ncp / D, n80 = lam80[k] / ncp, k = k) }
  aoc <- function(e, p) ifelse(e < -30, 0, 1 - log1p(exp(e)) / p)          # stable in the far lower tail
  B <- list(
    matched = dv / w, sym = eh * abs(eh), cub = eh^3,
    ao = aoc(eh, pis), aoM = aoc(-eh, 1 - pis),
    stukel = cbind(eh^2 * (eh >= 0), -eh^2 * (eh < 0)),
    lg23 = cbind(eh^2, eh^3), pi23 = cbind(pis^2, pis^3), lg234 = cbind(eh^2, eh^3, eh^4))
  res <- lapply(B, closed)
  ## HL at G = 10, as the omnibus reference
  g <- pmin(ceiling(rank(pis, ties.method = "first") / (N / 10)), 10)
  hl_ncp <- sum(tapply(dv, g, sum)^2 / tapply(w, g, sum)) / N
  res$HL10 <- c(rho = hl_ncp / D, n80 = lam80[8] / hl_ncp, k = 8)
  ## Monte Carlo powers in the Gaussian limit
  mc_n80 <- function(u, I, stat) {
    L <- tryCatch(t(chol(I)), error = function(e) NULL); if (is.null(L)) return(NA_real_)
    null <- stat(E0 %*% t(L)); crit <- quantile(null, .95)
    pw <- function(n) mean(stat(sweep(E0 %*% t(L), 2, sqrt(n) * u, "+")) > crit)
    if (pw(1e8) < .8) return(Inf)
    exp(uniroot(function(t) pw(exp(t)) - .8, c(log(5), log(1e8)))$root)
  }
  sg <- score(B$lg234)
  giviti_stat <- function(S) {                             # forward selection on eta^2, eta^3, eta^4
    Z <- t(solve(t(chol(sg$I)), t(S)))                      # whitened, nested order
    z2 <- Z^2; step3 <- z2[, 2] > qchisq(.95, 1); step4 <- step3 & z2[, 3] > qchisq(.95, 1)
    if (GIVITI_START == 2) z2[, 1] + z2[, 2] * step3 + z2[, 3] * step4
    else stop("only start degree 2 implemented") }
  res$GiViTI <- c(rho = NA, n80 = mc_n80(sg$u, sg$I, giviti_stat), k = NA)
  sp <- score(cbind(B$ao, B$aoM, B$sym)); sdv <- sqrt(diag(sp$I))
  res$maxprobe <- c(rho = NA, n80 = mc_n80(sp$u, sp$I, function(S) apply(abs(sweep(S, 2, sdv, "/")), 1, max)), k = NA)
  data.frame(s = s, c0 = c0, link = lk, event = mean(p), auc = auc(p, eta), D = D,
             test = names(res), rho = sapply(res, `[`, "rho"), n80 = sapply(res, `[`, "n80"), k = sapply(res, `[`, "k"),
             row.names = NULL)
}

t0 <- Sys.time(); out <- list()
for (i in seq_len(nrow(DESIGNS))) for (lk in names(LINKS)) {
  out[[length(out) + 1]] <- one_cell(DESIGNS$s[i], DESIGNS$c0[i], lk)
  cat(sprintf("  s=%.0f c0=%+.0f %-8s done (%s)\n", DESIGNS$s[i], DESIGNS$c0[i], lk,
              format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); flush.console()
}
R <- do.call(rbind, out); write.csv(R, "scout_envelope.csv", row.names = FALSE)

TESTS <- c("matched", "sym", "ao", "aoM", "maxprobe", "stukel", "GiViTI", "lg23", "pi23", "HL10")
cat("\n=== n80 (patients for 80% power, alpha 0.05, population limit) ===\n")
for (c0 in unique(R$c0)) for (s in unique(R$s)) {
  z <- R[R$s == s & R$c0 == c0, ]
  cat(sprintf("\n-- design s=%.0f c0=%+.0f --\n", s, c0))
  cat(sprintf("  %-8s %6s %5s %9s |%s\n", "link", "event", "AUC", "D", paste(sprintf("%9s", TESTS), collapse = "")))
  for (lk in names(LINKS)) { zz <- z[z$link == lk, ]
    cat(sprintf("  %-8s %6.3f %5.3f %9.2e |%s\n", lk, zz$event[1], zz$auc[1], zz$D[1],
      paste(sprintf("%9.0f", sapply(TESTS, function(t) zz$n80[zz$test == t])), collapse = ""))) }
}
cat("\n=== rho (share of the detectable misfit captured), design s=1 c0=0 ===\n")
z <- R[R$s == 1 & R$c0 == 0 & !is.na(R$rho), ]
for (lk in names(LINKS)) { zz <- z[z$link == lk, ]
  cat(sprintf("  %-8s %s\n", lk, paste(sprintf("%s=%.3f", zz$test, zz$rho), collapse = "  "))) }

## validation against the finite-sample simulations already run (base design, G = n/25)
cat("\n=== validation: predicted limit power vs simulated (base design) ===\n")
pwk <- function(ncp, k) pchisq(qchisq(.95, k), k, ncp = ncp, lower.tail = FALSE)
b <- R[R$s == 1 & R$c0 == 0, ]
ncp1 <- function(lk, t) lam80[b$k[b$link == lk & b$test == t]] / b$n80[b$link == lk & b$test == t]
cat(sprintf("  cloglog n=1000: Stukel pred %.3f (sim 0.860) | AO pred %.3f (sim 0.885) | lg23 pred %.3f (sim sc.poly 0.780)\n",
  pwk(1000 * ncp1("cloglog", "stukel"), 2), pwk(1000 * ncp1("cloglog", "ao"), 1), pwk(1000 * ncp1("cloglog", "lg23"), 2)))
cat(sprintf("  probit  n=5000: Stukel pred %.3f (sim 0.275-0.291) | lg23 pred %.3f (sim sc.poly 0.270) | AO pred %.3f (sim 0.206)\n",
  pwk(5000 * ncp1("probit", "stukel"), 2), pwk(5000 * ncp1("probit", "lg23"), 2), pwk(5000 * ncp1("probit", "ao"), 1)))
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
