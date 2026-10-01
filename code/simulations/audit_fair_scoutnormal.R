## audit_fair_scoutnormal.R -- fairness audit. The population envelope of scout_envelope.R (same D, rho,
## n80 = lambda80(k)/(rho D), same Gaussian-limit Monte Carlo for GiViTI and the max probe), recomputed
## on covariate distributions OTHER than x ~ U(-3,3):
##   unif   x ~ U(-3,3)                      (anchor: must reproduce scout_envelope.csv)
##   normal x ~ N(0, 1.5^2)                  (tails sparsely populated)
##   skew   x = 1.5 * standardised chi2_4    (right-skewed, clinical-style risk distribution)
## Writes audit_fair_scoutnormal.csv only. 6 worker processes.
suppressPackageStartupMessages(library(parallel))
lam80 <- sapply(1:10, function(k) uniroot(function(l)
  pchisq(qchisq(.95, k), k, ncp = l, lower.tail = FALSE) - .8, c(0.01, 200))$root)
CELLS <- expand.grid(xd = c("unif", "normal", "skew"), s = c(1, 2), c0 = c(0, -2),
                     link = c("probit", "cauchit", "t4", "cloglog", "loglog"), stringsAsFactors = FALSE)
CELLS <- CELLS[!(CELLS$s == 2 & CELLS$c0 == -2), ]

one_cell <- function(i) {
  xd <- CELLS$xd[i]; s <- CELLS$s[i]; c0 <- CELLS$c0[i]; lk <- CELLS$link[i]
  N <- 400000L; set.seed(20260920)
  x <- switch(xd, unif = runif(N, -3, 3), normal = rnorm(N, 0, 1.5), skew = 1.5 * as.numeric(scale(rchisq(N, 4))))
  dd <- rbinom(N, 1, 0.5); X <- cbind(1, x, dd); lin <- 0.6 * x + 0.5 * dd
  E0 <- matrix(rnorm(3 * 100000), ncol = 3)
  FL <- list(probit = pnorm, t4 = function(e) pt(e, 4), cauchit = pcauchy,
             cloglog = function(e) 1 - exp(-exp(e)), loglog = function(e) exp(-exp(-e)))
  eta <- c0 + s * lin; p <- FL[[lk]](eta)
  fit <- suppressWarnings(glm.fit(X, p, family = quasibinomial()))
  pis <- fit$fitted.values; eh <- drop(X %*% fit$coefficients); w <- pis * (1 - pis); dv <- p - pis
  D <- mean(dv^2 / w)
  XtW <- t(X * w); A <- solve(XtW %*% X)
  info <- function(Z) { ZW <- t(Z * w); (ZW %*% Z - ZW %*% X %*% A %*% (XtW %*% Z)) / N }
  score <- function(Z) { Z <- as.matrix(Z); u <- drop(crossprod(Z, dv)) / N; I <- info(Z); list(u = u, I = I) }
  keep_cols <- function(Z) { Z <- as.matrix(Z); ok <- apply(Z, 2, function(z) all(is.finite(z)) && sd(z) > 1e-10)
    Z[, ok, drop = FALSE] }
  closed <- function(Z) { Z <- keep_cols(Z); k <- ncol(Z)
    if (!k) return(c(rho = NA, n80 = NA, k = 0))
    sc <- score(Z); ncp <- tryCatch(drop(t(sc$u) %*% solve(sc$I, sc$u)), error = function(e) NA_real_)
    c(rho = ncp / D, n80 = lam80[k] / ncp, k = k) }
  aoc <- function(e, p) ifelse(e < -30, 0, 1 - log1p(exp(e)) / p)
  B <- list(matched = dv / w, sym = eh * abs(eh), ao = aoc(eh, pis), aoM = aoc(-eh, 1 - pis),
            stukel = cbind(eh^2 * (eh >= 0), -eh^2 * (eh < 0)), lg23 = cbind(eh^2, eh^3), lg234 = cbind(eh^2, eh^3, eh^4))
  res <- lapply(B, closed)
  g <- pmin(ceiling(rank(pis, ties.method = "first") / (N / 10)), 10)
  hl_ncp <- sum(tapply(dv, g, sum)^2 / tapply(w, g, sum)) / N
  res$HL10 <- c(rho = hl_ncp / D, n80 = lam80[8] / hl_ncp, k = 8)
  mc_n80 <- function(u, I, stat) {
    L <- tryCatch(t(chol(I)), error = function(e) NULL); if (is.null(L)) return(NA_real_)
    null <- stat(E0 %*% t(L)); crit <- quantile(null, .95)
    pw <- function(n) mean(stat(sweep(E0 %*% t(L), 2, sqrt(n) * u, "+")) > crit)
    if (pw(1e8) < .8) return(Inf)
    exp(uniroot(function(t) pw(exp(t)) - .8, c(log(5), log(1e8)))$root)
  }
  sg <- score(B$lg234)
  giviti_stat <- function(S) { Z <- t(solve(t(chol(sg$I)), t(S)))
    z2 <- Z^2; step3 <- z2[, 2] > qchisq(.95, 1); step4 <- step3 & z2[, 3] > qchisq(.95, 1)
    z2[, 1] + z2[, 2] * step3 + z2[, 3] * step4 }
  res$GiViTI <- c(rho = NA, n80 = tryCatch(mc_n80(sg$u, sg$I, giviti_stat), error = function(e) NA_real_), k = NA)
  sp <- score(cbind(B$ao, B$aoM, B$sym)); sdv <- sqrt(diag(sp$I))
  res$maxprobe <- c(rho = NA, n80 = tryCatch(mc_n80(sp$u, sp$I, function(S) apply(abs(sweep(S, 2, sdv, "/")), 1, max)),
                                            error = function(e) NA_real_), k = NA)
  o <- order(eta); pp <- p[o]; q <- 1 - pp; auc <- sum(pp * (cumsum(q) - q)) / (sum(pp) * sum(q))
  data.frame(xd = xd, s = s, c0 = c0, link = lk, event = mean(p), auc = auc, D = D,
             test = names(res), rho = sapply(res, `[`, "rho"), n80 = sapply(res, `[`, "n80"), row.names = NULL)
}
t0 <- Sys.time()
cl <- makeCluster(6); clusterExport(cl, c("CELLS", "lam80", "one_cell"))
R <- do.call(rbind, parLapply(cl, seq_len(nrow(CELLS)), one_cell)); stopCluster(cl)
write.csv(R, "audit_fair_scoutnormal.csv", row.names = FALSE)
TT <- c("matched", "sym", "ao", "aoM", "maxprobe", "stukel", "GiViTI", "HL10")
options(width = 200)
for (s in c(1, 2)) for (c0 in c(0, -2)) for (xd in c("unif", "normal", "skew")) {
  z <- R[R$s == s & R$c0 == c0 & R$xd == xd, ]; if (!nrow(z)) next
  cat(sprintf("\n-- x %s  s=%g c0=%+g --\n  %-8s %6s %5s |%s\n", xd, s, c0, "link", "event", "AUC", paste(sprintf("%9s", TT), collapse = "")))
  for (lk in unique(z$link)) { zz <- z[z$link == lk, ]
    cat(sprintf("  %-8s %6.3f %5.3f |%s\n", lk, zz$event[1], zz$auc[1], paste(sprintf("%9.0f", sapply(TT, function(t) zz$n80[zz$test == t])), collapse = ""))) }
}
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
