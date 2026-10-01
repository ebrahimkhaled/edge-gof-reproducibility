## audit_fair_giviti.R -- fairness audit: was GiViTI run at its best defensible setting?
## givitiR::givitiCalibrationTest(o, e, devel = "internal", thres = 0.95, maxDeg = 4) is the default.
## The only user-facing knobs are thres (forward-selection gate, LR > qchisq(thres, 1)) and maxDeg.
## Arms, all on the SAME samples (base harness design x ~ U(-3,3), d ~ Bern(.5)):
##   G95_4 default | G95_3 maxDeg 3 | G80_4 thres .80 | G50_4 thres .50   (all selection-aware nulls)
##   LRdeg3, LRdeg4  fixed-degree polynomial LR in logit(e) vs degree 1 (no selection; chi2_2, chi2_3)
##   Stk.LR, Stk.joint, g.sym.sc, g.max3.sc  (references, harness definitions)
## Cells: cauchit s=1 n=2500 (symmetric tail), probit s=2 n=4000 (high discrimination), cloglog s=1
## n=600 (bow, where GiViTI leads). Each with its own logistic null. 6 workers. Env AFG_B (default 500).
suppressPackageStartupMessages({library(parallel); library(givitiR); library(mvtnorm)})
source("_pstar_giviti_harness.R")
B <- as.integer(Sys.getenv("AFG_B", "500"))
CELLS <- data.frame(link = c("cauchit", "probit", "cloglog"), s = c(1, 2, 1), c0 = 0, n = c(2500L, 4000L, 600L))
FL <- list(logit = plogis, probit = pnorm, cauchit = pcauchy, cloglog = function(e) 1 - exp(-exp(e)))
VAR <- list(G95_4 = c(.95, 4), G95_3 = c(.95, 3), G80_4 = c(.80, 4), G50_4 = c(.50, 4))
TESTS <- c(names(VAR), "m_default", "m_G50", "LRdeg3", "LRdeg4", "Stk.LR", "Stk.joint", "g.sym.sc", "g.max3.sc")
aoc <- function(e, p) ifelse(e < -30, 0, 1 - log1p(exp(e)) / p)
maxz_p <- function(z, R) { T <- max(abs(z)); if (!all(is.finite(c(z, R)))) return(NA_real_)
  1 - as.numeric(pmvnorm(lower = rep(-T, length(z)), upper = rep(T, length(z)), corr = R)) }

one <- function(lk, s, c0, n) {
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); X <- cbind(1, x, d); et <- c0 + s * (0.6 * x + 0.5 * d)
  y <- rbinom(n, 1, FL[[lk]](et))
  out <- setNames(rep(NA_real_, length(TESTS)), TESTS)
  G <- max(10L, as.integer(round(n / 25)))
  if (sum(y) < G || sum(1 - y) < G) return(out)
  fit <- tryCatch(suppressWarnings(glm(y ~ x + d, family = binomial(), control = list(maxit = 50))), error = function(e) NULL)
  if (is.null(fit)) return(out)
  ph <- pmin(pmax(fitted(fit), 1e-10), 1 - 1e-10); eh <- qlogis(ph); w <- ph * (1 - ph)
  for (v in names(VAR)) {
    tt <- tryCatch(suppressWarnings(givitiCalibrationTest(y, ph, devel = "internal", thres = VAR[[v]][1], maxDeg = VAR[[v]][2])),
                   error = function(e) NULL)
    if (!is.null(tt)) { out[v] <- tt$p.value
      if (v == "G95_4") out["m_default"] <- length(tt$estimate) - 1
      if (v == "G50_4") out["m_G50"] <- length(tt$estimate) - 1 }
  }
  dev <- sapply(1:4, function(k) tryCatch(suppressWarnings(glm(y ~ poly(eh, k, raw = TRUE), family = binomial(),
                                                               control = list(maxit = 50)))$deviance, error = function(e) NA_real_))
  out["LRdeg3"] <- pchisq(dev[1] - dev[3], 2, lower.tail = FALSE)
  out["LRdeg4"] <- pchisq(dev[1] - dev[4], 3, lower.tail = FALSE)
  out["Stk.LR"] <- stukel_p(fit)
  Zs <- cbind(za = 0.5 * eh^2 * (eh >= 0), zb = -0.5 * eh^2 * (eh < 0))
  if (all(apply(Zs, 2, sd) > 0)) {
    us <- drop(crossprod(Zs, y - ph)); ZW <- crossprod(Zs, w * X)
    Is <- crossprod(Zs, w * Zs) - ZW %*% solve(crossprod(X, w * X), t(ZW))
    out["Stk.joint"] <- tryCatch(pchisq(drop(t(us) %*% solve(Is, us)), 2, lower.tail = FALSE), error = function(e) NA_real_)
  }
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  o <- as.numeric(tapply(y, g, sum)); e <- as.numeric(tapply(ph, g, sum)); V <- as.numeric(tapply(w, g, sum))
  pb <- as.numeric(tapply(ph, g, mean)); eb <- qlogis(pb); r <- (o - e) / sqrt(V)
  U <- rowsum(w * X, g) / sqrt(V); XWX <- crossprod(X, w * X)
  ZOZ <- function(Z) crossprod(Z) - crossprod(Z, U) %*% solve(XWX, crossprod(U, Z))
  Zt <- cbind(sym = eb * abs(eb), ao = aoc(eb, pb), aoM = aoc(-eb, 1 - pb)) * sqrt(V)
  Ot <- tryCatch(ZOZ(Zt), error = function(e) NULL)
  if (!is.null(Ot)) {
    zt <- drop(crossprod(Zt, r)) / sqrt(pmax(diag(Ot), 1e-300))
    out["g.sym.sc"] <- 2 * pnorm(-abs(zt["sym"]))
    Rt <- tryCatch(cov2cor(Ot), error = function(e) NULL)
    if (!is.null(Rt)) out["g.max3.sc"] <- maxz_p(zt, Rt)
  }
  out
}

## selection-aware critical values of the GiViTI statistic at alpha = .05, internal mode
cdf <- get("givitiStatCdf", asNamespace("givitiR"))
cat("GiViTI internal-mode 5% critical value of the statistic, by selected degree m\n")
for (th in c(.95, .80, .50)) cat(sprintf("  thres %.2f (gate %.2f): m=2 %.2f  m=3 %.2f  m=4 %.2f   [fixed-degree LR: 2 df 5.99, 3 df 7.81]\n", th, qchisq(th, 1),
  qchisq(.95, 1), uniroot(function(t) 1 - cdf(t, 3, "internal", th) - .05, c(qchisq(th, 1) + 1e-6, 60))$root,
  uniroot(function(t) 1 - cdf(t, 4, "internal", th) - .05, c(2 * qchisq(th, 1) + 1e-6, 60))$root))

## the lead's mechanism file, recomputed
for (lk in c("logit", "cauchit")) { M <- read.csv(sprintf("giviti_mechanism_%s.csv", lk))
  cat(sprintf("mechanism %s: B=%d  P(m=2) %.3f  P(lr3>3.84) %.3f  rej %.3f | rej|m=2 %.3f  rej|m=3 %.3f  (n m=3: %d) | fixed deg3 LR %.3f\n", lk, nrow(M),
    mean(M$m == 2), mean(M$lr3 > qchisq(.95, 1)), mean(M$p <= .05), mean(M$p[M$m == 2] <= .05), mean(M$p[M$m == 3] <= .05), sum(M$m == 3),
    mean(M$lr2 + M$lr3 > qchisq(.95, 2)))) }

t0 <- Sys.time()
cl <- makeCluster(6)
invisible(clusterEvalQ(cl, {suppressPackageStartupMessages({library(givitiR); library(mvtnorm)}); NULL}))
clusterExport(cl, c("one", "FL", "VAR", "TESTS", "aoc", "maxz_p", "stukel_p"))
P <- list(); k <- 0L
for (i in seq_len(nrow(CELLS))) for (lk in c("logit", CELLS$link[i])) {
  k <- k + 1L; clusterSetRNGStream(cl, 20260930L + 131L * k)
  M <- t(parSapply(cl, seq_len(B), function(b, lk, s, c0, n) one(lk, s, c0, n), lk = lk, s = CELLS$s[i], c0 = CELLS$c0[i], n = CELLS$n[i]))
  P[[k]] <- data.frame(cell = i, link = lk, n = CELLS$n[i], M, check.names = FALSE)
  cat(sprintf("  cell %d %-8s n=%d  %s\n", i, lk, CELLS$n[i], format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); flush.console()
}
stopCluster(cl)
pv <- do.call(rbind, P); write.csv(pv, "audit_fair_giviti_pvalues.csv", row.names = FALSE)

crit <- function(p) { p <- p[is.finite(p)]; as.numeric(quantile(p, .05, type = 1)) }
PT <- c(names(VAR), "LRdeg3", "LRdeg4", "Stk.LR", "Stk.joint", "g.sym.sc", "g.max3.sc")
for (i in seq_len(nrow(CELLS))) {
  nul <- pv[pv$cell == i & pv$link == "logit", ]; alt <- pv[pv$cell == i & pv$link != "logit", ]
  cat(sprintf("\n-- %s s=%g n=%d (B=%d) --  default m under alt: %s | G50 m: %s\n", CELLS$link[i], CELLS$s[i], CELLS$n[i], B,
    paste(names(table(alt$m_default)), table(alt$m_default), sep = ":", collapse = " "),
    paste(names(table(alt$m_G50)), table(alt$m_G50), sep = ":", collapse = " ")))
  rej <- lapply(setNames(PT, PT), function(t) { a <- alt[[t]]; ifelse(is.finite(a), a <= crit(nul[[t]]), FALSE) })
  for (t in PT) {
    dG <- mean(rej[[t]]) - mean(rej$G95_4); n10 <- sum(rej[[t]] & !rej$G95_4); n01 <- sum(!rej[[t]] & rej$G95_4)
    dS <- mean(rej[[t]]) - mean(rej$g.sym.sc)
    cat(sprintf("  %-9s size %.3f  na %.3f  power(adj) %.3f  raw %.3f | vs default GiViTI %+.3f (McNemar p %s) | vs g.sym.sc %+.3f\n", t,
      mean(nul[[t]] <= .05, na.rm = TRUE), mean(!is.finite(alt[[t]])), mean(rej[[t]]), mean(alt[[t]] <= .05, na.rm = TRUE), dG,
      ifelse(n10 + n01 > 0, format.pval(binom.test(n10, n10 + n01)$p.value, digits = 2), "NA"), dS))
  }
}
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
