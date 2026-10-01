## audit_fair_design.R -- FAIRNESS AUDIT copy of run_L2_probe_G.R. Changes, and only these:
##   (1) the design line: x ~ N(0, 1.5^2) instead of U(-3,3), in one() AND in the oracle table;
##   (2) makeCluster(6) (shared machine);
##   (3) output files audit_fair_design_<TAG>_*.csv (so analyse_runL.R never picks them up).
## Everything else (tests, grouping, nulls, seeds) is the run L2 harness verbatim.
## Usage:  RUNL_M=25 Rscript audit_fair_design.R LINKS DESIGN NS B TAG
suppressPackageStartupMessages({library(parallel); library(givitiR); library(mvtnorm)})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
LINKS <- strsplit(args[1], ",")[[1]]
dsg <- as.numeric(strsplit(args[2], ",")[[1]]); S <- dsg[1]; C0 <- dsg[2]
NS <- as.integer(strsplit(args[3], ",")[[1]]); B <- as.integer(args[4]); TAG <- args[5]
M_PER_G <- as.integer(Sys.getenv("RUNL_M", "25"))
source("_pstar_giviti_harness.R")                     # hl_stat, stukel_p (LR refit), size_adjusted, na_rate

lnconv <- function(e, s) { K <- 40L; w <- qnorm((seq_len(K) - 0.5) / K)
  out <- numeric(length(e)); for (k in seq_len(K)) out <- out + plogis(e + s * w[k]); out / K }
FL <- list(logit = plogis, probit = pnorm, t4 = function(e) pt(e, 4), cauchit = pcauchy,
           lnorm2 = function(e) lnconv(e, 2), cloglog = function(e) 1 - exp(-exp(e)),
           loglog = function(e) exp(-exp(-e)))
stopifnot(all(LINKS %in% names(FL)))
aoc <- function(e, p) ifelse(e < -30, 0, 1 - log1p(exp(e)) / p)

oracle_tab <- function(lk) {
  set.seed(1); Np <- 400000L
  x <- rnorm(Np, 0, 1.5); d <- rbinom(Np, 1, 0.5); X <- cbind(1, x, d); et <- C0 + S * (0.6 * x + 0.5 * d)   # AUDIT: design line
  p <- FL[[lk]](et); fit <- suppressWarnings(glm.fit(X, p, family = quasibinomial()))
  eh <- drop(X %*% fit$coefficients); ph <- fit$fitted.values; phi <- (p - ph) / (ph * (1 - ph))
  br <- unique(quantile(eh, seq(0, 1, length.out = 401))); g <- cut(eh, br, include.lowest = TRUE)
  list(x = as.numeric(tapply(eh, g, mean)), y = as.numeric(tapply(phi, g, mean)))
}
ORC <- lapply(setNames(LINKS, LINKS), oracle_tab)
orc_eval <- function(lk, e) approx(ORC[[lk]]$x, ORC[[lk]]$y, e, rule = 2)$y
orc_mat <- function(e) { M <- do.call(cbind, lapply(LINKS, function(l) orc_eval(l, e)))
  colnames(M) <- paste0("orc.", LINKS); M }

TESTS <- c("EDGE.poly3", "EDGE.stk", "EDGE.stk.sc", "g.sym", "g.ao", "g.aoM", "g.max3", "g.sym.sc", "g.ao.sc", "g.aoM.sc", "g.max3.sc", paste0("g.orc.", LINKS),
           "u.sym", "u.ao", "u.aoM", "u.cub", "u.max3", paste0("u.orc.", LINKS),
           "Stk.joint", "Stk.marg", "Stk.LR", "GiViTI", "HL10", "HLG")

postfit_info <- function(Z, w, X) { ZW <- crossprod(Z, w * X)
  crossprod(Z, w * Z) - ZW %*% solve(crossprod(X, w * X), t(ZW)) }
maxz_p <- function(z, R) { T <- max(abs(z)); if (!all(is.finite(c(z, R)))) return(NA_real_)
  1 - as.numeric(pmvnorm(lower = rep(-T, length(z)), upper = rep(T, length(z)), corr = R)) }

one <- function(lk, n) {
  G <- max(10L, as.integer(round(n / M_PER_G)))
  x <- rnorm(n, 0, 1.5); d <- rbinom(n, 1, 0.5); X <- cbind(1, x, d); et <- C0 + S * (0.6 * x + 0.5 * d)   # AUDIT: design line
  y <- rbinom(n, 1, FL[[lk]](et))
  out <- setNames(rep(NA_real_, length(TESTS)), TESTS)
  if (sum(y) < G || sum(1 - y) < G) return(out)
  fit <- tryCatch(suppressWarnings(glm(y ~ x + d, family = binomial(), control = list(maxit = 50))), error = function(e) NULL)
  if (is.null(fit)) return(out)
  ph <- pmin(pmax(fitted(fit), 1e-10), 1 - 1e-10); eh <- qlogis(ph); w <- ph * (1 - ph)

  Zu <- cbind(sym = eh * abs(eh), ao = aoc(eh, ph), aoM = aoc(-eh, 1 - ph), cub = eh^3,
              orc_mat(eh))
  uu <- drop(crossprod(Zu, y - ph)); Iu <- tryCatch(postfit_info(Zu, w, X), error = function(e) NULL)
  if (!is.null(Iu)) {
    zu <- uu / sqrt(pmax(diag(Iu), 1e-300))
    for (nm in c("sym", "ao", "aoM", "cub", paste0("orc.", LINKS))) out[paste0("u.", nm)] <- 2 * pnorm(-abs(zu[nm]))
    k3 <- c("sym", "ao", "aoM"); R3 <- tryCatch(cov2cor(Iu[k3, k3]), error = function(e) NULL)
    if (!is.null(R3)) out["u.max3"] <- maxz_p(zu[k3], R3)
  }
  Zs <- cbind(za = 0.5 * eh^2 * (eh >= 0), zb = -0.5 * eh^2 * (eh < 0))
  if (all(apply(Zs, 2, sd) > 0)) {
    us <- drop(crossprod(Zs, y - ph)); Is <- postfit_info(Zs, w, X)
    js <- tryCatch(drop(t(us) %*% solve(Is, us)), error = function(e) NA_real_)
    out["Stk.joint"] <- pchisq(js, 2, lower.tail = FALSE)
    out["Stk.marg"]  <- pchisq(sum(us^2 / diag(Is)), 2, lower.tail = FALSE)
  }
  out["Stk.LR"] <- stukel_p(fit)
  out["GiViTI"] <- tryCatch(suppressWarnings(givitiCalibrationTest(y, ph, devel = "internal")$p.value), error = function(e) NA_real_)
  out["HL10"] <- unname(hl_stat(y, ph, 10)["p"]); out["HLG"] <- unname(hl_stat(y, ph, G)["p"])

  g <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  o <- as.numeric(tapply(y, g, sum)); e <- as.numeric(tapply(ph, g, sum)); V <- as.numeric(tapply(w, g, sum))
  pb <- as.numeric(tapply(ph, g, mean)); eb <- qlogis(pb); r <- (o - e) / sqrt(V)
  U <- rowsum(w * X, g) / sqrt(V); XWX <- crossprod(X, w * X)
  ZOZ <- function(Z) crossprod(Z) - crossprod(Z, U) %*% solve(XWX, crossprod(U, Z))
  edge_eig <- function(Z) { Z <- Z[, apply(Z, 2, sd) > 1e-10, drop = FALSE]; if (!ncol(Z)) return(NA_real_)
    lam <- tryCatch(Re(eigen(solve(crossprod(Z), ZOZ(Z)), only.values = TRUE)$values), error = function(e) NULL)
    if (is.null(lam)) return(NA_real_); lam <- lam[lam > 1e-9]; if (!length(lam)) return(NA_real_)
    Sx <- sum(qr.fitted(qr(Z), r)^2); cc <- sum(lam^2) / sum(lam); nu <- sum(lam)^2 / sum(lam^2)
    pchisq(Sx / cc, nu, lower.tail = FALSE) }
  out["EDGE.poly3"] <- edge_eig(as.matrix(poly(pb, 3)))
  out["EDGE.stk"]   <- edge_eig(cbind(eb, eb^2 * (eb >= 0), -eb^2 * (eb < 0)))
  Zk <- cbind(za = eb^2 * (eb >= 0), zb = -eb^2 * (eb < 0)) * sqrt(V)
  if (all(apply(Zk, 2, sd) > 1e-10)) { sk <- drop(crossprod(Zk, r)); Ik <- tryCatch(ZOZ(Zk), error = function(e) NULL)
    if (!is.null(Ik)) out["EDGE.stk.sc"] <- tryCatch(pchisq(drop(t(sk) %*% solve(Ik, sk)), 2, lower.tail = FALSE), error = function(e) NA_real_) }
  Zt <- cbind(sym = eb * abs(eb), ao = aoc(eb, pb), aoM = aoc(-eb, 1 - pb)) * sqrt(V)
  Ot <- tryCatch(ZOZ(Zt), error = function(e) NULL)
  if (!is.null(Ot)) {
    zt <- drop(crossprod(Zt, r)) / sqrt(pmax(diag(Ot), 1e-300))
    for (nm in c("sym", "ao", "aoM")) out[paste0("g.", nm, ".sc")] <- 2 * pnorm(-abs(zt[nm]))
    Rt <- tryCatch(cov2cor(Ot), error = function(e) NULL)
    if (!is.null(Rt)) out["g.max3.sc"] <- maxz_p(zt, Rt)
  }
  Zg <- cbind(sym = eb * abs(eb), ao = aoc(eb, pb), aoM = aoc(-eb, 1 - pb),
              orc_mat(eb))
  Og <- tryCatch(ZOZ(Zg), error = function(e) NULL)
  if (!is.null(Og)) {
    zg <- drop(crossprod(Zg, r)) / sqrt(pmax(diag(Og), 1e-300))
    for (nm in c("sym", "ao", "aoM", paste0("orc.", LINKS))) out[paste0("g.", nm)] <- 2 * pnorm(-abs(zg[nm]))
    k3 <- c("sym", "ao", "aoM"); R3 <- tryCatch(cov2cor(Og[k3, k3]), error = function(e) NULL)
    if (!is.null(R3)) out["g.max3"] <- maxz_p(zg[k3], R3)
  }
  out
}

cl <- makeCluster(6)                                                                                     # AUDIT: 6 cores
invisible(clusterEvalQ(cl, {suppressPackageStartupMessages({library(givitiR); library(mvtnorm)}); NULL}))
clusterExport(cl, c("M_PER_G", "one", "FL", "lnconv", "aoc", "ORC", "orc_eval", "orc_mat", "TESTS", "LINKS", "S", "C0",
                    "postfit_info", "maxz_p", "hl_stat", "stukel_p"))
SEED0 <- 20260922L + as.integer(1000 * S) + as.integer(10 * abs(C0))
t0 <- Sys.time(); P <- list(); cellno <- 0L
for (n in NS) for (lk in c("logit", LINKS)) {
  cellno <- cellno + 1L; clusterSetRNGStream(cl, SEED0 + 7919L * cellno + n)
  M <- t(parSapply(cl, seq_len(B), function(b, lk, n) one(lk, n), lk = lk, n = n))
  P[[length(P) + 1]] <- data.frame(link = lk, s = S, c0 = C0, n = n, rep = seq_len(B), M, check.names = FALSE)
  cat(sprintf("  %-8s n=%6d  %s\n", lk, n, format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); flush.console()
}
stopCluster(cl)
pv <- do.call(rbind, P)
write.csv(pv, sprintf("audit_fair_design_%s_pvalues.csv", TAG), row.names = FALSE)                      # AUDIT: file name

out <- list()
for (n in NS) {
  nul <- pv[pv$link == "logit" & pv$n == n, ]
  for (lk in LINKS) {
    alt <- pv[pv$link == lk & pv$n == n, ]
    for (t in TESTS) {
      if (grepl("orc\\.", t) && !grepl(paste0("orc\\.", lk, "$"), t)) next
      out[[length(out) + 1]] <- data.frame(link = lk, s = S, c0 = C0, n = n, test = sub(paste0("orc\\.", lk), "orc", t),
        size = mean(nul[[t]] <= .05, na.rm = TRUE), power_adj = size_adjusted(alt[[t]], nul[[t]]),
        power_raw = mean(alt[[t]] <= .05, na.rm = TRUE), na_alt = na_rate(alt[[t]]), na_null = na_rate(nul[[t]]))
    }
  }
}
Sm <- do.call(rbind, out)
write.csv(Sm, sprintf("audit_fair_design_%s_summary.csv", TAG), row.names = FALSE)                      # AUDIT: file name
SHOW <- c("g.orc", "g.sym", "g.sym.sc", "g.ao.sc", "g.max3", "g.max3.sc", "EDGE.poly3", "EDGE.stk", "EDGE.stk.sc", "u.orc", "u.sym", "u.max3",
          "Stk.joint", "Stk.marg", "Stk.LR", "GiViTI", "HL10")
for (lk in LINKS) for (n in NS) {
  z <- Sm[Sm$link == lk & Sm$n == n, ]
  cat(sprintf("\n-- %s  s=%g c0=%g  n=%d  (x ~ N(0,1.5^2)) --\n", lk, S, C0, n))
  for (t in SHOW) if (t %in% z$test) with(z[z$test == t, ],
    cat(sprintf("  %-11s power %.3f  raw %.3f  size %.3f  na %.3f\n", t, power_adj, power_raw, size, na_alt)))
}
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
