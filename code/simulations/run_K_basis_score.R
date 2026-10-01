## run_K_basis_score.R -- "can the basis vector be enhanced to catch smooth link misfit better?"
##
## The bow (cloglog) is where GiViTI and Stukel lead EDGE (0.885 / 0.870 vs 0.827 at n=1000, G=40).
## Both rivals are score / likelihood-ratio tests of "add smooth functions of eta to the linear
## predictor". EDGE tests the SAME directions but with a different weighting: it projects the
## standardized residual r_g = (o_g-e_g)/sqrt(V_g) with unit weights and an eigenvalue null,
## whereas the efficient score of that alternative is Z'(o-e) -- the same contrast weighted by
## sqrt(V_g) -- referred to the Schur-complement information. In EDGE's own notation that is
##     S* = (Zt'r)' (Zt' Omega Zt)^{-1} (Zt'r),   Zt = D^{1/2} Z,   ~ chi^2_k exactly,
## i.e. (a) scale the basis by sqrt(V_g), (b) use the Omega-metric instead of (Z'Z)^{-1} with
## eigen-weights. No refit, no eigenvalues, cannot separate. Directions in the fitted column space
## (constant, eta) are Omega-null and must be left out of a score basis, so k = 2 here.
## Arms (every EDGE variant computed from the null fit only):
##   EDGE.poly3            the paper's default (reference)
##   EDGE.poly3.logit      same construction, polynomial in eta-bar instead of pi-bar
##   EDGE.sc.poly          score form, Z = [eta^2, eta^3]                 chi^2_2
##   EDGE.sc.stk           score form, Z = Stukel's two sign-split columns chi^2_2 (Stukel, no refit)
##   EDGE.sc.AO            score form, Z = Aranda-Ordaz asymmetric score direction  chi^2_1
##   EDGE.sc.AO2           score form, Z = [AO, eta^2]                    chi^2_2
##   EDGE.sc.poly4         score form, Z = [eta^2, eta^3, eta^4]          chi^2_3
##   HL, Stukel (refit), GiViTI      references, from the harness
## Design: the paper's base design at n=1000, G = n/25 = 40 (the rule), scenarios null / cloglog /
## probit / quad / binint / osc4 (+ cloglog at n=500 G=20 and probit at n=5000 G=200). B = 2000.
## Size-adjusted against the matched null in every cell.
suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR); library(parallel)})
source("_pstar_giviti_harness.R"); source("_dgp_library.R")

SEED0 <- 20260918L; B <- 2000L
CELLS <- rbind(data.frame(n = 1000L, scen = c("null", "cloglog", "probit", "quad", "binint", "osc4")),
               data.frame(n = 500L,  scen = c("null", "cloglog")),
               data.frame(n = 5000L, scen = c("null", "probit")))
CELLS$G <- as.integer(round(CELLS$n / 25))

gen <- function(scn, n) {
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); eta <- 0.6 * x + 0.5 * d
  if (scn == "osc4") return(list(d = data.frame(y = rbinom(n, 1, plogis(0.8 * x + 1.5 * sin(4 * x))), x = x), f = y ~ x))
  pr <- switch(scn,
    null    = plogis(eta), cloglog = 1 - exp(-exp(eta)), probit = pnorm(eta),
    quad    = {b <- sqb(0.5); plogis(b[1] + b[2] * x + b[3] * x^2)},
    binint  = {b <- sib(0.4); plogis(b[1] + b[2] * x + b[3] * d + b[4] * x * d)})
  list(d = data.frame(y = rbinom(n, 1, pr), x = x, d = d), f = y ~ x + d)
}

## grouped quantities from the null fit, the package's grouping
grp_stats <- function(fit, G) {
  ph <- pmin(pmax(fitted(fit), 1e-6), 1 - 1e-6); y <- fit$y; X <- model.matrix(fit); w <- ph * (1 - ph)
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (length(ph) / G)), G)
  o <- as.numeric(tapply(y, g, sum)); e <- as.numeric(tapply(ph, g, sum)); V <- as.numeric(tapply(w, g, sum))
  pb <- as.numeric(tapply(ph, g, mean)); eb <- qlogis(pb)
  U <- rowsum(w * X, g) / sqrt(V)                       # G x p
  XWX <- crossprod(X * sqrt(w))
  list(o = o, e = e, V = V, pb = pb, eb = eb, r = (o - e) / sqrt(V), U = U, XWX = XWX)
}
## the paper's form: unit-weighted projection, eigen-weighted null (Satterthwaite)
edge_eig <- function(gs, Z) {
  Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]; if (!ncol(Z)) return(NA_real_)
  ZOZ <- crossprod(Z) - crossprod(Z, gs$U) %*% solve(gs$XWX, crossprod(gs$U, Z))
  lam <- Re(eigen(solve(crossprod(Z), ZOZ), only.values = TRUE)$values); lam <- lam[lam > 1e-9]
  if (!length(lam)) return(NA_real_)
  S <- sum(qr.fitted(qr(Z), gs$r)^2); cc <- sum(lam^2) / sum(lam); nu <- sum(lam)^2 / sum(lam^2)
  pchisq(S / cc, nu, lower.tail = FALSE)
}
## the score form: sqrt(V)-scaled basis, Omega-metric, plain chi^2_k
edge_score <- function(gs, Z) {
  Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]; if (!ncol(Z)) return(NA_real_)
  Zt <- Z * sqrt(gs$V)
  s <- crossprod(Zt, gs$r)                                                   # = Z'(o - e)
  I <- crossprod(Zt) - crossprod(Zt, gs$U) %*% solve(gs$XWX, crossprod(gs$U, Zt))
  S <- tryCatch(as.numeric(t(s) %*% solve(I, s)), error = function(e) NA_real_)
  if (!is.finite(S)) return(NA_real_)
  pchisq(S, ncol(Z), lower.tail = FALSE)
}
ao_col <- function(eb) { pb <- plogis(eb); 1 - log1p(exp(eb)) / pb }        # Aranda-Ordaz asymmetric score dir at lambda=1
one <- function(dat, G) {
  fit <- tryCatch(glm(dat$f, data = dat$d, family = binomial(), control = list(maxit = 50)), error = function(e) NULL)
  nm <- c("EDGE.poly3", "EDGE.poly3.logit", "EDGE.sc.poly", "EDGE.sc.stk", "EDGE.sc.AO", "EDGE.sc.AO2",
          "EDGE.sc.poly4", "HL", "Stukel", "GiViTI")
  if (is.null(fit)) return(setNames(rep(NA_real_, length(nm)), nm))
  gs <- grp_stats(fit, G); eb <- gs$eb
  stk2 <- cbind(eb^2 * (eb >= 0), -eb^2 * (eb < 0))
  setNames(c(
    edge_eig(gs, as.matrix(poly(gs$pb, 3))),
    edge_eig(gs, as.matrix(poly(eb, 3))),
    edge_score(gs, cbind(eb^2, eb^3)),
    edge_score(gs, stk2),
    edge_score(gs, cbind(ao_col(eb))),
    edge_score(gs, cbind(ao_col(eb), eb^2)),
    edge_score(gs, cbind(eb^2, eb^3, eb^4)),
    unname(hl_stat(fit$y, fitted(fit), G)["p"]),
    stukel_p(fit), giviti_p(fit)), nm)
}

cl <- makeCluster(max(1L, min(20L, detectCores() - 2L)))
invisible(clusterEvalQ(cl, {suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)}); NULL}))
clusterExport(cl, c("gen", "grp_stats", "edge_eig", "edge_score", "ao_col", "one", "hl_stat", "stukel_p", "giviti_p", "sqb", "sib", "lg"))
t0 <- Sys.time(); P <- list()
for (i in seq_len(nrow(CELLS))) {
  ce <- CELLS[i, ]; clusterSetRNGStream(cl, SEED0 + i)
  M <- t(parSapply(cl, seq_len(B), function(b, scn, n, G) one(gen(scn, n), G), scn = ce$scen, n = ce$n, G = ce$G))
  P[[i]] <- data.frame(scen = ce$scen, n = ce$n, G = ce$G, rep = seq_len(B), M)
  cat(sprintf("  [%d/%d] %-8s n=%5d G=%3d  %s\n", i, nrow(CELLS), ce$scen, ce$n, ce$G,
              format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); flush.console()
}
stopCluster(cl)
pv <- do.call(rbind, P); write.csv(pv, "runK_basis_score_pvalues.csv", row.names = FALSE)
TESTS <- setdiff(names(pv), c("scen", "n", "G", "rep"))
out <- list()
for (n in unique(CELLS$n)) {
  nul <- pv[pv$scen == "null" & pv$n == n, ]
  for (sc in setdiff(unique(CELLS$scen[CELLS$n == n]), "null")) { alt <- pv[pv$scen == sc & pv$n == n, ]
    for (t in TESTS) out[[length(out) + 1]] <- data.frame(n = n, G = alt$G[1], scen = sc, test = t,
      size = mean(nul[[t]] <= .05, na.rm = TRUE), power_adj = size_adjusted(alt[[t]], nul[[t]]),
      power_raw = mean(alt[[t]] <= .05, na.rm = TRUE), na = na_rate(alt[[t]])) } }
S <- do.call(rbind, out); write.csv(S, "runK_basis_score_summary.csv", row.names = FALSE)
cat("\n=== size-adjusted power (size in parentheses) ===\n")
for (sc in c("cloglog", "probit", "quad", "binint", "osc4")) for (n in unique(S$n[S$scen == sc])) {
  z <- S[S$scen == sc & S$n == n, ]
  cat(sprintf("\n-- %s n=%d G=%d --\n", sc, n, z$G[1]))
  for (t in TESTS) cat(sprintf("  %-18s %.3f (%.3f)  na=%.3f\n", t, z$power_adj[z$test == t], z$size[z$test == t], z$na[z$test == t]))
}
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
