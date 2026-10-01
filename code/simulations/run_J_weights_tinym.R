## run_J_weights_tinym.R -- the growing-partition theorem's two testable predictions, checked
## on the base design at n = 2000:
##   (1) the plug-in weights lambda_j(G) = eig[(Z'Z)^{-1} Z' Omega Z] CONVERGE as G grows;
##   (2) the closed-form null holds at ANY partition, including m = n/G = 5 and 2 observations per
##       group, because the only normality the statistic uses is that of k weighted sums over
##       all n observations, not per-group normality.
## Size-adjusted power against the cloglog bow is recorded at the same G's, to see the plateau
## the theorem predicts once the finite-G discretisation loss has been paid off.
## The package does not return the weights, so they are recomputed here from the fitted model;
## the self-check at the top stops the run if this recomputation disagrees with edge.gof().
suppressPackageStartupMessages({library(ebrahim.gof); library(parallel)})

n <- 2000L; B <- 2000L; SEED <- 20260913L
GS <- c(10L, 20L, 50L, 100L, 200L, 400L, 1000L)          # m = 200 ... 2

gen <- function(n, scen) {
  x <- runif(n, -3, 3); d <- rbinom(n, 1, .5); eta <- 0.6 * x + 0.5 * d
  p <- switch(scen, null = plogis(eta), cloglog = 1 - exp(-exp(eta)))
  data.frame(y = rbinom(n, 1, p), x = x, d = d)
}
## EDGE-poly3 by hand: groups, residual, basis, Omega, weights, Satterthwaite p
edge_hand <- function(fit, G) {
  ph <- fitted(fit); y <- fit$y; X <- model.matrix(fit); w <- ph * (1 - ph)
  ph <- pmin(pmax(ph, 1e-6), 1 - 1e-6); w <- ph * (1 - ph)
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (length(ph) / G)), G)   # the package's grouping
  o <- tapply(y, g, sum); e <- tapply(ph, g, sum); V <- tapply(w, g, sum); pb <- tapply(ph, g, mean)
  r <- as.numeric((o - e) / sqrt(V))
  Z <- poly(as.numeric(pb), 3)
  U <- (rowsum(w * X, g) / sqrt(as.numeric(V)))
  XWX <- crossprod(X * sqrt(w))
  ZOZ <- crossprod(Z) - crossprod(Z, U) %*% solve(XWX, crossprod(U, Z))
  lam <- sort(Re(eigen(solve(crossprod(Z), ZOZ), only.values = TRUE)$values), decreasing = TRUE)
  S <- sum(qr.fitted(qr(Z), r)^2)
  cc <- sum(lam^2) / sum(lam); nu <- sum(lam)^2 / sum(lam^2)
  c(p = pchisq(S / cc, nu, lower.tail = FALSE), S = S, l1 = lam[1], l2 = lam[2], l3 = lam[3])
}
## ---- self-check against the package at G = 10 -------------------------------------------
set.seed(1); D <- gen(n, "null"); f0 <- glm(y ~ x + d, data = D, family = binomial())
a <- edge_hand(f0, 10); b <- edge.gof(f0, G = 10, basis = "poly3")
cat(sprintf("self-check G=10: hand S=%.6f p=%.6f | package S=%.6f p=%.6f\n", a["S"], a["p"], b$Test_Statistic, b$p_value))
if (abs(a["p"] - b$p_value) > 1e-4) stop("hand recomputation does not match edge.gof(); fix before trusting the weights")

cl <- makeCluster(max(1L, min(20L, detectCores() - 2L)))
clusterExport(cl, c("gen", "edge_hand", "GS", "n"))
P <- list()
for (scen in c("null", "cloglog")) {
  clusterSetRNGStream(cl, SEED + (scen == "cloglog"))
  M <- t(parSapply(cl, seq_len(B), function(b, scen) {
    D <- gen(n, scen); fit <- glm(y ~ x + d, data = D, family = binomial())
    unlist(lapply(GS, function(G) tryCatch(edge_hand(fit, G), error = function(e) c(p = NA, S = NA, l1 = NA, l2 = NA, l3 = NA))))
  }, scen = scen))
  P[[scen]] <- M; cat(scen, "done\n"); flush.console()
}
stopCluster(cl)
res <- list()
for (j in seq_along(GS)) {
  cp <- (j - 1) * 5 + 1; cl1 <- cp + 2
  pn <- P$null[, cp]; pa <- P$cloglog[, cp]
  crit <- quantile(pn[is.finite(pn)], 0.05, type = 1)
  res[[j]] <- data.frame(n = n, G = GS[j], m = round(n / GS[j]),
    size = mean(pn <= 0.05, na.rm = TRUE), na_null = mean(!is.finite(pn)),
    power_raw = mean(pa <= 0.05, na.rm = TRUE), power_adj = mean(pa[is.finite(pa)] <= crit),
    l1 = mean(P$null[, cl1], na.rm = TRUE), l2 = mean(P$null[, cl1 + 1], na.rm = TRUE),
    l3 = mean(P$null[, cl1 + 2], na.rm = TRUE), l3_sd = sd(P$null[, cl1 + 2], na.rm = TRUE))
}
R <- do.call(rbind, res)
write.csv(R, "runJ_weights_tinym.csv", row.names = FALSE)
saveRDS(P, "runJ_weights_tinym_pvalues.rds")
print(R, digits = 4)
