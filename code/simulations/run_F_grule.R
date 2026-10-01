## run_F_grule.R -- is there an automatic G rule, and does EDGE still lose once G is chosen well?
##
## THE DISTINCTION THAT MAKES THIS LEGITIMATE
## ------------------------------------------
## The paper forbids choosing the BASIS from the data, because that is selection: you look at the
## residuals, pick the basis that rejects, and the null no longer holds. Choosing G from n is a
## different act. n is known before any outcome is seen, so a rule G = f(n) is deterministic given
## the design and introduces NO multiplicity. It is pre-specifiable in a protocol.
##
## What it IS outside of is the paper's own asymptotics: assumption (A1) says "the number of groups
## G is fixed". A rule that grows G with n sits outside that theorem, so its validity is an
## EMPIRICAL question and has to be established here, not assumed.
##
## So this script asks, in order:
##   F1. As the groups get finer, where does SIZE break? That is the binding constraint, and it is
##       naturally expressed in observations per group, m = n/G, not in G.
##   F2. Inside the size-valid region, which m maximises power?
##   F3. With that rule applied, does EDGE still lose to GiViTI and to Stukel?
##
## Everything is size-adjusted against a matched null at the SAME (n, G), so no test is credited
## with power it bought by rejecting too often.

suppressPackageStartupMessages({library(parallel)})
source("_pstar_giviti_harness.R")
source("_dgp_library.R")

SEED0 <- 20260914L
NS <- c(500L, 1000L, 2000L, 5000L, 20000L)
MS <- c(15L, 25L, 50L, 100L, 200L)          # target observations per group
SCEN <- c("null", "cloglog", "probit")

gen <- function(scn, n) {
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); eta <- 0.6 * x + 0.5 * d
  pr <- switch(scn, null = plogis(eta), probit = pnorm(eta), cloglog = linkp("cloglog", eta), stop(scn))
  list(d = data.frame(y = rbinom(n, 1, pr), x = x, d = d), f = y ~ x + d, cat = "d")
}

cells <- do.call(rbind, lapply(NS, function(n) {
  G <- pmax(4L, pmin(as.integer(round(n / MS)), 400L))
  unique(data.frame(n = n, m = MS, G = G, stringsAsFactors = FALSE))
}))
cells <- do.call(rbind, lapply(SCEN, function(s) cbind(scen = s, cells)))
cells$B <- ifelse(cells$n >= 20000, 800L, 2000L)
cat(sprintf("run F: %d cells\n", nrow(cells)))

ncores <- max(1L, min(20L, parallel::detectCores() - 2L))
cl <- parallel::makeCluster(ncores); on.exit(parallel::stopCluster(cl), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)}); NULL }))
parallel::clusterExport(cl, c("hl_stat", "stukel_p", "giviti_p", "one_rep", "gen", "linkp", "inv_stukel", "lg", "logit"))

t0 <- Sys.time(); raw <- list()
for (i in seq_len(nrow(cells))) {
  ce <- cells[i, ]
  parallel::clusterSetRNGStream(cl, SEED0 + i * 1031L)
  m <- t(parallel::parSapply(cl, seq_len(ce$B), function(b, s, n, G) one_rep(gen(s, n), G = G),
                             s = ce$scen, n = ce$n, G = ce$G))
  raw[[i]] <- data.frame(scen = ce$scen, n = ce$n, m = ce$m, G = ce$G, B = ce$B,
                         rep = seq_len(ce$B), m, check.names = FALSE)
  cat(sprintf("  [%2d/%2d] %-7s n=%5d m=%3d G=%3d  (%s)\n", i, nrow(cells), ce$scen, ce$n, ce$m, ce$G,
              format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); utils::flush.console()
}
pv <- do.call(rbind, raw)
utils::write.csv(pv, "runF_grule_pvalues.csv", row.names = FALSE)

out <- list()
for (n in NS) for (mm in MS) {
  G <- unique(cells$G[cells$n == n & cells$m == mm]); if (!length(G)) next
  nb <- pv[pv$scen == "null" & pv$n == n & pv$m == mm, ]
  if (!nrow(nb)) next
  for (sc in c("cloglog", "probit")) {
    ab <- pv[pv$scen == sc & pv$n == n & pv$m == mm, ]
    for (t in TESTS) out[[length(out) + 1]] <- data.frame(
      n = n, m = mm, G = G[1], scen = sc, test = t, B = nb$B[1],
      size = raw_power(nb[[t]]), power_adj = size_adjusted(ab[[t]], nb[[t]]),
      na_null = na_rate(nb[[t]]))
  }
}
S <- do.call(rbind, out)
utils::write.csv(S, "runF_grule_summary.csv", row.names = FALSE)

cat("\n============================================================\n")
cat(" F1. SIZE of EDGE-poly3 as the groups get finer (nominal 0.05)\n")
cat("     rows = n, cols = observations per group. This is the constraint.\n")
cat("============================================================\n")
M <- matrix(NA_real_, length(NS), length(MS), dimnames = list(paste0("n=", NS), paste0("m=", MS)))
for (i in seq_along(NS)) for (j in seq_along(MS)) {
  v <- S$size[S$n == NS[i] & S$m == MS[j] & S$test == "EDGE.poly3" & S$scen == "cloglog"]
  if (length(v)) M[i, j] <- round(v[1], 3)
}
print(M)
cat("\n  G actually used:\n")
GM <- M; for (i in seq_along(NS)) for (j in seq_along(MS)) {
  v <- S$G[S$n == NS[i] & S$m == MS[j] & S$test == "EDGE.poly3" & S$scen == "cloglog"]
  GM[i, j] <- if (length(v)) v[1] else NA }
print(GM)

for (sc in c("cloglog", "probit")) {
  cat(sprintf("\n============================================================\n F2. POWER on %s, size-adjusted. rows = n, cols = obs/group\n============================================================\n", sc))
  for (t in c("EDGE.poly3", "EDGE.stk", "GiViTI", "Stukel", "HL")) {
    P <- M
    for (i in seq_along(NS)) for (j in seq_along(MS)) {
      v <- S$power_adj[S$n == NS[i] & S$m == MS[j] & S$test == t & S$scen == sc]
      P[i, j] <- if (length(v)) round(v[1], 3) else NA }
    cat("\n ", t, "\n"); print(P)
  }
}
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
