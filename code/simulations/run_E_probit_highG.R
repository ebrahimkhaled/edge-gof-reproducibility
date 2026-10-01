## run_E_probit_highG.R -- the one combination not yet measured: probit at LARGE n AND large G.
##
## Two separate advantages have been established for the directed test:
##   (i) on probit it beats GiViTI at every n (n=20,000: 0.703 vs 0.531), and
##  (ii) raising G buys it power on cloglog (0.667 -> 0.816 from G=5 to G=20) while the omnibus
##       tests collapse.
## Do those stack? GiViTI has no partition, so its column is constant by construction and serves
## as a fixed reference line across the sweep.

suppressPackageStartupMessages({library(parallel)})
source("_pstar_giviti_harness.R")
source("_dgp_library.R")

N <- 20000L
GS <- c(10L, 20L, 40L, 80L)
B <- 800L
SEED0 <- 20260913L

gen <- function(scn, n) {
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); eta <- 0.6 * x + 0.5 * d
  pr <- if (scn == "probit") pnorm(eta) else plogis(eta)
  list(d = data.frame(y = rbinom(n, 1, pr), x = x, d = d), f = y ~ x + d, cat = "d")
}

cells <- expand.grid(scen = c("null", "probit"), G = GS, stringsAsFactors = FALSE)
cat(sprintf("run E: %d cells, n=%d, B=%d\n", nrow(cells), N, B))

ncores <- max(1L, min(20L, parallel::detectCores() - 2L))
cl <- parallel::makeCluster(ncores); on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterEvalQ(cl, {suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)}); NULL})
parallel::clusterExport(cl, c("hl_stat", "stukel_p", "giviti_p", "one_rep", "gen", "N"))

t0 <- Sys.time(); raw <- list()
for (i in seq_len(nrow(cells))) {
  ce <- cells[i, ]
  parallel::clusterSetRNGStream(cl, SEED0 + i * 1013L)
  m <- t(parallel::parSapply(cl, seq_len(B), function(b, s, G, n) one_rep(gen(s, n), G = G),
                             s = ce$scen, G = ce$G, n = N))
  raw[[i]] <- data.frame(scen = ce$scen, G = ce$G, rep = seq_len(B), m, check.names = FALSE)
  cat(sprintf("  [%d/%d] %-7s G=%2d  (%s)\n", i, nrow(cells), ce$scen, ce$G,
              format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); utils::flush.console()
}
pv <- do.call(rbind, raw)
utils::write.csv(pv, "runE_probit_highG_pvalues.csv", row.names = FALSE)

out <- list()
for (G in GS) {
  nb <- pv[pv$scen == "null" & pv$G == G, ]; ab <- pv[pv$scen == "probit" & pv$G == G, ]
  for (t in TESTS) out[[length(out) + 1]] <- data.frame(
    G = G, test = t, size = raw_power(nb[[t]]), power_adj = size_adjusted(ab[[t]], nb[[t]]))
}
S <- do.call(rbind, out)
utils::write.csv(S, "runE_probit_highG_summary.csv", row.names = FALSE)

cat("\n=== PROBIT at n=20,000: size-adjusted power as the partition is refined ===\n")
cat(sprintf("  %-12s %8s %8s %8s %8s\n", "test", "G=10", "G=20", "G=40", "G=80"))
for (t in c("EDGE.poly3", "EDGE.stk", "EDGE.poly2", "GiViTI", "Stukel", "HL_F", "HL")) {
  v <- sapply(GS, function(g) S$power_adj[S$G == g & S$test == t])
  cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", t, v[1], v[2], v[3], v[4]))
}
cat(sprintf("\n  MCSE at B=%d is about %.3f\n", B, sqrt(.25 / B)))
cat("\n  GiViTI has no partition, so its row is one test measured four times:\n")
gv <- sapply(GS, function(g) S$power_adj[S$G == g & S$test == "GiViTI"])
cat(sprintf("    range %.3f against 2 MCSE of %.3f -> %s\n", max(gv) - min(gv), 2 * sqrt(.25 / B),
            ifelse(max(gv) - min(gv) < 3 * 2 * sqrt(.25 / B), "noise, as expected", "unexpected trend")))
cat("\n  best directed basis minus GiViTI, by G:\n")
for (g in GS) {
  be <- max(sapply(c("EDGE.poly3", "EDGE.stk", "EDGE.poly2"),
                   function(t) S$power_adj[S$G == g & S$test == t]))
  cat(sprintf("    G=%2d  %+.3f\n", g, be - S$power_adj[S$G == g & S$test == "GiViTI"]))
}
cat("\n=== size check (nominal 0.05) ===\n")
print(S[S$test %in% c("EDGE.poly3", "GiViTI", "HL"), c("G", "test", "size")], row.names = FALSE, digits = 3)
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
