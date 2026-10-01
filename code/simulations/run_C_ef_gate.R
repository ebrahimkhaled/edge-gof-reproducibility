## run_C_ef_gate.R -- THE GATE. Run this BEFORE deciding whether an EF paper exists.
##
## Settled already, and not re-litigated here: as an omnibus test, T_HL_F = T_HL - C cannot be
## first-order superior to HL against anything, because a linear term folded inside a quadratic
## form contributes to the non-centrality only at second order. That is why the "EF is a more
## powerful test" claim is dead and why the paper was rejected.
##
## The one question left open: the correction is O(G/sqrt(n)), so its effect is a function of
## BOTH G and n -- and nothing in this project has ever mapped that surface. If the HL_F - HL
## power gap is flat and near zero everywhere, there is no paper and EF is finished. If the gap
## has structure in (n, G) -- in particular if it REVERSES SIGN -- then there is a genuine,
## checkable, second-order finding: "the moment correction helps or hurts depending on G/sqrt(n),
## and here is where the boundary is". That would be a short honest note, not a new test.
##
## Design: cloglog (the one departure where HL_F has ever beaten HL) plus the matched null at
## every (n, G), so power is SIZE-ADJUSTED. A raw-power comparison between two tests with
## different realised size is not a comparison.

suppressPackageStartupMessages({library(parallel)})
source("_pstar_giviti_harness.R")
source("_dgp_library.R")

NS    <- c(200L, 500L, 1000L, 2000L, 5000L)
GS    <- c(6L, 8L, 10L, 12L, 16L, 20L)
B     <- 1500L
SEED0 <- 20260911L

gen_gate <- function(scenario, n) {
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
  eta <- 0.6 * x + 0.5 * d
  pr <- if (scenario == "null") plogis(eta) else linkp("cloglog", eta)
  list(d = data.frame(y = rbinom(n, 1, pr), x = x, d = d), f = y ~ x + d, cat = "d")
}

cells <- expand.grid(scenario = c("null", "cloglog"), n = NS, G = GS,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
cat(sprintf("run C (EF gate): %d cells x B=%d\n", nrow(cells), B))

ncores <- max(1L, min(20L, parallel::detectCores() - 2L))
cl <- parallel::makeCluster(ncores)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)}); NULL })
parallel::clusterExport(cl, c("hl_stat", "stukel_p", "giviti_p", "one_rep", "gen_gate",
                              "linkp", "inv_stukel", "lg", "logit"))

t0 <- Sys.time()
raw <- vector("list", nrow(cells))
for (i in seq_len(nrow(cells))) {
  ce <- cells[i, ]
  parallel::clusterSetRNGStream(cl, SEED0 + i * 1000L)
  m <- t(parallel::parSapply(cl, seq_len(B), function(b, sc, n, G) {
    one_rep(gen_gate(sc, n), G = G)
  }, sc = ce$scenario, n = ce$n, G = ce$G))
  raw[[i]] <- data.frame(scenario = ce$scenario, n = ce$n, G = ce$G,
                         rep = seq_len(B), m, check.names = FALSE)
  cat(sprintf("  [%2d/%2d] %-8s n=%5d G=%2d  (%s)\n", i, nrow(cells), ce$scenario, ce$n, ce$G,
              format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); utils::flush.console()
}
pv <- do.call(rbind, raw)
utils::write.csv(pv, "runC_ef_gate_pvalues.csv", row.names = FALSE)
cat("per-rep p-values written (", nrow(pv), "rows )\n")

out <- list()
for (n in NS) for (G in GS) {
  nb <- pv[pv$scenario == "null"    & pv$n == n & pv$G == G, ]
  ab <- pv[pv$scenario == "cloglog" & pv$n == n & pv$G == G, ]
  row <- data.frame(n = n, G = G, GoverSqrtN = G / sqrt(n))
  for (tst in c("HL_F", "HL", "EDGE.poly3", "Stukel", "GiViTI")) {
    row[[paste0("size_", tst)]]  <- raw_power(nb[[tst]])
    row[[paste0("adj_",  tst)]]  <- size_adjusted(ab[[tst]], nb[[tst]])
    row[[paste0("raw_",  tst)]]  <- raw_power(ab[[tst]])
  }
  row$gap_adj <- row$adj_HL_F - row$adj_HL          # THE quantity the gate is about
  row$gap_raw <- row$raw_HL_F - row$raw_HL
  row$mcse_gap <- sqrt(2 * 0.25 / B)                # conservative bound on the paired gap
  out[[length(out) + 1]] <- row
}
S <- do.call(rbind, out)
utils::write.csv(S, "runC_ef_gate_summary.csv", row.names = FALSE)

cat("\n=== SIZE-ADJUSTED HL_F - HL POWER GAP ON CLOGLOG (the gate) ===\n")
cat("conservative MCSE on each gap:", round(sqrt(2 * 0.25 / B), 4), "\n\n")
M <- matrix(NA_real_, length(NS), length(GS), dimnames = list(paste0("n=", NS), paste0("G=", GS)))
for (i in seq_along(NS)) for (j in seq_along(GS)) {
  v <- S$gap_adj[S$n == NS[i] & S$G == GS[j]]
  if (length(v)) M[i, j] <- round(v, 3)
}
print(M)
cat("\nrange of the gap:", sprintf("%.3f to %.3f", min(M, na.rm = TRUE), max(M, na.rm = TRUE)), "\n")
cat("cells where the gap exceeds 2 MCSE in MAGNITUDE:",
    sum(abs(M) > 2 * sqrt(2 * 0.25 / B), na.rm = TRUE), "of", sum(!is.na(M)), "\n")
cat("sign changes present:", (min(M, na.rm = TRUE) < -2 * sqrt(2 * 0.25 / B)) &&
                             (max(M, na.rm = TRUE) >  2 * sqrt(2 * 0.25 / B)), "\n")
cat("\nsummary written: runC_ef_gate_summary.csv\n")
cat("total wall time:", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
