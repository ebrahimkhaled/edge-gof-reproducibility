## run_D_probit_G_giviti.R -- three questions the author asked on 2026-09-12.
##
## D1. PROBIT. probe_probit.R shows the population signal is 74% cubic, i.e. ALREADY inside the
##     poly3 basis, and 4.4x smaller in norm than cloglog (=> ~20x smaller non-centrality). So the
##     prediction is: a probit-specific basis cannot help, but MORE DATA can, and the requirement
##     should be roughly 20x the cloglog sample size. Test it by sweeping n.
##
## D2. GROUPS. Does raising G past 10 buy power? EDGE spends a fixed k<=3 degrees of freedom
##     whatever G is, so unlike HL it should not be penalised for a finer partition -- but a finer
##     partition also means fewer observations per group and a noisier residual. Sweep G.
##
## D3. GiViTI in the two regimes where it was never tested: rough misfit, and sparse data with a
##     low event rate. Those are the only places an honest case for EDGE over GiViTI could live.

suppressPackageStartupMessages({library(parallel)})
source("_pstar_giviti_harness.R")
source("_dgp_library.R")

SEED0 <- 20260912L

gen <- function(scn, n) {
  if (scn == "sparse_cloglog") {              # low event rate, where refits separate
    x <- as.numeric(scale(rchisq(n, 4)))
    pr <- 1 - exp(-exp(-4.0 + 0.9 * x))
    return(list(d = data.frame(y = rbinom(n, 1, pr), x = x), f = y ~ x, cat = NA))
  }
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); eta <- 0.6 * x + 0.5 * d
  pr <- switch(scn,
    null     = plogis(eta),
    probit   = pnorm(eta),
    cloglog  = linkp("cloglog", eta),
    osc4     = plogis(0.8 * x + 1.5 * sin(4 * x)),
    sawtooth = plogis(0.8 * x + 1.2 * (2 * (x / 1.5 - floor(x / 1.5 + 0.5)))),
    stop(scn))
  if (scn %in% c("osc4", "sawtooth"))
    return(list(d = data.frame(y = rbinom(n, 1, pr), x = x), f = y ~ x, cat = NA))
  list(d = data.frame(y = rbinom(n, 1, pr), x = x, d = d), f = y ~ x + d, cat = "d")
}

## null partner for size adjustment: same design, correctly specified
gen_null_for <- function(scn, n) {
  if (scn == "sparse_cloglog") {
    x <- as.numeric(scale(rchisq(n, 4)))
    return(list(d = data.frame(y = rbinom(n, 1, plogis(-4.0 + 0.9 * x)), x = x), f = y ~ x, cat = NA))
  }
  if (scn %in% c("osc4", "sawtooth")) {
    x <- runif(n, -3, 3)
    return(list(d = data.frame(y = rbinom(n, 1, plogis(0.8 * x)), x = x), f = y ~ x, cat = NA))
  }
  gen("null", n)
}

CELLS <- rbind(
  ## D1 -- probit needs how much data?  cloglog at n=1000 is the yardstick
  data.frame(block = "D1_probit",  scen = "probit",  n = c(1000, 5000, 20000, 50000), G = 10, B = 2000),
  data.frame(block = "D1_probit",  scen = "cloglog", n = 1000,                        G = 10, B = 2000),
  ## D2 -- does a finer partition buy power?
  data.frame(block = "D2_groups",  scen = "cloglog", n = 1000, G = c(5, 10, 20, 40, 60), B = 2000),
  ## D3 -- GiViTI where it was never tested
  data.frame(block = "D3_giviti",  scen = "osc4",           n = 1000, G = 10, B = 2000),
  data.frame(block = "D3_giviti",  scen = "sawtooth",       n = 1000, G = 10, B = 2000),
  data.frame(block = "D3_giviti",  scen = "sparse_cloglog", n = c(500, 2000), G = 10, B = 2000)
)
CELLS$B[CELLS$n >= 20000] <- 800L     # the big-n cells are expensive; 800 still gives MCSE ~.017
cat(sprintf("run D: %d cells\n", nrow(CELLS)))

ncores <- max(1L, min(20L, parallel::detectCores() - 2L))
cl <- parallel::makeCluster(ncores); on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterEvalQ(cl, {suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)}); NULL})
parallel::clusterExport(cl, c("hl_stat", "stukel_p", "giviti_p", "one_rep", "gen", "gen_null_for",
                              "linkp", "inv_stukel", "lg", "logit"))

t0 <- Sys.time(); out <- list()
for (i in seq_len(nrow(CELLS))) {
  ce <- CELLS[i, ]
  parallel::clusterSetRNGStream(cl, SEED0 + i * 977L)
  alt <- t(parallel::parSapply(cl, seq_len(ce$B), function(b, s, n, G) one_rep(gen(s, n), G = G),
                               s = ce$scen, n = ce$n, G = ce$G))
  parallel::clusterSetRNGStream(cl, SEED0 + i * 977L + 500000L)
  nul <- t(parallel::parSapply(cl, seq_len(ce$B), function(b, s, n, G) one_rep(gen_null_for(s, n), G = G),
                               s = ce$scen, n = ce$n, G = ce$G))
  for (tst in TESTS) {
    out[[length(out) + 1]] <- data.frame(
      block = ce$block, scen = ce$scen, n = ce$n, G = ce$G, B = ce$B, test = tst,
      size = raw_power(nul[, tst]), power_raw = raw_power(alt[, tst]),
      power_adj = size_adjusted(alt[, tst], nul[, tst]),
      na_alt = na_rate(alt[, tst]), na_null = na_rate(nul[, tst]))
  }
  cat(sprintf("  [%2d/%2d] %-16s n=%5d G=%2d B=%4d  (%s)\n", i, nrow(CELLS), ce$scen, ce$n, ce$G, ce$B,
              format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); utils::flush.console()
}
R <- do.call(rbind, out)
utils::write.csv(R, "runD_probit_G_giviti.csv", row.names = FALSE)

sh <- function(blk, cap) {
  cat("\n==============================================================\n", cap, "\n")
  z <- R[R$block == blk & R$test %in% c("EDGE.poly2","EDGE.poly3","EDGE.stk","GiViTI","Stukel","HL","HL_F"), ]
  z$cell <- if (blk == "D2_groups") paste0("G=", z$G) else paste0(z$scen, " n=", z$n)
  w <- reshape(z[, c("test","cell","power_adj")], idvar = "test", timevar = "cell", direction = "wide")
  names(w) <- sub("power_adj.", "", names(w), fixed = TRUE)
  print(w, row.names = FALSE, digits = 3)
}
sh("D1_probit", "D1 -- PROBIT: size-adjusted power as n grows (cloglog n=1000 is the yardstick)")
sh("D2_groups", "D2 -- GROUPS: size-adjusted power on cloglog at n=1000")
sh("D3_giviti", "D3 -- GiViTI where it was never tested")

cat("\n-- D3 decline rates (a test that will not run is not a comparator) --\n")
z <- R[R$block == "D3_giviti" & R$test %in% c("GiViTI","Stukel","EDGE.poly3","HL"), c("scen","n","test","na_alt","size")]
print(z, row.names = FALSE, digits = 3)
cat("\nwritten: runD_probit_G_giviti.csv | wall", format(round(difftime(Sys.time(), t0, units="mins"),1)), "\n")
