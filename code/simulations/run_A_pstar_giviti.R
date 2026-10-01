## run_A_pstar_giviti.R -- Reviewer 1 points 2 and 3, answered in one sweep.
##
## POINT 3 (the number of covariates is never varied): p* enters the null BY CONSTRUCTION,
## since Omega = I_G - U (X'WX)^{-1} U' subtracts a rank-p* matrix from the G x G residual
## covariance. As p* grows toward G the correction removes more of the residual space, and
## when G <= p* the construction has no room left. That is a real boundary and it has never
## been measured. Design: hold the SIGNAL fixed and add pure-noise covariates to the fitted
## model, so the only thing changing is p*.
##
## POINT 2 (GiViTI is cited 10x as the motivation and never run): it is in the battery here,
## on identical samples, so the comparison is paired rather than across runs.
##
## Every replication's p-values are written to disk (standing project rule), not just the
## summaries.

suppressPackageStartupMessages({library(parallel)})
source("_pstar_giviti_harness.R")
source("_dgp_library.R")

N      <- 1000L
B      <- 2000L
PSTARS <- c(3L, 5L, 10L, 20L)
GRIDG  <- c(10L, 20L)
SEED0  <- 20260910L

## --- DGPs with the signal fixed and p* - p_base pure-noise covariates appended -----------
## scenario "null" is correctly specified; the others carry the departure on x1 only.
gen_pstar <- function(scenario, n, pstar) {
  x1 <- runif(n, -3, 3)
  d  <- rbinom(n, 1, 0.5)
  eta_base <- 0.6 * x1 + 0.5 * d
  pr <- switch(scenario,
    null    = plogis(eta_base),
    cloglog = linkp("cloglog", eta_base),
    quad    = { b <- sqb(0.5); plogis(b[1] + b[2] * x1 + b[3] * x1^2) },
    binint  = { b <- sib(0.4); plogis(b[1] + b[2] * x1 + b[3] * d + b[4] * x1 * d) },
    stop("scenario: ", scenario))
  df <- data.frame(y = rbinom(n, 1, pr), x1 = x1, d = d)
  ## p_base = 3 estimated coefficients (intercept, x1, d); pad with pure noise
  k <- pstar - 3L
  if (k > 0) {
    Z <- matrix(rnorm(n * k), n, k); colnames(Z) <- paste0("z", seq_len(k))
    df <- cbind(df, Z)
  }
  f <- stats::as.formula(paste("y ~ x1 + d", if (k > 0) paste("+", paste0("z", seq_len(k), collapse = " + ")) else ""))
  list(d = df, f = f, cat = "d")
}

SCEN <- c("null", "cloglog", "quad", "binint")
cells <- expand.grid(scenario = SCEN, pstar = PSTARS, G = GRIDG,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
cat(sprintf("run A: %d cells x B=%d  (n=%d)\n", nrow(cells), B, N))

ncores <- max(1L, min(20L, parallel::detectCores() - 2L))
cl <- parallel::makeCluster(ncores)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)})
  NULL
})
parallel::clusterExport(cl, c("hl_stat", "stukel_p", "giviti_p", "one_rep", "gen_pstar",
                              "linkp", "sqb", "sib", "inv_stukel", "lg", "logit", "N"))

t_start <- Sys.time()
raw <- vector("list", nrow(cells))
for (i in seq_len(nrow(cells))) {
  ce <- cells[i, ]
  parallel::clusterSetRNGStream(cl, SEED0 + i * 1000L)
  m <- parallel::parSapply(cl, seq_len(B), function(b, sc, ps, G, n) {
    one_rep(gen_pstar(sc, n, ps), G = G)
  }, sc = ce$scenario, ps = ce$pstar, G = ce$G, n = N)
  m <- t(m)
  raw[[i]] <- data.frame(scenario = ce$scenario, pstar = ce$pstar, G = ce$G,
                         rep = seq_len(B), m, check.names = FALSE)
  cat(sprintf("  [%2d/%2d] %-8s p*=%2d G=%2d  done  (%s)\n", i, nrow(cells),
              ce$scenario, ce$pstar, ce$G,
              format(round(difftime(Sys.time(), t_start, units = "mins"), 1))))
  utils::flush.console()
}
pv <- do.call(rbind, raw)
utils::write.csv(pv, "runA_pstar_giviti_pvalues.csv", row.names = FALSE)
cat("per-rep p-values written: runA_pstar_giviti_pvalues.csv (", nrow(pv), "rows )\n")

## --- summarise: raw size, raw power, size-adjusted power, and the NA (declined) rate -----
out <- list()
for (ps in PSTARS) for (G in GRIDG) {
  nullblk <- pv[pv$scenario == "null" & pv$pstar == ps & pv$G == G, ]
  for (sc in setdiff(SCEN, "null")) {
    altblk <- pv[pv$scenario == sc & pv$pstar == ps & pv$G == G, ]
    for (tst in TESTS) {
      out[[length(out) + 1]] <- data.frame(
        scenario = sc, pstar = ps, G = G, test = tst,
        size_raw   = raw_power(nullblk[[tst]]),
        power_raw  = raw_power(altblk[[tst]]),
        power_adj  = size_adjusted(altblk[[tst]], nullblk[[tst]]),
        na_null    = na_rate(nullblk[[tst]]),
        na_alt     = na_rate(altblk[[tst]]),
        B = B)
    }
  }
}
S <- do.call(rbind, out)
S$mcse_size <- mcse(S$size_raw, B)
utils::write.csv(S, "runA_pstar_giviti_summary.csv", row.names = FALSE)
cat("summary written: runA_pstar_giviti_summary.csv\n")
cat("\ntotal wall time:", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)), "\n")
