## run_G_validate_rule.R -- HOLD-OUT validation of whatever G rule run F suggests.
##
## Why this script exists at all. Run F sweeps (n, m) on cloglog and probit and will show which m
## gives the most power. Reporting that same maximum as the rule's performance would be circular:
## the rule was chosen to maximise it. So the rule is fixed from run F and then evaluated HERE, on
##   (a) departure families run F never saw -- omitted quadratic, omitted interaction, a rough
##       oscillation -- and
##   (b) a sparse, low-event-rate design, where the binding constraint is EVENTS per group rather
##       than observations per group and any rule fitted at a 50% event rate should be expected to
##       strain, and
##   (c) fresh random seeds throughout.
##
## The comparison reported is the one the author actually asked for: with G chosen by the rule
## rather than left at 10, does EDGE still lose to GiViTI and to Stukel?
##
## RULE_M is set from run F's output before this is run. It is a number, fixed in advance of
## seeing any of the scenarios below.

suppressPackageStartupMessages({library(parallel)})
source("_pstar_giviti_harness.R")
source("_dgp_library.R")

RULE_M <- as.integer(Sys.getenv("RULE_M", "25"))   # observations per group, set from run F
## GMAX was 200 in the first run (2026-09-12), which silently made the n=10,000 "rule" arm G=200
## (50 per group) instead of the rule's own 400 -- caught by the numbers audit of 2026-09-13.
## Re-run with GMAX=400 so the rule is validated at its own m at every n.
GMIN <- 6L; GMAX <- as.integer(Sys.getenv("GMAX", "400"))
g_of <- function(n) max(GMIN, min(GMAX, as.integer(round(n / RULE_M))))

SEED0 <- 20260915L                                  # fresh seeds, not run F's
NS <- c(500L, 2000L, 10000L)
SCEN <- c("null", "quad", "binint", "osc4", "null_sparse", "cloglog_sparse")

gen <- function(scn, n) {
  if (grepl("sparse", scn)) {                       # event rate ~1.5%: events per group binds
    x <- as.numeric(scale(rchisq(n, 4)))
    eta <- -4.0 + 0.9 * x
    pr <- if (scn == "null_sparse") plogis(eta) else 1 - exp(-exp(eta))
    return(list(d = data.frame(y = rbinom(n, 1, pr), x = x), f = y ~ x, cat = NA))
  }
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
  if (scn == "osc4") return(list(d = data.frame(y = rbinom(n, 1, plogis(0.8 * x + 1.5 * sin(4 * x))), x = x),
                                 f = y ~ x, cat = NA))
  pr <- switch(scn,
    null   = plogis(0.6 * x + 0.5 * d),
    quad   = {b <- sqb(0.5); plogis(b[1] + b[2] * x + b[3] * x^2)},
    binint = {b <- sib(0.4); plogis(b[1] + b[2] * x + b[3] * d + b[4] * x * d)},
    stop(scn))
  list(d = data.frame(y = rbinom(n, 1, pr), x = x, d = d), f = y ~ x + d, cat = "d")
}

## every scenario is run at BOTH the conventional G=10 and the rule's G, so the comparison is
## like-for-like and the rule has to earn its keep against the status quo.
cells <- expand.grid(scen = SCEN, n = NS, arm = c("G10", "rule"), stringsAsFactors = FALSE)
cells$G <- ifelse(cells$arm == "G10", 10L, sapply(cells$n, g_of))
cells$B <- ifelse(cells$n >= 10000, 1000L, 2000L)
cells <- cells[!duplicated(cells[, c("scen", "n", "G")]), ]
cat(sprintf("run G: RULE_M=%d obs/group -> G = %s for n = %s\n", RULE_M,
            paste(sapply(NS, g_of), collapse = "/"), paste(NS, collapse = "/")))
cat(sprintf("  %d cells\n", nrow(cells)))

ncores <- max(1L, min(20L, parallel::detectCores() - 2L))
cl <- parallel::makeCluster(ncores); on.exit(parallel::stopCluster(cl), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)}); NULL }))
parallel::clusterExport(cl, c("hl_stat", "stukel_p", "giviti_p", "one_rep", "gen",
                              "sqb", "sib", "lg", "logit", "inv_stukel", "linkp"))

t0 <- Sys.time(); raw <- list()
for (i in seq_len(nrow(cells))) {
  ce <- cells[i, ]
  parallel::clusterSetRNGStream(cl, SEED0 + i * 1049L)
  mm <- t(parallel::parSapply(cl, seq_len(ce$B), function(b, s, n, G) one_rep(gen(s, n), G = G),
                              s = ce$scen, n = ce$n, G = ce$G))
  raw[[i]] <- data.frame(scen = ce$scen, n = ce$n, G = ce$G, arm = ce$arm, B = ce$B,
                         rep = seq_len(ce$B), mm, check.names = FALSE)
  cat(sprintf("  [%2d/%2d] %-15s n=%5d G=%3d (%s)  %s\n", i, nrow(cells), ce$scen, ce$n, ce$G, ce$arm,
              format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); utils::flush.console()
}
pv <- do.call(rbind, raw)
utils::write.csv(pv, "runG_validate_pvalues.csv", row.names = FALSE)

nullof <- function(s) if (grepl("sparse", s)) "null_sparse" else "null"
out <- list()
for (n in NS) for (arm in c("G10", "rule")) {
  G <- unique(cells$G[cells$n == n & cells$arm == arm]); if (!length(G)) next
  for (sc in setdiff(SCEN, c("null", "null_sparse"))) {
    nb <- pv[pv$scen == nullof(sc) & pv$n == n & pv$G == G[1], ]
    ab <- pv[pv$scen == sc & pv$n == n & pv$G == G[1], ]
    if (!nrow(nb) || !nrow(ab)) next
    for (t in TESTS) out[[length(out) + 1]] <- data.frame(
      n = n, arm = arm, G = G[1], scen = sc, test = t,
      size = raw_power(nb[[t]]), power_adj = size_adjusted(ab[[t]], nb[[t]]),
      na_alt = na_rate(ab[[t]]))
  }
}
S <- do.call(rbind, out)
utils::write.csv(S, "runG_validate_summary.csv", row.names = FALSE)

cat("\n=================================================================\n")
cat(" HOLD-OUT: does the rule beat the conventional G=10?  (EDGE-poly3)\n")
cat("=================================================================\n")
for (sc in setdiff(SCEN, c("null", "null_sparse"))) {
  cat(sprintf("\n-- %s --\n", sc))
  for (n in NS) {
    a <- S$power_adj[S$scen == sc & S$n == n & S$arm == "G10"   & S$test == "EDGE.poly3"]
    b <- S$power_adj[S$scen == sc & S$n == n & S$arm == "rule"  & S$test == "EDGE.poly3"]
    sa <- S$size[S$scen == sc & S$n == n & S$arm == "G10"  & S$test == "EDGE.poly3"]
    sb <- S$size[S$scen == sc & S$n == n & S$arm == "rule" & S$test == "EDGE.poly3"]
    if (length(a) && length(b))
      cat(sprintf("   n=%5d  G=10: %.3f (size %.3f)   rule G=%3d: %.3f (size %.3f)   %+.3f\n",
                  n, a, sa, g_of(n), b, sb, b - a))
  }
}
cat("\n=================================================================\n")
cat(" THE QUESTION ASKED: with G by the rule, does EDGE still lose?\n")
cat("=================================================================\n")
for (sc in setdiff(SCEN, c("null", "null_sparse"))) {
  cat(sprintf("\n-- %s (rule arm) --\n", sc))
  for (n in NS) {
    z <- S[S$scen == sc & S$n == n & S$arm == "rule", ]
    if (!nrow(z)) next
    gv <- function(t) {v <- z$power_adj[z$test == t]; if (length(v)) v[1] else NA}
    cat(sprintf("   n=%5d  EDGE3 %.3f | EDGEstk %.3f | GiViTI %.3f | Stukel %.3f | HL %.3f\n",
                n, gv("EDGE.poly3"), gv("EDGE.stk"), gv("GiViTI"), gv("Stukel"), gv("HL")))
  }
}
cat("\n=================================================================\n")
cat(" SIZE under the rule -- the thing that can break (nominal 0.05)\n")
cat("=================================================================\n")
z <- S[S$arm == "rule" & S$test %in% c("EDGE.poly3", "EDGE.stk", "GiViTI", "HL"), ]
z$mcse_out <- round((z$size - .05) / sqrt(.05 * .95 / z$n * 0 + .05 * .95 / 2000), 1)
print(unique(z[, c("scen", "n", "G", "test", "size", "mcse_out")]), row.names = FALSE, digits = 3)
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
