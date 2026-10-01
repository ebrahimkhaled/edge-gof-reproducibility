## run_K2_ao_basis.R -- follow-up to run K: the Aranda-Ordaz (AO) link-family basis.
## Run K found that a single column phi_AO(eta) = 1 - log(1+e^eta)/pi, the score direction of the
## Aranda-Ordaz asymmetric family at the logistic (lambda = 1), beats GiViTI and Stukel on the
## cloglog bow (0.881 vs 0.872 / 0.857 at n = 1000) while being blind elsewhere. Here:
##   (i)  the same column in the PAPER'S construction (unit-weighted projection, eigen null) so the
##        paper can present EDGE-ao as a fourth basis, not a new statistic;
##   (ii) the score form again, for comparison;
##   (iii) n in {500, 1000, 2000}; cloglog, LOG-LOG (the mirror bow, which should load on the same
##        direction with the opposite sign), probit, binint, osc4; matched null per cell.
suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR); library(parallel)})
source("_pstar_giviti_harness.R"); source("_dgp_library.R")
L <- grep("^cl <- makeCluster", readLines("run_K_basis_score.R")) - 1
source(textConnection(paste(readLines("run_K_basis_score.R")[1:L], collapse = "\n")))   # gen, grp_stats, edge_eig, edge_score, ao_col

SEED0 <- 20260919L; B <- 2000L
gen2 <- function(scn, n) {
  if (scn == "loglog") { x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); eta <- 0.6 * x + 0.5 * d
    return(list(d = data.frame(y = rbinom(n, 1, exp(-exp(-eta))), x = x, d = d), f = y ~ x + d)) }
  gen(scn, n)
}
one2 <- function(dat, G) {
  fit <- tryCatch(glm(dat$f, data = dat$d, family = binomial(), control = list(maxit = 50)), error = function(e) NULL)
  nm <- c("EDGE.poly3", "EDGE.stk", "EDGE.ao", "EDGE.ao.score", "EDGE.ao.poly3", "HL", "Stukel", "GiViTI")
  if (is.null(fit)) return(setNames(rep(NA_real_, length(nm)), nm))
  gs <- grp_stats(fit, G); eb <- gs$eb; ao <- ao_col(eb)
  stk <- cbind(eb, eb^2 * (eb >= 0), -eb^2 * (eb < 0))
  safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)     # a singular k=4 system is NA, not a crash
  setNames(c(safe(edge_eig(gs, as.matrix(poly(gs$pb, 3)))), safe(edge_eig(gs, stk)),
             safe(edge_eig(gs, cbind(ao))), safe(edge_score(gs, cbind(ao))),
             safe(edge_eig(gs, cbind(ao, as.matrix(poly(gs$pb, 3))))),
             unname(hl_stat(fit$y, fitted(fit), G)["p"]), stukel_p(fit), giviti_p(fit)), nm)
}
CELLS <- rbind(expand.grid(n = c(500L, 1000L, 2000L), scen = c("null", "cloglog", "loglog"), stringsAsFactors = FALSE),
               data.frame(n = 1000L, scen = c("probit", "binint", "osc4")),
               data.frame(n = 5000L, scen = c("null", "probit")))
CELLS$G <- as.integer(round(CELLS$n / 25))
cl <- makeCluster(max(1L, min(20L, detectCores() - 2L)))
invisible(clusterEvalQ(cl, {suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)}); NULL}))
clusterExport(cl, c("gen", "gen2", "grp_stats", "edge_eig", "edge_score", "ao_col", "one2", "hl_stat", "stukel_p", "giviti_p", "sqb", "sib", "lg"))
t0 <- Sys.time(); P <- list()
for (i in seq_len(nrow(CELLS))) { ce <- CELLS[i, ]; clusterSetRNGStream(cl, SEED0 + i)
  M <- t(parSapply(cl, seq_len(B), function(b, scn, n, G) one2(gen2(scn, n), G), scn = ce$scen, n = ce$n, G = ce$G))
  P[[i]] <- data.frame(scen = ce$scen, n = ce$n, G = ce$G, rep = seq_len(B), M)
  cat(sprintf("  [%d/%d] %-8s n=%5d G=%3d  %s\n", i, nrow(CELLS), ce$scen, ce$n, ce$G, format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); flush.console() }
stopCluster(cl)
pv <- do.call(rbind, P); write.csv(pv, "runK2_ao_basis_pvalues.csv", row.names = FALSE)
TESTS <- setdiff(names(pv), c("scen", "n", "G", "rep")); out <- list()
for (n in unique(CELLS$n)) { nul <- pv[pv$scen == "null" & pv$n == n, ]
  for (sc in setdiff(unique(CELLS$scen[CELLS$n == n]), "null")) { alt <- pv[pv$scen == sc & pv$n == n, ]
    for (t in TESTS) out[[length(out) + 1]] <- data.frame(n = n, G = alt$G[1], scen = sc, test = t,
      size = mean(nul[[t]] <= .05, na.rm = TRUE), power_adj = size_adjusted(alt[[t]], nul[[t]]), na = na_rate(alt[[t]])) } }
S <- do.call(rbind, out); write.csv(S, "runK2_ao_basis_summary.csv", row.names = FALSE)
cat("\n=== size-adjusted power (size) ===\n")
for (sc in c("cloglog", "loglog", "probit", "binint", "osc4")) for (n in sort(unique(S$n[S$scen == sc]))) {
  z <- S[S$scen == sc & S$n == n, ]; cat(sprintf("\n-- %s n=%d G=%d --\n", sc, n, z$G[1]))
  for (t in TESTS) cat(sprintf("  %-14s %.3f (%.3f)\n", t, z$power_adj[z$test == t], z$size[z$test == t])) }
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
