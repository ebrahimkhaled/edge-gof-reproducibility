## run_I_bigdata_p3.R -- Section 7 of EDGE paper 3: run_I_bigdata.R on the corrected cohort.
##
## Two changes from run_I_bigdata.R, both from the paper-3 referee round, and nothing else:
##  1. the cohort is _cohort_p3.R's: one admission per patient, no death or hospice discharges, split by patient
##     (paper 2 split admissions, so patients recurred across the halves);
##  2. the external-mode statistic uses the paper's DEFAULT basis, a cubic polynomial in the group's mean risk,
##     plus the intercept column that external validation needs (4 df). Paper 2 used powers of the group logit,
##     a different space from the one the belt figure draws.
## The Stukel and slope bases, the Hosmer-Lemeshow statistic, the Cox/Miller test and GiViTI are as before.
##
##   Rscript run_I_bigdata_p3.R   -> runI_p3_{dev,val,calcurve}.csv
## ---- archive paths (inserted by make_archive.py; the convention is in README.md, "How to run") -------------
EDGE_ARCHIVE_ROOT <- local({
  r <- Sys.getenv("EDGE_ARCHIVE_ROOT")
  if (!nzchar(r)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    r <- if (length(f)) file.path(dirname(normalizePath(f[1], winslash = "/")), "..", "..") else getwd()
  }
  r <- normalizePath(r, winslash = "/", mustWork = FALSE)
  Sys.setenv(EDGE_ARCHIVE_ROOT = r)        # so worker processes started from here resolve the same root
  r
})
edge_path <- function(...) file.path(EDGE_ARCHIVE_ROOT, ...)
edge_battery <- function(...) {           # the author's battery/ folder: analysis/ -> results/analysis, rest -> results/blocks
  if (...length() == 0L) return(edge_path("results", "blocks"))
  p <- file.path(...)
  ifelse(p == "analysis" | startsWith(p, "analysis/"), edge_path("results", p), edge_path("results", "blocks", p))
}
edge_out <- function(...) { d <- edge_path("output", ...); dir.create(d, showWarnings = FALSE, recursive = TRUE); d }
## ---------------------------------------------------------------------------------------------------------------

suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)})
source("_pstar_giviti_harness.R")     # hl_stat, stukel_p, giviti_p (internal mode)

source("_cohort_p3.R")
CO <- cohort_p3(getwd())
Ddev <- CO$Ddev; Dval <- CO$Dval; f <- CO$f
cat(sprintf("admissions in file %d | cohort (one per patient, no death/hospice) %d | development %d | validation %d | event rate dev %.3f val %.3f
",
            CO$n_raw, CO$n, nrow(Ddev), nrow(Dval), mean(Ddev$y), mean(Dval$y)))

fit <- stats::glm(f, data = Ddev, family = stats::binomial())
cat(sprintf("model: %d coefficients, AUC-ish check: mean fitted %.4f vs mean y %.4f\n",
            length(coef(fit)), mean(fitted(fit)), mean(Ddev$y)))
pval <- as.numeric(stats::predict(fit, newdata = Dval, type = "response"))

## ---- external-mode tests (frozen predictions), as in run_H ----------------------------------
## equal-frequency groups by rank, the rule of def.gof()/edge.gof() in ebrahim.gof 2.9.0, so that
## edge.gof(y = y, predicted_probs = p, external = TRUE) reproduces these numbers
grp <- function(p, G) as.integer(pmin(ceiling(rank(p, ties.method = "first") / (length(p) / G)), G))
resid_ext <- function(y, p, G) { g <- grp(p, G)
  o <- tapply(y, g, sum); e <- tapply(p, g, sum); v <- tapply(p * (1 - p), g, sum)
  list(r = as.numeric((o - e) / sqrt(v)), eta = as.numeric(stats::qlogis(tapply(p, g, mean))), K = max(g)) }
edge_ext <- function(y, p, G, basis) { z <- resid_ext(y, p, G); eta <- z$eta; r <- z$r
  pb <- stats::plogis(eta)
  Z <- switch(basis, slope = cbind(1, eta), poly3 = cbind(1, pb, pb^2, pb^3),
              stukel = cbind(1, eta, eta^2 * (eta >= 0), -eta^2 * (eta < 0)))
  Q <- qr(Z); Z <- Z[, Q$pivot[seq_len(Q$rank)], drop = FALSE]
  S <- sum(qr.fitted(qr(Z), r)^2)
  c(stat = S, df = ncol(Z), p = stats::pchisq(S, ncol(Z), lower.tail = FALSE)) }
hl_ext <- function(y, p, G) { z <- resid_ext(y, p, G); S <- sum(z$r^2)
  c(stat = S, df = z$K, p = stats::pchisq(S, z$K, lower.tail = FALSE)) }
slope_lrt <- function(y, p) { lp <- stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  f1 <- stats::glm(y ~ lp, family = stats::binomial())
  ll0 <- sum(y * log(p) + (1 - y) * log(1 - p)); S <- 2 * (as.numeric(stats::logLik(f1)) - ll0)
  c(stat = S, df = 2, p = stats::pchisq(S, 2, lower.tail = FALSE), a = coef(f1)[1], b = coef(f1)[2]) }

GS <- c(10L, 20L, 50L, 100L, 200L, 500L, 1000L, 2000L)
cat("\n=================================================================\n")
cat(" DEVELOPMENT HALF, internal mode (fitting correction applied)\n")
cat(sprintf(" n=%d. m = n/G shown. p-values.\n", nrow(Ddev)))
cat("=================================================================\n")
cat(sprintf("  %5s %6s | %-9s %-9s %-9s | %-9s %-9s\n", "G", "m", "EDGE.poly3", "EDGE.stk", "EDGE.poly2", "HL", "HL_F"))
dev_rows <- list()
for (G in GS) {
  e3 <- tryCatch(edge.gof(fit, G = G, basis = "poly3")$p_value, error = function(e) NA)
  es <- tryCatch(edge.gof(fit, G = G, basis = "stukel")$p_value, error = function(e) NA)
  e2 <- tryCatch(edge.gof(fit, G = G, basis = "poly2")$p_value, error = function(e) NA)
  hl <- unname(hl_stat(Ddev$y, fitted(fit), G)["p"])
  ef <- tryCatch(ef.gof(fit, G = G)$p_value, error = function(e) NA)
  dev_rows[[length(dev_rows) + 1]] <- data.frame(half = "dev", G = G, m = round(nrow(Ddev) / G),
    EDGE.poly3 = e3, EDGE.stk = es, EDGE.poly2 = e2, HL = hl, HL_F = ef)
  cat(sprintf("  %5d %6d | %-9.2e %-9.2e %-9.2e | %-9.2e %-9.2e\n", G, round(nrow(Ddev) / G), e3, es, e2, hl, ef))
}
cat(sprintf("  once: Stukel p=%.2e | GiViTI p=%.2e\n", stukel_p(fit), giviti_p(fit)))

cat("\n=================================================================\n")
cat(" VALIDATION HALF, external mode (frozen predictions, Omega = I)\n")
cat(sprintf(" n=%d. Statistic (df) and p-value.\n", nrow(Dval)))
cat("=================================================================\n")
cat(sprintf("  %5s %6s | %-16s %-16s %-16s | %-18s\n", "G", "m", "EDGE.stk4", "EDGE.poly4", "EDGE.slope2", "HL.ext"))
val_rows <- list()
for (G in GS) {
  a <- edge_ext(Dval$y, pval, G, "stukel"); b <- edge_ext(Dval$y, pval, G, "poly3")
  s <- edge_ext(Dval$y, pval, G, "slope"); h <- hl_ext(Dval$y, pval, G)
  val_rows[[length(val_rows) + 1]] <- data.frame(half = "val", G = G, m = round(nrow(Dval) / G),
    stk_stat = a["stat"], stk_p = a["p"], poly_stat = b["stat"], poly_p = b["p"],
    slope_stat = s["stat"], slope_p = s["p"], hl_stat = h["stat"], hl_df = h["df"], hl_p = h["p"])
  cat(sprintf("  %5d %6d | %6.1f(%d) %-7.1e %6.1f(%d) %-7.1e %6.1f(%d) %-7.1e | %8.1f(%d) %-7.1e\n",
              G, round(nrow(Dval) / G), a["stat"], a["df"], a["p"], b["stat"], b["df"], b["p"],
              s["stat"], s["df"], s["p"], h["stat"], h["df"], h["p"]))
}
sl <- slope_lrt(Dval$y, pval)
cat(sprintf("  once: Cox/Miller slope LRT stat %.1f (2) p=%.2e  intercept a=%.3f slope b=%.3f\n",
            sl["stat"], sl["p"], sl["a"], sl["b"]))
gv <- tryCatch(givitiCalibrationTest(Dval$y, pval, devel = "external")$p.value, error = function(e) NA)
cat(sprintf("  once: GiViTI external p=%.2e\n", gv))

cat("\n=================================================================\n")
cat(" THE SHAPE: validation-half calibration curve, deciles of predicted risk\n")
cat("=================================================================\n")
z <- resid_ext(Dval$y, pval, 10); g <- grp(pval, 10)
cc <- data.frame(decile = 1:10, mean_pred = tapply(pval, g, mean), obs_rate = tapply(Dval$y, g, mean),
                 n = tapply(Dval$y, g, length), r_g = z$r)
cc$gap_pts <- 100 * (cc$obs_rate - cc$mean_pred)
print(cc, row.names = FALSE, digits = 4)
cat(sprintf("\n  gap sign pattern (obs - pred): %s\n", paste(ifelse(cc$gap_pts > 0, "+", "-"), collapse = " ")))
cat("  (an S / overconfidence reads -- at the low end and + at the high end, or the reverse)\n")

cat("\n=================================================================\n")
cat(" EDGE's diagnosis on the validation half, G=100: which direction carries it?\n")
cat("=================================================================\n")
zz <- resid_ext(Dval$y, pval, 100); eta <- zz$eta; r <- zz$r
Z <- cbind(c1 = 1, eta = eta, sq_pos = eta^2 * (eta >= 0), sq_neg = -eta^2 * (eta < 0))
Q <- qr(Z); R <- qr.R(Q); Qm <- qr.Q(Q)
comp <- as.numeric(crossprod(Qm, r))^2                 # squared projection onto successive orthonormal dirs
names(comp) <- colnames(Z)[Q$pivot[seq_len(Q$rank)]]
cat("  squared length of r along each successive orthogonalised basis direction:\n")
for (nm in names(comp)) cat(sprintf("    %-8s %8.2f  (%4.1f%% of the stukel-basis statistic)\n", nm, comp[nm], 100 * comp[nm] / sum(comp)))
cat(sprintf("  total ||r||^2 (HL, G=100) = %.1f over %d df; EDGE-stk4 captures %.1f of it in 4 df.\n",
            sum(r^2), length(r), sum(comp)))

write.csv(rbind(do.call(rbind, dev_rows)), edge_path("results/cohort/runI_p3_dev.csv"), row.names = FALSE)
write.csv(do.call(rbind, val_rows), edge_path("results/cohort/runI_p3_val.csv"), row.names = FALSE)
write.csv(cc, edge_path("results/cohort/runI_p3_calcurve.csv"), row.names = FALSE)
cat("\nwritten: runI_p3_{dev,val,calcurve}.csv\n")
