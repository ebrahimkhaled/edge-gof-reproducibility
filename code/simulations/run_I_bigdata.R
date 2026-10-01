## run_I_bigdata.R -- the G lever on REAL data at scale: UCI "Diabetes 130-US hospitals",
## n = 101,766 encounters, outcome = readmission within 30 days (rate 11.2%).
##
## Why this dataset. The whole point of the growing-G result is that it only bites when n is large
## enough that G can be pushed into the hundreds or thousands. GUSTO-I (n=40,830) is already in the
## family; this is 2.5x larger, public, clinical, and has a rare-ish outcome, so the events-per-group
## constraint is live too (11.2% of 50k validation cases ~ 5,600 events).
##
## Design. Random 50/50 split, seed fixed.
##   DEVELOPMENT half: fit a deliberately plain main-effects logistic model, then run the
##   internal-mode tests (fitting correction Omega applied) over a G sweep.
##   VALIDATION half: freeze that model's predictions and run the external-mode tests (Omega = I,
##   plain chi-squared_k) over the same G sweep. This is the frozen-model use case that clinical
##   validation actually is.
## Both halves report the raw calibration curve too, so the SHAPE of any misfit is visible and
## EDGE's per-direction breakdown can be read against it.
##
## Honesty note written in advance: at n ~ 50,000 every test will reject something -- that is the
## large-n over-rejection the paper already names as a limitation. The demonstration here is NOT
## "does it reject" but (a) whether the VERDICT is stable across G for EDGE and unstable for HL,
## (b) whether EDGE keeps returning a sensible statistic at G = 2000 (m = 25), and (c) what shape
## carries the departure.

suppressPackageStartupMessages({library(ebrahim.gof); library(givitiR)})
source("_pstar_giviti_harness.R")     # hl_stat, stukel_p, giviti_p (internal mode)

D <- read.csv("../data_large/diabetic_data.csv", stringsAsFactors = FALSE, na.strings = c("?", ""))
D <- D[D$gender %in% c("Female", "Male"), ]
D$y <- as.integer(D$readmitted == "<30")
age_mid <- c("[0-10)"=5,"[10-20)"=15,"[20-30)"=25,"[30-40)"=35,"[40-50)"=45,"[50-60)"=55,
             "[60-70)"=65,"[70-80)"=75,"[80-90)"=85,"[90-100)"=95)
D$age_n <- unname(age_mid[D$age])
D$female <- as.integer(D$gender == "Female")
D$insulin <- factor(D$insulin, levels = c("No", "Steady", "Up", "Down"))
D$change <- as.integer(D$change == "Ch"); D$diabetesMed <- as.integer(D$diabetesMed == "Yes")
f <- y ~ age_n + female + time_in_hospital + num_lab_procedures + num_procedures + num_medications +
         number_outpatient + number_emergency + number_inpatient + number_diagnoses + insulin +
         change + diabetesMed
D <- D[stats::complete.cases(D[, all.vars(f)]), ]
set.seed(20260917)
dev <- sample(nrow(D), floor(nrow(D) / 2))
Ddev <- D[dev, ]; Dval <- D[-dev, ]
cat(sprintf("n total %d | development %d | validation %d | event rate dev %.3f val %.3f\n",
            nrow(D), nrow(Ddev), nrow(Dval), mean(Ddev$y), mean(Dval$y)))

fit <- stats::glm(f, data = Ddev, family = stats::binomial())
cat(sprintf("model: %d coefficients, AUC-ish check: mean fitted %.4f vs mean y %.4f\n",
            length(coef(fit)), mean(fitted(fit)), mean(Ddev$y)))
pval <- as.numeric(stats::predict(fit, newdata = Dval, type = "response"))

## ---- external-mode tests (frozen predictions), as in run_H ----------------------------------
grp <- function(p, G) { br <- stats::quantile(p, seq(0, 1, length.out = G + 1), type = 1)
  br[1] <- -Inf; br[length(br)] <- Inf; as.integer(cut(p, unique(br), include.lowest = TRUE)) }
resid_ext <- function(y, p, G) { g <- grp(p, G)
  o <- tapply(y, g, sum); e <- tapply(p, g, sum); v <- tapply(p * (1 - p), g, sum)
  list(r = as.numeric((o - e) / sqrt(v)), eta = as.numeric(stats::qlogis(tapply(p, g, mean))), K = max(g)) }
edge_ext <- function(y, p, G, basis) { z <- resid_ext(y, p, G); eta <- z$eta; r <- z$r
  Z <- switch(basis, slope = cbind(1, eta), poly3 = cbind(1, eta, eta^2, eta^3),
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

write.csv(rbind(do.call(rbind, dev_rows)), "runI_bigdata_dev.csv", row.names = FALSE)
write.csv(do.call(rbind, val_rows), "runI_bigdata_val.csv", row.names = FALSE)
write.csv(cc, "runI_bigdata_calcurve.csv", row.names = FALSE)
cat("\nwritten: runI_bigdata_{dev,val,calcurve}.csv\n")
