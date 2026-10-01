## audit_numbers_raw.R -- (1) the claim families again at the NOMINAL cut-off p <= .05 (no estimated critical value,
## so McNemar is exact); (2) between-run spread of the same cell across independent run L / L2 replications
## (different seeds) against the binomial MC SE; (3) GiViTI mechanism detail for C4.
C <- read.csv("audit_numbers_paired_full.csv", stringsAsFactors = FALSE)
files <- setdiff(list.files(pattern = "^runL_.*_pvalues\\.csv$"), "runL_validate_base_pvalues.csv")
PL <- lapply(setNames(files, sub("^runL_(.*)_pvalues\\.csv$", "\\1", files)), function(f) read.csv(f, check.names = FALSE))
oc <- function(t, lk) if (grepl("orc$", t)) paste0(t, ".", lk) else t
raw <- function(v) is.finite(v) & v <= .05
out <- list()
for (i in seq_len(nrow(C))) { r <- C[i, ]; P <- PL[[r$tag]]; alt <- P[P$link == r$link & P$n == r$n, ]
  A <- raw(alt[[oc(r$probe, r$link)]]); B <- raw(alt[[oc(r$rival, r$link)]]); n10 <- sum(A & !B); n01 <- sum(!A & B)
  out[[i]] <- data.frame(raw_probe = mean(A), raw_rival = mean(B), raw_diff = mean(A) - mean(B), raw_n10 = n10, raw_n01 = n01,
    raw_z = (n10 - n01) / sqrt(max(1, (n10 + n01) - (n10 - n01)^2 / length(A))),
    raw_p = if (n10 + n01 > 0) binom.test(n10, n10 + n01)$p.value else NA_real_) }
C <- cbind(C, do.call(rbind, out))
ok <- is.finite(C$raw_p); C$raw_holm_all <- NA; C$raw_holm_all[ok] <- p.adjust(C$raw_p[ok], "holm")
C$raw_bh_all <- NA; C$raw_bh_all[ok] <- p.adjust(C$raw_p[ok], "BH")
cat(sprintf("=== nominal-cutoff McNemar over all %d pairs: raw p<.05 %d, Holm %d, BH %d ===\n", nrow(C), sum(C$raw_p < .05, na.rm = TRUE),
  sum(C$raw_holm_all < .05, na.rm = TRUE), sum(C$raw_bh_all < .05, na.rm = TRUE)))
C$raw_holm_fam <- NA
for (f in setdiff(unique(C$fam), c("", NA))) { k <- which(C$fam == f & is.finite(C$raw_p)); C$raw_holm_fam[k] <- p.adjust(C$raw_p[k], "holm")
  z <- C[C$fam == f, ]; z <- z[order(z$link, z$tag, z$n, z$probe, z$rival), ]
  cat(sprintf("\n--- family %s (nominal cut-off): p<.05 %d of %d; Holm-family %d; Holm-global %d ---\n", f, sum(z$raw_p < .05, na.rm = TRUE), nrow(z),
    sum(z$raw_holm_fam < .05, na.rm = TRUE), sum(z$raw_holm_all < .05, na.rm = TRUE)))
  for (j in seq_len(nrow(z))) with(z[j, ], cat(sprintf("  %-17s %-7s %-6d %-9s %-7s adj %+.3f (z_boot %+.1f) | raw %.3f vs %.3f diff %+.3f n10 %3d n01 %3d z %+5.1f p %.1e Holm_fam %.1e Holm_all %.1e\n",
    tag, link, n, probe, rival, diff, z_boot, raw_probe, raw_rival, raw_diff, raw_n10, raw_n01, raw_z, raw_p, raw_holm_fam, raw_holm_all))) }
write.csv(C, "audit_numbers_paired_full.csv", row.names = FALSE)

## (2) between-run spread for the same design cell
crit1 <- function(p) { p <- p[is.finite(p)]; as.numeric(quantile(p, .05, type = 1)) }
groups <- list(list(tags = c("base_asym", "base_asym_m10", "base_asym_m5"), links = c("cloglog", "loglog"), ns = c(600, 1000)),
               list(tags = c("base_cauchit", "base_cauchit_m10"), links = "cauchit", ns = 2500),
               list(tags = c("low_probit", "low_probit_m10"), links = "probit", ns = 5000),
               list(tags = c("s2_probit", "s2_probit_m25sc", "s2_probit_m10"), links = "probit", ns = c(4000, 6500)),
               list(tags = c("s2_t4", "s2_t4_m10"), links = "t4", ns = 6000))
cat("\n=== between-run spread of G-independent tests (independent seeds): adjusted and raw power per run; chi-square heterogeneity with binomial variance ===\n")
het <- list()
for (g in groups) for (lk in g$links) for (n in g$ns) for (t in c("GiViTI", "Stk.LR", "Stk.joint", "u.orc", "u.sym", "HL10")) {
  adj <- sapply(g$tags, function(tg) { P <- PL[[tg]]; nul <- P[P$link == "logit" & P$n == n, ]; alt <- P[P$link == lk & P$n == n, ]
    cn <- oc(t, lk); mean(is.finite(alt[[cn]]) & alt[[cn]] <= crit1(nul[[cn]])) })
  rw <- sapply(g$tags, function(tg) { P <- PL[[tg]]; alt <- P[P$link == lk & P$n == n, ]; mean(raw(alt[[oc(t, lk)]])) })
  B <- sapply(g$tags, function(tg) sum(PL[[tg]]$link == lk & PL[[tg]]$n == n))
  hp <- function(x) { pb <- sum(x * B) / sum(B); X2 <- sum((x - pb)^2 / (pb * (1 - pb) / B)); pchisq(X2, length(x) - 1, lower.tail = FALSE) }
  het[[length(het) + 1]] <- data.frame(link = lk, n = n, test = t, runs = length(adj), adj = paste(sprintf("%.3f", adj), collapse = "/"),
    adj_range = diff(range(adj)), adj_het_p = hp(adj), raw = paste(sprintf("%.3f", rw), collapse = "/"), raw_range = diff(range(rw)),
    raw_het_p = hp(rw), binom_se = sqrt(mean(adj) * (1 - mean(adj)) / B[1])) }
H <- do.call(rbind, het); print(H, row.names = FALSE, digits = 3)
cat(sprintf("\n  adjusted-power heterogeneity p < .05 in %d of %d rows; raw-power heterogeneity p < .05 in %d of %d rows\n",
  sum(H$adj_het_p < .05), nrow(H), sum(H$raw_het_p < .05), nrow(H)))
cat(sprintf("  KS of heterogeneity p-values vs U(0,1): adjusted p = %.3g, raw p = %.3g\n", ks.test(H$adj_het_p, "punif")$p.value, ks.test(H$raw_het_p, "punif")$p.value))

## (3) GiViTI mechanism: what does a fixed degree-3 LR do on the samples that pass the gate?
M <- read.csv("giviti_mechanism_cauchit.csv")
k3 <- M$m >= 3
cat(sprintf("\n=== C4 detail (cauchit n=2500, B=%d) ===\n  gate passed (m>=3): %d; GiViTI rejects %d (%.3f); fixed deg-3 LR (2 df) rejects %d (%.3f); LR3 alone > 3.84 by construction\n",
  nrow(M), sum(k3), sum(M$p[k3] <= .05), mean(M$p[k3] <= .05), sum(M$lr2[k3] + M$lr3[k3] > qchisq(.95, 2)), mean(M$lr2[k3] + M$lr3[k3] > qchisq(.95, 2))))
cat(sprintf("  gate failed (m=2): %d; GiViTI rejects %.3f; fixed deg-3 LR rejects %.3f; deg-2 LR (1 df) rejects %.3f\n",
  sum(!k3), mean(M$p[!k3] <= .05), mean(M$lr2[!k3] + M$lr3[!k3] > qchisq(.95, 2)), mean(M$lr2[!k3] > qchisq(.95, 1))))
cat(sprintf("  GiViTI statistic among m>=3: mean %.2f; deviance gain lr2+lr3 mean %.2f\n", mean(M$stat[k3]), mean(M$lr2[k3] + M$lr3[k3])))
