## audit_fair_analyse.R -- fairness audit: paired comparisons for the claims (reads audit_fairness_pairs.csv
## from audit_fairness_summary.R and the normal-x files from audit_fair_design.R). Prints only.
P <- read.csv("audit_fairness_pairs.csv", stringsAsFactors = FALSE)
options(width = 220)
sel <- function(tags, probes, rivals, link = NULL) {
  z <- P[P$tag %in% tags & P$probe %in% probes & P$rival %in% rivals, ]
  if (!is.null(link)) z <- z[z$link %in% link, ]
  z$p_mcnemar <- signif(z$p_mcnemar, 2)
  z[, c("tag", "link", "n", "probe", "rival", "pow_probe", "pow_rival", "diff", "se", "p_mcnemar")]
}
WIN <- c("base_probit", "base_cauchit", "base_cauchit_m10", "s2_cauchit", "s2_t4", "s2_t4_m10", "s2_probit", "s2_probit_m25sc",
         "s2_probit_m10", "low_probit", "low_probit_m10")
cat("=== C3: symmetric probe (unit and score form) vs each rival ===\n")
print(sel(WIN, c("g.sym", "g.sym.sc"), c("Stk.joint", "Stk.LR", "GiViTI", "EDGE.poly3")), digits = 3, row.names = FALSE)
cat("\n=== C5: bow -- matched probe and oracle vs GiViTI and Stk.LR ===\n")
ASY <- c("base_asym", "base_asym_m10", "base_asym_m5")
print(rbind(sel(ASY, c("g.orc", "u.orc", "g.ao", "g.ao.sc", "g.max3.sc"), c("GiViTI", "Stk.LR"), "cloglog"),
            sel(ASY, c("g.orc", "u.orc", "g.aoM", "g.aoM.sc", "g.max3.sc"), c("GiViTI", "Stk.LR"), "loglog")), digits = 3, row.names = FALSE)
cat("\n=== C7: the general-purpose max (score form, run L2; unit form, run L1) vs every rival ===\n")
L2 <- c("base_asym_m10", "base_asym_m5", "base_cauchit_m10", "low_probit_m10", "s2_probit_m25sc", "s2_probit_m10", "s2_t4_m10")
z <- sel(L2, "g.max3.sc", c("GiViTI", "Stk.LR", "Stk.joint", "EDGE.poly3")); print(z, digits = 3, row.names = FALSE)
cat(sprintf("  g.max3.sc - GiViTI range %+.3f .. %+.3f ; - Stk.LR range %+.3f .. %+.3f ; - best rival (min over cells) %+.3f\n",
  min(z$diff[z$rival == "GiViTI"]), max(z$diff[z$rival == "GiViTI"]), min(z$diff[z$rival == "Stk.LR"]), max(z$diff[z$rival == "Stk.LR"]),
  min(tapply(seq_len(nrow(z)), paste(z$tag, z$link, z$n), function(i) z$pow_probe[i[1]] - max(z$pow_rival[i])))))
L1 <- setdiff(unique(P$tag), L2)
z1 <- sel(L1, "g.max3", c("GiViTI", "Stk.LR", "Stk.joint", "EDGE.poly3"))
cat(sprintf("  unit g.max3 (run L1) - GiViTI range %+.3f .. %+.3f ; - Stk.LR range %+.3f .. %+.3f\n",
  min(z1$diff[z1$rival == "GiViTI"]), max(z1$diff[z1$rival == "GiViTI"]), min(z1$diff[z1$rival == "Stk.LR"]), max(z1$diff[z1$rival == "Stk.LR"])))

cat("\n=== AUDIT DESIGN CHECK: x ~ N(0, 1.5^2), paired, each test at its matched-null 5% quantile ===\n")
crit <- function(p) { p <- p[is.finite(p)]; as.numeric(quantile(p, .05, type = 1)) }
for (f in list.files(pattern = "^audit_fair_design_.*_pvalues\\.csv$")) {
  pv <- read.csv(f, check.names = FALSE)
  for (n in unique(pv$n)) { nul <- pv[pv$link == "logit" & pv$n == n, ]
    for (lk in setdiff(unique(pv$link), "logit")) { alt <- pv[pv$link == lk & pv$n == n, ]
      TT <- c("g.sym.sc", "g.sym", "u.sym", "g.max3.sc", "g.max3", "EDGE.poly3", paste0("g.orc.", lk), "Stk.joint", "Stk.LR", "GiViTI", "HL10")
      rej <- lapply(setNames(TT, TT), function(t) { a <- alt[[t]]; ifelse(is.finite(a), a <= crit(nul[[t]]), FALSE) })
      cat(sprintf("\n-- %s n=%d  B=%d --\n", lk, n, nrow(alt)))
      for (t in TT) cat(sprintf("  %-14s size %.3f  power(adj) %.3f  raw %.3f  na %.3f\n", t, mean(nul[[t]] <= .05, na.rm = TRUE),
                                mean(rej[[t]]), mean(alt[[t]] <= .05, na.rm = TRUE), mean(!is.finite(alt[[t]]))))
      for (a in c("g.sym.sc", "g.sym", "u.sym", "g.max3.sc", "EDGE.poly3")) for (b in c("Stk.joint", "Stk.LR", "GiViTI")) {
        A <- rej[[a]]; Bv <- rej[[b]]; n10 <- sum(A & !Bv); n01 <- sum(!A & Bv); d <- mean(A) - mean(Bv)
        se <- sqrt((n10 + n01) / length(A) - d^2) / sqrt(length(A))
        cat(sprintf("  %-10s vs %-9s diff %+.3f (se %.3f)  McNemar p %s\n", a, b, d, se,
                    ifelse(n10 + n01 > 0, format.pval(binom.test(n10, n10 + n01)$p.value, digits = 2), "NA")))
      }
    }
  }
}
