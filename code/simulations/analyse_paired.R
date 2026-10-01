## analyse_paired.R -- is the EDGE vs GiViTI / Stukel difference statistically real, once G is
## chosen well?
##
## Every test in this harness is computed on THE SAME fitted model within a replication, so the
## rejections are PAIRED. Comparing two rejection rates as if they were independent proportions
## throws that away and inflates the standard error. The right instrument is McNemar's test on the
## discordant pairs, which is what is used here, with an exact binomial CI on the difference.
##
## Every comparison is made at each test's own SIZE-ADJUSTED critical value, taken from the matched
## null run in the same cell, so no test is credited with power it bought by rejecting too often.

alpha <- 0.05

crit <- function(pnull, a = alpha) {
  p <- pnull[is.finite(pnull)]
  if (length(p) < 50) return(NA_real_)
  as.numeric(stats::quantile(p, probs = a, na.rm = TRUE, type = 1))
}

## McNemar on paired rejections + an exact CI for the difference in proportions
paired_cmp <- function(pa_alt, pb_alt, pa_null, pb_null, lab_a, lab_b, cell) {
  ca <- crit(pa_null); cb <- crit(pb_null)
  ok <- is.finite(pa_alt) & is.finite(pb_alt)
  A <- pa_alt[ok] <= ca; B <- pb_alt[ok] <= cb
  n <- sum(ok)
  b <- sum(A & !B)      # A rejects, B does not
  c <- sum(!A & B)      # B rejects, A does not
  d <- b + c            # discordant pairs -- all the information there is
  est <- mean(A) - mean(B)
  if (d == 0) return(data.frame(cell = cell, a = lab_a, b = lab_b, n = n,
                                pow_a = mean(A), pow_b = mean(B), diff = est,
                                discordant = 0, p = NA_real_, lo = NA_real_, hi = NA_real_,
                                verdict = "no discordant pairs"))
  ht <- stats::binom.test(b, d, p = 0.5)                  # exact McNemar
  ci <- ht$conf.int * d / n                               # back onto the difference scale
  lo <- (2 * ci[1] - d / n); hi <- (2 * ci[2] - d / n)
  data.frame(cell = cell, a = lab_a, b = lab_b, n = n,
             pow_a = mean(A), pow_b = mean(B), diff = est, discordant = d,
             p = ht$p.value, lo = lo, hi = hi,
             verdict = ifelse(ht$p.value < 0.05,
                              ifelse(est > 0, "A better", "B better"), "not significant"))
}

P <- read.csv("runF_grule_pvalues.csv", stringsAsFactors = FALSE)
P$m_actual <- round(P$n / P$G)

out <- list()
for (sc in c("cloglog", "probit")) {
  for (nn in sort(unique(P$n))) {
    cells <- unique(P[P$n == nn, c("m", "G")])
    for (i in seq_len(nrow(cells))) {
      G <- cells$G[i]
      nul <- P[P$scen == "null" & P$n == nn & P$G == G, ]
      alt <- P[P$scen == sc & P$n == nn & P$G == G, ]
      if (!nrow(nul) || !nrow(alt)) next
      lab <- sprintf("%s n=%d G=%d (m=%d)", sc, nn, G, round(nn / G))
      for (edge in c("EDGE.poly3", "EDGE.stk")) {
        for (rival in c("GiViTI", "Stukel")) {
          out[[length(out) + 1]] <- paired_cmp(alt[[edge]], alt[[rival]],
                                               nul[[edge]], nul[[rival]],
                                               edge, rival, lab)
        }
      }
    }
  }
}
R <- do.call(rbind, out)
R <- R[is.finite(R$diff), ]
write.csv(R, "paired_edge_vs_rivals.csv", row.names = FALSE)

cat("============================================================================\n")
cat(" AT EDGE'S BEST PARTITION: is the remaining gap statistically real?\n")
cat(" Paired McNemar, size-adjusted, same samples. n = replications.\n")
cat("============================================================================\n")
best <- c("cloglog n=500 G=33 (m=15)", "cloglog n=1000 G=67 (m=15)",
          "cloglog n=2000 G=133 (m=15)", "probit n=5000 G=333 (m=15)",
          "probit n=20000 G=400 (m=50)")
z <- R[R$cell %in% best & R$a == "EDGE.poly3", ]
for (i in seq_len(nrow(z))) with(z[i, ],
  cat(sprintf("  %-30s EDGE %.3f vs %-7s %.3f   diff %+.3f  [%+.3f, %+.3f]  p=%s  %s\n",
      cell, pow_a, b, pow_b, diff, lo, hi,
      ifelse(p < 1e-4, "<1e-4", sprintf("%.4f", p)), verdict)))

cat("\n============================================================================\n")
cat(" WHERE IS THE DIFFERENCE NOT SIGNIFICANT? (every unsaturated cell)\n")
cat("============================================================================\n")
u <- R[R$pow_a < 0.995 & R$pow_b < 0.995 & R$a == "EDGE.poly3", ]
ns <- u[u$verdict == "not significant", ]
cat(sprintf("  %d of %d unsaturated comparisons are NOT significant at 0.05\n", nrow(ns), nrow(u)))
if (nrow(ns)) for (i in seq_len(nrow(ns))) with(ns[i, ],
  cat(sprintf("    %-30s vs %-7s  diff %+.3f  [%+.3f, %+.3f]  p=%.3f\n", cell, b, diff, lo, hi, p)))

cat("\n============================================================================\n")
cat(" DOES A FINER PARTITION CLOSE THE GAP? gap at the conventional G=10 vs at EDGE's best\n")
cat("============================================================================\n")
for (sc in c("cloglog", "probit")) for (nn in c(1000, 2000, 5000, 20000)) {
  g10 <- R[grepl(sprintf("^%s n=%d G=10 ", sc, nn), R$cell) & R$a == "EDGE.poly3", ]
  fine <- R[grepl(sprintf("^%s n=%d ", sc, nn), R$cell) & R$a == "EDGE.poly3", ]
  if (!nrow(g10) || !nrow(fine)) next
  fine <- fine[fine$m_ok <- TRUE, ]
  bestrow <- do.call(rbind, lapply(split(fine, fine$b), function(d) d[which.max(d$pow_a), ]))
  for (rv in unique(g10$b)) {
    a <- g10[g10$b == rv, ]; bb <- bestrow[bestrow$b == rv, ]
    if (!nrow(a) || !nrow(bb)) next
    cat(sprintf("  %-8s n=%5d vs %-7s : G=10 diff %+.3f (p=%.3g)  ->  %s diff %+.3f (p=%.3g)\n",
        sc, nn, rv, a$diff, a$p, sub(".*G=", "G=", bb$cell), bb$diff, bb$p))
  }
}
cat("\nwritten: paired_edge_vs_rivals.csv\n")
