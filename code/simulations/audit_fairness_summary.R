## audit_fairness_summary.R -- fairness audit (independent recomputation, reads run L p-value files only)
## For every run L cell: raw size, raw power, size-adjusted power (NA = no rejection, crit from the
## matched null), NA rates, and paired differences for the probes against every rival, with McNemar p.
## Writes audit_fairness_cells.csv and audit_fairness_pairs.csv; prints the key tables.
files <- setdiff(list.files(pattern = "^runL_.*_pvalues\\.csv$"), "runL_validate_base_pvalues.csv")
PL <- lapply(setNames(files, sub("^runL_(.*)_pvalues\\.csv$", "\\1", files)), function(f) read.csv(f, check.names = FALSE))
crit <- function(p) { p <- p[is.finite(p)]; if (length(p) < 50) NA_real_ else as.numeric(quantile(p, .05, type = 1)) }
TESTS <- c("g.orc", "g.sym", "g.ao", "g.aoM", "g.max3", "g.sym.sc", "g.ao.sc", "g.aoM.sc", "g.max3.sc",
           "EDGE.poly3", "EDGE.stk", "EDGE.stk.sc", "u.orc", "u.sym", "u.ao", "u.aoM", "u.cub", "u.max3",
           "Stk.joint", "Stk.marg", "Stk.LR", "GiViTI", "HL10", "HLG")
cells <- list(); pairs <- list(); REJ <- list()
for (tg in names(PL)) for (n in unique(PL[[tg]]$n)) {
  Q <- PL[[tg]][PL[[tg]]$n == n, ]; nul <- Q[Q$link == "logit", ]
  for (lk in setdiff(unique(Q$link), "logit")) {
    alt <- Q[Q$link == lk, ]; rej <- list()
    for (t in TESTS) {
      cn <- if (grepl("orc$", t)) paste0(t, ".", lk) else t
      if (!cn %in% names(alt)) next
      cr <- crit(nul[[cn]]); a <- alt[[cn]]
      rej[[t]] <- ifelse(is.finite(a), a <= cr, FALSE)
      cells[[length(cells) + 1]] <- data.frame(tag = tg, link = lk, s = Q$s[1], c0 = Q$c0[1], n = n, B = nrow(alt), test = t,
        size_raw = mean(nul[[cn]] <= .05, na.rm = TRUE), na_null = mean(!is.finite(nul[[cn]])),
        pow_raw = mean(is.finite(a) & a <= .05), pow_adj = mean(rej[[t]]),
        pow_adj_finite = mean(a[is.finite(a)] <= cr), na_alt = mean(!is.finite(a)), crit = cr)
    }
    key <- paste(tg, lk, n); REJ[[key]] <- rej
    PR <- intersect(c("g.sym", "g.sym.sc", "g.ao", "g.ao.sc", "g.aoM", "g.aoM.sc", "g.max3", "g.max3.sc", "g.orc",
                      "u.sym", "u.ao", "u.aoM", "u.max3", "u.orc"), names(rej))
    RV <- intersect(c("GiViTI", "Stk.LR", "Stk.joint", "EDGE.poly3", "u.sym", "u.ao", "u.aoM", "u.max3"), names(rej))
    for (a in PR) for (b in RV) { if (a == b) next
      A <- rej[[a]]; Bv <- rej[[b]]; n10 <- sum(A & !Bv); n01 <- sum(!A & Bv)
      d <- mean(A) - mean(Bv); se <- sqrt((n10 + n01) / length(A) - d^2) / sqrt(length(A))
      p <- if (n10 + n01 > 0) binom.test(n10, n10 + n01)$p.value else NA_real_
      pairs[[length(pairs) + 1]] <- data.frame(tag = tg, link = lk, s = Q$s[1], c0 = Q$c0[1], n = n, probe = a, rival = b,
        pow_probe = mean(A), pow_rival = mean(Bv), diff = d, se = se, p_mcnemar = p)
    }
  }
}
C <- do.call(rbind, cells); P <- do.call(rbind, pairs)
write.csv(C, "audit_fairness_cells.csv", row.names = FALSE); write.csv(P, "audit_fairness_pairs.csv", row.names = FALSE)
options(width = 250)
cat("=== raw size of the headline tests in every cell (B nulls per cell) ===\n")
W <- reshape(C[C$test %in% c("g.sym", "g.sym.sc", "g.ao", "g.ao.sc", "g.aoM.sc", "g.max3", "g.max3.sc", "u.sym", "Stk.joint", "Stk.marg", "Stk.LR", "GiViTI", "EDGE.poly3"),
               c("tag", "link", "n", "test", "size_raw")], idvar = c("tag", "link", "n"), timevar = "test", direction = "wide")
names(W) <- sub("size_raw\\.", "", names(W)); print(W, digits = 3, row.names = FALSE)
cat("\n=== NA rates (declined) under null / alternative for Stk.LR and GiViTI ===\n")
print(C[C$test %in% c("Stk.LR", "GiViTI") & (C$na_null > 0 | C$na_alt > 0), c("tag", "link", "n", "test", "na_null", "na_alt", "pow_adj", "pow_adj_finite")], digits = 3, row.names = FALSE)
cat("\n=== size-adjusted power (NA = no rejection) and raw power ===\n")
SH <- c("g.orc", "g.sym", "g.sym.sc", "g.ao", "g.ao.sc", "g.aoM", "g.aoM.sc", "g.max3", "g.max3.sc", "u.sym", "u.ao", "u.aoM", "u.max3", "EDGE.poly3", "Stk.joint", "Stk.LR", "GiViTI", "HL10")
for (v in c("pow_adj", "pow_raw")) {
  W <- reshape(C[C$test %in% SH, c("tag", "link", "n", "test", v)], idvar = c("tag", "link", "n"), timevar = "test", direction = "wide")
  names(W) <- sub(paste0(v, "\\."), "", names(W)); cat("\n--", v, "--\n"); print(W, digits = 3, row.names = FALSE)
}
cat("\n=== best rival among {GiViTI, Stk.LR, Stk.joint, EDGE.poly3} vs g.max3.sc / g.max3 / g.sym.sc ===\n")
for (k in unique(paste(C$tag, C$link, C$n))) {
  z <- C[paste(C$tag, C$link, C$n) == k, ]; pw <- setNames(z$pow_adj, z$test)
  rv <- pw[intersect(c("GiViTI", "Stk.LR", "Stk.joint", "EDGE.poly3"), names(pw))]
  cat(sprintf("%-40s best rival %-10s %.3f | max3.sc %s | max3 %s | sym.sc %s | u.max3 %s | vsGiViTI(max3.sc) %s | vsStkLR(max3.sc) %s\n", k,
    names(which.max(rv)), max(rv),
    ifelse("g.max3.sc" %in% names(pw), sprintf("%.3f", pw["g.max3.sc"]), "  -  "),
    sprintf("%.3f", pw["g.max3"]), ifelse("g.sym.sc" %in% names(pw), sprintf("%.3f", pw["g.sym.sc"]), "  -  "), sprintf("%.3f", pw["u.max3"]),
    ifelse("g.max3.sc" %in% names(pw), sprintf("%+.3f", pw["g.max3.sc"] - pw["GiViTI"]), "  -  "),
    ifelse("g.max3.sc" %in% names(pw), sprintf("%+.3f", pw["g.max3.sc"] - pw["Stk.LR"]), "  -  ")))
}
cat("\n=== grouped score-form probe vs its UNGROUPED classical score test (the fair 1-df rival) ===\n")
print(P[P$rival %in% c("u.sym", "u.ao", "u.aoM", "u.max3") & sub("^g\\.", "", sub("\\.sc$", "", P$probe)) == sub("^u\\.", "", P$rival),
        c("tag", "link", "n", "probe", "rival", "pow_probe", "pow_rival", "diff", "se", "p_mcnemar")], digits = 3, row.names = FALSE)
