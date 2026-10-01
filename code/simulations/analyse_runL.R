## analyse_runL.R -- run L, paired: do the link probes beat GiViTI and Stukel on the same samples,
## and does finite-sample power match the population scout's prediction?
##
## Every test in run L is computed from the same fit in a replicate, so rejections are paired and
## compared by McNemar's exact test on discordant pairs, each test at its own size-adjusted critical
## value from the matched null (the 5% quantile of its null p-values at the same n and design).
## Predicted power for the one-df probes and for Stukel's two-df test is the local-limit value
## 1 - F_chi2_k(q_k; ncp), ncp = n * lambda80(k) / n80 from scout_envelope.csv.
files <- setdiff(list.files(pattern = "^runL_.*_pvalues\\.csv$"), "runL_validate_base_pvalues.csv")
stopifnot(length(files) > 0)
## one data frame per run: the oracle columns are named per link, so runs cannot be row-bound
PL <- lapply(setNames(files, sub("^runL_(.*)_pvalues\\.csv$", "\\1", files)), function(f) read.csv(f, check.names = FALSE))
SC <- read.csv("scout_envelope.csv", stringsAsFactors = FALSE)
lam80 <- sapply(1:3, function(k) uniroot(function(l) pchisq(qchisq(.95, k), k, ncp = l, lower.tail = FALSE) - .8, c(.01, 200))$root)

crit <- function(p) { p <- p[is.finite(p)]; if (length(p) < 50) NA_real_ else as.numeric(quantile(p, .05, type = 1)) }
## ".sc" = score form (shapes times sqrt(V_g), plain chi-square null); present only in run L2 files
PROBES <- c("g.orc", "g.ao", "g.aoM", "g.sym", "g.max3", "g.ao.sc", "g.aoM.sc", "g.sym.sc", "g.max3.sc",
            "EDGE.poly3", "EDGE.stk", "EDGE.stk.sc", "u.orc", "u.sym", "u.max3")
## GiViTI.t50 and Cubic.LR exist only in run L3 files; u.sym (Stukel's one-parameter symmetric score test) is a rival
## to the grouped symmetric probe, as the claims audit requires
RIVALS <- c("GiViTI", "GiViTI.t50", "Cubic.LR", "Stk.LR", "Stk.joint", "u.sym")
SCOUT_NAME <- c(g.orc = "matched", u.orc = "matched", g.ao = "ao", g.aoM = "aoM", g.sym = "sym", u.sym = "sym",
                g.ao.sc = "ao", g.aoM.sc = "aoM", g.sym.sc = "sym", Stk.joint = "stukel", Stk.LR = "stukel")
SCOUT_K <- c(matched = 1, ao = 1, aoM = 1, sym = 1, stukel = 2)

rows <- list(); cmp <- list()
for (tg in names(PL)) for (n in unique(PL[[tg]]$n)) {
  Q <- PL[[tg]][PL[[tg]]$n == n, ]; s <- Q$s[1]; c0 <- Q$c0[1]
  nul <- Q[Q$link == "logit", ]
  for (lk in setdiff(unique(Q$link), "logit")) {
    alt <- Q[Q$link == lk, ]
    col_of <- function(t) if (grepl("orc$", t)) paste0(t, ".", lk) else t
    rej <- list()
    for (t in c(PROBES, RIVALS)) { cn <- col_of(t); if (!cn %in% names(alt)) next
      cr <- crit(nul[[cn]]); a <- alt[[cn]]
      rej[[t]] <- ifelse(is.finite(a), a <= cr, FALSE)                       # a declined test does not reject
      sn <- SCOUT_NAME[t]; pred <- NA_real_
      if (!is.na(sn)) { z <- SC[SC$s == s & SC$c0 == c0 & SC$link == lk & SC$test == sn, ]
        if (nrow(z) && is.finite(z$n80)) { k <- SCOUT_K[[sn]]; pred <- pchisq(qchisq(.95, k), k, ncp = n * lam80[k] / z$n80, lower.tail = FALSE) } }
      rows[[length(rows) + 1]] <- data.frame(tag = tg, link = lk, s = s, c0 = c0, n = n, test = t,
        size = mean(nul[[cn]] <= .05, na.rm = TRUE), power_adj = mean(rej[[t]]), declined = mean(!is.finite(a)),
        predicted = pred)
    }
    for (a in intersect(PROBES, names(rej))) for (b in intersect(RIVALS, names(rej))) {
      A <- rej[[a]]; Bv <- rej[[b]]; n10 <- sum(A & !Bv); n01 <- sum(!A & Bv)
      p <- if (n10 + n01 > 0) binom.test(n10, n10 + n01)$p.value else NA_real_
      cmp[[length(cmp) + 1]] <- data.frame(tag = tg, link = lk, s = s, c0 = c0, n = n, probe = a, rival = b,
        pow_probe = mean(A), pow_rival = mean(Bv), diff = mean(A) - mean(Bv), discordant = n10 + n01, p = p,
        verdict = ifelse(is.na(p) | p >= .05, "not significant", ifelse(mean(A) > mean(Bv), "probe better", "rival better")))
    }
  }
}
R <- do.call(rbind, rows); C <- do.call(rbind, cmp)
write.csv(R, "runL_power_vs_scout.csv", row.names = FALSE); write.csv(C, "runL_paired.csv", row.names = FALSE)

cat("=== size-adjusted power, with the population-limit prediction where it exists ===\n")
for (i in seq_len(nrow(unique(R[, c("tag", "link", "n")])))) {
  key <- unique(R[, c("tag", "link", "n")])[i, ]; z <- R[R$tag == key$tag & R$link == key$link & R$n == key$n, ]
  cat(sprintf("\n-- %s  s=%g c0=%g  n=%d --\n", key$link, z$s[1], z$c0[1], key$n))
  for (j in seq_len(nrow(z))) cat(sprintf("  %-10s power %.3f  pred %5s  size %.3f  declined %.3f\n", z$test[j], z$power_adj[j],
    ifelse(is.na(z$predicted[j]), "  -  ", sprintf("%.3f", z$predicted[j])), z$size[j], z$declined[j]))
}
cat("\n=== paired: best general-purpose probe (g.max3) and matched probes vs each rival ===\n")
k <- C[C$probe %in% c("g.max3", "g.ao", "g.aoM", "g.sym"), ]
for (j in seq_len(nrow(k))) with(k[j, ], cat(sprintf("  %-8s s=%g c0=%+g n=%-6d %-6s vs %-9s %.3f vs %.3f  diff %+.3f  p=%s  %s\n",
  link, s, c0, n, probe, rival, pow_probe, pow_rival, diff, ifelse(is.na(p), "NA", format.pval(p, digits = 2)), verdict)))
