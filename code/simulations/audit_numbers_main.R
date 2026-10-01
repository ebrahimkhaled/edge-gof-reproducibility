## audit_numbers_main.R -- NUMBERS / MONTE CARLO ERROR / MULTIPLICITY audit of run L and run L2 (claims C1-C7).
## Reads only; writes audit_numbers_*.csv. Single core.
suppressPackageStartupMessages(library(stats))
set.seed(20260913)
files <- setdiff(list.files(pattern = "^runL_.*_pvalues\\.csv$"), "runL_validate_base_pvalues.csv")
PL <- lapply(setNames(files, sub("^runL_(.*)_pvalues\\.csv$", "\\1", files)), function(f) read.csv(f, check.names = FALSE))
SC <- read.csv("scout_envelope.csv", stringsAsFactors = FALSE)
lam80 <- sapply(1:3, function(k) uniroot(function(l) pchisq(qchisq(.95, k), k, ncp = l, lower.tail = FALSE) - .8, c(.01, 200))$root)
crit1 <- function(p) { p <- p[is.finite(p)]; if (length(p) < 50) NA_real_ else as.numeric(quantile(p, .05, type = 1)) }
oc <- function(t, lk) if (grepl("orc$", t)) paste0(t, ".", lk) else t
rejv <- function(a, cr) ifelse(is.finite(a), a <= cr, FALSE)
ALLT <- c("g.orc", "g.ao", "g.aoM", "g.sym", "g.max3", "g.ao.sc", "g.aoM.sc", "g.sym.sc", "g.max3.sc", "EDGE.poly3",
          "EDGE.stk", "EDGE.stk.sc", "u.orc", "u.sym", "u.ao", "u.aoM", "u.cub", "u.max3", "Stk.joint", "Stk.marg",
          "Stk.LR", "GiViTI", "HL10", "HLG")

## ============ A. realised size of every test in every null cell ============
sz <- list()
for (tg in names(PL)) { P <- PL[[tg]]; alts <- setdiff(unique(P$link), "logit")
  for (n in unique(P$n)) { nul <- P[P$link == "logit" & P$n == n, ]
    for (lk in alts) for (t in ALLT) { if (lk != alts[1] && !grepl("orc$", t)) next
      cn <- oc(t, lk); if (!cn %in% names(nul)) next
      v <- nul[[cn]]; Bf <- sum(is.finite(v)); s <- mean(v[is.finite(v)] <= .05)
      sz[[length(sz) + 1]] <- data.frame(tag = tg, s = nul$s[1], c0 = nul$c0[1], n = n,
        test = if (grepl("orc$", t)) paste0(t, ".", lk) else t, Bnull = nrow(nul), Bfinite = Bf, size = s,
        na = mean(!is.finite(v)), z = (s - .05) / sqrt(.05 * .95 / Bf)) } } }
SZ <- do.call(rbind, sz); write.csv(SZ, "audit_numbers_size.csv", row.names = FALSE)
cat("=== A. realised size at alpha = .05 (null p <= .05 among finite), range over cells; |z| > 3 = outside 3 MC SE ===\n")
for (t in unique(sub("\\.(probit|cauchit|t4|cloglog|loglog)$", "", SZ$test))) {
  z <- SZ[sub("\\.(probit|cauchit|t4|cloglog|loglog)$", "", SZ$test) == t, ]
  cat(sprintf("  %-11s cells %2d  size %.4f .. %.4f  mean %.4f  |z|>3: %d  (max NA %.3f)  worst: %s\n", t, nrow(z),
    min(z$size), max(z$size), mean(z$size), sum(abs(z$z) > 3), max(z$na),
    paste(sprintf("%s/n=%d %.3f", z$tag, z$n, z$size)[order(-abs(z$z))][1:min(2, nrow(z))], collapse = "; ")))
}
cat("\n  Stk.marg by cell:\n"); z <- SZ[SZ$test == "Stk.marg", ]
cat(paste(sprintf("    %-18s n=%-6d B=%d size %.4f (z %+.1f)", z$tag, z$n, z$Bfinite, z$size, z$z), collapse = "\n"), "\n")
cat("\n  Stk.LR NA rate under the null by cell:", paste(sprintf("%s/%d %.3f", SZ$tag[SZ$test == "Stk.LR"], SZ$n[SZ$test == "Stk.LR"], SZ$na[SZ$test == "Stk.LR"]), collapse = "; "), "\n")

## ============ B. size-adjusted power in every alternative cell + paired comparisons ============
PROBES <- c("g.orc", "g.ao", "g.aoM", "g.sym", "g.max3", "g.ao.sc", "g.aoM.sc", "g.sym.sc", "g.max3.sc",
            "EDGE.poly3", "EDGE.stk", "EDGE.stk.sc", "u.orc", "u.sym", "u.max3")
RIVALS <- c("GiViTI", "Stk.LR", "Stk.joint")
POWT <- unique(c(PROBES, RIVALS, "u.ao", "u.aoM", "HL10", "Stk.marg"))
pw <- list(); cmp <- list(); cells <- list()
for (tg in names(PL)) { P <- PL[[tg]]
  for (n in unique(P$n)) { nul <- P[P$link == "logit" & P$n == n, ]
    for (lk in setdiff(unique(P$link), "logit")) { alt <- P[P$link == lk & P$n == n, ]
      rej <- list(); cr <- c()
      for (t in POWT) { cn <- oc(t, lk); if (!cn %in% names(alt)) next
        cr[t] <- crit1(nul[[cn]]); rej[[t]] <- rejv(alt[[cn]], cr[t])
        pw[[length(pw) + 1]] <- data.frame(tag = tg, link = lk, s = alt$s[1], c0 = alt$c0[1], n = n, test = t, B = nrow(alt),
          power = mean(rej[[t]]), se = sqrt(mean(rej[[t]]) * (1 - mean(rej[[t]])) / nrow(alt)), crit = cr[t]) }
      cells[[length(cells) + 1]] <- list(tag = tg, link = lk, n = n, nul = nul, alt = alt)
      for (a in intersect(PROBES, names(rej))) for (b in intersect(RIVALS, names(rej))) {
        A <- rej[[a]]; Bv <- rej[[b]]; n10 <- sum(A & !Bv); n01 <- sum(!A & Bv); B <- length(A)
        p <- if (n10 + n01 > 0) binom.test(n10, n10 + n01)$p.value else NA_real_
        d <- (n10 - n01) / B; se <- sqrt((n10 + n01) - (n10 - n01)^2 / B) / B
        cmp[[length(cmp) + 1]] <- data.frame(tag = tg, link = lk, s = alt$s[1], c0 = alt$c0[1], n = n, probe = a, rival = b, B = B,
          pow_probe = mean(A), pow_rival = mean(Bv), diff = d, n10 = n10, n01 = n01, discordant = n10 + n01,
          se_pair = se, zratio = d / se, p = p) } } } }
PW <- do.call(rbind, pw); C <- do.call(rbind, cmp)
## reproduce the lead's runL_paired.csv
old <- read.csv("runL_paired.csv", stringsAsFactors = FALSE)
mm <- merge(old, C, by = c("tag", "link", "n", "probe", "rival"), suffixes = c(".old", ".new"))
cat(sprintf("\n=== B. reproduction of runL_paired.csv: %d rows old, %d new, %d matched; max|diff| %.2e, max|disc| %d, max|log p ratio| %.2e ===\n",
  nrow(old), nrow(C), nrow(mm), max(abs(mm$diff.old - mm$diff.new)), max(abs(mm$discordant.old - mm$discordant.new)),
  max(abs(log(mm$p.old) - log(mm$p.new)), na.rm = TRUE)))
pvs <- read.csv("runL_power_vs_scout.csv", stringsAsFactors = FALSE)
m2 <- merge(pvs, PW, by = c("tag", "link", "n", "test"))
cat(sprintf("    reproduction of runL_power_vs_scout.csv power_adj: %d rows matched, max|diff| %.2e\n", nrow(m2), max(abs(m2$power_adj - m2$power))))
## compare with the harness summary rule (NA alternatives dropped from the denominator)
sm <- do.call(rbind, lapply(names(PL), function(tg) { z <- read.csv(sprintf("runL_%s_summary.csv", tg), stringsAsFactors = FALSE); z$tag <- tg; z }))
m3 <- merge(sm, PW, by = c("tag", "link", "n", "test"))
m3$d <- m3$power_adj - m3$power
cat(sprintf("    summary.csv power_adj (NA dropped) vs NA-as-no-reject: %d rows, max|diff| %.4f; rows with |diff|>.005: %s\n", nrow(m3),
  max(abs(m3$d), na.rm = TRUE), paste(sprintf("%s/%s/%d/%s %+.3f", m3$tag, m3$link, m3$n, m3$test, m3$d)[which(abs(m3$d) > .005)], collapse = "; ")))

## multiplicity: global and within claim families
C$p_holm_all <- NA; C$p_bh_all <- NA; ok <- is.finite(C$p)
C$p_holm_all[ok] <- p.adjust(C$p[ok], "holm"); C$p_bh_all[ok] <- p.adjust(C$p[ok], "BH")
matched <- function(link, probe) (link == "cloglog" & probe %in% c("g.ao", "g.ao.sc")) | (link == "loglog" & probe %in% c("g.aoM", "g.aoM.sc"))
C$fam <- ""
C$fam[C$link %in% c("probit", "cauchit", "t4") & C$c0 == 0 & C$probe %in% c("g.sym", "g.sym.sc") & C$rival %in% c("GiViTI", "Stk.LR")] <- "C3"
C$fam[C$link %in% c("cloglog", "loglog") & C$rival == "GiViTI" & (C$probe %in% c("g.orc", "u.orc") | matched(C$link, C$probe))] <- "C5"
C$fam[C$probe == "g.max3.sc" & C$rival %in% c("Stk.LR", "GiViTI")] <- "C7"
C$fam[C$link == "probit" & C$c0 == -2 & C$probe %in% c("g.sym", "g.sym.sc", "u.sym") & C$rival %in% c("GiViTI", "Stk.LR")] <- "LOWP"
C$p_holm_fam <- NA; C$p_bh_fam <- NA
for (f in setdiff(unique(C$fam), "")) { k <- C$fam == f & is.finite(C$p)
  C$p_holm_fam[k] <- p.adjust(C$p[k], "holm"); C$p_bh_fam[k] <- p.adjust(C$p[k], "BH") }
cat(sprintf("\n=== multiplicity over ALL %d McNemar tests (finite p %d) ===\n", nrow(C), sum(ok)))
cat(sprintf("  raw p<.05: %d | Holm<.05: %d | BH q<.05: %d | Bonferroni<.05: %d\n", sum(C$p < .05, na.rm = TRUE),
  sum(C$p_holm_all < .05, na.rm = TRUE), sum(C$p_bh_all < .05, na.rm = TRUE), sum(p.adjust(C$p[ok], "bonferroni") < .05)))
cat(sprintf("  raw p<.05 that fail global Holm: %d ; fail global BH: %d\n", sum(C$p < .05 & C$p_holm_all >= .05, na.rm = TRUE),
  sum(C$p < .05 & C$p_bh_all >= .05, na.rm = TRUE)))
cat(sprintf("  |zratio| < 2.8 among raw p<.05: %d\n", sum(C$p < .05 & abs(C$zratio) < 2.8, na.rm = TRUE)))

## bootstrap SD of the paired difference INCLUDING the estimated size-adjusted critical values
NBOOT <- 1000
bootsd <- rep(NA_real_, nrow(C))
famrows <- which(C$fam != "")
for (cc in cells) { k <- famrows[C$tag[famrows] == cc$tag & C$link[famrows] == cc$link & C$n[famrows] == cc$n]
  if (!length(k)) next
  tt <- unique(c(C$probe[k], C$rival[k])); Bn <- nrow(cc$nul); Ba <- nrow(cc$alt)
  PN <- sapply(tt, function(t) cc$nul[[oc(t, cc$link)]]); PA <- sapply(tt, function(t) cc$alt[[oc(t, cc$link)]])
  BP <- matrix(NA_real_, NBOOT, length(tt), dimnames = list(NULL, tt))
  for (bb in seq_len(NBOOT)) { i0 <- sample.int(Bn, Bn, TRUE); i1 <- sample.int(Ba, Ba, TRUE)
    for (j in seq_along(tt)) { cr <- crit1(PN[i0, j]); a <- PA[i1, j]; BP[bb, j] <- mean(is.finite(a) & a <= cr) } }
  ## paired: same alternative resample for all tests; the difference SD uses the common resample
  for (r in k) bootsd[r] <- sd(BP[, C$probe[r]] - BP[, C$rival[r]]) }
C$boot_sd <- bootsd; C$z_boot <- C$diff / C$boot_sd
write.csv(C, "audit_numbers_paired_full.csv", row.names = FALSE)
show <- function(f) { z <- C[C$fam == f, ]; z <- z[order(z$link, z$tag, z$n, z$probe, z$rival), ]
  cat(sprintf("\n=== family %s: %d comparisons; raw p<.05 %d; family Holm<.05 %d; family BH<.05 %d; global Holm<.05 %d ===\n", f, nrow(z),
    sum(z$p < .05, na.rm = TRUE), sum(z$p_holm_fam < .05, na.rm = TRUE), sum(z$p_bh_fam < .05, na.rm = TRUE), sum(z$p_holm_all < .05, na.rm = TRUE)))
  cat("  tag                 link     n      probe      rival     B     probe  rival   diff   n10 n01  se_pair z_pair boot_sd z_boot   p        Holm_fam BH_fam  Holm_all\n")
  for (j in seq_len(nrow(z))) with(z[j, ], cat(sprintf("  %-19s %-8s %-6d %-10s %-9s %-5d %.3f  %.3f  %+.3f %4d %4d  %.4f  %+5.1f  %.4f  %+5.1f  %.1e  %.1e  %.1e  %.1e\n",
    tag, link, n, probe, rival, B, pow_probe, pow_rival, diff, n10, n01, se_pair, zratio, boot_sd, z_boot, p, p_holm_fam, p_bh_fam, p_holm_all))) }
for (f in c("C3", "C5", "C7", "LOWP")) show(f)

## C7 ranges
cat("\n=== C7: g.max3.sc minus rival, range over every run L2 cell ===\n")
for (rv in c("Stk.LR", "GiViTI")) { z <- C[C$probe == "g.max3.sc" & C$rival == rv, ]
  cat(sprintf("  vs %-7s diff %+.3f .. %+.3f over %d cells\n", rv, min(z$diff), max(z$diff), nrow(z))) }
z <- PW[PW$test %in% c("g.max3.sc", "Stk.LR", "GiViTI", "Stk.joint"), ]
zz <- reshape(z[, c("tag", "link", "n", "test", "power")], idvar = c("tag", "link", "n"), timevar = "test", direction = "wide")
zz <- zz[is.finite(zz$power.g.max3.sc), ]
zz$best_rival <- pmax(zz$power.Stk.LR, zz$power.GiViTI, zz$power.Stk.joint)
zz$gap <- zz$power.g.max3.sc - zz$best_rival
print(zz[order(zz$gap), ], row.names = FALSE, digits = 3)
cat("\n  g.max3 (UNIT weights) minus best of Stk.LR/GiViTI/Stk.joint in batch-1 cells (no .sc there):\n")
z <- PW[PW$test %in% c("g.max3", "Stk.LR", "GiViTI", "Stk.joint"), ]
zz <- reshape(z[, c("tag", "link", "n", "test", "power")], idvar = c("tag", "link", "n"), timevar = "test", direction = "wide")
zz$gap <- zz$power.g.max3 - pmax(zz$power.Stk.LR, zz$power.GiViTI, zz$power.Stk.joint)
print(zz[order(zz$gap), ], row.names = FALSE, digits = 3)

## ============ C. key power numbers for every cell ============
cat("\n=== C. size-adjusted power (NA = no rejection), key tests ===\n")
KT <- c("g.orc", "u.orc", "g.sym", "g.sym.sc", "u.sym", "g.ao", "g.ao.sc", "g.aoM", "g.aoM.sc", "g.max3", "g.max3.sc", "EDGE.poly3", "EDGE.stk.sc", "Stk.LR", "Stk.joint", "GiViTI", "HL10")
for (cc in cells) { z <- PW[PW$tag == cc$tag & PW$link == cc$link & PW$n == cc$n, ]
  cat(sprintf("  %-18s %-7s n=%-6d B=%d | %s\n", cc$tag, cc$link, cc$n, z$B[1],
    paste(sprintf("%s %.3f", z$test[match(intersect(KT, z$test), z$test)], z$power[match(intersect(KT, z$test), z$test)]), collapse = "  "))) }
write.csv(PW, "audit_numbers_power.csv", row.names = FALSE)

## ============ D. n for 80% power by interpolation between the two simulated n, with bootstrap intervals ============
KDF <- c(g.orc = 1, u.orc = 1, g.sym = 1, g.sym.sc = 1, u.sym = 1, g.ao = 1, g.ao.sc = 1, g.aoM = 1, g.aoM.sc = 1, Stk.LR = 2, Stk.joint = 2, EDGE.stk.sc = 2)
SCN <- c(g.orc = "matched", u.orc = "matched", g.sym = "sym", g.sym.sc = "sym", u.sym = "sym", g.ao = "ao", g.ao.sc = "ao",
         g.aoM = "aoM", g.aoM.sc = "aoM", g.max3 = "maxprobe", g.max3.sc = "maxprobe", u.max3 = "maxprobe",
         Stk.LR = "stukel", Stk.joint = "stukel", GiViTI = "GiViTI", HL10 = "HL10")
DT <- c("g.orc", "u.orc", "g.sym", "g.sym.sc", "u.sym", "g.ao", "g.ao.sc", "g.aoM", "g.aoM.sc", "g.max3", "g.max3.sc", "u.max3",
        "EDGE.poly3", "EDGE.stk.sc", "Stk.LR", "Stk.joint", "GiViTI", "HL10")
ncp_of <- function(p, k) { if (!is.finite(p) || p <= .05) return(0); if (p >= .9999) return(NA_real_)
  uniroot(function(l) pchisq(qchisq(.95, k), k, ncp = l, lower.tail = FALSE) - p, c(0, 400))$root }
interp <- function(n1, n2, p1, p2, fx) { p1 <- min(max(p1, 1e-3), 1 - 1e-3); p2 <- min(max(p2, 1e-3), 1 - 1e-3)
  y1 <- qnorm(p1); y2 <- qnorm(p2); x1 <- fx(n1); x2 <- fx(n2); b <- (y2 - y1) / (x2 - x1); if (b <= 0) return(Inf)
  (qnorm(.8) - (y1 - b * x1)) / b }
dres <- list()
for (tg in names(PL)) { P <- PL[[tg]]; ns <- sort(unique(P$n)); if (length(ns) != 2) next
  for (lk in setdiff(unique(P$link), "logit")) {
    s <- P$s[1]; c0 <- P$c0[1]
    nul <- lapply(ns, function(n) P[P$link == "logit" & P$n == n, ]); alt <- lapply(ns, function(n) P[P$link == lk & P$n == n, ])
    mt <- SC[SC$s == s & SC$c0 == c0 & SC$link == lk & SC$test == "matched", "n80"]
    for (t in DT) { cn <- oc(t, lk); if (!cn %in% names(alt[[1]])) next
      pw2 <- sapply(1:2, function(i) mean(rejv(alt[[i]][[cn]], crit1(nul[[i]][[cn]]))))
      n_sqrt <- interp(ns[1], ns[2], pw2[1], pw2[2], sqrt)^2
      n_lin <- interp(ns[1], ns[2], pw2[1], pw2[2], identity)
      n_log <- exp(interp(ns[1], ns[2], pw2[1], pw2[2], log))
      n_ncp <- NA_real_; n_ncp0 <- c(NA, NA)
      if (t %in% names(KDF)) { k <- KDF[[t]]; nc <- sapply(pw2, ncp_of, k = k)
        if (all(is.finite(nc))) { b <- (nc[2] - nc[1]) / (ns[2] - ns[1]); a <- nc[1] - b * ns[1]
          n_ncp <- if (b > 0) (lam80[k] - a) / b else Inf; n_ncp0 <- ns * lam80[k] / nc } }
      ## bootstrap of the sqrt-probit interpolation
      bs <- replicate(400, { pb <- sapply(1:2, function(i) { i0 <- sample.int(nrow(nul[[i]]), replace = TRUE); i1 <- sample.int(nrow(alt[[i]]), replace = TRUE)
        mean(rejv(alt[[i]][[cn]][i1], crit1(nul[[i]][[cn]][i0]))) }); interp(ns[1], ns[2], pb[1], pb[2], sqrt)^2 })
      sn <- if (t %in% names(SCN)) SCN[[t]] else NA
      scout <- if (!is.na(sn)) { v <- SC[SC$s == s & SC$c0 == c0 & SC$link == lk & SC$test == sn, "n80"]; if (length(v)) v else NA } else NA
      ceil <- sapply(ns, function(n) pchisq(qchisq(.95, 1), 1, ncp = n * lam80[1] / mt, lower.tail = FALSE))
      dres[[length(dres) + 1]] <- data.frame(tag = tg, link = lk, s = s, c0 = c0, test = t, n1 = ns[1], n2 = ns[2],
        pow1 = pw2[1], pow2 = pw2[2], ceil1 = ceil[1], ceil2 = ceil[2], n80_sqrtprobit = n_sqrt, lo95 = quantile(bs, .025, na.rm = TRUE),
        hi95 = quantile(bs, .975, na.rm = TRUE), frac_boot_inf = mean(!is.finite(bs)), n80_linprobit = n_lin, n80_logprobit = n_log,
        n80_ncp_line = n_ncp, n80_ncp_origin_n1 = n_ncp0[1], n80_ncp_origin_n2 = n_ncp0[2], scout_n80 = scout, scout_matched = mt) } } }
D <- do.call(rbind, dres); rownames(D) <- NULL
write.csv(D, "audit_numbers_n80.csv", row.names = FALSE)
cat("\n=== D. n for 80% power: probit(power) linear in sqrt(n) between the two simulated n; bootstrap 95% (400 resamples, null and alternative); other interpolations; scout ===\n")
for (j in seq_len(nrow(D))) with(D[j, ], cat(sprintf("  %-16s %-7s %-11s n=%d/%d pow %.3f/%.3f ceil %.3f/%.3f | n80 sqrt %6.0f [%6.0f, %6.0f] inf %.2f | lin %6.0f log %6.0f | ncp line %6.0f origin %6.0f/%6.0f | scout %7.0f matched %6.0f\n",
  tag, link, test, n1, n2, pow1, pow2, ceil1, ceil2, n80_sqrtprobit, lo95, hi95, frac_boot_inf, n80_linprobit, n80_logprobit, n80_ncp_line,
  n80_ncp_origin_n1, n80_ncp_origin_n2, scout_n80, scout_matched)))
## single-n cells (m10 files): through-origin ncp only
cat("\n  single-n cells, through-origin ncp estimate n*lambda80/ncp:\n")
for (tg in names(PL)) { P <- PL[[tg]]; ns <- sort(unique(P$n)); if (length(ns) != 1) next
  for (lk in setdiff(unique(P$link), "logit")) { nul <- P[P$link == "logit", ]; alt <- P[P$link == lk, ]
    for (t in intersect(names(KDF), DT)) { cn <- oc(t, lk); if (!cn %in% names(alt)) next
      p <- mean(rejv(alt[[cn]], crit1(nul[[cn]]))); k <- KDF[[t]]; nc <- ncp_of(p, k)
      cat(sprintf("    %-16s %-7s %-10s n=%d power %.3f -> n80 %6.0f\n", tg, lk, t, ns, p, ns * lam80[k] / nc)) } } }

## ============ E. GiViTI mechanism numbers (C4) ============
cat("\n=== E. giviti_mechanism_{logit,cauchit}.csv ===\n")
ci <- function(k, m) { q <- binom.test(k, m)$conf.int; sprintf("%.3f [%.3f, %.3f] (%d/%d)", k / m, q[1], q[2], k, m) }
for (lk in c("logit", "cauchit")) { M <- read.csv(sprintf("giviti_mechanism_%s.csv", lk)); B <- nrow(M)
  cat(sprintf("  %s: B=%d, GiViTI NA %d; m table: %s\n", lk, B, sum(is.na(M$p)), paste(names(table(M$m)), table(M$m), sep = ":", collapse = " ")))
  cat("    P(m==2)", ci(sum(M$m == 2, na.rm = TRUE), sum(!is.na(M$m))), " | P(lr3>3.84)", ci(sum(M$lr3 > 3.84), B), "\n")
  cat("    cross-tab m vs lr3>3.84:\n"); print(table(m = M$m, lr3_pass = M$lr3 > 3.84))
  cat("    P(GiViTI p<=.05)", ci(sum(M$p <= .05, na.rm = TRUE), sum(!is.na(M$p))), "\n")
  for (k in sort(unique(M$m))) cat(sprintf("      m=%d: rej %s\n", k, ci(sum(M$p[M$m == k] <= .05, na.rm = TRUE), sum(M$m == k, na.rm = TRUE))))
  cat("    fixed-degree LR vs degree 1: deg2", ci(sum(M$lr2 > qchisq(.95, 1)), B), " deg3", ci(sum(M$lr2 + M$lr3 > qchisq(.95, 2)), B),
      " deg4", ci(sum(M$lr2 + M$lr3 + M$lr4 > qchisq(.95, 3)), B), "\n")
  ## paired: GiViTI vs fixed deg-3 LR on the same samples
  A <- M$lr2 + M$lr3 > qchisq(.95, 2); G <- M$p <= .05; G[is.na(G)] <- FALSE
  cat(sprintf("    McNemar fixed-deg3 vs GiViTI: n10 %d n01 %d p %.2e\n", sum(A & !G), sum(!A & G), binom.test(sum(A & !G), sum(A & !G) + sum(!A & G))$p.value))
}
cat("\nDONE\n")
