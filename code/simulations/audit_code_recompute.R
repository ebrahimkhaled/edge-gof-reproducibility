## audit_code_recompute.R -- code-lens audit (read-only): recompute the headline numbers of run L / L2 directly from
## the per-replicate p-value files, independently of analyse_runL.R. Writes detail tables to the scratchpad only.
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
setwd(edge_path("code/simulations"))
OUT <- file.path(tempdir())
options(width = 250)
files <- setdiff(list.files(pattern = "^runL_.*_pvalues\\.csv$"), "runL_validate_base_pvalues.csv")
PL <- lapply(setNames(files, sub("^runL_(.*)_pvalues\\.csv$", "\\1", files)), function(f) read.csv(f, check.names = FALSE))
cat("tags:", names(PL), "\n")

## ---------------------------------------------------------------- A. null calibration
cat("\n==== A. null calibration (logit cells) ====\n")
A <- list()
for (tg in names(PL)) for (n in unique(PL[[tg]]$n)) {
  nul <- PL[[tg]][PL[[tg]]$link == "logit" & PL[[tg]]$n == n, ]
  cols <- setdiff(names(nul), c("link", "s", "c0", "n", "rep"))
  for (cn in cols) { p <- nul[[cn]]; f <- p[is.finite(p)]; if (length(f) < 50) next
    ks <- suppressWarnings(ks.test(f, "punif")$p.value)
    A[[length(A) + 1]] <- data.frame(tag = tg, n = n, test = cn, B = nrow(nul), na = mean(!is.finite(p)),
      s01 = mean(f <= .01), s05 = mean(f <= .05), s10 = mean(f <= .10), ks_p = ks, pmin = min(f), pmax = max(f)) }
}
A <- do.call(rbind, A); write.csv(A, file.path(OUT, "audit_code_null.csv"), row.names = FALSE)
cat("any p outside [0,1]:", any(A$pmin < 0 | A$pmax > 1), "\n")
foc <- c("g.sym", "g.max3", "g.sym.sc", "g.ao.sc", "g.aoM.sc", "g.max3.sc", "EDGE.stk.sc", "u.sym", "u.max3", "Stk.joint", "Stk.marg", "Stk.LR", "GiViTI")
z <- A[A$test %in% foc, ]; z$flag05 <- ifelse(abs(z$s05 - .05) > 3 * sqrt(.05 * .95 / z$B), "*", "")
print(z[order(z$test, z$tag, z$n), c("test", "tag", "n", "B", "na", "s01", "s05", "s10", "ks_p", "flag05")], digits = 3, row.names = FALSE)
cat("\n-- pooled over all run L2 null cells (score form) --\n")
for (cn in c("g.sym.sc", "g.ao.sc", "g.aoM.sc", "g.max3.sc", "EDGE.stk.sc")) {
  p <- unlist(lapply(names(PL), function(tg) { Q <- PL[[tg]]; if (!cn %in% names(Q)) return(NULL); Q[Q$link == "logit", cn] }))
  f <- p[is.finite(p)]
  cat(sprintf("  %-12s N=%5d na=%.4f  size .01 %.4f  .05 %.4f  .10 %.4f  (SE at .05 %.4f)  KS p=%.3g\n", cn, length(p), mean(!is.finite(p)),
    mean(f <= .01), mean(f <= .05), mean(f <= .10), sqrt(.05 * .95 / length(f)), suppressWarnings(ks.test(f, "punif")$p.value)))
}
cat("\n-- Stukel variants: realised size at 0.05 across every null cell --\n")
for (cn in c("Stk.marg", "Stk.joint", "Stk.LR")) { zz <- A[A$test == cn, ]
  cat(sprintf("  %-9s range %.3f-%.3f  (cells %d; na range %.3f-%.3f)\n", cn, min(zz$s05), max(zz$s05), nrow(zz), min(zz$na), max(zz$na))) }

## ---------------------------------------------------------------- B. size-adjusted power, two conventions
crit <- function(p) { p <- p[is.finite(p)]; if (length(p) < 50) NA_real_ else as.numeric(quantile(p, .05, type = 1)) }
cell <- function(tg, lk, n, tests, show = TRUE) {
  Q <- PL[[tg]]; nul <- Q[Q$link == "logit" & Q$n == n, ]; alt <- Q[Q$link == lk & Q$n == n, ]
  res <- lapply(tests, function(t) { if (!t %in% names(alt)) return(NULL)
    cr <- crit(nul[[t]]); a <- alt[[t]]
    data.frame(tag = tg, link = lk, n = n, test = t, B = nrow(alt), crit = cr, ties_at_crit = sum(nul[[t]] == cr, na.rm = TRUE),
      pow_NAasNo = mean(ifelse(is.finite(a), a <= cr, FALSE)), pow_NAexcl = mean(a[is.finite(a)] <= cr),
      pow_raw = mean(a <= .05, na.rm = TRUE), size = mean(nul[[t]] <= .05, na.rm = TRUE), na_alt = mean(!is.finite(a)), na_null = mean(!is.finite(nul[[t]]))) })
  r <- do.call(rbind, res); if (show) print(r, digits = 3, row.names = FALSE); invisible(list(tab = r, alt = alt, nul = nul))
}
mcnemar <- function(L, a, b) { cra <- crit(L$nul[[a]]); crb <- crit(L$nul[[b]])
  A <- ifelse(is.finite(L$alt[[a]]), L$alt[[a]] <= cra, FALSE); Bv <- ifelse(is.finite(L$alt[[b]]), L$alt[[b]] <= crb, FALSE)
  n10 <- sum(A & !Bv); n01 <- sum(!A & Bv)
  cat(sprintf("   McNemar %s vs %s: %.3f vs %.3f diff %+.3f  (n10=%d n01=%d, exact p=%.2g)\n", a, b, mean(A), mean(Bv), mean(A) - mean(Bv), n10, n01,
    binom.test(n10, n10 + n01)$p.value)) }

cat("\n==== B1. HEADLINE: cauchit base n=4000 (and 2500) ====\n")
T1 <- c("g.sym", "g.max3", "u.sym", "g.orc.cauchit", "u.orc.cauchit", "Stk.LR", "Stk.joint", "Stk.marg", "GiViTI", "HL10", "EDGE.poly3", "EDGE.stk")
for (n in c(2500, 4000)) { L <- cell("base_cauchit", "cauchit", n, T1); mcnemar(L, "g.sym", "GiViTI"); mcnemar(L, "g.sym", "Stk.LR") }
L <- cell("base_cauchit_m10", "cauchit", 2500, c(T1, "g.sym.sc", "g.max3.sc")); mcnemar(L, "g.sym", "GiViTI"); mcnemar(L, "g.sym.sc", "Stk.LR")

cat("\n==== B2. HEADLINE: probit base n=16000 (and 10000) ====\n")
T2 <- c("g.sym", "g.max3", "u.sym", "g.orc.probit", "u.orc.probit", "Stk.LR", "Stk.joint", "GiViTI", "HL10", "EDGE.poly3", "EDGE.stk")
for (n in c(10000, 16000)) { L <- cell("base_probit", "probit", n, T2); mcnemar(L, "g.sym", "Stk.LR"); mcnemar(L, "g.sym", "GiViTI") }

cat("\n==== B3. HEADLINE: probit s=2 n=6500 (and 4000), score form ====\n")
T3 <- c("g.sym.sc", "g.max3.sc", "EDGE.stk.sc", "g.sym", "g.max3", "u.sym", "u.max3", "g.orc.probit", "u.orc.probit", "Stk.LR", "Stk.joint", "GiViTI", "EDGE.poly3")
for (tg in c("s2_probit_m25sc", "s2_probit_m10")) for (n in c(4000, 6500)) { cat("--", tg, "\n"); L <- cell(tg, "probit", n, T3)
  mcnemar(L, "g.sym.sc", "Stk.LR"); mcnemar(L, "g.sym.sc", "GiViTI"); mcnemar(L, "g.max3.sc", "Stk.LR") }
for (n in c(4000, 6500)) { cat("-- s2_probit (batch 1)\n"); L <- cell("s2_probit", "probit", n, T3) }

cat("\n==== B4. other symmetric cells (C3) ====\n")
T4 <- function(lk) c("g.sym", "g.sym.sc", "g.max3", "g.max3.sc", "u.sym", paste0("g.orc.", lk), paste0("u.orc.", lk), "Stk.LR", "Stk.joint", "GiViTI", "EDGE.poly3")
for (x in list(c("low_probit", "probit", 5000), c("low_probit", "probit", 8000), c("low_probit_m10", "probit", 5000),
               c("s2_cauchit", "cauchit", 600), c("s2_cauchit", "cauchit", 900), c("s2_t4", "t4", 6000), c("s2_t4", "t4", 9000), c("s2_t4_m10", "t4", 6000))) {
  cat("--", x[1], "\n"); cell(x[1], x[2], as.integer(x[3]), T4(x[2])) }

cat("\n==== B5. bow cells (C5) ====\n")
for (tg in c("base_asym", "base_asym_m10", "base_asym_m5")) for (lk in c("cloglog", "loglog")) for (n in c(600, 1000)) {
  cat("--", tg, "\n"); L <- cell(tg, lk, n, c("g.ao", "g.aoM", "g.ao.sc", "g.aoM.sc", "g.max3", "g.max3.sc", "u.ao", "u.aoM", paste0("g.orc.", lk), paste0("u.orc.", lk), "GiViTI", "Stk.LR"))
  mt <- if (lk == "cloglog") "g.ao" else "g.aoM"; mcnemar(L, mt, "GiViTI"); mcnemar(L, paste0("g.orc.", lk), "GiViTI"); mcnemar(L, paste0("u.orc.", lk), "GiViTI") }

## ---------------------------------------------------------------- C. every cell: max3 vs rivals, and the NA convention
cat("\n==== C. g.max3(.sc) minus rivals, every cell; and |NA-as-no minus NA-excluded| ====\n")
allc <- list()
for (tg in names(PL)) { Q <- PL[[tg]]; for (n in unique(Q$n)) for (lk in setdiff(unique(Q$link), "logit")) {
  tests <- intersect(c("g.max3", "g.max3.sc", "g.sym", "g.sym.sc", "Stk.LR", "Stk.joint", "GiViTI", paste0("g.orc.", lk), paste0("u.orc.", lk)), names(Q))
  allc[[length(allc) + 1]] <- cell(tg, lk, n, tests, show = FALSE)$tab } }
allc <- do.call(rbind, allc); write.csv(allc, file.path(OUT, "audit_code_power.csv"), row.names = FALSE)
cat("max |pow_NAasNo - pow_NAexcl| over all:", max(abs(allc$pow_NAasNo - allc$pow_NAexcl)), "\n")
print(allc[abs(allc$pow_NAasNo - allc$pow_NAexcl) > .003, c("tag", "link", "n", "test", "pow_NAasNo", "pow_NAexcl", "na_alt", "na_null")], digits = 3, row.names = FALSE)
cat("max ties at the critical value:", max(allc$ties_at_crit), "\n")
W <- reshape(allc[, c("tag", "link", "n", "test", "pow_NAasNo")], idvar = c("tag", "link", "n"), timevar = "test", direction = "wide")
names(W) <- sub("pow_NAasNo\\.", "", names(W))
W$max3sc_vs_LR <- W$g.max3.sc - W$Stk.LR; W$max3sc_vs_Gi <- W$g.max3.sc - W$GiViTI
W$max3_vs_LR <- W$g.max3 - W$Stk.LR; W$max3_vs_Gi <- W$g.max3 - W$GiViTI
W$sym_vs_Gi <- W$g.sym - W$GiViTI; W$sym_vs_LR <- W$g.sym - W$Stk.LR; W$symsc_vs_Gi <- W$g.sym.sc - W$GiViTI; W$symsc_vs_LR <- W$g.sym.sc - W$Stk.LR
keep <- c("tag", "link", "n", "g.sym", "g.sym.sc", "g.max3", "g.max3.sc", "Stk.LR", "GiViTI", "max3sc_vs_LR", "max3sc_vs_Gi", "max3_vs_LR", "max3_vs_Gi", "sym_vs_LR", "sym_vs_Gi", "symsc_vs_LR", "symsc_vs_Gi")
print(W[order(W$link, W$tag, W$n), keep], digits = 3, row.names = FALSE)
cat("\nrange max3.sc - Stk.LR:", range(W$max3sc_vs_LR, na.rm = TRUE), " | max3.sc - GiViTI:", range(W$max3sc_vs_Gi, na.rm = TRUE), "\n")
cat("range max3 (unit) - Stk.LR:", range(W$max3_vs_LR, na.rm = TRUE), " | max3 (unit) - GiViTI:", range(W$max3_vs_Gi, na.rm = TRUE), "\n")
P <- read.csv("runL_paired.csv")
chk <- merge(P[P$probe %in% c("g.sym", "g.sym.sc", "g.max3.sc") & P$rival %in% c("Stk.LR", "GiViTI"), c("tag", "link", "n", "probe", "rival", "pow_probe", "pow_rival")],
  allc[, c("tag", "link", "n", "test", "pow_NAasNo")], by.x = c("tag", "link", "n", "probe"), by.y = c("tag", "link", "n", "test"))
cat("analyse_runL.R runL_paired.csv pow_probe vs my recomputation, max abs diff:", max(abs(chk$pow_probe - chk$pow_NAasNo)), "over", nrow(chk), "rows\n")

## ---------------------------------------------------------------- D. seeds: batch-2 cells that reuse batch-1 seeds
cat("\n==== D. seed reuse: identical-seed cells in batch 1 vs batch 2 ====\n")
cmpseed <- function(t1, t2, n, cols) for (lk in unique(PL[[t1]]$link)) {
  a <- PL[[t1]][PL[[t1]]$link == lk & PL[[t1]]$n == n, ]; b <- PL[[t2]][PL[[t2]]$link == lk & PL[[t2]]$n == n, ]
  if (!nrow(a) || !nrow(b)) next
  for (cn in intersect(cols, intersect(names(a), names(b)))) { eq <- which(abs(a[[cn]] - b[[cn]]) < 1e-10)
    cat(sprintf("  %-15s vs %-16s %-7s n=%-5d %-10s identical reps %4d / %4d  first: %s | rej@.05 %.3f vs %.3f\n", t1, t2, lk, n, cn, length(eq), nrow(a),
      paste(head(eq, 8), collapse = ","), mean(a[[cn]] <= .05, na.rm = TRUE), mean(b[[cn]] <= .05, na.rm = TRUE))) } }
cmpseed("s2_probit", "s2_probit_m25sc", 4000, c("g.sym", "g.max3", "u.sym", "Stk.LR", "GiViTI"))
cmpseed("s2_probit", "s2_probit_m25sc", 6500, c("g.sym", "u.sym", "Stk.LR", "GiViTI"))
cmpseed("s2_probit", "s2_probit_m10", 6500, c("u.sym", "Stk.LR"))
cmpseed("base_asym", "base_asym_m10", 1000, c("u.ao", "Stk.LR", "GiViTI"))
cmpseed("base_cauchit", "base_cauchit_m10", 2500, c("u.sym", "Stk.LR"))
cmpseed("s2_probit_m25sc", "s2_probit_m10", 6500, c("u.sym", "Stk.LR"))
cat("-- within s2_probit_m25sc (same seeds and same G as s2_probit): paired power difference g.sym, size-adjusted --\n")
for (n in c(4000, 6500)) { a <- cell("s2_probit", "probit", n, c("g.sym", "u.sym", "Stk.LR", "GiViTI"), FALSE)$tab; b <- cell("s2_probit_m25sc", "probit", n, c("g.sym", "u.sym", "Stk.LR", "GiViTI"), FALSE)$tab
  cat(sprintf("  n=%d  batch1 %s | batch2 %s  (MC SE of a difference ~ %.3f)\n", n, paste(sprintf("%s %.3f", a$test, a$pow_NAasNo), collapse = " "),
    paste(sprintf("%.3f", b$pow_NAasNo), collapse = " "), sqrt(2 * .25 / 1000))) }
cat("-- independence null vs alternative within a file: cor of u.sym p by rep index --\n")
for (tg in c("base_probit", "s2_probit_m25sc", "base_cauchit")) { Q <- PL[[tg]]; n <- max(Q$n); lk <- setdiff(unique(Q$link), "logit")[1]
  cat(sprintf("  %s n=%d cor(null, alt) = %.3f ; identical p-values across the two cells: %d\n", tg, n,
    cor(Q[Q$link == "logit" & Q$n == n, "u.sym"], Q[Q$link == lk & Q$n == n, "u.sym"], use = "complete"),
    length(intersect(Q[Q$link == "logit" & Q$n == n, "Stk.LR"], Q[Q$link == lk & Q$n == n, "Stk.LR"])))) }

## ---------------------------------------------------------------- E. scout envelope
cat("\n==== E. scout_envelope.csv ====\n")
SC <- read.csv("scout_envelope.csv", stringsAsFactors = FALSE)
cat("max rho over all finite rows:", max(SC$rho, na.rm = TRUE), " rows with rho > 1 + 1e-6:", sum(SC$rho > 1 + 1e-6, na.rm = TRUE), "\n")
lam80 <- sapply(1:3, function(k) uniroot(function(l) pchisq(qchisq(.95, k), k, ncp = l, lower.tail = FALSE) - .8, c(.01, 200), tol = 1e-10)$root)
cat("check matched: n80*D vs lambda80(1):", range(SC$n80[SC$test == "matched"] * SC$D[SC$test == "matched"]), lam80[1], "\n")
for (x in list(c(1, 0, "probit"), c(1, -2, "probit"), c(2, 0, "probit"), c(1, 0, "cauchit"), c(2, 0, "cauchit"), c(2, 0, "t4"), c(1, 0, "cloglog"), c(1, 0, "loglog"))) {
  z <- SC[SC$s == as.numeric(x[1]) & SC$c0 == as.numeric(x[2]) & SC$link == x[3], ]
  cat(sprintf("  s=%s c0=%s %-8s event %.3f auc %.3f | %s\n", x[1], x[2], x[3], z$event[1], z$auc[1],
    paste(sprintf("%s n80=%.0f rho=%.3f", z$test, z$n80, z$rho)[z$test %in% c("matched", "sym", "ao", "aoM", "maxprobe", "stukel", "GiViTI", "HL10")], collapse = "; "))) }
pw1 <- function(n, n80, k) pchisq(qchisq(.95, k), k, ncp = n * lam80[k] / n80, lower.tail = FALSE)
b <- SC[SC$s == 1 & SC$c0 == 0 & SC$link == "probit", ]
cat(sprintf("  probit base, local-limit power at n=16000: matched %.3f  sym %.3f  stukel %.3f\n",
  pw1(16000, b$n80[b$test == "matched"], 1), pw1(16000, b$n80[b$test == "sym"], 1), pw1(16000, b$n80[b$test == "stukel"], 2)))

## ---------------------------------------------------------------- F. GiViTI mechanism csv
cat("\n==== F. giviti_mechanism_*.csv ====\n")
for (lk in c("logit", "cauchit")) { M <- read.csv(sprintf("giviti_mechanism_%s.csv", lk))
  cat(sprintf("  %s B=%d | degree table: %s | P(lr3>3.84) %.3f | P(m=2) %.3f | rej %.3f | rej|m=2 %.3f rej|m>=3 %.3f | fixed deg3 LR %.3f | na %d\n", lk, nrow(M),
    paste(names(table(M$m)), table(M$m), sep = ":", collapse = " "), mean(M$lr3 > 3.84), mean(M$m == 2, na.rm = TRUE), mean(M$p <= .05, na.rm = TRUE),
    mean(M$p[M$m == 2] <= .05, na.rm = TRUE), mean(M$p[M$m >= 3] <= .05, na.rm = TRUE), mean(M$lr2 + M$lr3 > qchisq(.95, 2)), sum(is.na(M$m))))
  cat(sprintf("     consistency: m>=3 iff lr3>3.84 holds in %.3f of reps; m=4 iff (lr3>3.84 & lr4>3.84) in %.3f\n",
    mean((M$m >= 3) == (M$lr3 > qchisq(.95, 1)), na.rm = TRUE), mean((M$m == 4) == (M$lr3 > qchisq(.95, 1) & M$lr4 > qchisq(.95, 1)), na.rm = TRUE))) }
