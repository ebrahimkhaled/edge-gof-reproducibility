## analyse_block9d.R -- block 9d read exactly as its pre-declaration fixed in advance.
## Contract: paper_EDGE/theory/PREDECLARATION_block9d_robust.md (sha256 799a85aa...), sections 2-4.
## Rejection at p <= 0.05; no p-value counts as no rejection; 3 nominal SE at B = 1000 is 0.0207.
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
suppressMessages(library(data.table))
SIM <- edge_path("code/simulations")
D9d <- edge_battery("9d"); D9 <- edge_battery("9")
ALPHA <- 0.05; SE3 <- 3 * sqrt(ALPHA * (1 - ALPHA) / 1000)
TESTS <- c("EDGE.poly3.u.Grule", "EDGE.sym.u.Grule", "EDGE.poly3.sc.Grule", "EDGE.sym.sc.Grule",
           "Stk.joint", "Stk.sym1", "GiViTI", "Cubic.LR", "HL.Grule", "HLF.Grule", "Stk.LR")
CLAIM2 <- c("Stk.joint", "Stk.sym1", "GiViTI", "Cubic.LR")      # D9d.2 names these four
rej <- function(p) mean(is.finite(p) & p <= ALPHA)

R <- rbindlist(lapply(c("logit", "probit"), function(tr) rbindlist(lapply(c(0, 5, 10, 25), function(k)
  rbindlist(lapply(c("ML", "robust"), function(est) {
    x <- fread(file.path(D9d, sprintf("%s_k%02d_%s_pvalues.csv.gz", tr, k, est)))
    data.table(truth = tr, k = k, est = est, test = TESTS,
               rate = vapply(TESTS, function(t) rej(x[[t]]), numeric(1)),
               no_p = vapply(TESTS, function(t) mean(!is.finite(x[[t]])), numeric(1)),
               converged = if (est == "robust") mean(x$robust_converged %in% 1) else NA_real_,
               b_x = mean(x$b.x, na.rm = TRUE))
  }))))))
fwrite(R, file.path(D9d, "analysis_9d.csv"))
cat(sprintf("band: 0.05 +/- 3 SE = %.4f to %.4f\n", 0.05 - SE3, 0.05 + SE3))

cat("\n===== rule 4: the ML arm reproduces block 9 (different seeds, so a consistency check) =====\n")
b9 <- fread(file.path(D9, "analysis", "cell_test.csv"))
map <- data.table(k = c(0, 5, 10), cell = c("logit_clean_n1000", "logit_C1_r005_n1000", "logit_C1_r010_n1000"))
chk <- rbindlist(lapply(seq_len(nrow(map)), function(i) {
  a <- R[truth == "logit" & est == "ML" & k == map$k[i] & test %in% CLAIM2]
  b <- b9[cell == map$cell[i] & test %in% CLAIM2, .(test, b9 = rejection)]
  m <- merge(a[, .(test, d9d = rate)], b, by = "test"); m[, k := map$k[i]]; m
}))
chk[, diff := d9d - b9]
print(chk[, .(k, test, block9d = round(d9d, 3), block9 = round(b9, 3), diff = round(diff, 3))], row.names = FALSE)
cat("largest |difference|:", sprintf("%.3f", max(abs(chk$diff))), "  (3 nominal SE on a difference of two independent rates:",
    sprintf("%.3f)\n", 3 * sqrt(2 * ALPHA * (1 - ALPHA) / 1000)))

cat("\n===== D9d.1 (GATE): under a clean logistic truth, every test holds its level after a robust fit =====\n")
g <- R[truth == "logit" & k == 0 & test != "Stk.LR"]
print(dcast(g, test ~ est, value.var = "rate")[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)],
      row.names = FALSE)
d1 <- g[est == "robust", .(test, rate, inside = abs(rate - 0.05) <= SE3)]
cat("\nrobust-fit tests outside the band:", if (all(d1$inside)) "none" else paste(d1[inside == FALSE, test], collapse = ", "), "\n")
D1 <- all(d1$inside)
cat("D9d.1 VERDICT:", if (D1) "HOLDS" else "FAILS", "\n")
if (!D1) cat("(per section 4, the comparisons below then measure the reference mismatch, not the contamination)\n")

cat("\n===== D9d.2: at k = 10, the robust fit repairs the four ungrouped tests =====\n")
t2 <- dcast(R[truth == "logit" & k == 10 & test %in% c(CLAIM2, "EDGE.poly3.u.Grule", "EDGE.sym.u.Grule")],
            test ~ est, value.var = "rate")
t2[, repaired := abs(robust - 0.05) <= SE3]
print(t2[, .(test, ML = round(ML, 3), robust = round(robust, 3), within_band_after_robust = repaired)], row.names = FALSE)
D2 <- all(t2[test %in% CLAIM2, repaired])
cat("D9d.2 VERDICT:", if (D2) "HOLDS -- robust fitting repairs all four" else "FAILS", "\n")

cat("\n===== the whole k grid under the logistic truth, both estimators (false alarms) =====\n")
w <- dcast(R[truth == "logit" & test %in% c(CLAIM2, "EDGE.poly3.u.Grule")], test + est ~ k, value.var = "rate")
print(w[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)], row.names = FALSE)

cat("\n===== D9d.3 (reported): what the robust fit costs in power, probit truth, k = 0 =====\n")
t3 <- dcast(R[truth == "probit" & k == 0 & test != "Stk.LR"], test ~ est, value.var = "rate")
t3[, change := robust - ML]
print(t3[, .(test, ML = round(ML, 3), robust = round(robust, 3), change = round(change, 3))][order(change)],
      row.names = FALSE)

cat("\n===== D9d.4 (reported): EDGE under the robust fit, every k, logistic truth =====\n")
print(dcast(R[truth == "logit" & grepl("^EDGE", test)], test + est ~ k, value.var = "rate")[
      , lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)], row.names = FALSE)

cat("\n===== D9d.5 (reported) =====\n")
cat("robust convergence:", sprintf("%.4f", min(R[est == "robust", converged])), "(minimum over the eight configurations)\n")
print(unique(R[, .(truth, k, est, mean_b_x = round(b_x, 4))])[order(truth, k, est)], row.names = FALSE)
cat("\nDISCLOSED OMISSION: the declaration asks for the Bianco-Yohai estimator (method = 'BY') to be computed\n")
cat("where it converges and reported without a claim. It was NOT computed in this run -- only Mqle was.\n")
cat("No claim depends on it (D9d.1-D9d.4 are Mqle); it is recorded here rather than left silent.\n")
cat("\nStk.LR under the robust arm compares an ML augmented fit with a robust base and is NOT a\n")
cat("likelihood-ratio test; it is excluded from D9d.1 and carries no claim.\n")
