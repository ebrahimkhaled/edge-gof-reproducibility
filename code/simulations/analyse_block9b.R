## analyse_block9b.R -- block 9b read exactly as its pre-declaration fixed in advance.
## Contract: paper_EDGE/theory/PREDECLARATION_block9b_slow_rivals_robustness.md (sha256 72d703a9...), sections 2-4.
##   Rule 1: the projection test and BAGofT reject at p < 0.05 (their p-values lie on grids), every other test at p <= 0.05.
##   Rule 2: a test with no p-value counts as no rejection, and that rate is reported.
##   Rule 3: each rival is paired with EDGE-poly3 unit (rule G) on the same replicates, exact McNemar, Holm within 9b.
##   Bound of B9b.1-3: 0.05 + 3 nominal SE at B = 100 = 0.1154.
## BAGofT's p-value is its p.value column ("BAGofT"); p.value3 is stored beside it and carries no claim.
##
##   Rscript analyse_block9b.R              the finished cells only; refuses if any is missing
##   Rscript analyse_block9b.R --provisional  reads .part files too, and says so on every table (a code check, not a result)
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
D <- edge_battery("9b")
suppressMessages({ source(file.path(SIM, "_battery_tests.R")); source(file.path(SIM, "_block9b_rivals.R")) })
PROV <- "--provisional" %in% commandArgs(TRUE)

ALPHA <- 0.05; BOUND <- ALPHA + 3 * sqrt(ALPHA * (1 - ALPHA) / 100)
GRID  <- c("proj", "BAGofT")
EDGE  <- "EDGE.poly3.u.Grule"
REPORTED <- c("EDGE.sym.u.Grule", "Stk.joint", "Stk.sym1", "Stk.LR", "GiViTI", "Cubic.LR", "HL.G10")
TESTS <- c(EDGE, GRID, REPORTED)
rej <- function(p, t) if (t %in% GRID) is.finite(p) & p < ALPHA - 1e-9 else is.finite(p) & p <= ALPHA

C <- b9b_cells()
X <- list()
for (i in seq_len(nrow(C))) {
  f <- file.path(D, sprintf("%s_pvalues.csv.gz", C$cell[i]))
  if (!file.exists(f) && PROV && file.exists(paste0(f, ".part"))) f <- paste0(f, ".part")
  if (!file.exists(f)) { if (PROV) next else stop("block 9b: ", C$cell[i], " is not finished; nothing is read before all three are") }
  X[[C$cell[i]]] <- fread(f)
}
tag <- if (PROV) "  [PROVISIONAL: .part files, not a result]" else ""
cat(sprintf("block 9b%s\ncells read: %s\n", tag, paste(sprintf("%s (%d reps)", names(X), vapply(X, nrow, 1L)), collapse = ", ")))
cat(sprintf("bound for B9b.1-3: 0.05 + 3 SE(B = 100) = %.4f\n", BOUND))

## ---- the rates (rules 1 and 2) --------------------------------------------------------------------------------------
R <- rbindlist(lapply(names(X), function(ce) {
  x <- X[[ce]]
  rbindlist(lapply(TESTS, function(t) {
    r <- rej(x[[t]], t)
    data.table(cell = ce, test = t, reps = nrow(x), rate = mean(r), mcse = sqrt(mean(r) * (1 - mean(r)) / nrow(x)),
               no_p = mean(!is.finite(x[[t]])), rule = if (t %in% GRID) "p < 0.05" else "p <= 0.05")
  }))
}))
R[, cell := factor(cell, levels = C$cell)]
out_dir <- file.path(D, if (PROV) "analysis_provisional" else "analysis")
dir.create(out_dir, showWarnings = FALSE)
fwrite(R, file.path(out_dir, "rates_9b.csv"))
cat("\n===== false-alarm rates (the model is correct for every uncorrupted record) =====\n")
w <- dcast(R, test ~ cell, value.var = "rate")
w <- w[match(TESTS, test)]
print(w[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)], row.names = FALSE)
if (any(R$no_p > 0)) { cat("\nno p-value (counted as no rejection):\n"); print(R[no_p > 0, .(cell, test, no_p)], row.names = FALSE) }

rate <- function(t, ce) R[test == t & cell == ce, rate]
have <- function(ce) ce %in% names(X)
verdict <- function(ok) if (is.na(ok)) "NOT DECIDABLE (cells missing)" else if (ok) "HOLDS" else "FAILS"

cat("\n===== B9b.1: BAGofT's rate is at most", sprintf("%.4f", BOUND), "in all three cells =====\n")
b1 <- if (all(have(C$cell))) all(vapply(C$cell, function(ce) rate("BAGofT", ce) <= BOUND, TRUE)) else NA
cat(sprintf("  %s: %.3f\n", C$cell[have(C$cell)], vapply(C$cell[have(C$cell)], function(ce) rate("BAGofT", ce), 1)), sep = "")
cat("B9b.1 VERDICT:", verdict(b1), tag, "\n")

cat("\n===== B9b.2: the projection test is at most", sprintf("%.4f", BOUND), "in clean and x4, and above it in x8 =====\n")
b2 <- if (all(have(C$cell))) rate("proj", "clean") <= BOUND && rate("proj", "x4") <= BOUND && rate("proj", "x8") > BOUND else NA
cat(sprintf("  %s: %.3f\n", C$cell[have(C$cell)], vapply(C$cell[have(C$cell)], function(ce) rate("proj", ce), 1)), sep = "")
cat("B9b.2 VERDICT:", verdict(b2), tag, "\n")
if (isFALSE(b2)) cat("  which half fails:", if (rate("proj", "clean") > BOUND || rate("proj", "x4") > BOUND) "the robust half (clean/x4)" else "",
                     if (rate("proj", "x8") <= BOUND) "the breaking half (x8 stays within the bound: section 4 puts the projection test in the robust group)" else "", "\n")

cat("\n===== B9b.3: EDGE-poly3 unit (rule G) is at most", sprintf("%.4f", BOUND), "in all three cells =====\n")
b3 <- if (all(have(C$cell))) all(vapply(C$cell, function(ce) rate(EDGE, ce) <= BOUND, TRUE)) else NA
cat(sprintf("  %s: %.3f\n", C$cell[have(C$cell)], vapply(C$cell[have(C$cell)], function(ce) rate(EDGE, ce), 1)), sep = "")
cat("B9b.3 VERDICT:", verdict(b3), tag, "\n")

## ---- rule 3: the pairs ----------------------------------------------------------------------------------------------
P <- rbindlist(lapply(names(X), function(ce) rbindlist(lapply(GRID, function(t) {
  x <- X[[ce]]; r1 <- rej(x[[t]], t); r2 <- rej(x[[EDGE]], EDGE); nb <- sum(r1 & !r2); nc <- sum(!r1 & r2)
  data.table(cell = ce, rival = t, reps = nrow(x), rival_rate = mean(r1), edge_rate = mean(r2), both = sum(r1 & r2),
             rival_only = nb, edge_only = nc, neither = sum(!r1 & !r2),
             mcnemar_p = if (nb + nc == 0) 1 else stats::binom.test(nb, nb + nc, 0.5)$p.value)
}))))
P[, holm_p := stats::p.adjust(mcnemar_p, method = "holm")]
fwrite(P, file.path(out_dir, "paired_9b.csv"))
cat("\n===== rule 3: each rival against EDGE-poly3 unit on the same replicates (exact McNemar, Holm within 9b) =====\n")
print(P[, .(cell, rival, rival_rate = round(rival_rate, 3), edge_rate = round(edge_rate, 3), rival_only, edge_only,
            mcnemar_p = signif(mcnemar_p, 3), holm_p = signif(holm_p, 3))], row.names = FALSE)

## ---- B9b.4 (reported): the fits and the cost ---------------------------------------------------------------------------
## The fitted coefficients were not stored per replicate; each data set is regenerated from its seed (the replicate's
## own path, b9b_gen after set.seed) and kept only if EDGE-poly3 unit at G = 10 recomputed on it matches the stored value
## to 1e-8, as block 8's identity gate. Beside the coefficients: the largest |fitted linear predictor| among the corrupted
## records, the quantity the declaration's cells are described by.
cat("\n===== B9b.4 (reported): fitted coefficients, corrupted records' linear predictors, seconds =====\n")
RNGkind("L'Ecuyer-CMRG")
F <- rbindlist(lapply(names(X), function(ce) {
  cc <- as.list(C[C$cell == ce, ]); x <- X[[ce]]
  rbindlist(lapply(seq_len(nrow(x)), function(j) {
    set.seed(cc$seed_base + x$rep[j])
    g <- b9b_gen(cc$mult)
    fq <- bt_fit(list(d = g$d, f = g$f))
    pe <- tryCatch(bt_edge_arm(fq, 10L, cc)[["EDGE.poly3.u"]], error = function(e) NA_real_)
    eta <- as.numeric(fq$fit$linear.predictors)
    data.table(cell = ce, rep = x$rep[j], id_ok = isTRUE(abs(pe - x$EDGE.poly3.u.G10[j]) <= 1e-8),
               b0 = coef(fq$fit)[[1]], bx = coef(fq$fit)[["x"]], bd = coef(fq$fit)[["d"]],
               max_eta_corrupt = if (length(g$corrupt)) max(abs(eta[g$corrupt])) else NA_real_)
  }))
}))
fwrite(F, file.path(out_dir, "fits_9b.csv"))
cat(sprintf("identity gate: %d of %d regenerated data sets match their stored EDGE p-value to 1e-8\n", sum(F$id_ok), nrow(F)))
print(F[id_ok == TRUE, .(reps = .N, mean_b0 = round(mean(b0), 3), mean_bx = round(mean(bx), 3), mean_bd = round(mean(bd), 3),
                        median_max_eta_corrupt = round(stats::median(max_eta_corrupt), 1),
                        max_max_eta_corrupt = round(max(max_eta_corrupt), 1)), by = cell], row.names = FALSE)
cat("(truth: intercept 0, slope on x 0.6, slope on d 0.5)\n\n")
S <- rbindlist(lapply(names(X), function(ce) {
  x <- X[[ce]]
  data.table(cell = ce, test = c("projection", "BAGofT"),
             median_sec = c(stats::median(x$sec_proj, na.rm = TRUE), stats::median(x$sec_bagoft, na.rm = TRUE)),
             q25 = c(stats::quantile(x$sec_proj, 0.25, na.rm = TRUE), stats::quantile(x$sec_bagoft, 0.25, na.rm = TRUE)),
             q75 = c(stats::quantile(x$sec_proj, 0.75, na.rm = TRUE), stats::quantile(x$sec_bagoft, 0.75, na.rm = TRUE)))
}))
fwrite(S, file.path(out_dir, "seconds_9b.csv"))
print(S[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 1) else z)], row.names = FALSE)
cat("(seconds per data set, one core each, 20 data sets in parallel)\n")
cat("\nwritten:", out_dir, tag, "\n")
