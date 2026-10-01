## analyse_block8L.R -- block 8L read exactly as its pre-declaration fixed in advance.
## Contract: paper_EDGE/theory/PREDECLARATION_block8L_lecessie.md (sha256 0c1a1731...), sections 3 and 4.
## Reads _summary.csv and _paired.csv written by run_M_block8L.R --summary; computes nothing new about the tests.
##
##   Rscript analyse_block8L.R            battery/8L
##   Rscript analyse_block8L.R --test     battery/8L_test (a code check on the test run, not a result)
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
TEST <- "--test" %in% commandArgs(TRUE)
D <- edge_battery(if (TEST) "8L_test" else "8L")
tag <- if (TEST) "  [TEST RUN: a code check, not a result]" else ""
S  <- fread(file.path(D, "_summary.csv"))
PR <- fread(file.path(D, "_paired.csv"))
MEM <- unique(fread(edge_battery("analysis", "rule_A_membership.csv"))[, .(cell, family_name)], by = "cell")
S <- merge(S, MEM, by = "cell", all.x = TRUE, sort = FALSE)
r3 <- function(z) round(z, 3)
verdict <- function(ok) if (is.na(ok)) "NOT DECIDABLE (cells missing)" else if (ok) "HOLDS" else "FAILS"
cat(sprintf("block 8L%s: %d cells in the summary (%d Part A, %d Part B)\n", tag, nrow(S), sum(S$part == "A"), sum(S$part == "B")))

cat("\n===== data: the identity gate and the p-values =====\n")
cat(sprintf("identity failures: %d; largest |EDGE p regenerated - stored| among kept: %.2e; no p-value: %d; degenerate: %d\n",
            sum(S$identity_failures), max(S$max_identity_absdiff, na.rm = TRUE), sum(S$no_pvalue), sum(S$degenerate)))

## ---- L1: level in the ten Part A null cells ----------------------------------------------------------------------
cat("\n===== L1 (level): each Part A null cell within 0.05 +/- 3 SE (B = 500: 0.021 to 0.079) =====\n")
N <- S[part == "A" & role == "null"]
print(N[, .(cell, reps = reps_kept, rate = r3(lecessie_rejection), band = sprintf("%.3f-%.3f", band_lo, band_hi), in_band)], row.names = FALSE)
L1 <- if (nrow(N) == 10L) all(N$in_band) else NA
cat("L1 VERDICT:", verdict(L1), tag, "\n")
if (isFALSE(L1)) cat("  outside the band:", paste(N[in_band == FALSE, cell], collapse = ", "),
                     "-- per L1, the power comparisons of their alternatives carry this beside them\n")

## ---- L2: power, reported ---------------------------------------------------------------------------------------------
cat("\n===== L2 (reported): size-adjusted power against the twenty Part A alternatives, and the pairing with EDGE-poly3 =====\n")
A <- S[part == "A" & role == "alternative"]
A <- merge(A, PR[, .(cell, lecessie_only, edge_only, mcnemar_p, holm_p)], by = "cell", sort = FALSE)
nullok <- setNames(N$in_band, N$cell)
A[, `:=`(lead_or_tie = edge_size_adj > lecessie_size_adj - 0.01,                   # Table 1's lead-or-tie margin
         null_in_band = nullok[null_cell])]
print(A[order(family_name, cell), .(cell, family = family_name, lecessie = r3(lecessie_size_adj), EDGE = r3(edge_size_adj),
        diff = r3(edge_size_adj - lecessie_size_adj), lead_or_tie, lc_only = lecessie_only, edge_only,
        holm_p = signif(holm_p, 2), null_in_band)], row.names = FALSE)
cat(sprintf("\nEDGE-poly3 (unit, rule G) leads or ties le Cessie in %d of %d cells (margin 0.01, size-adjusted)\n",
            sum(A$lead_or_tie), nrow(A)))
cat(sprintf("significant after Holm (raw decisions at 0.05, paired): EDGE ahead in %d, le Cessie ahead in %d, neither in %d\n",
            sum(A$holm_p <= 0.05 & A$edge_only > A$lecessie_only), sum(A$holm_p <= 0.05 & A$lecessie_only > A$edge_only),
            sum(A$holm_p > 0.05)))
fam <- A[, .(cells = .N, edge_leads_or_ties = sum(lead_or_tie), mean_lecessie = r3(mean(lecessie_size_adj)),
             mean_EDGE = r3(mean(edge_size_adj))), by = family_name]
print(fam, row.names = FALSE)

## beside the resampling rivals of block 8, same cells (block 8's own summary; the projection test on the same replicates 1-500)
S8 <- fread(edge_battery("8", "_summary.csv"))[role == "alternative" & alpha == 0.05]
rv <- dcast(S8[test %in% c("proj", "BAGofT")], cell ~ test, value.var = "size_adj_power")
cmp <- merge(A[, .(cell, family_name, lecessie = lecessie_size_adj, EDGE = edge_size_adj)], rv, by = "cell", all.x = TRUE)
cat("\nbeside the two resampling rivals (block 8, size-adjusted; BAGofT on replicates 1-200):\n")
print(cmp[order(family_name, cell), lapply(.SD, function(z) if (is.numeric(z)) r3(z) else z)], row.names = FALSE)
cat(sprintf("mean over the %d cells: le Cessie %.3f, EDGE %.3f, projection %.3f, BAGofT %.3f\n", nrow(cmp),
            mean(cmp$lecessie), mean(cmp$EDGE), mean(cmp$proj, na.rm = TRUE), mean(cmp$BAGofT, na.rm = TRUE)))

## ---- L3: corrupted covariates ------------------------------------------------------------------------------------------
cat("\n===== L3: Part B false-alarm rate of le Cessie above 0.10 at k = 5 and at k = 10 =====\n")
B <- S[part == "B"][order(k)]
print(B[, .(cell, k, reps = reps_kept, lecessie = r3(lecessie_rejection), EDGE = r3(edge_rejection))], row.names = FALSE)
k5 <- B[k == 5, lecessie_rejection]; k10 <- B[k == 10, lecessie_rejection]
L3 <- if (length(k5) == 1L && length(k10) == 1L) k5 > 0.10 && k10 > 0.10 else NA
cat("L3 VERDICT:", verdict(L3), tag, "\n")
if (isFALSE(L3)) cat("  per L3, le Cessie is reported as robust to this corruption beside the grouped tests, in the same place and type\n")

## ---- L4: reported -------------------------------------------------------------------------------------------------------
cat("\n===== L4 (reported) =====\n")
cl <- B[k == 0]
if (nrow(cl)) cat(sprintf("clean cell: %.3f (band at B = 1000: %.3f-%.3f, %s)\n", cl$lecessie_rejection, cl$band_lo, cl$band_hi,
                          if (isTRUE(cl$in_band)) "inside" else "OUTSIDE"))
cat(sprintf("k = 1: %s   k = 2: %s\n", paste(r3(B[k == 1, lecessie_rejection])), paste(r3(B[k == 2, lecessie_rejection]))))
cat(sprintf("median seconds per data set: %.2f (Part A, n = 380-1000), %.2f (Part B, n = 1000)\n",
            stats::median(S[part == "A", median_sec]), stats::median(S[part == "B", median_sec])))
fwrite(A, file.path(D, "analysis_8L_partA.csv")); fwrite(B, file.path(D, "analysis_8L_partB.csv"))
cat("\nwritten:", D, tag, "\n")
