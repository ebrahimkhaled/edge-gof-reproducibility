## analysis_verify2_edge.R -- the four edge cases of the E13 readings, each on its own synthetic root, with analyse_M_battery.R
## and the second implementation compared output by output, and then the expected answer asserted directly:
##   A  a family entirely missing (family 3, no summary row at all; and family 5, to see the note that lists only ten cells)
##   B  one H1 cell without a value for GiViTI only (and the same cell also failing size, so E13.3 may not overwrite "no value")
##   C  a matched null that declined in exactly 0.95 of its replicates (with 0.9499 and 0.99 as controls)
##   D  a size exactly at the nominal 3-MCSE limit, in the size gate and in a paired test
## No real battery result is read.
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

source(file.path(edge_path("code/simulations"), "analysis_verify2_impl.R"))
sink(file.path(V2, "verify2_edge.log"), split = TRUE)
cat("analysis_verify2_edge.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
shift_mat <- function(v) matrix(v, nrow = 4, ncol = 6, byrow = TRUE, dimnames = list(V2_BASES4, NULL))
set.seed(8820)
## the cells the assertions need per-replicate files for (every fread costs seconds while the battery runs): the H1 cell of
## case B, the family 1 cells that share the declined null of case C, and the two family 2 cells of case D
f2D <- FAMS[family == 2]; nf2 <- v2_null(f2D$block, f2D$cell)
kf2 <- paste(nf2$null_block, nf2$null_cell); iD2 <- 1L; iD3 <- which(kf2 != kf2[1])[1]
cD2 <- f2D$cell[iD2]; cD3 <- f2D$cell[iD3]; nD2 <- nf2[iD2]; nD3 <- nf2[iD3]
PFE <- unique(rbind(data.table(block = "2", cell = c("probit_auc_n3500", "probit_base_n9400")),
                    f2D[c(iD2, iD3), .(block, cell)]))
## the edge roots carry blocks 1a-4, where every cell of Section A and H1-H5 lives (both implementations read the same root)
P0 <- function(extra = list()) c(list(noise = 0.02, rho = 0.4, nofile = integer(0), dropcol = integer(0),
                                      blocks = c("1a", "1b", "2", "3", "4"), reps = 120L, pfile = PFE,
                                      k = function(tests, ce) runif(length(tests), 0.5, 0.8),
                                      off = function(tests, ce) 0,
                                      shift = shift_mat(rep(0.02, 24))), extra)
flat <- function(m) gsub("\n", " | ", m)
as_real <- function(dir, extra = character(0)) {                 # the same root read as if it were the battery
  v2_finished(dir)
  v2_run_builder(c("--root", dir, "--out", file.path(dir, "out_real"), extra), bat = dir)
}

## ---- A: a family with no summary at all -------------------------------------------------------------------------------------------
cat("\n---- A: a family entirely missing ----\n")
dA <- file.path(V2, "edgeA_family3")
v2_make(dA, 8820011, P0(list(dropcells = FAMS[family == 3, .(block, cell)])))
oA <- v2_compare(dA, "edgeA_family3", "rule")
SET <- "edgeA_family3 expected"
fmA <- oA$builder$families; hyA <- oA$builder$hypotheses; dcA <- oA$builder$decision
check("A: family 3 has no mean and no cell with a value, in every basis and form",
      all(is.na(fmA[family == "3"]$mean_power)) && all(fmA[family == "3"]$with_value == 0) && all(fmA[family == "3"]$cells == 10),
      sprintf("%d rows", nrow(fmA[family == "3"])))
check("A: the macro-average has no value in every basis and form", all(is.na(fmA[family == "macro"]$mean_power)))
check("A: no basis is decided while both forms pass the size gate",
      all(dcA[eligible_unit & eligible_score]$choice == "no decision") && all(grepl("family 3", dcA[eligible_unit & eligible_score]$reason)),
      paste(dcA$basis, dcA$choice, collapse = "; "))
check("A: H4 and H4.ruleG have no value (two census cells are in family 3), H1, H2, H3 and H5 keep theirs",
      all(hyA[hypothesis %in% c("H4", "H4.ruleG")]$verdict == "no value") &&
        all(hyA[hypothesis %in% c("H4", "H4.ruleG")]$final_verdict == "no value") &&
        all(hyA[hypothesis %in% c("H1", "H2", "H3", "H5")]$verdict != "no value"),
      paste(sprintf("%s %s", hyA$hypothesis, hyA$verdict), collapse = "; "))
check("A: every missing cell-test pair is listed (10 cells x 6 variants of Section A, 2 census cells x 12 tests of H4)",
      nrow(oA$builder$missing_values[set == "Section A families"]) == 60 && nrow(oA$builder$missing_values[set == "H4"]) == 24 &&
        setequal(oA$builder$missing_values[set == "Section A families"]$cell, FAMS[family == 3]$cell),
      sprintf("%d pairs in all", nrow(oA$builder$missing_values)))
rA <- as_real(dA)
check("A: read as the real root it refuses, naming both patterns, every one of the ten cells and no \"...\"",
      grepl("refusing the real root", rA$msg) && grepl("pattern with no cell: family 3", rA$msg) &&
        all(vapply(FAMS[family == 3]$cell, grepl, logical(1), x = rA$msg, fixed = TRUE)) &&
        grepl("Section A families: 60 cell-test pairs in 10 cells", rA$msg, fixed = TRUE), substr(flat(rA$msg), 1, 300))
mvA <- fread(file.path(dA, "out_real", "missing_values.csv"))
check("A: missing_values.csv is written even when it refuses", nrow(mvA) == 84, sprintf("%d rows", nrow(mvA)))
rA2 <- as_real(dA, "--force")
check("A: with --force it runs and every item that needs family 3 has no value", rA2$msg == "ran" &&
        all(is.na(rA2$res$families[family %in% c("3", "macro")]$mean_power)) &&
        all(rA2$res$hypotheses[hypothesis %in% c("H4", "H4.ruleG")]$verdict == "no value"), rA2$msg)

dA5 <- file.path(V2, "edgeA_family5")
v2_make(dA5, 8820012, P0(list(dropcells = FAMS[family == 5, .(block, cell)])))
SET <- "edgeA_family5 expected"
oA5 <- list(builder = v2_run_builder(c("--root", dA5, "--out", file.path(dA5, "out_rule")))$res)   # a light root: assertions only
rA5 <- as_real(dA5)
check("A5: 13 missing family cells -- the refusal lists every pair, but the 'family cells have no summary' note stops at ten",
      grepl("13 family cells have no summary", rA5$msg, fixed = TRUE) && grepl(", ...", rA5$msg, fixed = TRUE) &&
        grepl("Section A families: 78 cell-test pairs in 13 cells", rA5$msg, fixed = TRUE) &&
        all(vapply(FAMS[family == 5]$cell, grepl, logical(1), x = rA5$msg, fixed = TRUE)),
      "finding 3: the note truncates, the pair list does not")
check("A5: family 5 and the macro-average have no mean; H4 keeps its value (family 5 is not in the census)",
      all(is.na(oA5$builder$families[family %in% c("5", "macro")]$mean_power)) &&
        all(oA5$builder$hypotheses[hypothesis == "H4"]$verdict != "no value"))

## ---- B: one H1 cell without GiViTI ---------------------------------------------------------------------------------------------------
cat("\n---- B: one family 1 cell has no value for GiViTI ----\n")
dB <- file.path(V2, "edgeB_h1_givit")
v2_make(dB, 8820021, P0(list(blank = data.table(block = "2", cell = "^probit_auc_n3500$", test = "^GiViTI$"))))
oB <- v2_compare(dB, "edgeB_h1_givit", "rule")
oB10 <- v2_compare(dB, "edgeB_h1_givit", "10")
SET <- "edgeB expected"
hB <- oB$builder$hypotheses
check("B: H1 has no value in both forms and names the cell; H2, H3, H5 keep their values",
      all(hB[hypothesis == "H1"]$verdict == "no value") && all(hB[hypothesis == "H1"]$final_verdict == "no value") &&
        all(hB[hypothesis == "H1"]$cells_without_value == "2/probit_auc_n3500") &&
        all(hB[hypothesis %in% c("H2", "H3", "H5")]$verdict != "no value"),
      paste(sprintf("%s %s", hB$hypothesis, hB$verdict), collapse = "; "))
check("B: the auc design mean has no value while base and e12 keep theirs, and auc is not counted as a design below -0.014",
      all(is.na(hB[hypothesis == "H1"]$mean_auc)) && all(is.finite(hB[hypothesis == "H1"]$mean_base)) &&
        all(is.finite(hB[hypothesis == "H1"]$mean_e12)) && !any(grepl("auc", hB[hypothesis == "H1"]$designs_below_margin)),
      sprintf("base %.4f, auc %s, e12 %.4f", hB[hypothesis == "H1"]$mean_base[1], hB[hypothesis == "H1"]$mean_auc[1],
              hB[hypothesis == "H1"]$mean_e12[1]))
check("B: Section A is untouched -- every family mean and the macro-average have a value, and the bases are decided",
      all(is.finite(oB$builder$families$mean_power)) && all(oB$builder$decision$choice != "no decision"),
      paste(oB$builder$decision$basis, oB$builder$decision$choice, collapse = "; "))
check("B: the same at G = 10 (GiViTI has no arm): H1 has no value there too, Section A still complete",
      all(oB10$builder$hypotheses[hypothesis == "H1"]$verdict == "no value") && all(is.finite(oB10$builder$families$mean_power)))
check("B: exactly one missing pair, in H1", nrow(oB$builder$missing_values) == 1 && oB$builder$missing_values$set == "H1" &&
        oB$builder$missing_values$test == "GiViTI")
rB <- as_real(dB)
check("B: read as the real root it refuses and names H1 only", grepl("refusing the real root", rB$msg) &&
        grepl("H1: 1 cell-test pairs in 1 cells", rB$msg, fixed = TRUE) && !grepl("Section A families", rB$msg, fixed = TRUE) &&
        !grepl("H3:", rB$msg, fixed = TRUE), flat(rB$msg))
mcB <- oB$builder$mcnemar[hypothesis == "H1" & cell == "probit_auc_n3500"]
check("B: the paired test of that cell still runs (E12.6 is per cell; only the mean has no value)", nrow(mcB) == 2 && all(mcB$used),
      sprintf("used %d of %d", sum(mcB$used), nrow(mcB)))

dB2 <- file.path(V2, "edgeB_h1_givit_sizefail")
nB <- v2_null("2", "probit_auc_n3500")
v2_make(dB2, 8820022, P0(list(blank = data.table(block = "2", cell = "^probit_auc_n3500$", test = "^GiViTI$"),
                              post = function(S) v2_place_null(S, nB$null_block, nB$null_cell, "GiViTI", 0.01))))
SET <- "edgeB2 expected"
oB2 <- list(builder = v2_run_builder(c("--root", dB2, "--out", file.path(dB2, "out_rule")))$res)   # a light root: assertions only
hB2 <- oB2$builder$hypotheses[hypothesis == "H1"]
check("B2: the cell also fails size -- E13.3 lists it, but the row keeps \"no value\" as its final verdict (it is not \"unresolved\")",
      all(hB2$verdict == "no value") && all(hB2$final_verdict == "no value") && all(hB2$size_fail_cells >= 1) &&
        all(grepl("2/probit_auc_n3500", hB2$size_fail_cell_names, fixed = TRUE)),
      sprintf("size-failing cells %d, statistic without them %s", hB2$size_fail_cells[1], format(hB2$statistic_without_size_failures[1])))

## ---- C: a matched null that declined in exactly 0.95 of its replicates -----------------------------------------------------------------
cat("\n---- C: a matched null declined at exactly 0.95 ----\n")
nC1 <- v2_null("2", "probit_base_n9400")          # family 1: GiViTI declined at exactly 0.95
nC2 <- v2_null("3", "cauchit_n2000")              # family 1: Stk.joint declined at 0.9499, the control
nC3 <- v2_null("3", "loglog_n500")                # family 2: EDGE.poly3.u.Grule declined at exactly 0.99
dC <- file.path(V2, "edgeC_declined")
v2_make(dC, 8820031, P0(list(declined = rbind(data.table(block = nC1$null_block, cell = nC1$null_cell, test = "GiViTI", value = 0.95),
                                              data.table(block = nC2$null_block, cell = nC2$null_cell, test = "Stk.joint", value = 0.9499),
                                              data.table(block = nC3$null_block, cell = nC3$null_cell, test = "EDGE.poly3.u.Grule", value = 0.99)))))
oC <- v2_compare(dC, "edgeC_declined", "rule")
SET <- "edgeC expected"
VC <- oC$builder$V
shareC <- CT[null_block == nC1$null_block & null_cell == nC1$null_cell & role == "alternative"]
vC <- VC[block == "2" & cell == "probit_base_n9400" & test == "GiViTI"]
check("C: declined exactly 0.95 removes the size-adjusted power at 0.05 and at 0.10, and keeps it at 0.01 (1 - alpha = 0.99)",
      is.na(vC[abs(alpha - 0.05) < 1e-9]$size_adj_power) && is.na(vC[abs(alpha - 0.10) < 1e-9]$size_adj_power) &&
        is.finite(vC[abs(alpha - 0.01) < 1e-9]$size_adj_power) && all(abs(vC$null_declined - 0.95) < 1e-12) &&
        vC[abs(alpha - 0.05) < 1e-9]$status == BE$AN_NULL_DECLINED,
      sprintf("null %s/%s, shared by %d alternative cells", nC1$null_block, nC1$null_cell, nrow(shareC)))
check("C: every alternative cell of that null loses the value, and H1 has no value and names them",
      all(is.na(VC[paste(block, cell) %in% ckey(shareC) & test == "GiViTI" & abs(alpha - 0.05) < 1e-9]$size_adj_power)) &&
        all(oC$builder$hypotheses[hypothesis == "H1"]$verdict == "no value") &&
        all(grepl("probit_base_n9400", oC$builder$hypotheses[hypothesis == "H1"]$cells_without_value, fixed = TRUE)),
      oC$builder$hypotheses[hypothesis == "H1"]$cells_without_value[1])
vC2 <- VC[block == "3" & cell == "cauchit_n2000" & test == "Stk.joint" & abs(alpha - 0.05) < 1e-9]
check("C: 0.9499 is below 1 - alpha and keeps the value (H2 keeps its verdict)", is.finite(vC2$size_adj_power) &&
        abs(vC2$null_declined - 0.9499) < 1e-12 && all(oC$builder$hypotheses[hypothesis == "H2"]$verdict != "no value"))
vC3 <- VC[block == "3" & cell == "loglog_n500" & test == "EDGE.poly3.u.Grule"]
check("C: declined exactly 0.99 removes the value at alpha 0.01 as well as at 0.05 and 0.10",
      all(is.na(vC3$size_adj_power)), sprintf("%d rows, declined %.4f", nrow(vC3), vC3$null_declined[1]))
mcC <- oC$builder$mcnemar[hypothesis == "H1" & cell %in% shareC$cell & form == "unit"]
check("C: NOTE -- the paired tests of those cells still run, although the matched null gave no p-value in 95% of replicates",
      TRUE, sprintf("%d of %d rows compared (the null's own size is near zero, so it 'holds size')", sum(mcC$used), nrow(mcC)))
check("C: the refusal names the declined-null reason for those pairs",
      any(grepl("matched null gave no p-value", oC$builder$missing_values$reason)),
      sprintf("%d of %d missing pairs are from the declined rule",
              sum(grepl("matched null gave no p-value", oC$builder$missing_values$reason)), nrow(oC$builder$missing_values)))

## ---- D: a size exactly at the nominal limit --------------------------------------------------------------------------------------------
cat("\n---- D: a size exactly at the nominal 3-MCSE limit ----\n")
nD <- CT[block == "1b" & role == "null" & B == 2000]$cell[1:3]
lim05 <- 0.05 + 3 * sqrt(0.05 * 0.95 / 2000); lim01 <- 0.01 + 3 * sqrt(0.01 * 0.99 / 2000)
cat(sprintf("D: paired-test nulls %s/%s (exactly at the limit, cell %s) and %s/%s (1e-6 above, cell %s)\n",
            nD2$null_block, nD2$null_cell, cD2, nD3$null_block, nD3$null_cell, cD3))
dD <- file.path(V2, "edgeD_bound")
v2_make(dD, 8820041, P0(list(
  size = rbind(data.table(block = "1b", cell = nD[1], test = "EDGE.poly3.u.Grule", alpha = 0.05, rejection = lim05),
               data.table(block = "1b", cell = nD[2], test = "EDGE.poly3.u.Grule", alpha = 0.05, rejection = lim05 + 1e-6),
               data.table(block = "1b", cell = nD[3], test = "EDGE.sym.u.Grule", alpha = 0.01, rejection = lim01),
               data.table(block = "1b", cell = nD[1], test = "EDGE.stk.u.Grule", alpha = 0.05, rejection = floor(lim05 * 2000) / 2000),
               data.table(block = "1b", cell = nD[2], test = "EDGE.stk.u.Grule", alpha = 0.05, rejection = ceiling(lim05 * 2000) / 2000)),
  post = function(S) { S <- v2_place_null(S, nD2$null_block, nD2$null_cell, "GiViTI", 0)
                       v2_place_null(S, nD3$null_block, nD3$null_cell, "GiViTI", 1e-6) })))
oD <- v2_compare(dD, "edgeD_bound", "rule")
SET <- "edgeD expected"
gD <- oD$builder$gate
g1 <- gD[cell == nD[1] & test == "EDGE.poly3.u.Grule"]; g2 <- gD[cell == nD[2] & test == "EDGE.poly3.u.Grule"]
g3 <- gD[cell == nD[3] & test == "EDGE.sym.u.Grule"]
check("D: a size exactly at 0.05 + 3 sqrt(0.05 x 0.95 / B) passes the gate; 1e-6 above it fails",
      isTRUE(g1$pass_05) && isFALSE(g2$pass_05) && abs(g1$z_05 - 3) < 1e-9 && g2$z_05 > 3,
      sprintf("limit %.10f: z %.9f passes, z %.9f fails", lim05, g1$z_05, g2$z_05))
check("D: the same at alpha 0.01", isTRUE(g3$pass_01) && abs(g3$z_01 - 3) < 1e-9, sprintf("limit %.10f, z %.9f", lim01, g3$z_01))
gs1 <- gD[cell == nD[1] & test == "EDGE.stk.u.Grule"]; gs2 <- gD[cell == nD[2] & test == "EDGE.stk.u.Grule"]
check("D: on the attainable grid (multiples of 1/B) the largest size below the limit passes and the smallest above fails",
      isTRUE(gs1$pass_05) && isFALSE(gs2$pass_05),
      sprintf("%.4f passes, %.4f fails (limit %.6f)", gs1$size_05, gs2$size_05, lim05))
mcD <- oD$builder$mcnemar
d2 <- mcD[hypothesis == "H3" & cell == cD2 & form == "unit"]
d3 <- mcD[hypothesis == "H3" & cell == cD3 & form == "unit"]
check("D: a matched null exactly at the limit holds size in a paired test; 1e-6 above it does not",
      nrow(d2) == 1 && nrow(d3) == 1 && isTRUE(d2$holds_b) && isFALSE(d3$holds_b) && isTRUE(d2$used) && isFALSE(d3$used) &&
        abs(d2$z_nominal_b - 3) < 1e-9, sprintf("z at the limit %.9f, just above %.9f", d2$z_nominal_b, d3$z_nominal_b))
hD <- oD$builder$hypotheses[hypothesis == "H3"]
check("D: E13.3 counts the cell just above the limit as a size failure and the cell at the limit as holding size",
      all(hD$size_fail_cells >= 1) && all(grepl(paste0(f2D$block[iD3], "/", cD3), hD$size_fail_cell_names, fixed = TRUE)) &&
        !any(grepl(paste0(f2D$block[iD2], "/", cD2), hD$size_fail_cell_names, fixed = TRUE)),
      sprintf("size-failing cells: %s", hD$size_fail_cell_names[1]))

v2_write_checks(file.path(V2, "verify2_edge_checks.csv"))
sink()
