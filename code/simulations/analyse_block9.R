## analyse_block9.R -- block 9 read exactly as its pre-declaration fixed in advance.
## Contract: paper_EDGE/theory/PREDECLARATION_block9_contamination.md (sha256 99eb4c27...), sections 2 and 3.
##
## Rules applied here, and nowhere departed from:
##   2.1 under the logistic truth a rejection is a FALSE ALARM, reported with the nominal SE sqrt(.05*.95/B);
##   2.2 under probit and cloglog, power is compared with the clean cell of the same truth and n;
##   2.3 a test with no p-value counts as NO rejection, and its rate is reported;
##   2.4 EDGE (unit) against each ungrouped test, paired inside the cell, exact McNemar, Holm within block 9;
##   2.5 the claims are judged on C1 and C3 only; C2 and C4 are reported without a claim.
## Every test rejects at p <= alpha: block 9 has no Monte Carlo rival whose p-values lie on a grid.
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
SIMDIR <- edge_path("code/simulations")
DIR    <- edge_battery("9")
OUT    <- file.path(DIR, "analysis")
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_block9_contam.R"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ALPHA <- 0.05
EDGE_U <- c("EDGE.poly3.u.Grule", "EDGE.sym.u.Grule")
UNGROUPED <- c("Stk.joint", "Stk.sym1", "Stk.LR", "GiViTI", "Cubic.LR")
DECLARED <- c(EDGE_U, "EDGE.poly3.sc.Grule", "EDGE.sym.sc.Grule",
              "EDGE.poly3.u.G10", "EDGE.sym.u.G10", "EDGE.poly3.sc.G10", "EDGE.sym.sc.G10",
              UNGROUPED, "HL.G10", "HL.Grule", "HLF.G10", "HLF.Grule")

C <- b9_cells()
have <- file.exists(file.path(DIR, paste0(C$cell, "_pvalues.csv.gz")))
if (!all(have)) {
  cat("MISSING", sum(!have), "cells:\n"); print(C$cell[!have]); stop("block 9 is not finished")
}
read_cell <- function(cl) as.data.frame(fread(file.path(DIR, paste0(cl, "_pvalues.csv.gz"))))

## rule 2.3: no p-value counts as no rejection, so NA is FALSE, never dropped
rej <- function(p) is.finite(p) & p <= ALPHA

## ---- 1. one row per cell x test ---------------------------------------------------------------------
rows <- list()
for (i in seq_len(nrow(C))) {
  ce <- C[i, ]; D <- read_cell(ce$cell)
  B  <- nrow(D)
  for (t in DECLARED) {
    if (!t %in% names(D)) next
    r <- rej(D[[t]])
    rows[[length(rows) + 1]] <- data.table(
      cell = ce$cell, truth = ce$truth, n = ce$n, corruption = ce$corruption, rate = ce$rate, k = ce$k,
      test = t, B = B, rejection = mean(r), no_pvalue = mean(!is.finite(D[[t]])),
      mcse_nominal = sqrt(ALPHA * (1 - ALPHA) / B))
  }
  rows[[length(rows) + 1]] <- data.table(
    cell = ce$cell, truth = ce$truth, n = ce$n, corruption = ce$corruption, rate = ce$rate, k = ce$k,
    test = "_fit", B = B, rejection = NA_real_, no_pvalue = NA_real_, mcse_nominal = NA_real_)
}
S <- rbindlist(rows)

## the fitted coefficients and the corrupted-in-top-group counts, one row per cell (B9.4)
fitrows <- rbindlist(lapply(seq_len(nrow(C)), function(i) {
  ce <- C[i, ]; D <- read_cell(ce$cell)
  data.table(cell = ce$cell, truth = ce$truth, n = ce$n, corruption = ce$corruption, rate = ce$rate, k = ce$k,
             b_intercept = mean(D$b.intercept, na.rm = TRUE), b_x = mean(D$b.x, na.rm = TRUE),
             b_d = mean(D$b.d, na.rm = TRUE),
             corrupt_in_top_Grule = mean(D[["corrupt_in_top.Grule"]], na.rm = TRUE),
             corrupt_in_top_G10 = mean(D[["corrupt_in_top.G10"]], na.rm = TRUE),
             errors = sum(D[["flag.b9_error"]] %in% 1))
}))
fwrite(S[test != "_fit"], file.path(OUT, "cell_test.csv"))
fwrite(fitrows, file.path(OUT, "fits.csv"))

clean_of <- function(truth, n) sprintf("%s_clean_n%d", truth, n)
get_rate <- function(cl, t) S[cell == cl & test == t, rejection][1]

## ---- 2. B9.1: false alarms under the logistic truth, C1 at 0.1% to 0.5% ------------------------------
cat("\n=====================  B9.1  false alarms, logistic truth, corrupted covariates  =====================\n")
cat("Claim: EDGE-poly3 (unit) and EDGE-sym (unit) stay within 3 nominal SE of 0.05 in every C1 cell at\n")
cat("rate 0.1-0.5%; at least two of Stk.joint, Stk.sym1, GiViTI, Cubic.LR exceed 0.10 in at least half.\n\n")
B91_RATES <- c(0.001, 0.002, 0.005)
b91 <- S[truth == "logit" & corruption == "C1" & rate %in% B91_RATES & test %in% c(EDGE_U, UNGROUPED)]
lim <- 0.05 + 3 * sqrt(ALPHA * (1 - ALPHA) / 1000)
cat(sprintf("3 nominal SE band at B = 1000: %.4f to %.4f\n\n", 0.05 - 3 * sqrt(ALPHA * .95 / 1000), lim))
W <- dcast(b91, n + rate ~ test, value.var = "rejection")
print(W, row.names = FALSE)
edge_ok <- all(b91[test %in% EDGE_U, abs(rejection - 0.05) <= 3 * mcse_nominal])
ncells  <- nrow(unique(b91[, .(n, rate)]))
over <- b91[test %in% c("Stk.joint", "Stk.sym1", "GiViTI", "Cubic.LR"),
            .(cells_over_0.10 = sum(rejection > 0.10)), by = test]
cat("\nEDGE unit forms inside the 3-SE band in every such cell:", edge_ok, "\n")
cat("cells in the B9.1 set:", ncells, "; half is", ncells / 2, "\n")
print(over, row.names = FALSE)
b91_second <- sum(over$cells_over_0.10 >= ncells / 2) >= 2
cat("at least two ungrouped tests over 0.10 in at least half:", b91_second, "\n")
cat("\nB9.1 VERDICT:", if (edge_ok && b91_second) "HOLDS" else "FAILS", "\n")

## ---- 3. B9.6: C3 moves no test -----------------------------------------------------------------------
cat("\n=====================  B9.6  the unit error (C3) moves no test  =====================\n")
cat("Against the claim: any test whose rate changes by more than 3 nominal SE from its clean cell.\n\n")
b96 <- rbindlist(lapply(seq_len(nrow(C)), function(i) {
  ce <- C[i, ]; if (ce$corruption != "C3") return(NULL)
  cl <- clean_of(ce$truth, ce$n)
  rbindlist(lapply(c(EDGE_U, UNGROUPED), function(t) {
    a <- get_rate(ce$cell, t); b <- get_rate(cl, t)
    data.table(cell = ce$cell, truth = ce$truth, n = ce$n, rate = ce$rate, test = t,
               contaminated = a, clean = b, change = a - b,
               beyond_3se = abs(a - b) > 3 * sqrt(ALPHA * .95 / 1000))
  }))
}))
fwrite(b96, file.path(OUT, "B9_6_c3.csv"))
cat("cells x tests examined:", nrow(b96), "; moved by more than 3 nominal SE:", sum(b96$beyond_3se), "\n")
if (any(b96$beyond_3se)) print(b96[beyond_3se == TRUE][order(-abs(change))], row.names = FALSE)
cat("\nB9.6 VERDICT:", if (!any(b96$beyond_3se)) "HOLDS" else "FAILS", "\n")

## ---- 4. B9.3: power retention at 0.1%, measured above it ---------------------------------------------
cat("\n=====================  B9.3  power against a real misfit  =====================\n")
cat("Claim: at 0.1% corrupted covariates EDGE's unit forms lose at most 0.05 of their clean power\n")
cat("against the probit and cloglog truths. Above 0.1% the loss is measured, without a claim.\n\n")
b93 <- rbindlist(lapply(seq_len(nrow(C)), function(i) {
  ce <- C[i, ]; if (!(ce$truth %in% c("probit", "cloglog") && ce$corruption %in% c("C1", "C3"))) return(NULL)
  cl <- clean_of(ce$truth, ce$n)
  rbindlist(lapply(EDGE_U, function(t) {
    a <- get_rate(ce$cell, t); b <- get_rate(cl, t)
    data.table(truth = ce$truth, n = ce$n, corruption = ce$corruption, rate = ce$rate, test = t,
               clean = b, contaminated = a, loss = b - a)
  }))
}))
fwrite(b93, file.path(OUT, "B9_3_power.csv"))
print(dcast(b93[corruption == "C1"], truth + n + rate ~ test, value.var = "loss"), row.names = FALSE)
claim <- b93[corruption == "C1" & rate == 0.001]
cat("\nat rate 0.001, largest loss:", sprintf("%.3f", max(claim$loss)), "\n")
cat("B9.3 VERDICT:", if (max(claim$loss) <= 0.05) "HOLDS" else "FAILS", "\n")

## ---- 5. B9.4: is the corruption influential for the fit? ----------------------------------------------
cat("\n=====================  B9.4  the corruption is not influential for the fit  =====================\n")
cat("Reported: mean fitted coefficients differ from the clean cell by less than 10%.\n\n")
b94 <- rbindlist(lapply(seq_len(nrow(fitrows)), function(i) {
  f <- fitrows[i]; if (!f$corruption %in% c("C1", "C3")) return(NULL)
  g <- fitrows[cell == clean_of(f$truth, f$n)]
  data.table(cell = f$cell, truth = f$truth, n = f$n, corruption = f$corruption, rate = f$rate,
             d_intercept = (f$b_intercept - g$b_intercept) / abs(g$b_intercept),
             d_x = (f$b_x - g$b_x) / abs(g$b_x), d_d = (f$b_d - g$b_d) / abs(g$b_d))
}))
b94[, worst := pmax(abs(d_intercept), abs(d_x), abs(d_d))]
fwrite(b94, file.path(OUT, "B9_4_fit.csv"))
print(b94[order(-worst)][1:10, .(cell, corruption, rate, d_intercept = round(d_intercept, 3),
                                 d_x = round(d_x, 3), d_d = round(d_d, 3), worst = round(worst, 3))],
      row.names = FALSE)
cat("\ncells where a coefficient moves 10% or more:", sum(b94$worst >= 0.10), "of", nrow(b94), "\n")
cat("B9.4:", if (all(b94$worst < 0.10)) "the corruption is not influential anywhere" else
        "influential in the cells listed above; B9.1 is judged without them (section 3)", "\n")

## ---- 6. rule 2.4: paired tests, EDGE unit against each ungrouped test ---------------------------------
cat("\n=====================  rule 2.4  paired comparisons, Holm within block 9  =====================\n")
pair <- rbindlist(lapply(seq_len(nrow(C)), function(i) {
  ce <- C[i, ]; if (!ce$corruption %in% c("C1", "C3")) return(NULL)
  D <- read_cell(ce$cell)
  rbindlist(lapply(EDGE_U, function(e) rbindlist(lapply(UNGROUPED, function(u) {
    a <- rej(D[[e]]); b <- rej(D[[u]])
    nb <- sum(a & !b); nc <- sum(!a & b); m <- nb + nc
    data.table(cell = ce$cell, truth = ce$truth, n = ce$n, corruption = ce$corruption, rate = ce$rate,
               edge = e, ungrouped = u, edge_rate = mean(a), ungrouped_rate = mean(b),
               edge_only = nb, ungrouped_only = nc,
               p = if (m == 0) 1 else min(1, 2 * pbinom(min(nb, nc), m, 0.5)))
  }))))
}))
pair[, holm_p := p.adjust(p, "holm")]
pair[, holm_reject := holm_p <= ALPHA]
fwrite(pair, file.path(OUT, "paired.csv"))
cat("comparisons in the Holm family:", nrow(pair), "; significant after Holm:", sum(pair$holm_reject), "\n")
cat("\nthe ten largest gaps in favour of the ungrouped test (logistic truth = false alarms):\n")
print(pair[truth == "logit"][order(-(ungrouped_rate - edge_rate))][1:10,
      .(cell, edge, ungrouped, edge_rate, ungrouped_rate, holm_p = signif(holm_p, 3), holm_reject)],
      row.names = FALSE)

## ---- 7. B9.2 and B9.5: reported without a claim -------------------------------------------------------
cat("\n=====================  B9.2  missed events (C2), reported without a claim  =====================\n")
c2 <- S[corruption == "C2" & truth == "logit" & test %in% c(EDGE_U, UNGROUPED)]
print(dcast(c2, n + rate ~ test, value.var = "rejection"), row.names = FALSE)

cat("\n=====================  C4  duplicated high-risk records, reported without a claim  ===============\n")
c4 <- S[corruption == "C4" & test %in% c(EDGE_U, UNGROUPED)]
print(dcast(c4, truth + n ~ test, value.var = "rejection"), row.names = FALSE)

cat("\n=====================  B9.5  the score form, reported  =====================\n")
sc <- S[truth == "logit" & corruption == "C1" &
        test %in% c("EDGE.poly3.u.Grule", "EDGE.poly3.sc.Grule", "EDGE.sym.u.Grule", "EDGE.sym.sc.Grule")]
print(dcast(sc, n + rate ~ test, value.var = "rejection"), row.names = FALSE)

cat("\n=====================  how often a corrupted record reaches the top group  =====================\n")
print(fitrows[corruption %in% c("C1", "C3"),
              .(cell, k, mean_in_top_Grule = round(corrupt_in_top_Grule, 3),
                mean_in_top_G10 = round(corrupt_in_top_G10, 3))][order(cell)], row.names = FALSE)

cat("\nwritten:", OUT, "\n")
