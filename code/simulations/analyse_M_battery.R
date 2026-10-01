## analyse_M_battery.R -- analysis of the EDGE restructure battery as fixed before any result was read in
## paper_EDGE/theory/PREDECLARATION_restructure_battery.md: Section A (unit or score weighting), D (H1-H6) and C3 (paired
## tests), in the reading of E12 and E13.
##
## Usage
##   Rscript analyse_M_battery.R                        root battery, out battery/analysis, rule-G arm, forms chosen by Section A
##   Rscript analyse_M_battery.R --arm 10               the same computations at G = 10 (Supporting Information), suffix _G10
##   Rscript analyse_M_battery.R --forms unit           EDGE forms of H1-H5 fixed (unit or score) instead of the Section A choice
##   Rscript analyse_M_battery.R --mcse realised        sensitivity only: size checks against sqrt(r(1-r)/B) (not E13.1; test roots,
##                                                      or the real root with --force)
##   Rscript analyse_M_battery.R --root battery/dryrun --out battery/_review/analysis/dryrun_out      (a test root)
## On the real root it refuses until battery/launch.log records "launch finished" after its last launch; it refuses a named
## cell pattern of E12.4/E12.5 that matches no cell with a summary (or a family cell with no summary); it refuses when a cell
## of Section A or H1-H5 has no size-adjusted power for a test the run uses (e.g. "matched null not run yet"), listing every
## such cell and test (also in missing_values.csv); and it refuses --mcse realised. --force overrides all four, with a note.
##
## Size checks (E13.1): the size gate of A (E12.3), "holds size" in C1 and E12.6, and the size failures of E13.3 use the nominal
## standard error sqrt(alpha (1 - alpha) / B), B the replicates of the null cell (as run_M_rivals.R does for block 8). The
## standard error of the realised size (the summary's mcse) and its z are written beside it and decide nothing.
## Missing values (E13.2): a mean over a declared set of cells (a family of A, H1-H3, H5) has no value unless every cell has a
## value for every test in it; the macro-average has no value unless all six family means have one, and A then makes no
## decision when both forms are eligible; H4 has no value when a census cell lacks a value of a compared test. The cells are
## listed, and the mean (or count) over the cells with a value is written beside it. A size-adjusted power whose matched null
## gave no p-value in at least 1 - alpha of its replicates counts as no value.
## Size failures (E13.3): each hypothesis row lists the cells where either test fails size in its matched null (H4: EDGE-poly3
## or the best of HL_F, HL, HL_w, PH in that cell, the pair of E13.4), the statistic and verdict without those cells, and a
## final verdict: the E12.5 verdict when the two agree, otherwise "unresolved (size)".
## H4 paired tests (E13.4): EDGE-poly3 against the best partition test, over the detectable, non-saturated cells only; Holm
## over those cells, the rule-G version as its own family.
##
## Reads <root>/<block>/_summary.csv for blocks 0-7, the cell table (<root>/cells.csv, else battery/cells.csv) and, for the
## paired tests, <root>/<block>/<cell>_pvalues.csv.gz.
##
## Outputs in --out (suffix _G10 with --arm 10)
##   rule_A_gate.csv        each EDGE variant (poly3, sym, stk x unit, score) at the arm in every null cell whose summary has it:
##                          size at 0.05 and 0.01, both MCSEs, z and pass under both readings, and the pass in use (E12.3, E13.1)
##   rule_A_membership.csv  the alternative cells of the six families (E12.4)
##   rule_A_families.csv    family means of size-adjusted power per variant and the macro-average (and the means over the cells
##                          with a value, and the cells without one)
##   rule_A_decision.txt    per basis: eligibility, macro-averages, difference, worst family loss, choice, reason
##   hypotheses.csv         H1-H5 (E12.5) in both forms, the forms used marked; H4 at G = 10 with the rule-G version beside it;
##                          size failures, the statistic and verdict without them, and the final verdict (E13.3)
##   h4_census.csv          the H4 census cell by cell (the rule of make_headline_recount.R)
##   declined.csv           H6: declined rate of every test in every cell, overall and within events < G and events >= G
##   mcnemar.csv            exact McNemar per cell of each hypothesis on the raw decisions, Holm within each hypothesis (E12.6, E13.4)
##   missing_values.csv     every cell and test of Section A and H1-H5 without a size-adjusted power, with the reason (E13.2)
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

SIMDIR <- edge_path("code/simulations")
suppressPackageStartupMessages(library(data.table))
BAT <- edge_battery()

AN_BLOCKS <- c("0", "1a", "1b", "2", "3", "4", "5", "6", "7")     # the battery of run_M_all.R; block 8 (E9, E11) is analysed on its own
AN_BASES  <- c("poly3", "sym", "stk")                              # the bases of Section A; poly2 is a sensitivity row (B)
AN_FORM   <- c(unit = "u", score = "sc")
ALPHA <- 0.05
TIE   <- 0.014      # the paper's tie margin (A, C3)
LOSS  <- 0.10       # A rule 1: the variant with the higher macro-average may not lose more than this in any family
ZMAX  <- 3          # size gate (A) and C1: size at most alpha + 3 MCSE
EPS   <- 1e-9       # floating-point slack at a declared threshold, so "at least 0.10" holds at 0.10
H1_GAIN <- 0.10; H2_MARGIN <- -0.014; H3_MARGIN <- -0.05; H5_COST <- 0.02
H4_DETECT <- 0.15; H4_SAT <- 0.97; H4_HOLDS <- 19L                 # make_headline_recount.R compares these without slack; so does H4
G10_ONLY <- c("HL_w", "PH", "Tsiatis", "Xie", "PR")                # rivals run at G = 10 only (E0.6)
AN_MCSE <- c("null", "realised")   # MCSE of the size checks: nominal sqrt(alpha (1 - alpha) / B) (E13.1) or, as a sensitivity, sqrt(r (1 - r) / B)
AN_NULL_DECLINED <- "no value: the matched null gave no p-value in at least 1 - alpha of its replicates (E13.2)"

## E12.4: the alternative cells of the six families, by block and cell name
AN_FAM <- data.frame(stringsAsFactors = FALSE,
  family  = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6),
  block   = c("2", "2", "2", "3", "3", "3", "2", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "4", "4", "4", "4"),
  pattern = c("^probit_(base|auc|e12)_n[0-9]+$", "^cauchit_(base|auc|e12)_n[0-9]+$", "^t4_(base|auc|e12)_n[0-9]+$",
              "^link_probit_n[0-9]+$", "^cauchit_n[0-9]+$", "^t4_n[0-9]+$",
              "^loglog_(base|auc|e12)_n[0-9]+$", "^link_cloglog_n[0-9]+$", "^loglog_n[0-9]+$", "^stk_asym_n[0-9]+$",
              "^stk_long_n[0-9]+$", "^stk_short_n[0-9]+$",
              "^quad_[0-9.]+_n[0-9]+$", "^binint_[0-9.]+_n[0-9]+$", "^contint_[0-9.]+_n[0-9]+$",
              "^rough_osc2_n[0-9]+$", "^rough_osc4_n[0-9]+$", "^rough_sawtooth_n[0-9]+$", "^osc4_n[0-9]+$", "^sawtooth_n[0-9]+$",
              "^crossover_n1000$", "^crossover_n2000$"))
AN_FAM_NAME <- c("symmetric tails", "asymmetric links", "Stukel symmetric family", "omitted terms", "rough misfit", "off-index")
## E12.5: the 8 census cells of H4 (with the block 3 cells at n = 1000 of families 2-4) and the cells of H5
AN_H4_CENSUS <- data.frame(block = "3", pattern = "^census_[a-z0-9_]+_n1000$", stringsAsFactors = FALSE)
AN_H5 <- data.frame(set = c("base", "skew"), block = c("2", "4"),
                    pattern = c("^(probit|cauchit)_base_n[0-9]+$", "^(probit|cauchit)_skew_n[0-9]+$"), stringsAsFactors = FALSE)
AN_H4_PART <- function(ver) c(paste0("HLF.", ver), paste0("HL.", ver), "HL_w", "PH")
AN_H4_ALL6 <- function(ver) c(AN_H4_PART(ver), "Tsiatis", "Xie")

## ---- paths, options, launch guard ---------------------------------------------------------------------------------------
an_path <- function(p) {
  if (!grepl("^([A-Za-z]:|/)", p)) p <- file.path(SIMDIR, p)
  sub("/+$", "", normalizePath(p, winslash = "/", mustWork = FALSE))
}
an_is_real <- function(root, bat = BAT) tolower(an_path(root)) == tolower(an_path(bat))

an_opts <- function(a) {
  opt <- list(root = "battery", out = file.path("battery", "analysis"), arm = "rule", forms = "auto", mcse = AN_MCSE[1], force = FALSE)
  i <- 1L
  while (i <= length(a)) {
    key <- sub("^--", "", a[i]); val <- NULL
    if (grepl("=", key)) { val <- sub("^[^=]*=", "", key); key <- sub("=.*$", "", key) }
    if (!key %in% names(opt)) stop("unknown option --", key, call. = FALSE)
    if (key == "force") { opt$force <- TRUE; i <- i + 1L; next }
    if (is.null(val)) { i <- i + 1L; val <- a[i] }
    if (is.na(val)) stop("--", key, " needs a value", call. = FALSE)
    opt[[key]] <- val; i <- i + 1L
  }
  if (!opt$arm %in% c("rule", "10")) stop("--arm must be rule or 10", call. = FALSE)
  if (!opt$forms %in% c("auto", "unit", "score")) stop("--forms must be auto, unit or score", call. = FALSE)
  if (!opt$mcse %in% AN_MCSE) stop("--mcse must be null or realised", call. = FALSE)
  opt
}

## the reading of run_M_rivals.R (E11.4): a "launch finished" line after the last "launch from step" line
an_launch_finished <- function(lg) {
  L <- if (file.exists(lg)) readLines(lg, warn = FALSE) else character(0)
  last_start <- max(c(0L, grep("launch from step", L, fixed = TRUE)))
  any(grepl("launch finished", L[seq_along(L) > last_start], fixed = TRUE))
}
an_guard <- function(root, force = FALSE, bat = BAT) {
  if (!an_is_real(root, bat)) return(invisible(FALSE))
  lg <- file.path(an_path(bat), "launch.log")
  if (!an_launch_finished(lg)) {
    if (!force) stop("refusing the real root: ", lg, " does not record \"launch finished\" after its last launch ",
                     "(pass --force to override)", call. = FALSE)
    cat("note: --force, the real root is analysed before launch.log records \"launch finished\"\n")
  }
  invisible(TRUE)
}

## ---- inputs ----------------------------------------------------------------------------------------------------------------
an_read_summaries <- function(root) {
  L <- list()
  for (b in AN_BLOCKS) {
    f <- file.path(root, b, "_summary.csv")
    if (!file.exists(f)) next
    L[[b]] <- fread(f, colClasses = list(character = c("block", "cell", "role", "null_cell", "subset", "test", "status"),
                                         numeric = c("n", "B", "alpha", "rejection", "mcse", "size_adj_power", "null_size",
                                                     "declined", "rejection_given_p")))
  }
  if (!length(L)) stop("no block summary under ", root, call. = FALSE)
  S <- rbindlist(L, use.names = TRUE, fill = TRUE)
  S[is.na(status), status := ""]
  S[is.na(null_cell), null_cell := ""]
  S
}
## E12.1: rows with subset = all; a rival that was not run (n >= 10,000) has no value rather than a rate of 0
an_values <- function(S) {
  V <- S[subset == "all" & !grepl("^flag\\.", test) & is.finite(alpha)]
  V[, alpha := round(alpha, 6)]
  V[startsWith(status, "not run"), `:=`(rejection = NA_real_, mcse = NA_real_, size_adj_power = NA_real_)]
  V
}
an_get <- function(V, block, cell, test, col, alpha = ALPHA) {
  if (!length(cell)) return(numeric(0))
  q <- data.table(block = as.character(block), cell = as.character(cell), test = as.character(test), alpha = round(alpha, 6))
  V[q, on = c("block", "cell", "test", "alpha"), mult = "first"][[col]]
}
## E14: a test that gives no p-value in at least 1 - alpha of the replicates of BOTH the alternative and its matched null
## cannot run in that cell at all (a covariate-space test on a model with only binary covariates, say). It is left out of
## that cell's comparisons, as make_headline_recount.R does, instead of taking the cell's value away.
an_applicable <- function(V, block, cell, test, alpha = ALPHA) {
  if (!length(cell)) return(logical(0))
  d  <- an_get(V, block, cell, test, "declined", alpha)
  nd <- an_get(V, block, cell, test, "null_declined", alpha)
  st <- an_get(V, block, cell, test, "status", alpha)
  gone <- (is.finite(d) & d >= 1 - alpha - EPS) | (!is.na(st) & st %in% c("no p-value in any replicate", "not run (n >= 10,000)"))
  !(gone & (!is.finite(nd) | nd >= 1 - alpha - EPS))
}

an_cell_table <- function(root) {
  f <- file.path(root, "cells.csv")
  if (!file.exists(f)) f <- edge_battery("cells.csv")
  C <- fread(f, colClasses = list(character = c("block", "cell", "role", "null_block", "null_cell", "design", "link")))
  C <- C[, .(block, cell, role, n, design, link, null_block, null_cell)]
  C[is.na(null_block), null_block := ""]; C[is.na(null_cell), null_cell := ""]
  C
}
an_mean <- function(v) if (any(is.finite(v))) mean(v[is.finite(v)]) else NA_real_
an_max  <- function(v) if (any(is.finite(v))) max(v[is.finite(v)]) else NA_real_
## a declared mean over a set of cells (E12.4 "the plain mean over its cells"): no value unless every cell has one
an_mean_all <- function(v) if (length(v) && all(is.finite(v))) mean(v) else NA_real_
an_max_all  <- function(v) if (length(v) && all(is.finite(v))) max(v) else NA_real_
an_z <- function(r, a, se) ifelse(se > 0, (r - a) / se, ifelse(r > a, Inf, ifelse(r < a, -Inf, 0)))
## E13.2: a size-adjusted power whose matched null (the null of the cell table) gave no p-value in at least 1 - alpha of its
## replicates, at the same test and alpha, has no value; 'null_declined' keeps that rate for every alternative row
an_null_declined <- function(V, C) {
  V <- copy(V)
  m <- match(paste(V$block, V$cell), paste(C$block, C$cell))
  nb <- C$null_block[m]; nc <- C$null_cell[m]
  V[, null_declined := NA_real_]
  i <- which(!V$role %in% "null" & !is.na(nc) & nzchar(nc))
  if (length(i)) set(V, i, "null_declined", an_get(V, nb[i], nc[i], V$test[i], "declined", V$alpha[i]))
  j <- which(is.finite(V$size_adj_power) & is.finite(V$null_declined) & V$null_declined >= 1 - V$alpha - EPS)
  if (length(j)) set(V, j, c("size_adj_power", "status"), list(NA_real_, AN_NULL_DECLINED))
  V
}
## E13.1, C1, E12.6: a test holds size in a cell's matched null when its size at 0.05 is at most 0.05 + 3 sqrt(0.05 x 0.95 / B),
## B the null cell's replicates (the realised MCSE only with the sensitivity reading); no matched null or no size = does not hold
an_holds <- function(V, C, block, cell, test, reading = AN_MCSE[1]) {
  if (!length(cell)) return(logical(0))
  m <- match(paste(block, cell), paste(C$block, C$cell))
  nb <- C$null_block[m]; nc <- C$null_cell[m]
  nb[is.na(nb)] <- ""; nc[is.na(nc)] <- ""
  size <- an_get(V, nb, nc, test, "rejection")
  se <- if (reading == "null") sqrt(ALPHA * (1 - ALPHA) / an_get(V, nb, nc, test, "B")) else an_get(V, nb, nc, test, "mcse")
  nzchar(nc) & is.finite(size) & is.finite(se) & size <= ALPHA + ZMAX * se + EPS
}
an_cells <- function(block, cell) paste(sprintf("%s/%s", block, cell), collapse = " ")

## ---- Section A: families, size gate, decision -------------------------------------------------------------------------------
an_membership <- function(C, have) {
  M <- rbindlist(lapply(seq_len(nrow(AN_FAM)), function(i) {
    z <- C[block == AN_FAM$block[i] & grepl(AN_FAM$pattern[i], cell), .(block, cell, role, n, design)]
    z[, `:=`(family = AN_FAM$family[i], family_name = AN_FAM_NAME[AN_FAM$family[i]], pattern = AN_FAM$pattern[i])]
    z
  }))
  if (!nrow(M)) M <- data.table(block = character(0), cell = character(0), role = character(0), n = numeric(0), design = character(0),
                                family = numeric(0), family_name = character(0), pattern = character(0))
  M[, in_summary := paste(block, cell) %in% have]
  M
}
## every named pattern of E12.4 and E12.5 with the cells it matches in the cell table and among the cells with a summary
an_pattern_check <- function(C, have) {
  P <- rbind(data.table(set = paste("family", AN_FAM$family), block = AN_FAM$block, pattern = AN_FAM$pattern),
             data.table(set = "H4 census", block = AN_H4_CENSUS$block, pattern = AN_H4_CENSUS$pattern),
             data.table(set = paste("H5", AN_H5$set), block = AN_H5$block, pattern = AN_H5$pattern))
  hit <- lapply(seq_len(nrow(P)), function(i) C$block == P$block[i] & grepl(P$pattern[i], C$cell))
  P[, declared := vapply(hit, sum, integer(1))]
  P[, with_summary := vapply(hit, function(h) sum(h & paste(C$block, C$cell) %in% have), integer(1))]
  P
}
## the size-adjusted powers a run needs: the six EDGE variants in every family cell (at the arm, and at rule G, whose choice fixes
## the forms of H1-H5), and the tests of each hypothesis in both forms (H4: EDGE-poly3 and the six rivals, both versions)
an_need <- function(set, cells, tests)
  data.table(set = rep(set, nrow(cells) * length(tests)), block = rep(cells$block, each = length(tests)),
             cell = rep(cells$cell, each = length(tests)), test = rep(tests, nrow(cells)))
an_value_check <- function(V, M, H4c, H5c, arm) {
  ef <- function(b, g) paste("EDGE", b, AN_FORM, g, sep = ".")
  fam <- unlist(lapply(unique(c(arm, "Grule")), function(g) unlist(lapply(AN_BASES, ef, g = g))))
  h4 <- unique(c(ef("poly3", "G10"), ef("poly3", "Grule"), AN_H4_ALL6("G10"), AN_H4_ALL6("Grule")))
  N <- rbind(an_need("Section A families", M, fam),
             an_need("H1", M[family == 1], c(ef("poly3", arm), "GiViTI")),
             an_need("H2", M[family == 1], c(ef("sym", arm), "Stk.joint")),
             an_need("H3", M[family == 2], c(ef("poly3", arm), "GiViTI")),
             an_need("H4", H4c, h4),
             an_need("H5", H5c[set == "base"], c("Stk.sym1", ef("sym", arm))),
             an_need("H5.skew", H5c[set == "skew"], c("Stk.sym1", ef("sym", arm))))
  N[, `:=`(value = an_get(V, block, cell, test, "size_adj_power"), B = an_get(V, block, cell, test, "B"),
           status = an_get(V, block, cell, test, "status"))]
  X <- N[!is.finite(value)]
  X[, reason := ifelse(!is.finite(B), "no summary row", ifelse(!is.na(status) & nzchar(status), status, "no size-adjusted power"))]
  X[, applicable := an_applicable(V, block, cell, test)]                                            # E14
  X[applicable %in% FALSE, reason := paste(reason, "| not applicable: the test cannot run in this cell (E14)")]
  X[, .(set, block, cell, test, reason, applicable)]
}
an_value_notes <- function(X) {                                                                     # E14: reported, never a refusal
  X <- X[applicable %in% FALSE]
  if (!nrow(X)) return(character(0))
  G <- X[, .(k = .N, cells = uniqueN(paste(block, cell)), ex = paste(sprintf("%s %s %s", block, cell, test), collapse = ", ")), by = set]
  sprintf("%s: %d cell-test pairs in %d cells are not applicable and are left out of that cell's comparisons (E14): %s",
          G$set, G$k, G$cells, G$ex)
}
an_value_messages <- function(X) {
  X <- X[applicable %in% TRUE]                                                                      # E14: only the tests that can run
  if (!nrow(X)) return(character(0))
  G <- X[, .(k = .N, cells = uniqueN(paste(block, cell)), ex = paste(sprintf("%s %s %s (%s)", block, cell, test, reason), collapse = ", ")),
         by = set]
  sprintf("%s: %d cell-test pairs in %d cells have no size-adjusted power: %s", G$set, G$k, G$cells, G$ex)     # every pair (E13.2)
}

an_print_membership <- function(M) {
  cat("\nSection A families (E12.4), alternative cells by block and name:\n")
  for (f in 1:6) {
    z <- M[family == f]
    cat(sprintf("  family %d  %-24s %3d cells, %3d with a summary\n", f, AN_FAM_NAME[f], nrow(z), sum(z$in_summary)))
    for (j in which(AN_FAM$family == f)) {
      zz <- z[block == AN_FAM$block[j] & pattern == AN_FAM$pattern[j]]
      cat(sprintf("    block %s  %-34s %2d  %s\n", AN_FAM$block[j], AN_FAM$pattern[j], nrow(zz), paste(zz$cell, collapse = " ")))
    }
  }
}

## E12.3: every null cell, in any block, whose summary holds the variant at the arm; size within 3 MCSE of alpha or below.
## Both readings of the MCSE are computed; 'reading' picks the one that decides (z_05, pass_05, z_01, pass_01, pass).
an_gate <- function(V, arm, reading = AN_MCSE[1]) {
  tests <- paste("EDGE", rep(AN_BASES, each = 2), rep(AN_FORM, 3), arm, sep = ".")
  g <- V[role == "null" & test %in% tests & alpha %in% c(0.01, 0.05) & is.finite(rejection)]
  g[, mcse_null := sqrt(alpha * (1 - alpha) / B)]
  g[, `:=`(z_realised = an_z(rejection, alpha, mcse), z_null = an_z(rejection, alpha, mcse_null),
           ok_realised = is.finite(mcse) & rejection <= alpha + ZMAX * mcse + EPS,
           ok_null = is.finite(mcse_null) & rejection <= alpha + ZMAX * mcse_null + EPS)]
  if (reading == "null") g[, `:=`(z = z_null, ok = ok_null)] else g[, `:=`(z = z_realised, ok = ok_realised)]
  a5 <- g[alpha == 0.05, .(block, cell, n, B, test, size_05 = rejection, z_05 = z, pass_05 = ok, mcse_05 = mcse, mcse_null_05 = mcse_null,
                           z_realised_05 = z_realised, z_null_05 = z_null, pass_realised_05 = ok_realised, pass_null_05 = ok_null)]
  a1 <- g[alpha == 0.01, .(block, cell, test, size_01 = rejection, z_01 = z, pass_01 = ok, mcse_01 = mcse, mcse_null_01 = mcse_null,
                           z_realised_01 = z_realised, z_null_01 = z_null, pass_realised_01 = ok_realised, pass_null_01 = ok_null)]
  W <- merge(a5, a1, by = c("block", "cell", "test"), all = TRUE)
  W[, pass := !(pass_05 %in% FALSE) & !(pass_01 %in% FALSE)]
  W[, `:=`(pass_realised = !(pass_realised_05 %in% FALSE) & !(pass_realised_01 %in% FALSE),
           pass_null = !(pass_null_05 %in% FALSE) & !(pass_null_01 %in% FALSE), mcse_reading = reading)]
  W[, c("basis", "form") := .(sub("^EDGE\\.([a-z0-9]+)\\..*$", "\\1", test), sub("^EDGE\\.[a-z0-9]+\\.([a-z]+)\\..*$", "\\1", test))]
  setcolorder(W, c("basis", "form", "test", "block", "cell", "n", "B", "mcse_reading", "pass"))
  W[order(basis, form, !pass, block, cell)]
}
an_eligibility <- function(W) {
  E <- CJ(basis = AN_BASES, form = unname(AN_FORM))
  E[, test := sprintf("EDGE.%s.%s", basis, form)]
  X <- W[, .(null_cells = .N, failing = sum(!pass), eligible = all(pass),
             failing_cells = paste(sprintf("%s/%s (z05 %.2f, z01 %.2f)", block, cell, z_05, z_01)[!pass], collapse = "; ")),
         by = .(basis, form)]
  E <- merge(E[, .(basis, form)], X, by = c("basis", "form"), all.x = TRUE)
  E[is.na(null_cells), `:=`(null_cells = 0L, failing = 0L, eligible = FALSE, failing_cells = "no null cell holds this variant")]
  E
}
## E13.2: a mean that is published beside the deciding one is still a mean, so it is NA unless every cell has a value
an_avail <- function(v) if (!length(v) || any(!is.finite(v))) NA_real_ else an_mean(v)

an_family_means <- function(V, M, arm) {
  rbindlist(lapply(AN_BASES, function(b) rbindlist(lapply(AN_FORM, function(f) {
    t <- paste("EDGE", b, f, arm, sep = ".")
    p <- an_get(V, M$block, M$cell, rep(t, nrow(M)), "size_adj_power")
    fm <- rbindlist(lapply(1:6, function(k) { ik <- M$family == k; pk <- p[ik]
      data.table(family = as.character(k), family_name = AN_FAM_NAME[k], cells = length(pk), with_value = sum(is.finite(pk)),
                 mean_power = an_mean_all(pk), mean_power_available = an_avail(pk),
                 cells_without_value = an_cells(M$block[ik & !is.finite(p)], M$cell[ik & !is.finite(p)])) }))
    rbind(fm, data.table(family = "macro", family_name = "mean of the six family means", cells = sum(fm$cells), with_value = sum(fm$with_value),
                         mean_power = an_mean_all(fm$mean_power), mean_power_available = an_avail(fm$mean_power_available),
                         cells_without_value = an_cells(M$block[M$family %in% 1:6 & !is.finite(p)], M$cell[M$family %in% 1:6 & !is.finite(p)])))[, `:=`(basis = b, form = f, test = t)][]
  }))))[, .(basis, form, test, family, family_name, cells, with_value, mean_power, mean_power_available, cells_without_value)]
}
## A and E12.3, per basis
an_decide <- function(E, Fm) {
  rbindlist(lapply(AN_BASES, function(b) {
    fu <- Fm[basis == b & form == "u" & family != "macro"][order(as.integer(family))]$mean_power
    fs <- Fm[basis == b & form == "sc" & family != "macro"][order(as.integer(family))]$mean_power
    eu <- E[basis == b & form == "u"]; es <- E[basis == b & form == "sc"]
    ok_u <- isTRUE(eu$eligible); ok_s <- isTRUE(es$eligible)
    mu <- an_mean_all(fu); ms <- an_mean_all(fs); d <- ms - mu   # A: the six family means with equal weight, no value otherwise
    loss_s <- fu - fs; loss_u <- fs - fu                         # what each form loses against the other, per family
    wl_s <- an_max_all(loss_s); wl_u <- an_max_all(loss_u)
    ff_s <- if (is.finite(wl_s)) which.max(loss_s) else NA_integer_
    ff_u <- if (is.finite(wl_u)) which.max(loss_u) else NA_integer_
    inc_u <- which(!is.finite(fu)); inc_s <- which(!is.finite(fs))
    cu <- Fm[basis == b & form == "u" & family == "macro"]$cells_without_value; cs <- Fm[basis == b & form == "sc" & family == "macro"]$cells_without_value
    cells_note <- paste(c(if (length(cu) && nzchar(cu)) paste("unit", cu), if (length(cs) && nzchar(cs)) paste("score", cs)), collapse = "; ")
    if (!ok_u && !ok_s) {
      ch <- "u"; why <- sprintf("neither form passes the size gate (unit fails in %d null cells, score in %d); the unit form is kept (E12.3)", eu$failing, es$failing)
    } else if (ok_u != ok_s) {
      ch <- if (ok_u) "u" else "sc"
      why <- sprintf("only the %s form passes the size gate (the %s form fails in %d null cells); it is chosen (A eligibility, E12.3)",
                     if (ok_u) "unit" else "score", if (ok_u) "score" else "unit", if (ok_u) es$failing else eu$failing)
    } else if (!is.finite(d)) {
      ch <- NA_character_
      why <- sprintf("no decision: both forms pass the size gate and A needs all six family means in both forms, but %s (a family cell has no size-adjusted power)",
                     paste(c(if (length(inc_u)) sprintf("the unit form has no mean in family %s", paste(inc_u, collapse = ", ")),
                             if (length(inc_s)) sprintf("the score form has no mean in family %s", paste(inc_s, collapse = ", "))), collapse = " and "))
      why <- sprintf("%s; cells without a size-adjusted power: %s (E13.2)", why, cells_note)
    } else if (d >= TIE - EPS) {
      if (wl_s > LOSS + EPS) {
        ch <- "u"; why <- sprintf("score higher by %.4f (at least %.3f) but loses %.4f (more than %.2f) in family %d (%s); the unit form is chosen (A rule 1)",
                                  d, TIE, wl_s, LOSS, ff_s, AN_FAM_NAME[ff_s])
      } else {
        ch <- "sc"; why <- sprintf("score higher by %.4f (at least %.3f) and loses no family by more than %.2f (largest loss %.4f); the score form is chosen (A rule 1)",
                                   d, TIE, LOSS, wl_s)
      }
    } else if (-d >= TIE - EPS) {
      if (wl_u > LOSS + EPS) {
        ch <- "sc"; why <- sprintf("unit higher by %.4f (at least %.3f) but loses %.4f (more than %.2f) in family %d (%s); the score form is chosen (A rule 1)",
                                   -d, TIE, wl_u, LOSS, ff_u, AN_FAM_NAME[ff_u])
      } else {
        ch <- "u"; why <- sprintf("unit higher by %.4f (at least %.3f) and loses no family by more than %.2f (largest loss %.4f); the unit form is chosen (A rule 1)",
                                  -d, TIE, LOSS, wl_u)
      }
    } else {
      ch <- "u"; why <- sprintf("difference %+.4f is below the tie margin %.3f; the unit form is kept (A rule 2)", d, TIE)
    }
    miss <- which(!is.finite(fu) | !is.finite(fs))
    if (length(miss) && !is.na(ch)) why <- paste0(why, sprintf(" [families with a cell without a size-adjusted power: %s; cells: %s]", paste(miss, collapse = ", "), cells_note))
    data.table(basis = b, eligible_unit = ok_u, eligible_score = ok_s, failing_unit = eu$failing, failing_score = es$failing,
               macro_unit = mu, macro_score = ms, difference = d, worst_loss_score = wl_s, worst_loss_score_family = ff_s,
               worst_loss_unit = wl_u, worst_loss_unit_family = ff_u, choice = if (is.na(ch)) "no decision" else if (ch == "u") "unit" else "score",
               choice_code = ch, reason = why)
  }))
}
an_write_decision <- function(D, E, Fm, file, arm, root, forms_note, reading = AN_MCSE[1]) {
  num <- function(x) if (is.finite(x)) sprintf("%.4f", x) else "NA"
  L <- c(sprintf("Section A weighting decision, arm %s", arm), sprintf("root: %s", root),
         sprintf("written: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
         "rule: PREDECLARATION_restructure_battery.md A and E12.3 (size gate alpha + 3 MCSE at 0.05 and 0.01 in every null cell;",
         "      tie margin 0.014 on the macro-average; no family loss above 0.10; one eligible form is chosen; none -> unit;",
         "      a family with a cell without a size-adjusted power has no mean, and then no decision when both forms are eligible (E13.2))",
         sprintf("MCSE of the size gate: %s", if (reading == "null") "null, sqrt(alpha (1 - alpha) / B), the nominal standard error of E13.1 (the realised one is in rule_A_gate.csv and decides nothing)"
                 else "realised, sqrt(r (1 - r) / B) (the summary's mcse): a SENSITIVITY reading; E13.1 decides with sqrt(alpha (1 - alpha) / B)"), "")
  for (i in seq_len(nrow(D))) {
    b <- D$basis[i]
    fu <- Fm[basis == b & form == "u"]; fs <- Fm[basis == b & form == "sc"]
    eu <- E[basis == b & form == "u"]; es <- E[basis == b & form == "sc"]
    L <- c(L, sprintf("basis %s  (EDGE.%s.u.%s against EDGE.%s.sc.%s)", b, b, arm, b, arm),
           sprintf("  eligible        unit %s (%d of %d null cells fail)   score %s (%d of %d)",
                   if (D$eligible_unit[i]) "yes" else "no", eu$failing, eu$null_cells,
                   if (D$eligible_score[i]) "yes" else "no", es$failing, es$null_cells))
    if (nzchar(eu$failing_cells) && eu$failing > 0) L <- c(L, paste("    unit fails:  ", eu$failing_cells))
    if (nzchar(es$failing_cells) && es$failing > 0) L <- c(L, paste("    score fails: ", es$failing_cells))
    for (k in as.character(1:6))
      L <- c(L, sprintf("  family %s %-24s unit %s  score %s  (cells %d, with a value %d / %d)", k, AN_FAM_NAME[as.integer(k)],
                        num(fu[family == k]$mean_power), num(fs[family == k]$mean_power), fu[family == k]$cells,
                        fu[family == k]$with_value, fs[family == k]$with_value),
             if (nzchar(fu[family == k]$cells_without_value) || nzchar(fs[family == k]$cells_without_value))
               sprintf("    without a size-adjusted power: unit [%s]  score [%s]", fu[family == k]$cells_without_value, fs[family == k]$cells_without_value))
    L <- c(L, sprintf("  macro-average   unit %s  score %s  difference (score - unit) %s", num(D$macro_unit[i]), num(D$macro_score[i]),
                      if (is.finite(D$difference[i])) sprintf("%+.4f", D$difference[i]) else "NA"),
           sprintf("  worst family loss  score against unit %s (family %s)   unit against score %s (family %s)",
                   num(D$worst_loss_score[i]), D$worst_loss_score_family[i], num(D$worst_loss_unit[i]), D$worst_loss_unit_family[i]),
           sprintf("  choice          %s", D$choice[i]), sprintf("  reason          %s", D$reason[i]), "")
  }
  L <- c(L, forms_note, "", sprintf("choice %s %s", D$basis, D$choice))
  writeLines(L, file)
  invisible(L)
}

## ---- H1-H5 (E12.5) --------------------------------------------------------------------------------------------------------
an_h1_design <- function(block, cell) ifelse(block == "2", sub("^[a-z0-9]+_(base|auc|e12)_n[0-9]+$", "\\1", cell), "base")
an_pair <- function(V, cells, test_a, test_b) {
  k <- nrow(cells)
  pa <- an_get(V, cells$block, cells$cell, rep(test_a, k), "size_adj_power")
  pb <- an_get(V, cells$block, cells$cell, rep(test_b, k), "size_adj_power")
  data.table(block = cells$block, cell = cells$cell, n = cells$n, test_a = rep(test_a, k), test_b = rep(test_b, k),
             power_a = pa, power_b = pb, diff = pa - pb)
}
## the census rule of make_headline_recount.R, expression for expression: detectable = some compared test reaches 0.15;
## saturated = EDGE-poly3 and the four partition tests all reach 0.97; EDGE-poly3 leads or ties when it exceeds the best
## partition test minus 0.014 (strictly, without slack, as there). Rival columns are named by test.
an_census <- function(W, edge, part, all6) {
  W <- as.data.frame(W, stringsAsFactors = FALSE)
  if (!nrow(W)) return(data.table(best_partition = numeric(0), best_partition_name = character(0), best_all6 = numeric(0),
                                  best_all6_name = character(0), detectable = logical(0), saturated = logical(0),
                                  beats_partition = logical(0), beats_all6 = logical(0), in_headline = logical(0)))
  top <- function(v) if (any(!is.na(v))) names(which.max(v)) else NA_character_
  big <- function(v) if (any(!is.na(v))) max(v, na.rm = TRUE) else NA_real_
  low <- function(v) if (any(!is.na(v))) min(v, na.rm = TRUE) else NA_real_
  K <- rbindlist(lapply(seq_len(nrow(W)), function(i) {
    p3 <- W[[edge]][i]; rp <- unlist(W[i, part]); r6 <- unlist(W[i, all6]); cmp <- c(p3, rp)
    data.table(best_partition = big(rp), best_partition_name = top(rp), best_all6 = big(r6), best_all6_name = top(r6),
               detectable = big(c(p3, r6)) >= H4_DETECT, saturated = low(cmp) >= H4_SAT,
               beats_partition = p3 > big(rp) - TIE, beats_all6 = p3 > big(r6) - TIE)
  }))
  K[, in_headline := detectable & !saturated]
  K
}
an_h4 <- function(V, cells, f, ver) {
  edge <- paste("EDGE.poly3", f, ver, sep = ".")
  W <- data.frame(block = cells$block, cell = cells$cell, n = cells$n, stringsAsFactors = FALSE)
  for (t in c(edge, AN_H4_ALL6(ver))) W[[t]] <- an_get(V, W$block, W$cell, rep(t, nrow(W)), "size_adj_power")
  K <- an_census(W, edge, AN_H4_PART(ver), AN_H4_ALL6(ver))
  R <- as.data.table(W[, AN_H4_ALL6(ver), drop = FALSE])
  setnames(R, c("HLF", "HL", "HL_w", "PH", "Tsiatis", "Xie"))                 # the arm is in 'version'; best_*_name keep the test names
  tn <- c(edge, AN_H4_ALL6(ver)); sn <- c("EDGE.poly3", "HLF", "HL", "HL_w", "PH", "Tsiatis", "Xie")
  AP <- as.data.table(lapply(tn, function(t) an_applicable(V, W$block, W$cell, rep(t, nrow(W)))))    # E14
  setnames(AP, paste0("ap_", sn))
  na_names <- vapply(seq_len(nrow(AP)), function(i) paste(sn[!unlist(AP[i, ])], collapse = " "), "")
  miss_ap <- vapply(seq_len(nrow(W)), function(i) any(unlist(AP[i, ]) & !is.finite(unlist(W[i, tn]))), logical(1))
  cbind(data.table(block = W$block, cell = W$cell, n = W$n, test_edge = rep(edge, nrow(W)), edge_poly3 = W[[edge]]), R, K,
        data.table(not_applicable = na_names, incomplete = miss_ap))
}
an_hrow <- function(hypothesis, arm, form, used, test_a, test_b, cells, with_value, statistic, rule, threshold, verdict,
                    mean_base = NA_real_, mean_auc = NA_real_, mean_e12 = NA_real_, designs_below = "",
                    headline_cells = NA_integer_, led_or_tied = NA_integer_, share = NA_real_, note = "", statistic_available = NA_real_,
                    cells_without_value = "", size = NULL) {
  if (with_value < cells) note <- paste0(note, if (nzchar(note)) "; ", sprintf("%d of %d cells have no value", cells - with_value, cells))
  if (is.null(size)) size <- list(k = NA_integer_, cells = "", statistic = NA_real_, verdict = NA_character_)
  data.table(hypothesis, arm, form, forms_used = used, test_a, test_b, cells, cells_with_value = with_value, statistic,
             statistic_available, rule, threshold, verdict, size_fail_cells = size$k, statistic_without_size_failures = size$statistic,
             verdict_without_size_failures = size$verdict, final_verdict = an_final(verdict, size$verdict),
             mean_base, mean_auc, mean_e12, designs_below_margin = designs_below, headline_cells, led_or_tied, share, note,
             cells_without_value, size_fail_cell_names = size$cells)
}
## E13.3: the E12.5 verdict when the verdict without the size-failing cells agrees with it, otherwise "unresolved (size)"; an item
## with no value (E13.2) and a reported row keep their verdict
an_final <- function(v, vw) if (is.na(vw) || v %in% c("no value", "reported") || identical(v, vw)) v else "unresolved (size)"

an_hypotheses <- function(V, M, H4c, H5c, arm, forms_used, C = an_cell_table(BAT), reading = AN_MCSE[1]) {
  H <- list(); K <- list(); PR <- list()
  pairs <- function(x, hyp, a, fn, used, in_head = NA)
    PR[[length(PR) + 1L]] <<- cbind(data.table(hypothesis = hyp, arm = a, form = fn, forms_used = used),
                                    x[, .(block, cell, n, test_a, test_b, power_a, power_b, diff)], in_headline = in_head)
  f1 <- unique(M[family == 1, .(block, cell, n)]); f2 <- unique(M[family == 2, .(block, cell, n)])
  hb <- H5c[set == "base"]; hs <- H5c[set == "skew"]
  vd <- function(st, ok) if (!is.finite(st)) "no value" else if (ok) "holds" else "against"
  ## E13.3: the cells where either test of a pair fails size in its matched null; E13.2: the cells without a value
  sfail <- function(x) !(an_holds(V, C, x$block, x$cell, x$test_a, reading) & an_holds(V, C, x$block, x$cell, x$test_b, reading))
  nov <- function(x) an_cells(x$block[!is.finite(x$diff)], x$cell[!is.finite(x$diff)])
  ## a mean rule (H2, H3, H5): the verdict over all cells and without the size-failing cells
  mrule <- function(x, ok) {
    sf <- sfail(x); st <- an_mean_all(x$diff); sw <- an_mean_all(x$diff[!sf])
    list(st = st, v = vd(st, is.finite(st) && ok(st)),
         size = list(k = sum(sf), cells = an_cells(x$block[sf], x$cell[sf]), statistic = sw, verdict = vd(sw, is.finite(sw) && ok(sw))))
  }
  for (fn in names(AN_FORM)) {
    f <- AN_FORM[[fn]]
    e3 <- paste("EDGE.poly3", f, arm, sep = "."); es <- paste("EDGE.sym", f, arm, sep = ".")
    u3 <- identical(unname(forms_used[["poly3"]]), f); us <- identical(unname(forms_used[["sym"]]), f)

    ## H1: family 1, EDGE-poly3 minus GiViTI; per design (block 3 cells count as base)
    x <- an_pair(V, f1, e3, "GiViTI"); dsg <- an_h1_design(x$block, x$cell); sf <- sfail(x)
    h1 <- function(keep) {
      st <- an_mean_all(x$diff[keep])
      dm <- vapply(c("base", "auc", "e12"), function(d) an_mean_all(x$diff[keep & dsg == d]), numeric(1))
      below <- names(dm)[is.finite(dm) & dm < -TIE - EPS]
      list(st = st, dm = dm, below = below, v = vd(st, is.finite(st) && st >= H1_GAIN - EPS && !length(below)))
    }
    a <- h1(rep(TRUE, nrow(x))); w <- h1(!sf); st <- a$st
    H[[length(H) + 1L]] <- an_hrow("H1", arm, fn, u3, e3, "GiViTI", nrow(x), sum(is.finite(x$diff)), st,
      "mean(EDGE-poly3 - GiViTI) >= 0.10 over family 1, and no design mean below -0.014", H1_GAIN, a$v,
      mean_base = a$dm[["base"]], mean_auc = a$dm[["auc"]], mean_e12 = a$dm[["e12"]], designs_below = paste(a$below, collapse = " "),
      note = paste0(sprintf("mean %s 0.10; cells per design: base %d, auc %d, e12 %d",
                            if (!is.finite(st)) "has no value, so no verdict on" else if (st >= H1_GAIN - EPS) "reaches" else "misses",
                            sum(dsg == "base"), sum(dsg == "auc"), sum(dsg == "e12")),
                    if (any(sf)) sprintf("; without the %d size-failing cells the designs below -0.014 are: %s", sum(sf),
                                         if (length(w$below)) paste(w$below, collapse = " ") else "none")),
      statistic_available = an_avail(x$diff), cells_without_value = nov(x),
      size = list(k = sum(sf), cells = an_cells(x$block[sf], x$cell[sf]), statistic = w$st, verdict = w$v))
    pairs(x, "H1", arm, fn, u3)

    ## H2: family 1, EDGE-sym minus Stukel's joint score
    x <- an_pair(V, f1, es, "Stk.joint"); r <- mrule(x, function(s) s >= H2_MARGIN - EPS)
    H[[length(H) + 1L]] <- an_hrow("H2", arm, fn, us, es, "Stk.joint", nrow(x), sum(is.finite(x$diff)), r$st,
      "mean(EDGE-sym - Stk.joint) >= -0.014 over family 1", H2_MARGIN, r$v,
      statistic_available = an_avail(x$diff), cells_without_value = nov(x), size = r$size)
    pairs(x, "H2", arm, fn, us)

    ## H3: family 2, EDGE-poly3 minus GiViTI
    x <- an_pair(V, f2, e3, "GiViTI"); r <- mrule(x, function(s) s >= H3_MARGIN - EPS)
    H[[length(H) + 1L]] <- an_hrow("H3", arm, fn, u3, e3, "GiViTI", nrow(x), sum(is.finite(x$diff)), r$st,
      "mean(EDGE-poly3 - GiViTI) >= -0.05 over family 2", H3_MARGIN, r$v,
      statistic_available = an_avail(x$diff), cells_without_value = nov(x), size = r$size)
    pairs(x, "H3", arm, fn, u3)

    ## H4: the census at G = 10, the rule-G version beside it. E13.2: no value when a cell lacks a value of a compared test.
    ## E13.3: a size failure = EDGE-poly3 or the best of HL_F, HL, HL_w, PH in that cell (the pair of E13.4) fails size.
    for (ver in c("G10", "Grule")) {
      k <- an_h4(V, H4c, f, ver)
      hyp <- if (ver == "G10") "H4" else "H4.ruleG"
      inc <- k$incomplete                                                    # E14: a test that cannot run does not make the cell incomplete
      sf <- !(an_holds(V, C, k$block, k$cell, k$test_edge, reading) & an_holds(V, C, k$block, k$cell, k$best_partition_name, reading))
      nh <- sum(k$in_headline %in% TRUE); nl <- sum(k$in_headline %in% TRUE & k$beats_partition %in% TRUE)
      nlw <- sum(k$in_headline %in% TRUE & k$beats_partition %in% TRUE & !sf)
      ninc <- sum(inc); novalue <- !nrow(k) || any(inc)
      na6 <- sum(k$in_headline %in% TRUE & k$beats_all6 %in% TRUE)
      v  <- if (novalue) "no value" else if (ver == "G10") (if (nl >= H4_HOLDS) "holds" else "against") else "reported"
      vw <- if (novalue) "no value" else if (ver == "G10") (if (nlw >= H4_HOLDS) "holds" else "against") else "reported"
      H[[length(H) + 1L]] <- an_hrow(hyp, ver, fn, u3, paste("EDGE.poly3", f, ver, sep = "."), paste(AN_H4_PART(ver), collapse = " "),
        nrow(k), sum(!inc), if (novalue) NA_real_ else nl, "at least 19 detectable, non-saturated cells led or tied (census of make_headline_recount.R)",
        H4_HOLDS, v, headline_cells = nh, led_or_tied = nl, share = if (nh > 0) nl / nh else NA_real_,
        note = paste0(sprintf("%sleads or ties all six (with Tsiatis, Xie) in %d", if (ver == "Grule") sprintf("beside H4; 19 reached: %s; ", if (nl >= H4_HOLDS) "yes" else "no") else "", na6),
                      if (ninc > 0) sprintf("; %d of %d cells lack a value of a compared test (the census rule would leave a missing test out; no verdict, E13.2)", ninc, nrow(k)),
                      if (any(sf)) sprintf("; led or tied without the %d size-failing cells: %d", sum(sf), nlw)),
        statistic_available = if (novalue) NA_real_ else nl, cells_without_value = an_cells(k$block[inc], k$cell[inc]),
        size = list(k = sum(sf), cells = an_cells(k$block[sf], k$cell[sf]), statistic = if (novalue) NA_real_ else nlw, verdict = vw))
      K[[length(K) + 1L]] <- cbind(data.table(version = ver, form = fn, forms_used = u3), k)
      px <- k[, .(block, cell, n, test_a = test_edge, test_b = best_partition_name, power_a = edge_poly3, power_b = best_partition,
                  diff = edge_poly3 - best_partition)]
      pairs(px, hyp, ver, fn, u3, if (novalue) rep(FALSE, nrow(k)) else k$in_headline)   # E13.2: no paired tests when H4 has no value
    }

    ## H5: base-design probit and cauchit cells of block 2, Stukel's one-parameter score minus EDGE-sym; skew cells reported
    x <- an_pair(V, hb, "Stk.sym1", es); r <- mrule(x, function(s) s <= H5_COST + EPS); st <- r$st
    H[[length(H) + 1L]] <- an_hrow("H5", arm, fn, us, "Stk.sym1", es, nrow(x), sum(is.finite(x$diff)), st,
      "mean(Stk.sym1 - EDGE-sym) <= 0.02 over the block 2 base probit and cauchit cells", H5_COST, r$v,
      statistic_available = an_avail(x$diff), cells_without_value = nov(x), size = r$size)
    pairs(x, "H5", arm, fn, us)
    y <- an_pair(V, hs, "Stk.sym1", es); sk <- an_mean_all(y$diff); sfy <- sfail(y); skw <- an_mean_all(y$diff[!sfy])
    H[[length(H) + 1L]] <- an_hrow("H5.skew", arm, fn, us, "Stk.sym1", es, nrow(y), sum(is.finite(y$diff)), sk,
      "reported: the same mean over the block 4 skew cells (skewed-covariate cost)", NA_real_, if (is.finite(sk)) "reported" else "no value",
      note = sprintf("skew cost %s the base cost; skew cell n %s the base cell n",
                     if (is.finite(sk) && is.finite(st) && sk > st) "larger than" else "not larger than",
                     if (setequal(y$n, x$n)) "equal to" else "differ from"), statistic_available = an_avail(y$diff), cells_without_value = nov(y),
      size = list(k = sum(sfy), cells = an_cells(y$block[sfy], y$cell[sfy]), statistic = skw, verdict = if (is.finite(skw)) "reported" else "no value"))
    pairs(y, "H5.skew", arm, fn, us)
  }
  list(hyp = rbindlist(H), census = rbindlist(K, fill = TRUE), pairs = rbindlist(PR, fill = TRUE))
}

## ---- C3 / E12.6: exact McNemar on the raw decisions over the shared replicates ------------------------------------------------
an_mcnemar_p <- function(b, c) { m <- b + c; if (m == 0) 1 else min(1, 2 * stats::pbinom(min(b, c), m, 0.5)) }

an_mcnemar <- function(root, V, C, pairs, reading = AN_MCSE[1]) {
  if (!nrow(pairs)) return(data.table())
  X <- merge(pairs, C[, .(block, cell, null_block, null_cell)], by = c("block", "cell"), all.x = TRUE, sort = FALSE)
  X[is.na(null_cell), `:=`(null_block = "", null_cell = "")]
  X[, `:=`(size_a = an_get(V, null_block, null_cell, test_a, "rejection"), mcse_a = an_get(V, null_block, null_cell, test_a, "mcse"),
           null_B_a = an_get(V, null_block, null_cell, test_a, "B"),
           size_b = an_get(V, null_block, null_cell, test_b, "rejection"), mcse_b = an_get(V, null_block, null_cell, test_b, "mcse"),
           null_B_b = an_get(V, null_block, null_cell, test_b, "B"))]
  se_of <- function(mc, Bn) if (reading == "null") sqrt(ALPHA * (1 - ALPHA) / Bn) else mc        # E13.1: nominal (the realised = sensitivity)
  X[, `:=`(limit_a = ALPHA + ZMAX * se_of(mcse_a, null_B_a), limit_b = ALPHA + ZMAX * se_of(mcse_b, null_B_b), mcse_reading = reading)]
  X[, `:=`(holds_a = is.finite(size_a) & is.finite(limit_a) & size_a <= limit_a + EPS,
           holds_b = is.finite(size_b) & is.finite(limit_b) & size_b <= limit_b + EPS)]
  ## reported beside the check and deciding nothing: the nominal and the realised standard error of each null size, with its z
  X[, `:=`(mcse_nominal_a = sqrt(ALPHA * (1 - ALPHA) / null_B_a), mcse_nominal_b = sqrt(ALPHA * (1 - ALPHA) / null_B_b))]
  X[, `:=`(z_nominal_a = an_z(size_a, ALPHA, mcse_nominal_a), z_realised_a = an_z(size_a, ALPHA, mcse_a),
           z_nominal_b = an_z(size_b, ALPHA, mcse_nominal_b), z_realised_b = an_z(size_b, ALPHA, mcse_b))]
  ## E13.4: the H4 paired tests (both versions) run over the detectable, non-saturated cells only
  X[, in_family := !startsWith(hypothesis, "H4") | in_headline %in% TRUE]
  X[, `:=`(B = NA_integer_, raw_rejection_a = NA_real_, raw_rejection_b = NA_real_, both = NA_integer_, a_only = NA_integer_,
           b_only = NA_integer_, neither = NA_integer_, p_mcnemar = NA_real_, note = "")]
  X[, key := paste(block, cell)]
  for (kk in unique(X$key)) {
    ii <- which(X$key == kk)
    f <- file.path(root, X$block[ii[1]], paste0(X$cell[ii[1]], "_pvalues.csv.gz"))
    cols <- unique(stats::na.omit(c(X$test_a[ii], X$test_b[ii])))
    P <- if (file.exists(f)) suppressWarnings(fread(f, select = c("rep", cols))) else NULL
    for (i in ii) {
      why <- character(0)
      if (!nzchar(X$null_cell[i])) why <- c(why, "no matched null")
      else {
        if (!is.finite(X$size_a[i])) why <- c(why, sprintf("no size of %s in %s", X$test_a[i], X$null_cell[i]))
        else if (!X$holds_a[i]) why <- c(why, sprintf("%s fails size (%.4f above %.4f, %s MCSE)", X$test_a[i], X$size_a[i], X$limit_a[i], reading))
        if (is.na(X$test_b[i]) || !is.finite(X$size_b[i])) why <- c(why, sprintf("no size of %s in %s", X$test_b[i], X$null_cell[i]))
        else if (!X$holds_b[i]) why <- c(why, sprintf("%s fails size (%.4f above %.4f, %s MCSE)", X$test_b[i], X$size_b[i], X$limit_b[i], reading))
      }
      if (!X$in_family[i]) why <- c(why, "not detectable or saturated: outside the H4 family of paired tests (E13.4)")
      if (is.null(P)) why <- c(why, "no per-replicate file")
      else if (is.na(X$test_b[i]) || !all(c(X$test_a[i], X$test_b[i]) %in% names(P))) why <- c(why, "test column missing in the per-replicate file")
      else {
        pa <- P[[X$test_a[i]]]; pb <- P[[X$test_b[i]]]
        da <- is.finite(pa) & pa <= ALPHA; db <- is.finite(pb) & pb <= ALPHA          # declined = no rejection (E0.5)
        set(X, i, c("B", "both", "a_only", "b_only", "neither"),
            list(length(pa), sum(da & db), sum(da & !db), sum(!da & db), sum(!da & !db)))
        set(X, i, c("raw_rejection_a", "raw_rejection_b"), list(mean(da), mean(db)))
        if (!length(why)) set(X, i, "p_mcnemar", an_mcnemar_p(sum(da & !db), sum(!da & db)))
      }
      set(X, i, "note", paste(why, collapse = "; "))
    }
  }
  X[, used := is.finite(p_mcnemar)]
  X[, p_holm := NA_real_]
  X[used == TRUE, p_holm := stats::p.adjust(p_mcnemar, method = "holm"), by = .(hypothesis, arm, form)]
  X[, holm_reject := p_holm <= ALPHA]
  X[, key := NULL]
  setcolorder(X, c("hypothesis", "arm", "form", "forms_used", "block", "cell", "n", "test_a", "test_b", "power_a", "power_b", "diff",
                   "in_headline", "in_family", "null_block", "null_cell", "mcse_reading", "size_a", "mcse_nominal_a", "z_nominal_a", "mcse_a", "z_realised_a",
                   "null_B_a", "limit_a", "holds_a", "size_b", "mcse_nominal_b", "z_nominal_b", "mcse_b", "z_realised_b", "null_B_b", "limit_b", "holds_b", "used",
                   "B", "raw_rejection_a", "raw_rejection_b", "both", "a_only", "b_only", "neither", "p_mcnemar", "p_holm", "holm_reject", "note"))
  X
}

## ---- H6 --------------------------------------------------------------------------------------------------------------------
## declined rate of every test in every cell: overall and within events < G / events >= G. The summary splits a cell only when
## both kinds of sample occur; a flag rate of 1 (or 0) means every sample has events < G (or >= G), and the 'all' row is then
## that subset while the other subset holds no sample.
an_declined <- function(S) {
  D <- S[!grepl("^flag\\.", test) & is.finite(alpha) & abs(alpha - ALPHA) < 1e-9,
         .(block, cell, role, n, test, subset, B, declined, rejection_05 = rejection, rejection_given_p_05 = rejection_given_p, status)]
  D[, source := "summary"]
  Fl <- S[grepl("^flag\\.evlt\\.", test) & subset == "all" & rejection %in% c(0, 1),
          .(block, cell, G = sub("flag.evlt.", "", test, fixed = TRUE), rate = rejection)]
  if (nrow(Fl)) {
    have <- unique(D[subset != "all", .(block, cell, subset)])
    one <- Fl[, .(block, cell, subset = paste0(ifelse(rate == 1, "events<", "events>="), G))][!have, on = c("block", "cell", "subset")]
    none <- Fl[, .(block, cell, subset = paste0(ifelse(rate == 1, "events>=", "events<"), G))][!have, on = c("block", "cell", "subset")]
    base <- D[subset == "all"][, subset := NULL]
    a1 <- base[one, on = c("block", "cell"), nomatch = NULL, allow.cartesian = TRUE]
    a1[, source := "all samples are in this subset (flag rate 0 or 1)"]
    a0 <- base[none, on = c("block", "cell"), nomatch = NULL, allow.cartesian = TRUE]
    a0[, `:=`(B = 0, declined = NA_real_, rejection_05 = NA_real_, rejection_given_p_05 = NA_real_, source = "no sample in this subset")]
    D <- rbind(D, a1, a0, use.names = TRUE)
  }
  D[, test_arm := ifelse(grepl("\\.G([0-9]+|rule)$", test), sub("^.*\\.(G([0-9]+|rule))$", "\\1", test), ifelse(test %in% G10_ONLY, "G10", "none"))]
  D[, subset_G := ifelse(subset == "all", "", sub("^events(<|>=)", "", subset))]
  D[, subset_matches_arm := subset == "all" | test_arm == "none" | subset_G == test_arm]
  setcolorder(D, c("block", "cell", "role", "n", "test", "test_arm", "subset", "subset_G", "subset_matches_arm", "B", "declined",
                   "rejection_05", "rejection_given_p_05", "status", "source"))
  D[order(match(block, AN_BLOCKS), cell, test, subset != "all", subset)]
}

## ---- main --------------------------------------------------------------------------------------------------------------------
an_main <- function(args = commandArgs(trailingOnly = TRUE), bat = BAT) {
  opt <- an_opts(args)
  root <- an_path(opt$root); out <- an_path(opt$out)
  real <- an_is_real(root, bat)
  an_guard(root, opt$force, bat)
  if (real && opt$mcse != AN_MCSE[1]) {
    if (!opt$force) stop("refusing the real root: --mcse realised is a sensitivity reading; E13.1 decides every size check with the nominal ",
                         "standard error sqrt(alpha (1 - alpha) / B) (pass --force to override)", call. = FALSE)
    cat("note: --force, the real root is analysed with the realised MCSE (a sensitivity reading, not E13.1)\n")
  }
  arm <- if (opt$arm == "rule") "Grule" else "G10"
  suf <- paste0(if (arm == "Grule") "" else "_G10", if (opt$mcse != AN_MCSE[1]) "_mcse_realised" else "")   # E13.1 files are never overwritten by the sensitivity reading
  S <- an_read_summaries(root); C <- an_cell_table(root); V <- an_null_declined(an_values(S), C)
  have <- unique(paste(S$block, S$cell))
  cat(sprintf("analyse_M_battery.R  %s\nroot %s%s\nout  %s\narm %s, forms %s, MCSE of the size checks %s; blocks read: %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), root, if (real) " (the real root)" else " (a test root)", out, arm, opt$forms, opt$mcse,
              paste(unique(S$block), collapse = " ")))

  M <- an_membership(C, have); PC <- an_pattern_check(C, have)
  an_print_membership(M)
  bad <- PC[with_summary == 0L]; miss <- M[in_summary == FALSE]; notalt <- M[role != "alternative"]
  msg <- c(if (nrow(bad)) sprintf("pattern with no cell: %s, block %s, %s (%d in the cell table)", bad$set, bad$block, bad$pattern, bad$declared),
           if (nrow(miss)) sprintf("%d family cells have no summary: %s", nrow(miss), paste(paste(miss$block, miss$cell), collapse = ", ")),
           if (nrow(notalt)) sprintf("family cells that are not alternatives: %s", paste(paste(notalt$block, notalt$cell), collapse = ", ")))
  H4c <- unique(rbind(M[family %in% 2:4 & block == "3" & n == 1000, .(block, cell, n)],
                      C[block == AN_H4_CENSUS$block & grepl(AN_H4_CENSUS$pattern, cell), .(block, cell, n)]))
  H5c <- rbindlist(lapply(seq_len(nrow(AN_H5)), function(i)
    C[block == AN_H5$block[i] & grepl(AN_H5$pattern[i], cell), .(set = AN_H5$set[i], block, cell, n)]))
  MV <- an_value_check(V, M, H4c, H5c, arm)
  msg <- c(msg, an_value_messages(MV))
  note14 <- an_value_notes(MV)                                                  # E14
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  fwrite(MV, file.path(out, paste0("missing_values", suf, ".csv")))            # E13.2: every cell and test without a value, also on a refusal
  if (length(note14)) cat(paste0("note: ", note14, "
"), sep = "")
  if (length(msg)) {
    if (real && !opt$force) stop("refusing the real root:\n  ", paste(msg, collapse = "\n  "), "\n(pass --force to override)", call. = FALSE)
    cat(paste0("note: ", msg, "\n"), sep = "")
  }

  ## Section A at the arm
  W <- an_gate(V, arm, opt$mcse); E <- an_eligibility(W); Fm <- an_family_means(V, M, arm); D <- an_decide(E, Fm)
  ## the forms of H1-H5: the headline decision (rule G, E12.2) unless fixed on the command line
  Dh <- if (arm == "Grule") D else an_decide(an_eligibility(an_gate(V, "Grule", opt$mcse)), an_family_means(V, M, "Grule"))
  forms_used <- if (opt$forms == "auto") setNames(Dh$choice_code, Dh$basis) else setNames(rep(AN_FORM[[opt$forms]], 3), AN_BASES)
  fname <- function(code) if (is.na(code)) "no decision" else names(AN_FORM)[AN_FORM == code]
  forms_note <- sprintf("forms used for H1-H5: poly3 %s, sym %s (%s%s)", fname(forms_used[["poly3"]]), fname(forms_used[["sym"]]),
                        if (opt$forms != "auto") paste("fixed by --forms", opt$forms)
                        else if (arm == "Grule") "Section A at the rule-G arm"
                        else sprintf("Section A at the rule-G arm, the headline; at G = 10 it would choose poly3 %s, sym %s", D$choice[D$basis == "poly3"], D$choice[D$basis == "sym"]),
                        if (anyNA(forms_used[c("poly3", "sym")])) "; with no decision both forms are reported and neither is marked as used" else "")

  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  fwrite(W, file.path(out, paste0("rule_A_gate", suf, ".csv")))
  fwrite(M, file.path(out, paste0("rule_A_membership", suf, ".csv")))
  fwrite(Fm, file.path(out, paste0("rule_A_families", suf, ".csv")))
  dec <- an_write_decision(D, E, Fm, file.path(out, paste0("rule_A_decision", suf, ".txt")), arm, root, forms_note, opt$mcse)
  cat("\n"); cat(dec, sep = "\n")

  R <- an_hypotheses(V, M, H4c, H5c, arm, forms_used, C, opt$mcse)
  fwrite(R$hyp, file.path(out, paste0("hypotheses", suf, ".csv")))
  fwrite(R$census, file.path(out, paste0("h4_census", suf, ".csv")))
  DC <- an_declined(S)
  fwrite(DC, file.path(out, paste0("declined", suf, ".csv")))
  MC <- an_mcnemar(root, V, C, R$pairs, opt$mcse)
  fwrite(MC, file.path(out, paste0("mcnemar", suf, ".csv")))

  cat("\nHypotheses (forms used):\n")
  hx <- R$hyp[forms_used == TRUE, .(hypothesis, arm, form, cells, cells_with_value, statistic = round(statistic, 4), verdict,
                                    size_fail_cells, without = round(statistic_without_size_failures, 4), verdict_without_size_failures, final_verdict)]
  print(hx, row.names = FALSE)
  if (nrow(MC)) {
    cat("\nPaired tests (forms used): cells, cells in the family (H4: detectable, non-saturated), compared, significant after Holm\n")
    print(MC[forms_used == TRUE, .(cells = .N, in_family = sum(in_family), compared = sum(used), holm_significant = sum(holm_reject %in% TRUE)), by = .(hypothesis, arm)],
          row.names = FALSE)
  }
  cat(sprintf("\nwritten in %s: %s\n", out, paste(c(paste0(c("rule_A_gate", "rule_A_membership", "rule_A_families"), suf, ".csv"),
                                                   paste0("rule_A_decision", suf, ".txt"),
                                                   paste0(c("hypotheses", "h4_census", "declined", "mcnemar", "missing_values"), suf, ".csv")), collapse = ", ")))
  invisible(list(S = S, V = V, C = C, M = M, patterns = PC, H4c = H4c, H5c = H5c, gate = W, eligibility = E, families = Fm,
                 decision = D, decision_rule = Dh, forms_used = forms_used, mcse_reading = opt$mcse, missing_values = MV, hypotheses = R$hyp, census = R$census,
                 pairs = R$pairs, declined = DC, mcnemar = MC))
}

if (sys.nframe() == 0L) an_main()
