## analyse_selftest.R -- known-answer checks of analyse_M_battery.R. Synthetic block summaries are built from the real cell
## table (battery/cells.csv) with controlled numbers; the paired tests are also checked on the dry run (battery/dryrun, B = 10).
## No value of the real battery is read; on the real root only the launch-log guard is called.
##   1.     the six family memberships (E12.4) = an independent assignment from the generator, link, family and param columns
##   (i)    score +0.02 in every family, both forms eligible -> score
##   (ii)   score +0.02 on average but -0.11 in family 6 -> unit
##   (iii)  difference 0.010 -> unit
##   (iv)   unit fails the gate in one null cell (z = 3.2 against the null MCSE), score eligible -> score, although its power is lower
##   (v)    both forms fail the gate -> unit
##   (vi)   H1, H2, H3, H5 exactly at their thresholds and just beyond; H1 against through one design below -0.014
##   (vii)  H4 census = make_headline_recount.R on the July sim_power_broad rows at n = 1000 and its 8 tab_extra rows
##   (viii) exact McNemar and Holm against a hand computation on the dry-run per-replicate files (both MCSE readings) and on
##          21 synthetic cells
##   (ix)   refusal on the real root before "launch finished" (and on a named pattern with no cell)
##   (x)    MCSE reading: a size between the two limits fails with the null MCSE and passes with --mcse realised (gate and E12.6)
##   (xi)   cells with a summary but no size-adjusted power: refusal on a finished real root, no mean / no decision / no value
##          with --force, per basis, per hypothesis and per arm
##   (xii)  E13: (a) the nominal MCSE decides the gate and the paired-test size check (0.024 at alpha 0.01, B = 500 fails it and
##          passes the realised one; 0.023 passes both), the real root refuses --mcse realised; (b) one family cell without a
##          size-adjusted power -> family, macro-average, H3 and H4 have no value, no basis is decided, every missing pair is
##          listed on a synthetic copy treated as the real root; (c) a matched null declined at 0.96 (and 0.95) -> no value;
##          (d) H2 and H4 verdicts that flip without the size-failing cells -> "unresolved (size)"; (e) the H4 Holm family holds
##          only detectable, non-saturated cells
##   plus   H6 subsets, --arm 10 outputs and --forms
## Writes battery/_review/analysis/selftest.log and the roots battery/_review/analysis/synth_*; exits 1 if a check fails.
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
OUTDIR <- edge_battery("_review", "analysis")
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)                                            # the battery is using the machine
source(file.path(SIMDIR, "analyse_M_battery.R"))           # functions only: the main call runs under Rscript alone
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUTDIR, "selftest.log"), split = TRUE)
cat("analyse_selftest.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
RES <- list()
check <- function(label, ok, detail = "") {
  RES[[length(RES) + 1]] <<- data.frame(check = label, ok = isTRUE(ok), detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "ok" else "FAIL", label, detail))
}
guarded <- function(label, expr) tryCatch(expr, error = function(e) check(label, FALSE, paste("error:", conditionMessage(e))))
near <- function(a, b, tol = 1e-9) length(a) == length(b) && all(is.finite(a) & is.finite(b)) && all(abs(a - b) <= tol)
run_quiet <- function(args, bat = BAT, log = NULL) {
  res <- NULL
  txt <- capture.output(res <- an_main(args, bat = bat))
  if (!is.null(log)) writeLines(txt, log)
  res
}

## ---- synthetic summaries from the real cell table ----------------------------------------------------------------------------
CT <- fread(edge_battery("cells.csv"),
            colClasses = list(character = c("block", "cell", "role", "null_block", "null_cell", "param", "design", "link", "family", "generator")))
## E12.4 assigned from the cell attributes, not from the names the analysis matches
st_family <- function(C) {
  b <- C$block; g <- C$generator; l <- C$link; fa <- C$family; pa <- C$param; alt <- C$role %in% "alternative"
  f <- rep(NA_integer_, nrow(C))
  f[alt & ((b %in% "2" & g %in% "design" & l %in% c("probit", "cauchit", "t4")) |
           (b %in% "3" & g %in% "dgp_alt" & fa %in% "link" & l %in% "probit") |
           (b %in% "3" & g %in% "design" & l %in% c("cauchit", "t4")))] <- 1L
  f[alt & ((b %in% "2" & g %in% "design" & l %in% "loglog") |
           (b %in% "3" & g %in% "dgp_alt" & fa %in% "link" & l %in% "cloglog") |
           (b %in% "3" & g %in% "design" & l %in% c("loglog", "stk_asym")))] <- 2L
  f[alt & b %in% "3" & g %in% "design" & l %in% c("stk_long", "stk_short")] <- 3L
  f[alt & b %in% "3" & g %in% "dgp_alt" & fa %in% c("quad", "binint", "contint")] <- 4L
  f[alt & ((b %in% "3" & g %in% "dgp_alt" & fa %in% "rough" & pa %in% c("osc2", "osc4", "sawtooth")) |
           (b %in% "4" & g %in% "dgp_alt" & fa %in% "rough"))] <- 5L
  f[alt & b %in% "4" & g %in% "crossover"] <- 6L
  f
}
C0 <- CT[block %in% c("1a", "1b", "2", "3", "4", "5", "6", "7") & role %in% c("null", "alternative")]
C0[, fam := st_family(C0)]
ST_EDGE  <- as.vector(outer(paste0("EDGE.", as.vector(outer(c("poly3", "poly2", "stk", "sym"), c("u", "sc"), paste, sep = "."))),
                            c("G10", "Grule"), paste, sep = "."))
ST_TESTS <- c(ST_EDGE, "HL.G10", "HL.Grule", "HLF.G10", "HLF.Grule", "HL_w", "PH", "Tsiatis", "Xie", "PR",
              "Stk.joint", "Stk.LR", "Stk.sym1", "Stk.marg", "GiViTI", "GiViTI.t50", "Cubic.LR")
set.seed(20260914)
BASE <- C0[, .(test = ST_TESTS), by = .(block, cell, role, null_cell, n, B, fam)]
BASE[, u := runif(1), by = .(block, cell)]
BASE[, `:=`(noise = runif(.N, -0.02, 0.02), dec = runif(.N, 0, 0.02), sz = runif(.N, 0.6, 0.9))]
BASE[, form := ifelse(startsWith(test, "EDGE."), sub("^EDGE\\.[a-z0-9]+\\.([a-z]+)\\..*$", "\\1", test), "")]
BASE[, unit_power := ifelse(startsWith(test, "EDGE."), 0.3 + 0.4 * u, 0.25 + 0.4 * u + noise)]
BASE[, flag_k := (match(paste(block, cell), unique(paste(block, cell))) %% 3L) + 1L]

## the size with z = (r - a) / sqrt(r (1 - r) / B) exactly (realised MCSE), and with z = (r - a) / sqrt(a (1 - a) / B) (null MCSE)
st_size_for_z <- function(z, a, B) { k <- z^2 / B; (2 * a + k + sqrt((2 * a + k)^2 - 4 * (1 + k) * a^2)) / (2 * (1 + k)) }
st_size_for_z_null <- function(z, a, B) a + z * sqrt(a * (1 - a) / B)

## delta(fam): score minus unit power of every EDGE basis at both arms in the family cells; gate: overrides of null sizes, by z
## against the null MCSE, or by a 'size' column
st_write_root <- function(dir, delta, gate = NULL) {
  unlink(dir, recursive = TRUE); dir.create(dir, recursive = TRUE)
  file.copy(edge_battery("cells.csv"), file.path(dir, "cells.csv"))
  X <- copy(BASE)
  X[, pw := unit_power + ifelse(form == "sc" & !is.na(fam), delta(fam), 0)]
  R <- rbindlist(lapply(c(0.01, 0.05, 0.10), function(a) {
    k <- if (a == 0.05) 1 else if (a == 0.01) 0.6 else 1.2
    X[, .(block, cell, role, null_cell, n, subset = "all", B, test, alpha = a,
          rejection = ifelse(role == "null", a * sz, pmin(1, pw * k + 0.01)),
          size_adj_power = ifelse(role == "null", NA_real_, pmin(1, pw * k)),
          null_size = ifelse(role == "null", NA_real_, a * 0.9), declined = dec)]
  }))
  if (!is.null(gate)) for (i in seq_len(nrow(gate))) {
    j <- which(R$block == gate$block[i] & R$cell == gate$cell[i] & R$test == gate$test[i] & abs(R$alpha - gate$alpha[i]) < 1e-12)
    if (length(j) != 1) stop("gate override not found: ", gate$cell[i], " ", gate$test[i])
    set(R, j, "rejection", if ("size" %in% names(gate)) gate$size[i] else st_size_for_z_null(gate$z[i], gate$alpha[i], R$B[j]))
  }
  R[, mcse := sqrt(rejection * (1 - rejection) / B)]
  R[, `:=`(rejection_given_p = pmin(1, rejection / (1 - declined)), status = "")]
  FL <- unique(X[, .(block, cell, role, null_cell, n, B, flag_k)])
  F2 <- rbindlist(lapply(c("G10", "Grule"), function(g)
    FL[, .(block, cell, role, null_cell, n, subset = "all", B, test = paste0("flag.evlt.", g), alpha = NA_real_,
           rejection = c(0, 0.25, 1)[flag_k], mcse = NA_real_, size_adj_power = NA_real_, null_size = NA_real_, declined = 0,
           rejection_given_p = NA_real_, status = "flag rate: share of replicates with the flag = 1")]))
  split_cells <- FL[flag_k == 2L & block %in% c("1b", "4"), paste(block, cell)]
  SP <- rbindlist(lapply(c("G10", "Grule"), function(g) {
    z <- R[paste(block, cell) %in% split_cells]
    rbind(copy(z)[, `:=`(subset = paste0("events<", g), B = round(B / 4), declined = pmin(1, declined + 0.3))],
          copy(z)[, `:=`(subset = paste0("events>=", g), B = B - round(B / 4))])
  }))
  S <- rbind(R, F2, SP, use.names = TRUE)
  setcolorder(S, c("block", "cell", "role", "null_cell", "n", "subset", "B", "test", "alpha", "rejection", "mcse", "size_adj_power",
                   "null_size", "declined", "rejection_given_p", "status"))
  for (b in unique(S$block)) { dir.create(file.path(dir, b), showWarnings = FALSE); fwrite(S[block == b], file.path(dir, b, "_summary.csv")) }
  invisible(S)
}
bases <- c("poly3", "sym", "stk")
nb <- function(block, cell) CT$B[CT$block == block & CT$cell == cell]
decision_is <- function(res, want, pattern) {
  D <- res$decision
  all(D$choice[match(bases, D$basis)] == want) && all(grepl(pattern, D$reason))
}

## ---- 1 and (i): memberships; score +0.02 everywhere ------------------------------------------------------------------------------
cat("1 / (i) families and score +0.02 in every family\n")
D1 <- file.path(OUTDIR, "synth_i")
guarded("(i)", {
  st_write_root(D1, function(f) rep(0.02, length(f)))
  r1 <- run_quiet(c("--root", D1, "--out", file.path(D1, "out")), log = file.path(D1, "run.out"))
  want <- C0[!is.na(fam), .(block, cell, fam)]
  got <- r1$M[, .(block, cell, fam = as.integer(family))]
  check("1: family memberships equal the independent assignment", nrow(fsetdiff(want, got)) == 0 && nrow(fsetdiff(got, want)) == 0,
        paste(sprintf("family %d: %d", 1:6, tabulate(got$fam, 6)), collapse = ", "))
  check("1: family sizes 32, 21, 10, 80, 13, 2 as E12.4 lists them", identical(tabulate(got$fam, 6), c(32L, 21L, 10L, 80L, 13L, 2L)))
  check("1: every named pattern matches a cell with a summary", all(r1$patterns$with_summary > 0),
        sprintf("%d patterns", nrow(r1$patterns)))
  nnull <- nrow(C0[role == "null"])
  check("(i) the gate holds every null cell of blocks 1a-7 for each variant", all(r1$gate[, .N, by = test]$N == nnull) && nrow(r1$gate[, .N, by = test]) == 6,
        sprintf("%d null cells x 6 variants", nnull))
  check("(i) both forms eligible for every basis", all(r1$decision$eligible_unit & r1$decision$eligible_score))
  check("(i) difference 0.02 for every basis", near(r1$decision$difference, rep(0.02, 3), 1e-12), paste(sprintf("%.6f", r1$decision$difference), collapse = " "))
  check("(i) choice score for poly3, sym, stk", decision_is(r1, "score", "A rule 1"), paste(r1$decision$reason[1]))
  check("(i) forms used for H1-H5 = score", identical(unname(r1$forms_used[c("poly3", "sym")]), c("sc", "sc")))
  fs <- c("rule_A_gate.csv", "rule_A_membership.csv", "rule_A_families.csv", "rule_A_decision.txt", "hypotheses.csv", "h4_census.csv", "declined.csv", "mcnemar.csv")
  check("(i) the eight outputs are written", all(file.exists(file.path(D1, "out", fs))))
  dl <- readLines(file.path(D1, "out", "rule_A_decision.txt"))
  check("(i) rule_A_decision.txt ends with the three choices", all(c("choice poly3 score", "choice sym score", "choice stk score") %in% dl))
  check("H4 and H5 cell sets: 29 census cells (21 of families 2-4 at n = 1000 + 8), 4 base and 4 skew",
        nrow(r1$H4c) == 29 && sum(grepl("^census_", r1$H4c$cell)) == 8 && nrow(r1$H5c[set == "base"]) == 4 && nrow(r1$H5c[set == "skew"]) == 4 &&
          setequal(r1$H5c[set == "base"]$n, r1$H5c[set == "skew"]$n))
  dsg <- table(an_h1_design(r1$M[family == 1]$block, r1$M[family == 1]$cell))
  check("H1 designs: base 20 (block 3 counted as base), auc 6, e12 6", identical(as.integer(dsg[c("base", "auc", "e12")]), c(20L, 6L, 6L)))
  ## H6
  DC <- r1$declined
  k1 <- unique(BASE[flag_k == 3L & block == "2", cell])[1]           # flag rate 1: every sample has events < G
  a <- DC[block == "2" & cell == k1 & test == "EDGE.poly3.u.Grule"]
  check("H6: rate-1 cell, events<Grule is the 'all' row and events>=Grule holds no sample",
        near(a[subset == "events<Grule"]$declined, a[subset == "all"]$declined) && a[subset == "events>=Grule"]$B == 0 &&
          grepl("all samples", a[subset == "events<Grule"]$source), k1)
  k2 <- unique(BASE[flag_k == 2L & block == "4", cell])[1]
  a2 <- DC[block == "4" & cell == k2 & test == "Stk.joint"]
  b2 <- nb("4", k2)
  check("H6: split cell reads the summary subsets (B/4 with events < G, declined + 0.3)",
        all(a2[subset != "all"]$source == "summary") && a2[subset == "events<G10"]$B == round(b2 / 4) &&
          near(a2[subset == "events<G10"]$declined, a2[subset == "all"]$declined + 0.3, 1e-12), k2)
  check("H6: every test of every synthetic cell has an 'all' row", nrow(DC[subset == "all"]) == nrow(BASE))
  check("H6: arms of tests and subsets (HL_w is G10; Stk.joint has no arm)",
        all(DC[test == "HL_w"]$test_arm == "G10") && all(DC[test == "Stk.joint"]$subset_matches_arm) &&
          all(!DC[test == "EDGE.poly3.u.G10" & subset_G == "Grule"]$subset_matches_arm))
  ## --arm 10 and --forms
  r1g <- run_quiet(c("--root", D1, "--out", file.path(D1, "out"), "--arm", "10"), log = file.path(D1, "run_G10.out"))
  fs10 <- sub("(\\.csv|\\.txt)$", "_G10\\1", fs)
  check("--arm 10 writes the eight outputs with suffix _G10", all(file.exists(file.path(D1, "out", fs10))))
  check("--arm 10 gate and families use the G10 variants", all(grepl("\\.G10$", r1g$gate$test)) && all(grepl("\\.G10$", r1g$families$test)))
  check("--arm 10: hypotheses H1-H3, H5 at G10, H4 at G10 with the rule-G version", all(r1g$hypotheses[hypothesis %in% c("H1", "H2", "H3", "H5")]$arm == "G10") &&
          setequal(r1g$hypotheses[grepl("^H4", hypothesis)]$arm, c("G10", "Grule")))
  ru <- run_quiet(c("--root", D1, "--out", file.path(D1, "out_unit"), "--forms", "unit"))
  check("--forms unit fixes the forms of H1-H5", identical(unname(ru$forms_used), rep("u", 3)) && all(ru$hypotheses[forms_used == TRUE]$form == "unit"))
})

## ---- (ii) score +0.02 on average, -0.11 in family 6 --------------------------------------------------------------------------------
cat("\n(ii) score +0.02 on average but -0.11 in family 6\n")
guarded("(ii)", {
  D2 <- file.path(OUTDIR, "synth_ii")
  st_write_root(D2, function(f) ifelse(f == 6L, -0.11, 0.046))            # (5 x 0.046 - 0.11) / 6 = 0.02
  r2 <- run_quiet(c("--root", D2, "--out", file.path(D2, "out")), log = file.path(D2, "run.out"))
  check("(ii) difference 0.02, worst score loss 0.11 in family 6", near(r2$decision$difference, rep(0.02, 3), 1e-12) &&
          near(r2$decision$worst_loss_score, rep(0.11, 3), 1e-12) && all(r2$decision$worst_loss_score_family == 6),
        paste(sprintf("%.6f / %.6f", r2$decision$difference, r2$decision$worst_loss_score), collapse = "; "))
  check("(ii) choice unit (score loses more than 0.10 in one family)", decision_is(r2, "unit", "loses 0.1100 .*family 6"), r2$decision$reason[1])
})

## ---- (iii) difference 0.010 ---------------------------------------------------------------------------------------------------
cat("\n(iii) difference 0.010\n")
guarded("(iii)", {
  D3 <- file.path(OUTDIR, "synth_iii")
  st_write_root(D3, function(f) rep(0.010, length(f)))
  r3 <- run_quiet(c("--root", D3, "--out", file.path(D3, "out")), log = file.path(D3, "run.out"))
  check("(iii) difference 0.010", near(r3$decision$difference, rep(0.010, 3), 1e-12))
  check("(iii) choice unit (below the tie margin)", decision_is(r3, "unit", "below the tie margin"), r3$decision$reason[1])
})

## ---- (iv) unit fails the gate once (z = 3.2) ----------------------------------------------------------------------------------------
cat("\n(iv) unit fails the gate in one null cell, score eligible\n")
guarded("(iv)", {
  D4 <- file.path(OUTDIR, "synth_iv")
  gate <- rbind(data.table(block = "1a", cell = "null_link_n1000", test = paste0("EDGE.", bases, ".u.Grule"), alpha = 0.05, z = 3.2),
                data.table(block = "1b", cell = "null_base_n610", test = paste0("EDGE.", bases, ".u.Grule"), alpha = 0.05, z = 2.95))
  st_write_root(D4, function(f) rep(-0.03, length(f)), gate)              # score power lower: the gate alone decides
  r4 <- run_quiet(c("--root", D4, "--out", file.path(D4, "out")), log = file.path(D4, "run.out"))
  G <- r4$gate
  fz <- G[form == "u" & pass == FALSE]
  check("(iv) unit fails in exactly one null cell per basis, z = 3.2 at 0.05", nrow(fz) == 3 && all(fz$cell == "null_link_n1000") &&
          near(fz$z_05, rep(3.2, 3), 1e-8) && all(!fz$pass_05) && all(fz$pass_01),
        paste(sprintf("%s z %.6f", fz$test, fz$z_05), collapse = "; "))
  pz <- G[form == "u" & cell == "null_base_n610"]
  check("(iv) z = 2.95 passes", nrow(pz) == 3 && all(pz$pass) && near(pz$z_05, rep(2.95, 3), 1e-8))
  check("(iv) score eligible, unit not", all(!r4$decision$eligible_unit & r4$decision$eligible_score))
  check("(iv) choice score although its macro-average is 0.03 lower", decision_is(r4, "score", "only the score form") &&
          near(r4$decision$difference, rep(-0.03, 3), 1e-12), r4$decision$reason[1])
  dl <- readLines(file.path(D4, "out", "rule_A_decision.txt"))
  check("(iv) the failing cell is listed with its z in rule_A_decision.txt", any(grepl("unit fails: .*1a/null_link_n1000 \\(z05 3.20", dl)))
  r4g <- run_quiet(c("--root", D4, "--out", file.path(D4, "out"), "--arm", "10"))
  check("(iv) --arm 10: at G10 both eligible and unit chosen, but H1-H5 keep the rule-G choice (score)",
        all(r4g$decision$choice == "unit") && identical(unname(r4g$forms_used[c("poly3", "sym")]), c("sc", "sc")))
})

## ---- (v) both fail ----------------------------------------------------------------------------------------------------------------
cat("\n(v) both forms fail the gate\n")
guarded("(v)", {
  D5 <- file.path(OUTDIR, "synth_v")
  gate <- rbind(data.table(block = "1a", cell = "null_link_n1000", test = paste0("EDGE.", bases, ".u.Grule"), alpha = 0.05, z = 3.2),
                data.table(block = "1b", cell = "null_base_n610", test = paste0("EDGE.", bases, ".sc.Grule"), alpha = 0.01, z = 3.5))
  st_write_root(D5, function(f) rep(0.03, length(f)), gate)
  r5 <- run_quiet(c("--root", D5, "--out", file.path(D5, "out")), log = file.path(D5, "run.out"))
  fs <- r5$gate[form == "sc" & pass == FALSE]
  check("(v) score fails at 0.01 (z = 3.5) in one cell", nrow(fs) == 3 && all(!fs$pass_01) && all(fs$pass_05) && near(fs$z_01, rep(3.5, 3), 1e-8))
  check("(v) neither eligible -> unit, although score is 0.03 higher", decision_is(r5, "unit", "neither form passes"), r5$decision$reason[1])
})

## ---- (vi) H1, H2, H3, H5 at their thresholds ----------------------------------------------------------------------------------------
cat("\n(vi) hypotheses at their thresholds and just beyond\n")
guarded("(vi)", {
  fu <- c(poly3 = "u", sym = "u", stk = "u")
  setp <- function(V, cells, test, value) {
    q <- data.table(block = cells$block, cell = cells$cell, test = test, alpha = 0.05, val = value)
    V[q, on = c("block", "cell", "test", "alpha"), size_adj_power := i.val]
    invisible(V)
  }
  getp <- function(V, cells, test) an_get(V, cells$block, cells$cell, rep(test, nrow(cells)), "size_adj_power")
  M <- r1$M; f1 <- M[family == 1]; f2 <- M[family == 2]; hb <- r1$H5c[set == "base"]; hs <- r1$H5c[set == "skew"]
  verdict <- function(V, hyp) { H <- an_hypotheses(V, M, r1$H4c, r1$H5c, "Grule", fu)$hyp; H[hypothesis == hyp & form == "unit"] }
  p3 <- getp(r1$V, f1, "EDGE.poly3.u.Grule"); ps <- getp(r1$V, f1, "EDGE.sym.u.Grule")
  for (cs in list(list(g = 0.10, want = "holds"), list(g = 0.0999, want = "against"))) {
    V <- copy(r1$V); setp(V, f1, "GiViTI", p3 - cs$g); h <- verdict(V, "H1")
    check(sprintf("(vi) H1 mean gain %.4f -> %s", cs$g, cs$want), h$verdict == cs$want && near(h$statistic, cs$g, 1e-12) && h$designs_below_margin == "",
          sprintf("statistic %.15f", h$statistic))
  }
  dsg <- an_h1_design(f1$block, f1$cell)
  for (cs in list(list(a = -0.014, want = "holds", below = ""), list(a = -0.0141, want = "against", below = "auc"))) {
    V <- copy(r1$V); setp(V, f1, "GiViTI", p3 - ifelse(dsg == "auc", cs$a, 0.13)); h <- verdict(V, "H1")
    check(sprintf("(vi) H1 auc design mean %.4f (overall mean %.4f) -> %s", cs$a, h$statistic, cs$want),
          h$verdict == cs$want && h$designs_below_margin == cs$below && h$statistic >= 0.10 && near(h$mean_auc, cs$a, 1e-12))
  }
  for (cs in list(list(g = -0.014, want = "holds"), list(g = -0.0141, want = "against"))) {
    V <- copy(r1$V); setp(V, f1, "Stk.joint", ps - cs$g); h <- verdict(V, "H2")
    check(sprintf("(vi) H2 mean %.4f -> %s", cs$g, cs$want), h$verdict == cs$want && near(h$statistic, cs$g, 1e-12))
  }
  p32 <- getp(r1$V, f2, "EDGE.poly3.u.Grule")
  for (cs in list(list(g = -0.05, want = "holds"), list(g = -0.0501, want = "against"))) {
    V <- copy(r1$V); setp(V, f2, "GiViTI", p32 - cs$g); h <- verdict(V, "H3")
    check(sprintf("(vi) H3 mean %.4f -> %s", cs$g, cs$want), h$verdict == cs$want && near(h$statistic, cs$g, 1e-12))
  }
  psb <- getp(r1$V, hb, "EDGE.sym.u.Grule"); pss <- getp(r1$V, hs, "EDGE.sym.u.Grule")
  for (cs in list(list(g = 0.02, want = "holds"), list(g = 0.0201, want = "against"))) {
    V <- copy(r1$V); setp(V, hb, "Stk.sym1", psb + cs$g); setp(V, hs, "Stk.sym1", pss + 0.05)
    h <- verdict(V, "H5"); k <- verdict(V, "H5.skew")
    check(sprintf("(vi) H5 base cost %.4f -> %s; skew cost 0.05 reported as larger", cs$g, cs$want),
          h$verdict == cs$want && near(h$statistic, cs$g, 1e-12) && k$verdict == "reported" && near(k$statistic, 0.05, 1e-12) && grepl("larger than", k$note))
  }
})

## ---- (vii) H4 census against make_headline_recount.R ----------------------------------------------------------------------------------
cat("\n(vii) H4 census against make_headline_recount.R\n")
guarded("(vii)", {
  src <- readLines(file.path(SIMDIR, "make_headline_recount.R"), warn = FALSE)
  dest <- file.path(OUTDIR, "headline_recount_rerun.csv")
  iw <- grep('write.csv(out, "headline_recount.csv"', src, fixed = TRUE)
  src[iw] <- sub('"headline_recount.csv"', sprintf('"%s"', dest), src[iw], fixed = TRUE)
  if (length(iw) != 1 || any(grepl('"headline_recount.csv"', src, fixed = TRUE))) stop("could not redirect the script's output file")
  env <- new.env(); owd <- getwd()
  txt <- capture.output(eval(parse(text = src), envir = env)); setwd(owd)
  cat("    make_headline_recount.R (output redirected): ", txt[1], "\n", sep = "")
  pb <- env$pb; pb1 <- pb[!duplicated(pb[, c("family", "param", "test")]), ]      # the script's g() takes the first matching row
  map <- c("DEF.poly3" = "EDGE.poly3.u.G10", EF = "HLF.G10", HL = "HL.G10", "HL-equalwidth" = "HL_w", "Pigeon-Heyse" = "PH", Tsiatis = "Tsiatis", Xie = "Xie")
  cells <- unique(pb[pb$family != "null", c("family", "param")])
  rows <- list()
  for (i in seq_len(nrow(cells))) for (old in names(map)) {
    v <- pb1$power_size_adj[pb1$family == cells$family[i] & pb1$param == cells$param[i] & pb1$test == old]
    rows[[length(rows) + 1]] <- data.table(block = "3", cell = sprintf("h4c%02d", i), test = map[[old]], alpha = 0.05, size_adj_power = if (length(v)) v[1] else NA_real_)
  }
  ex <- env$ex; exmap <- c(edge_poly3 = "EDGE.poly3.u.G10", EF = "HLF.G10", HL = "HL.G10", HLw = "HL_w", PH = "PH", Tsi = "Tsiatis", Xie = "Xie")
  for (i in seq_len(nrow(ex))) for (old in names(exmap))
    rows[[length(rows) + 1]] <- data.table(block = "3", cell = sprintf("h4c%02d", nrow(cells) + i), test = exmap[[old]], alpha = 0.05, size_adj_power = as.numeric(ex[i, old]))
  Vf <- rbindlist(rows)
  hc <- data.table(block = "3", cell = sprintf("h4c%02d", seq_len(nrow(cells) + nrow(ex))), n = 1000)
  K <- an_h4(Vf, hc, "u", "G10")
  O <- env$out
  same <- function(col) identical(as.logical(K[[col]]), as.logical(O[[if (col == "in_headline") "in_headline_22" else col]]))   # the script's name
  check("(vii) per-scenario detectable, saturated, leads/ties partition, leads/ties all six, in headline = the script",
        nrow(K) == nrow(O) && all(vapply(c("detectable", "saturated", "beats_partition", "beats_all6", "in_headline"), same, logical(1))),
        sprintf("%d scenarios (%d July cells + %d tab_extra)", nrow(K), nrow(cells), nrow(ex)))
  nm <- c(EF = "HLF.G10", HL = "HL.G10", "HL-equalwidth" = "HL_w", "Pigeon-Heyse" = "PH", HLw = "HL_w", PH = "PH", Tsiatis = "Tsiatis", Tsi = "Tsiatis", Xie = "Xie")
  check("(vii) best partition and best-of-six tests and values = the script",
        identical(unname(nm[O$best_partition_name]), K$best_partition_name) && identical(unname(nm[O$best_all6_name]), K$best_all6_name) &&
          all(abs(K$best_partition - O$best_partition) <= 5e-4) && all(abs(K$best_all6 - O$best_all6) <= 5e-4))
  inh <- env$inh
  cnt <- c(headline = sum(K$in_headline), led = sum(K$in_headline & K$beats_partition), led6 = sum(K$in_headline & K$beats_all6))
  check("(vii) census counts = the script (cells in headline, led or tied vs partition, vs all six)",
        identical(unname(cnt), c(nrow(inh), sum(inh$beats_partition), sum(inh$beats_all6))),
        sprintf("analysis %d / %d / %d; script %d / %d / %d", cnt[1], cnt[2], cnt[3], nrow(inh), sum(inh$beats_partition), sum(inh$beats_all6)))
  check("(vii) the July census holds H4's bar (at least 19)", cnt[["led"]] >= H4_HOLDS, sprintf("%d of %d", cnt[["led"]], cnt[["headline"]]))
})

## ---- (viii) McNemar and Holm by hand ----------------------------------------------------------------------------------------------------
cat("\n(viii) exact McNemar and Holm against a hand computation\n")
holm_hand <- function(p) { m <- length(p); o <- order(p); adj <- pmin(1, cummax((m - seq_len(m) + 1) * p[o])); out <- numeric(m); out[o] <- adj; out }
mcn_hand <- function(b10, b01) { m <- b10 + b01; if (m == 0) 1 else min(1, 2 * sum(dbinom(0:min(b10, b01), m, 0.5))) }
guarded("(viii) formula", {
  check("(viii) exact McNemar p: (0,5) 0.0625; (1,9) 22/1024; (3,3) 1; (0,0) 1; = binom.test on (2,12)",
        near(an_mcnemar_p(0, 5), 0.0625, 1e-15) && near(an_mcnemar_p(1, 9), 22 / 1024, 1e-15) && an_mcnemar_p(3, 3) == 1 && an_mcnemar_p(0, 0) == 1 &&
          near(an_mcnemar_p(2, 12), binom.test(2, 14, 0.5)$p.value, 1e-12))
})
guarded("(viii) dry run", {
  DR <- edge_battery("dryrun")
  compared <- c()
  for (rdg in AN_MCSE) {
    od <- file.path(OUTDIR, "dryrun_out", rdg)
    rd <- run_quiet(c("--root", "battery/dryrun", "--out", od, "--mcse", rdg), log = file.path(od, "run.out"))
    MC <- rd$mcnemar
    todo <- MC[!is.na(B)]
    bad <- 0L; nused <- 0L
    lim <- function(s) 0.05 + 3 * (if (rdg == "null") sqrt(0.05 * 0.95 / s$B) else s$mcse)
    for (i in seq_len(nrow(todo))) {
      r <- todo[i]
      P <- fread(file.path(DR, r$block, paste0(r$cell, "_pvalues.csv.gz")))
      ra <- !is.na(P[[r$test_a]]) & P[[r$test_a]] <= 0.05; rb <- !is.na(P[[r$test_b]]) & P[[r$test_b]] <= 0.05
      tab <- table(factor(ra, c(FALSE, TRUE)), factor(rb, c(FALSE, TRUE)))
      NS <- fread(file.path(DR, r$null_block, "_summary.csv"))
      sa <- NS[cell == r$null_cell & test == r$test_a & subset == "all" & abs(alpha - 0.05) < 1e-9]
      sb <- NS[cell == r$null_cell & test == r$test_b & subset == "all" & abs(alpha - 0.05) < 1e-9]
      hold <- nrow(sa) == 1 && nrow(sb) == 1 && sa$rejection <= lim(sa) + 1e-9 && sb$rejection <= lim(sb) + 1e-9
      infam <- !startsWith(r$hypothesis, "H4") || isTRUE(r$in_headline)      # E13.4: H4 pairs only in detectable, non-saturated cells
      use <- hold && infam
      ok <- r$both == tab["TRUE", "TRUE"] && r$a_only == tab["TRUE", "FALSE"] && r$b_only == tab["FALSE", "TRUE"] && r$neither == tab["FALSE", "FALSE"] &&
        r$used == use && r$in_family == infam && r$mcse_reading == rdg && (!use || near(r$p_mcnemar, mcn_hand(tab["TRUE", "FALSE"], tab["FALSE", "TRUE"]), 1e-12))
      nused <- nused + use
      if (!ok) { bad <- bad + 1L; cat("    mismatch:", rdg, r$hypothesis, r$form, r$cell, r$test_a, r$test_b, "\n") }
    }
    compared[rdg] <- nused
    check(sprintf("(viii) dry run, %s MCSE: counts, size holding and McNemar p = hand computation", rdg), bad == 0L && nrow(todo) > 0,
          sprintf("%d pairs with a per-replicate file, %d compared, %d mismatches", nrow(todo), nused, bad))
    U <- MC[used == TRUE]
    hh <- U[, .(cell, test_b, p_holm, hand = holm_hand(p_mcnemar)), by = .(hypothesis, arm, form)]
    check(sprintf("(viii) dry run, %s MCSE: Holm = hand computation within each hypothesis", rdg), nrow(U) > 0 && near(hh$p_holm, hh$hand, 1e-12),
          sprintf("%d families", nrow(unique(U[, .(hypothesis, arm, form)]))))
    check(sprintf("(viii) dry run, %s MCSE: no comparison without a matched null", rdg), nrow(MC[used == TRUE & !nzchar(null_cell)]) == 0 &&
            all(MC[cell == "cauchit_auc_n460"]$used == FALSE))
  }
  check("(viii) dry run: the null reading compares no more pairs than the realised one", compared[["null"]] <= compared[["realised"]],
        sprintf("compared: null %d, realised %d", compared[["null"]], compared[["realised"]]))
})
guarded("(viii) synthetic", {
  f2 <- C0[fam == 2L, .(block, cell)]
  exp <- list()
  for (k in seq_len(nrow(f2))) {
    Bk <- 200L; ea <- rep(0.5, Bk); gv <- rep(0.5, Bk)
    ea[1:20] <- 0.01; gv[1:20] <- 0.01                                   # both reject
    ia <- 20 + seq_len(k); ea[ia] <- 0.01; ea[ia[1]] <- 0.05             # EDGE alone; p = 0.05 exactly rejects (p <= alpha)
    gv[ia[length(ia)]] <- NA                                             # a declined GiViTI counts as no rejection
    if (k %% 2 == 1) { ea[21 + k] <- 0.06; gv[21 + k] <- 0.03 }          # GiViTI alone
    P <- data.table(rep = seq_len(Bk), EDGE.poly3.u.Grule = ea, EDGE.poly3.sc.Grule = ea, GiViTI = gv)
    fwrite(P, file.path(D1, f2$block[k], paste0(f2$cell[k], "_pvalues.csv.gz")), compress = "gzip")
    exp[[k]] <- data.table(block = f2$block[k], cell = f2$cell[k], both = 20L, a_only = k, b_only = k %% 2L, p = mcn_hand(k, k %% 2L))
  }
  E <- rbindlist(exp); E[, holm := holm_hand(p)]
  rs <- run_quiet(c("--root", D1, "--out", file.path(D1, "out_mcnemar")))
  X <- merge(rs$mcnemar[hypothesis == "H3" & form == "unit"], E, by = c("block", "cell"))
  check("(viii) synthetic 21 cells of H3: counts, exact p and Holm = hand", nrow(X) == 21 && all(X$used) && all(X$both.x == X$both.y) &&
          all(X$a_only.x == X$a_only.y) && all(X$b_only.x == X$b_only.y) && near(X$p_mcnemar, X$p, 1e-12) && near(X$p_holm, X$holm, 1e-12),
        sprintf("Holm-significant %d of 21 (hand %d)", sum(X$holm_reject), sum(E$holm <= 0.05)))
  unlink(file.path(D1, f2$block, paste0(f2$cell, "_pvalues.csv.gz")))
})

## ---- (ix) refusal ---------------------------------------------------------------------------------------------------------------
cat("\n(ix) refusal before \"launch finished\"\n")
guarded("(ix)", {
  lg <- file.path(OUTDIR, "synth_launch.log")
  st <- function(lines) { if (is.null(lines)) unlink(lg) else writeLines(lines, lg); an_launch_finished(lg) }
  s1 <- "2026-09-14 08:28:48  launch from step 1 of 19"; s7 <- "2026-09-14 12:00:00  launch from step 7 of 19"
  fin <- "2026-09-14 14:00:00  launch finished"; mid <- "2026-09-14 10:12:30  step 7: run block 3"
  check("(ix) launch.log reading: none / running / finished / restarted after finishing / finished again",
        !st(NULL) && !st(c(s1, mid)) && st(c(s1, mid, fin)) && !st(c(s1, fin, s7)) && st(c(s1, fin, s7, mid, fin)))
  DRr <- file.path(OUTDIR, "synth_realroot")
  st_write_root(DRr, function(f) rep(0.02, length(f)))
  out <- file.path(DRr, "out")
  writeLines(c(s1, mid), file.path(DRr, "launch.log"))
  e1 <- tryCatch({ run_quiet(c("--root", DRr, "--out", out), bat = DRr); "ran" }, error = function(e) conditionMessage(e))
  check("(ix) a root that is the battery refuses while the launch is not finished", grepl("refusing the real root", e1) && !file.exists(file.path(out, "hypotheses.csv")), e1)
  e2 <- tryCatch({ run_quiet(c("--root", DRr, "--out", out, "--force"), bat = DRr); "ran" }, error = function(e) conditionMessage(e))
  check("(ix) --force runs it", e2 == "ran" && file.exists(file.path(out, "hypotheses.csv")))
  unlink(out, recursive = TRUE)
  writeLines(c(s1, mid, fin), file.path(DRr, "launch.log"))
  e3 <- tryCatch({ run_quiet(c("--root", DRr, "--out", out), bat = DRr); "ran" }, error = function(e) conditionMessage(e))
  check("(ix) after \"launch finished\" it runs", e3 == "ran" && file.exists(file.path(out, "hypotheses.csv")), e3)
  S4 <- fread(file.path(DRr, "4", "_summary.csv"), colClasses = list(character = c("block", "cell")))
  fwrite(S4[!grepl("^crossover_", cell)], file.path(DRr, "4", "_summary.csv"))
  e4 <- tryCatch({ run_quiet(c("--root", DRr, "--out", out), bat = DRr); "ran" }, error = function(e) conditionMessage(e))
  check("(ix) finished real root with no crossover cell in its summaries refuses (named pattern with no cell)",
        grepl("refusing the real root", e4) && grepl("crossover_n1000", e4, fixed = TRUE), gsub("\n", " | ", e4))
  e5 <- tryCatch({ run_quiet(c("--root", DRr, "--out", out, "--force"), bat = DRr); "ran" }, error = function(e) conditionMessage(e))
  check("(ix) --force runs past the missing pattern", e5 == "ran")
  check("(ix) real-root test: battery (and battery/) is the real root, battery/dryrun is not; default root is battery",
        an_is_real("battery") && an_is_real(edge_battery()) && !an_is_real("battery/dryrun") && an_opts(character(0))$root == "battery" &&
          an_opts(character(0))$out == file.path("battery", "analysis"))
  check("(ix) bad options stop", inherits(tryCatch(an_opts(c("--arm", "5")), error = function(e) e), "error") &&
          inherits(tryCatch(an_opts(c("--forms", "both")), error = function(e) e), "error"))
  ## the real battery: only the guard, and only while its own reading of launch.log says the launch is not finished
  real_lg <- edge_battery("launch.log")
  L <- if (file.exists(real_lg)) readLines(real_lg, warn = FALSE) else character(0)
  starts <- grep("launch from step", L, fixed = TRUE); fins <- grep("launch finished", L, fixed = TRUE)
  finished <- length(fins) > 0 && (length(starts) == 0 || max(fins) > max(starts))
  if (!finished) {
    e6 <- tryCatch({ an_guard(edge_battery(), force = FALSE); "no refusal" }, error = function(e) conditionMessage(e))
    check("(ix) the real root refuses now (launch.log has no \"launch finished\" after its last launch)", grepl("refusing the real root", e6), e6)
  } else {
    check("(ix) the real root guard was not called: launch.log records \"launch finished\" (the results stay blind)", TRUE)
  }
})

## ---- (x) the MCSE reading of the size checks ------------------------------------------------------------------------------------
cat("\n(x) the MCSE of the size checks: null sqrt(a(1-a)/B) by default, realised sqrt(r(1-r)/B) with --mcse realised\n")
guarded("(x)", {
  DX <- file.path(OUTDIR, "synth_x")
  n500 <- CT[block == "4" & role == "null" & B == 500, cell][1]
  gate <- data.table(block = "4", cell = n500, test = paste0("EDGE.", bases, ".sc.Grule"), alpha = 0.01, size = 0.030)
  st_write_root(DX, function(f) rep(0.02, length(f)), gate)
  rn <- run_quiet(c("--root", DX, "--out", file.path(DX, "out")), log = file.path(DX, "run.out"))
  rr <- run_quiet(c("--root", DX, "--out", file.path(DX, "out_realised"), "--mcse", "realised"), log = file.path(DX, "run_realised.out"))
  g <- rn$gate[cell == n500 & form == "sc"]
  zn <- 0.02 / sqrt(0.01 * 0.99 / 500); zr <- 0.02 / sqrt(0.03 * 0.97 / 500)
  check("(x) size 0.030 at 0.01 with B = 500: z 4.49 against the null MCSE (fails), 2.62 against the realised MCSE (passes)",
        nrow(g) == 3 && near(g$size_01, rep(0.03, 3), 1e-12) && near(g$z_null_01, rep(zn, 3), 1e-9) && near(g$z_realised_01, rep(zr, 3), 1e-9) &&
          all(!g$pass_null_01) && all(g$pass_realised_01) && all(!g$pass) && near(g$z_01, rep(zn, 3), 1e-9) && all(g$mcse_reading == "null"),
        sprintf("%s: z %.4f / %.4f", n500, zn, zr))
  check("(x) largest passing size at 0.01 with B = 500: 0.0233 (null MCSE), 0.0345 (realised MCSE)",
        abs(st_size_for_z_null(3, 0.01, 500) - 0.0233) < 5e-5 && abs(st_size_for_z(3, 0.01, 500) - 0.0345) < 5e-5)
  check("(x) default (null MCSE): the score form fails the gate and unit is chosen", decision_is(rn, "unit", "only the unit form") &&
          all(!rn$decision$eligible_score) && all(rn$decision$eligible_unit), rn$decision$reason[1])
  check("(x) --mcse realised: both forms eligible and score chosen (+0.02)", decision_is(rr, "score", "A rule 1") &&
          all(rr$gate[cell == n500 & form == "sc"]$pass) && all(rr$gate$mcse_reading == "realised"), rr$decision$reason[1])
  dl <- readLines(file.path(DX, "out", "rule_A_decision.txt")); dr <- readLines(file.path(DX, "out_realised", "rule_A_decision_mcse_realised.txt"))
  check("(x) rule_A_decision.txt names the reading", any(grepl("^MCSE of the size gate: null", dl)) && any(grepl("^MCSE of the size gate: realised", dr)))
  check("(x) the sensitivity reading writes its own file names, so an E13.1 file is never overwritten",
        file.exists(file.path(DX, "out_realised", "hypotheses_mcse_realised.csv")) &&
          !file.exists(file.path(DX, "out_realised", "hypotheses.csv")) &&
          file.exists(file.path(DX, "out", "hypotheses.csv")))
  ## E12.6 uses the same reading: GiViTI in the matched null of a family 2 cell, size with realised z = 2.9
  pr <- rn$pairs[hypothesis == "H3" & form == "unit"][1]
  nc <- CT[block == pr$block & cell == pr$cell]
  Bn <- CT[block == nc$null_block & cell == nc$null_cell]$B
  sz <- st_size_for_z(2.9, 0.05, Bn)
  Vm <- copy(rn$V)
  Vm[block == nc$null_block & cell == nc$null_cell & test == "GiViTI" & alpha == 0.05, `:=`(rejection = sz, mcse = sqrt(sz * (1 - sz) / B))]
  mn <- an_mcnemar(DX, Vm, rn$C, pr, "null"); mr <- an_mcnemar(DX, Vm, rn$C, pr, "realised")
  check("(x) E12.6 follows the reading: a GiViTI null size with realised z 2.9 fails with the null MCSE and holds with the realised one",
        nrow(mn) == 1 && isFALSE(mn$holds_b) && isTRUE(mr$holds_b) && isTRUE(mn$holds_a) && grepl("GiViTI fails size .*null MCSE", mn$note) &&
          near(mn$limit_b, 0.05 + 3 * sqrt(0.0475 / Bn), 1e-12) && near(mr$limit_b, 0.05 + 3 * sqrt(sz * (1 - sz) / Bn), 1e-12),
        sprintf("%s/%s, null %s/%s, B %d, size %.5f", pr$block, pr$cell, nc$null_block, nc$null_cell, Bn, sz))
})

## ---- (xi) cells with a summary but no size-adjusted power --------------------------------------------------------------------------
cat("\n(xi) cells with a summary but no size-adjusted power (\"matched null not run yet\")\n")
guarded("(xi)", {
  DM <- file.path(OUTDIR, "synth_missing"); out <- file.path(DM, "out")
  S0 <- st_write_root(DM, function(f) rep(0.02, length(f)))
  writeLines(c("2026-09-14 08:28:48  launch from step 1 of 19", "2026-09-15 02:00:00  launch finished"), file.path(DM, "launch.log"))
  st_case <- function(what) {          # what: block, cell, test (a regular expression) -> size_adj_power NA, as run_M_battery.R writes it
    S <- copy(S0)
    for (k in seq_len(nrow(what)))
      S[block == what$block[k] & cell == what$cell[k] & grepl(what$test[k], test) & !startsWith(test, "flag."),
        `:=`(size_adj_power = NA_real_, status = "matched null not run yet")]
    for (b in unique(S$block)) fwrite(S[block == b], file.path(DM, b, "_summary.csv"))
    unlink(out, recursive = TRUE)
  }
  try_run <- function(extra = character(0), bat = DM) {
    res <- NULL
    msg <- tryCatch({ capture.output(res <- an_main(c("--root", DM, "--out", out, extra), bat = bat)); "ran" }, error = function(e) conditionMessage(e))
    list(msg = msg, res = res)
  }
  flat <- function(m) gsub("\n", " | ", m)

  ## A: the EDGE values of one family 1 cell (the review's case 1)
  st_case(data.table(block = "3", cell = "link_probit_n1000", test = "^EDGE\\."))
  a1 <- try_run()
  check("(xi) A: a finished real root refuses when one family 1 cell has no EDGE size-adjusted power (Section A, H1, H2 named; H3 not)",
        grepl("refusing the real root", a1$msg) && grepl("Section A families: 6 cell-test pairs in 1 cells", a1$msg, fixed = TRUE) &&
          grepl("3 link_probit_n1000 EDGE.poly3.u.Grule (matched null not run yet)", a1$msg, fixed = TRUE) &&
          grepl("H1: 2 cell-test pairs", a1$msg, fixed = TRUE) && grepl("H2: 2 cell-test pairs", a1$msg, fixed = TRUE) && !grepl("H3:", a1$msg, fixed = TRUE),
        flat(a1$msg))
  a2 <- try_run("--force")
  fm <- a2$res$families
  check("(xi) A --force: family 1 (31 of 32 cells) and the macro-average have no mean, and no available mean either (E13.2); no basis is decided",
        a2$msg == "ran" && all(fm[family == "1"]$with_value == 31) && all(is.na(fm[family %in% c("1", "macro")]$mean_power)) &&
          all(is.na(fm[family %in% c("1", "macro")]$mean_power_available)) && all(a2$res$decision$choice == "no decision") &&
          all(is.na(a2$res$decision$choice_code)) && all(grepl("^no decision: .*family 1", a2$res$decision$reason)), a2$res$decision$reason[1])
  h <- a2$res$hypotheses
  check("(xi) A --force: H1 and H2 have no value, with no available mean either (E13.2), H3 keeps its value, no form is marked as used",
        all(h[hypothesis %in% c("H1", "H2")]$verdict == "no value") && all(is.na(h[hypothesis %in% c("H1", "H2")]$statistic_available)) &&
          all(grepl("1 of 32 cells have no value", h[hypothesis %in% c("H1", "H2")]$note)) && all(is.finite(h[hypothesis == "H3"]$statistic)) &&
          !any(h$forms_used) && is.na(a2$res$forms_used[["poly3"]]))
  dl <- readLines(file.path(out, "rule_A_decision.txt"))
  check("(xi) A --force: rule_A_decision.txt ends with 'choice <basis> no decision' and says that neither form is marked",
        all(c("choice poly3 no decision", "choice sym no decision", "choice stk no decision") %in% dl) && any(grepl("neither is marked as used", dl)))
  a3 <- tryCatch({ run_quiet(c("--root", DM, "--out", out)); "ran" }, error = function(e) conditionMessage(e))
  check("(xi) A on a test root: runs without --force (notes only)", a3 == "ran", a3)

  ## B: both family 6 cells (the review's case 2)
  st_case(data.table(block = "4", cell = c("crossover_n1000", "crossover_n2000"), test = "^EDGE\\."))
  b1 <- try_run(); b2 <- try_run("--force")
  check("(xi) B: both family 6 cells without EDGE power: refusal; with --force family 6 and the macro-average have no mean and no basis is decided",
        grepl("refusing the real root", b1$msg) && grepl("Section A families: 12 cell-test pairs in 2 cells", b1$msg, fixed = TRUE) &&
          b2$msg == "ran" && all(is.na(b2$res$families[family %in% c("6", "macro")]$mean_power)) && all(b2$res$families[family == "6"]$with_value == 0) &&
          all(b2$res$decision$choice == "no decision") && all(grepl("family 6", b2$res$decision$reason)), flat(b1$msg))

  ## C: one variant in one family 6 cell
  st_case(data.table(block = "4", cell = "crossover_n2000", test = "^EDGE\\.stk\\.sc\\.Grule$"))
  c2 <- try_run("--force"); D <- c2$res$decision
  check("(xi) C: EDGE.stk.sc.Grule missing in one family 6 cell: stk has no decision, poly3 and sym are decided as before (score)",
        c2$msg == "ran" && D[basis == "stk"]$choice == "no decision" && grepl("the score form has no mean in family 6", D[basis == "stk"]$reason) &&
          !grepl("unit form has no mean", D[basis == "stk"]$reason) && all(D[basis != "stk"]$choice == "score") &&
          identical(unname(c2$res$forms_used[c("poly3", "sym")]), c("sc", "sc")), D[basis == "stk"]$reason)

  ## D: GiViTI in one family 2 cell
  st_case(data.table(block = "3", cell = "loglog_n1000", test = "^GiViTI$"))
  d1 <- try_run(); d2 <- try_run("--force")
  check("(xi) D: GiViTI missing in one family 2 cell: refusal names H3 only; with --force Section A is unchanged and H3 has no value",
        grepl("H3: 1 cell-test pairs in 1 cells", d1$msg, fixed = TRUE) && !grepl("Section A families", d1$msg, fixed = TRUE) &&
          !grepl("H1:", d1$msg, fixed = TRUE) && d2$msg == "ran" && all(d2$res$decision$choice == "score") &&
          all(d2$res$hypotheses[hypothesis == "H3"]$verdict == "no value") && all(d2$res$hypotheses[hypothesis == "H1"]$verdict != "no value"), flat(d1$msg))

  ## E: a rival of the census
  st_case(data.table(block = "3", cell = "census_logx_n1000", test = "^HL_w$"))
  e1 <- try_run(); e2 <- try_run("--force")
  check("(xi) E: HL_w missing in one census cell: refusal names H4; with --force the census rule runs and the note counts the cell",
        grepl("H4: 1 cell-test pairs in 1 cells", e1$msg, fixed = TRUE) && e2$msg == "ran" &&
          all(grepl("1 of 29 cells lack a value of a compared test", e2$res$hypotheses[hypothesis %in% c("H4", "H4.ruleG")]$note)), flat(e1$msg))

  ## F: a skew cell of H5
  st_case(data.table(block = "4", cell = "cauchit_skew_n2400", test = "^Stk\\.sym1$"))
  f1 <- try_run()
  check("(xi) F: Stk.sym1 missing in one skew cell: refusal names H5.skew", grepl("H5.skew: 1 cell-test pairs", f1$msg, fixed = TRUE) &&
          !grepl("H5:", f1$msg, fixed = TRUE), flat(f1$msg))

  ## G, H: the arms a run needs
  st_case(data.table(block = "3", cell = "link_probit_n1000", test = "\\.G10$"))
  g1 <- try_run(); g2 <- try_run(c("--arm", "10"))
  check("(xi) G: G = 10 values missing in a family 1 cell: the rule-G run is not refused, the --arm 10 run is",
        g1$msg == "ran" && nrow(g1$res$missing_values) == 0 && grepl("refusing the real root", g2$msg) && grepl("EDGE.poly3.u.G10", g2$msg, fixed = TRUE),
        flat(g2$msg))
  st_case(data.table(block = "3", cell = "census_corr_n1000", test = "^EDGE\\.poly3\\.sc\\.G10$"))
  h1 <- try_run()
  check("(xi) H: H4 uses G = 10 at every arm: EDGE.poly3.sc.G10 missing in a census cell refuses the rule-G run",
        grepl("H4: 1 cell-test pairs", h1$msg, fixed = TRUE) && !grepl("Section A families", h1$msg, fixed = TRUE), flat(h1$msg))

  st_case(data.table(block = character(0), cell = character(0), test = character(0)))
  z1 <- try_run()
  check("(xi) the complete finished root runs with no refusal and no missing value", z1$msg == "ran" && nrow(z1$res$missing_values) == 0, z1$msg)
})

## ---- (xii) E13 ------------------------------------------------------------------------------------------------------------------
cat("\n(xii) E13: nominal MCSE (a), no mean over fewer cells (b), null declined (c), size failures inside a hypothesis (d), H4 family (e)\n")
st_nominal_limit <- function(a, B) a + 3 * sqrt(a * (1 - a) / B)
st_realised_limit <- function(r, a, B) a + 3 * sqrt(r * (1 - r) / B)
st_flat <- function(m) gsub("\n", " | ", m)
st_getp <- function(V, cells, tt) an_get(V, cells$block, cells$cell, rep(tt, nrow(cells)), "size_adj_power")
st_setp <- function(V, cells, tt, value) {
  V[data.table(block = cells$block, cell = cells$cell, test = tt, alpha = 0.05, val = value), on = c("block", "cell", "test", "alpha"), size_adj_power := i.val]
  invisible(V)
}
st_failnull <- function(V, nb, nc, tt) {           # the test's size at 0.05 in that null just above the nominal limit
  Bq <- CT[block == nb & cell == nc]$B; s <- st_nominal_limit(0.05, Bq) + 0.002
  V[block == nb & cell == nc & test == tt & alpha == 0.05, `:=`(rejection = s, mcse = sqrt(s * (1 - s) / B))]
  invisible(V)
}

## (a) E13.1: the size gate and the paired-test size check use sqrt(alpha (1 - alpha) / B); the realised MCSE is reported only.
## The requested case (0.023 at alpha 0.01, B = 500) lies below the nominal limit 0.023349, so it passes both readings; 0.024 is
## the size just above that limit and is the one that separates them.
guarded("(xii a)", {
  DA <- file.path(OUTDIR, "synth_e13a")
  c500 <- CT[block == "4" & role == "null" & B == 500, cell][2:3]
  gate <- rbind(data.table(block = "4", cell = c500[1], test = paste0("EDGE.", bases, ".sc.Grule"), alpha = 0.01, size = 0.024),
                data.table(block = "4", cell = c500[2], test = paste0("EDGE.", bases, ".sc.Grule"), alpha = 0.01, size = 0.023))
  st_write_root(DA, function(f) rep(0.02, length(f)), gate)
  ra <- run_quiet(c("--root", DA, "--out", file.path(DA, "out")), log = file.path(DA, "run.out"))
  L <- st_nominal_limit(0.01, 500); Lr <- st_realised_limit(0.024, 0.01, 500)
  check("(xii a) at alpha 0.01 with B = 500 the nominal limit is 0.023349 (0.023 below, 0.024 above); the realised limit of 0.024 is 0.03053",
        abs(L - 0.0233494) < 1e-6 && 0.023 < L && 0.024 > L && abs(Lr - 0.03053) < 2e-5, sprintf("nominal %.6f, realised %.6f", L, Lr))
  g1 <- ra$gate[cell == c500[1] & form == "sc"]; g2 <- ra$gate[cell == c500[2] & form == "sc"]
  check("(xii a) size 0.024 at 0.01 (B = 500) fails the gate with the nominal MCSE though it passes with the realised one; realised z reported",
        nrow(g1) == 3 && all(!g1$pass) && all(!g1$pass_01) && all(!g1$pass_null_01) && all(g1$pass_realised_01) && all(g1$mcse_reading == "null") &&
          near(g1$z_01, g1$z_null_01, 1e-12) && near(g1$z_null_01, rep(0.014 / sqrt(0.01 * 0.99 / 500), 3), 1e-9) &&
          near(g1$z_realised_01, rep(0.014 / sqrt(0.024 * 0.976 / 500), 3), 1e-9),
        sprintf("%s: z nominal %.4f, z realised %.4f", c500[1], g1$z_null_01[1], g1$z_realised_01[1]))
  check("(xii a) size 0.023 at 0.01 (B = 500) passes under both readings", nrow(g2) == 3 && all(g2$pass) && all(g2$pass_null_01) && all(g2$pass_realised_01),
        sprintf("%s: z nominal %.4f", c500[2], g2$z_null_01[1]))
  check("(xii a) the score form fails in that one null cell and the unit form is chosen", decision_is(ra, "unit", "only the unit form") &&
          all(ra$eligibility[form == "sc"]$failing == 1), ra$decision$reason[1])
  pr <- ra$pairs[hypothesis == "H3" & form == "unit"][1]
  nc <- CT[block == pr$block & cell == pr$cell]
  Bn <- CT[block == nc$null_block & cell == nc$null_cell]$B
  Ln <- st_nominal_limit(0.05, Bn)
  res <- lapply(c(Ln - 5e-4, Ln + 5e-4), function(s) {
    Vm <- copy(ra$V)
    Vm[block == nc$null_block & cell == nc$null_cell & test == "GiViTI" & alpha == 0.05, `:=`(rejection = s, mcse = sqrt(s * (1 - s) / B))]
    list(s = s, m = an_mcnemar(DA, Vm, ra$C, pr), h = an_holds(Vm, ra$C, pr$block, pr$cell, "GiViTI"))
  })
  lo <- res[[1]]; hi <- res[[2]]
  check("(xii a) paired test by default: GiViTI holds size just below 0.05 + 3 sqrt(0.0475 / B) and fails just above, under the realised limit; both z reported",
        isTRUE(lo$m$holds_b) && isFALSE(hi$m$holds_b) && hi$s < st_realised_limit(hi$s, 0.05, Bn) && all(hi$m$mcse_reading == "null") &&
          near(hi$m$z_nominal_b, (hi$s - 0.05) / sqrt(0.0475 / Bn), 1e-9) && near(hi$m$z_realised_b, (hi$s - 0.05) / sqrt(hi$s * (1 - hi$s) / Bn), 1e-9) &&
          isTRUE(lo$h) && isFALSE(hi$h) && grepl("GiViTI fails size .*null MCSE", hi$m$note),
        sprintf("%s/%s, null %s/%s, B %d, limit %.5f", pr$block, pr$cell, nc$null_block, nc$null_cell, Bn, Ln))
})

## (b) E13.2 on a synthetic copy treated as the real root: one family cell without a size-adjusted power
guarded("(xii b)", {
  DB <- file.path(OUTDIR, "synth_e13b"); outb <- file.path(DB, "out")
  SB <- st_write_root(DB, function(f) rep(0.02, length(f)))
  writeLines(c("2026-09-14 08:28:48  launch from step 1 of 19", "2026-09-15 02:00:00  launch finished"), file.path(DB, "launch.log"))
  runb <- function(extra = character(0)) {
    res <- NULL
    msg <- tryCatch({ capture.output(res <- an_main(c("--root", DB, "--out", outb, extra), bat = DB)); "ran" }, error = function(e) conditionMessage(e))
    list(msg = msg, res = res)
  }
  putb <- function(S) { for (b in unique(S$block)) fwrite(S[block == b], file.path(DB, b, "_summary.csv")); unlink(outb, recursive = TRUE) }
  m1 <- runb(c("--mcse", "realised")); m2 <- runb(c("--mcse", "realised", "--force"))
  check("(xii a) the real root refuses --mcse realised (a sensitivity reading, E13.1); --force runs it",
        grepl("refusing the real root", m1$msg) && grepl("--mcse realised", m1$msg, fixed = TRUE) && m2$msg == "ran", st_flat(m1$msg))
  S <- copy(SB)
  S[block == "3" & cell == "link_cloglog_n1000" & startsWith(test, "EDGE."), `:=`(size_adj_power = NA_real_, status = "matched null not run yet")]
  putb(S)
  b1 <- runb(); b2 <- runb("--force")
  want <- sprintf("3 link_cloglog_n1000 EDGE.%s.%s.Grule (matched null not run yet)", rep(bases, each = 2), rep(c("u", "sc"), 3))
  check("(xii b) the refusal lists every missing pair (all 6 of Section A, no '...') and names H3 and H4, not H1",
        grepl("refusing the real root", b1$msg) && all(vapply(want, grepl, logical(1), x = b1$msg, fixed = TRUE)) && !grepl("...", b1$msg, fixed = TRUE) &&
          grepl("Section A families: 6 cell-test pairs in 1 cells", b1$msg, fixed = TRUE) && grepl("H3: 2 cell-test pairs in 1 cells", b1$msg, fixed = TRUE) &&
          grepl("H4: 4 cell-test pairs in 1 cells", b1$msg, fixed = TRUE) && !grepl("H1:", b1$msg, fixed = TRUE), st_flat(b1$msg))
  mv <- fread(file.path(outb, "missing_values.csv"))
  check("(xii b) missing_values.csv holds the 12 pairs (6 Section A, 2 H3, 4 H4)", nrow(mv) == 12 && all(mv$cell == "link_cloglog_n1000") &&
          identical(as.integer(table(mv$set)[c("Section A families", "H3", "H4")]), c(6L, 2L, 4L)))
  fm <- b2$res$families; D <- b2$res$decision; h <- b2$res$hypotheses
  check("(xii b) --force: family 2 (20 of 21 cells) and the macro-average have no mean in every basis and form, the cell is named, no basis is decided",
        b2$msg == "ran" && all(fm[family == "2"]$with_value == 20) && all(is.na(fm[family %in% c("2", "macro")]$mean_power)) &&
          all(fm[family %in% c("2", "macro")]$cells_without_value == "3/link_cloglog_n1000") && all(fm[family %in% c("1", "3", "4", "5", "6")]$cells_without_value == "") &&
          all(is.finite(fm[family %in% c("1", "3", "4", "5", "6")]$mean_power)) && all(D$choice == "no decision") && all(grepl("family 2", D$reason)) &&
          all(grepl("cells without a size-adjusted power: unit 3/link_cloglog_n1000; score 3/link_cloglog_n1000", D$reason, fixed = TRUE)), D$reason[1])
  hx <- h[hypothesis %in% c("H3", "H4", "H4.ruleG")]; hk <- h[hypothesis %in% c("H1", "H2", "H5")]
  check("(xii b) --force: H3, H4 and H4.ruleG have no value and name the cell, with no census count beside H4 either (E13.2); H1, H2, H5 keep their values",
        nrow(hx) == 6 && all(hx$verdict == "no value") && all(hx$final_verdict == "no value") && all(hx$cells_without_value == "3/link_cloglog_n1000") &&
          all(is.na(hx$statistic)) && all(is.na(hx$statistic_available)) &&
          all(hk$verdict != "no value") && all(hk$cells_without_value == ""))
  dl <- readLines(file.path(outb, "rule_A_decision.txt"))
  check("(xii b) rule_A_decision.txt names the cell under family 2",
        any(grepl("without a size-adjusted power: unit [3/link_cloglog_n1000]  score [3/link_cloglog_n1000]", dl, fixed = TRUE)))
})

## (c) E13.2: a matched null that gave no p-value in at least 1 - alpha of its replicates
guarded("(xii c)", {
  Vt <- data.table(block = c(rep("1b", 4), rep("2", 4)), cell = c(rep("nul", 4), rep("alt", 4)), role = rep(c("null", "alternative"), each = 4),
                   test = rep(c("T1", "T2", "T3", "T1"), 2), alpha = rep(c(0.05, 0.05, 0.05, 0.01), 2),
                   declined = c(0.96, 0.95, 0.9499, 0.96, 0, 0, 0, 0), size_adj_power = c(NA, NA, NA, NA, 0.5, 0.5, 0.5, 0.3), status = "")
  Ct <- data.table(block = c("1b", "2"), cell = c("nul", "alt"), null_block = c("", "1b"), null_cell = c("", "nul"))
  Vo <- an_null_declined(Vt, Ct); a <- Vo[block == "2"]
  check("(xii c) null declined 0.96 or 0.95 at alpha 0.05 -> no value; 0.9499 keeps it; 0.96 at alpha 0.01 (below 0.99) keeps it; input unchanged",
        is.na(a$size_adj_power[1]) && is.na(a$size_adj_power[2]) && near(a$size_adj_power[3:4], c(0.5, 0.3)) &&
          all(a$status[1:2] == AN_NULL_DECLINED) && all(a$status[3:4] == "") && near(a$null_declined, c(0.96, 0.95, 0.9499, 0.96)) &&
          all(is.na(Vo[block == "1b"]$null_declined)) && identical(Vt$size_adj_power, c(NA, NA, NA, NA, 0.5, 0.5, 0.5, 0.3)))
  S <- copy(SB)
  S[block == "1b" & cell == "null_crossover_n2000" & test == "EDGE.poly3.sc.Grule", declined := 0.96]
  putb(S)
  c1 <- runb(); c2 <- runb("--force")
  vv <- c2$res$V[block == "4" & cell == "crossover_n2000" & test == "EDGE.poly3.sc.Grule"]
  check("(xii c) a null declined at 0.96: the alternative's power at 0.05 has no value, the refusal lists it, and with --force poly3 has no decision (sym, stk score)",
        grepl("refusing the real root", c1$msg) && grepl(sprintf("4 crossover_n2000 EDGE.poly3.sc.Grule (%s)", AN_NULL_DECLINED), c1$msg, fixed = TRUE) &&
          c2$msg == "ran" && is.na(vv[alpha == 0.05]$size_adj_power) && vv[alpha == 0.05]$status == AN_NULL_DECLINED && is.finite(vv[alpha == 0.01]$size_adj_power) &&
          near(vv$null_declined, rep(0.96, 3)) && c2$res$decision[basis == "poly3"]$choice == "no decision" &&
          grepl("score 4/crossover_n2000", c2$res$decision[basis == "poly3"]$reason, fixed = TRUE) && all(c2$res$decision[basis != "poly3"]$choice == "score"),
        st_flat(c1$msg))
})

## (d) E13.3: a verdict that flips when the size-failing cells are removed
guarded("(xii d)", {
  fu <- c(poly3 = "u", sym = "u", stk = "u")
  f1 <- r1$M[family == 1, .(block, cell, n)]
  tgt <- CT[block == "2" & cell == "probit_base_n9400"]
  nul <- CT[match(paste(f1$block, f1$cell), paste(CT$block, CT$cell)), paste(null_block, null_cell)]
  shared <- nul == paste(tgt$null_block, tgt$null_cell); k <- sum(shared); m <- nrow(f1) - k
  x <- (-0.010 * nrow(f1) + 0.02 * m) / k                 # mean -0.010 over all cells, -0.020 over the others
  h2 <- function(dif, fail) {
    V <- copy(r1$V); ps <- st_getp(V, f1, "EDGE.sym.u.Grule"); st_setp(V, f1, "Stk.joint", ps - dif)
    if (fail) st_failnull(V, tgt$null_block, tgt$null_cell, "Stk.joint")
    an_hypotheses(V, r1$M, r1$H4c, r1$H5c, "Grule", fu, r1$C)$hyp[hypothesis == "H2" & form == "unit"]
  }
  hf <- h2(ifelse(shared, x, -0.02), TRUE)
  check("(xii d) H2 mean -0.010 holds over all cells but -0.020 without the size-failing cell -> final verdict unresolved (size)",
        hf$verdict == "holds" && near(hf$statistic, -0.010, 1e-12) && hf$size_fail_cells == k && grepl("2/probit_base_n9400", hf$size_fail_cell_names, fixed = TRUE) &&
          near(hf$statistic_without_size_failures, -0.02, 1e-12) && hf$verdict_without_size_failures == "against" && hf$final_verdict == "unresolved (size)",
        sprintf("%d cell(s) with Stk.joint failing size in %s/%s: %s", k, tgt$null_block, tgt$null_cell, hf$size_fail_cell_names))
  hn <- h2(rep(-0.010, nrow(f1)), TRUE)
  check("(xii d) the same size failure without a flip keeps the E12.5 verdict (holds) and still lists the cell",
        hn$verdict == "holds" && hn$size_fail_cells == k && near(hn$statistic_without_size_failures, -0.010, 1e-12) && hn$final_verdict == "holds")
  h0 <- h2(ifelse(shared, x, -0.02), FALSE)
  check("(xii d) no size failure: no cell listed, the mean without failures = the mean, final = verdict",
        h0$size_fail_cells == 0 && h0$size_fail_cell_names == "" && near(h0$statistic_without_size_failures, h0$statistic, 1e-12) && h0$final_verdict == "holds")
  ## H4: 19 of 29 cells led; the best partition test of one led census cell fails size -> 18 -> unresolved (size)
  h4c <- r1$H4c; census <- grepl("^census_", h4c$cell)
  led <- census | seq_len(nrow(h4c)) %in% which(!census)[1:11]
  tests7 <- c("EDGE.poly3.u.G10", "HLF.G10", "HL.G10", "HL_w", "PH", "Tsiatis", "Xie")
  V <- copy(r1$V)
  for (j in seq_along(tests7))
    st_setp(V, h4c, tests7[j], ifelse(led, c(0.60, 0.55, 0.50, 0.50, 0.50, 0.40, 0.40)[j], c(0.40, 0.60, 0.50, 0.50, 0.50, 0.40, 0.40)[j]))
  h4 <- function(V) an_hypotheses(V, r1$M, h4c, r1$H5c, "Grule", fu, r1$C)$hyp[hypothesis == "H4" & form == "unit"]
  q0 <- h4(V)
  V1 <- copy(V); st_failnull(V1, "1b", "null_census_logx_n1000", "HLF.G10"); st_failnull(V1, "1b", "null_census_skew_n1000", "PH")
  q1 <- h4(V1)
  check("(xii d) H4 with 19 of 29 cells led: holds, no size failure", sum(led) == 19 && q0$verdict == "holds" && q0$statistic == 19 &&
          q0$headline_cells == 29 && q0$size_fail_cells == 0 && q0$final_verdict == "holds")
  check("(xii d) H4: the best partition test fails size in a led cell -> 18 -> unresolved (size); a failing PH that is not the best test is no size failure",
        q1$verdict == "holds" && q1$statistic == 19 && q1$size_fail_cells == 1 && q1$size_fail_cell_names == "3/census_logx_n1000" &&
          q1$statistic_without_size_failures == 18 && q1$verdict_without_size_failures == "against" && q1$final_verdict == "unresolved (size)", q1$note)
})

## (e) E13.4: the H4 paired tests run over the detectable, non-saturated cells only, Holm over those cells
guarded("(xii e)", {
  DE <- file.path(OUTDIR, "synth_e13e")
  SE <- st_write_root(DE, function(f) rep(0.02, length(f)))
  h4c <- r1$H4c
  nodet <- c("quad_0.01_n1000", "binint_0.1_n1000", "contint_0.1_n1000"); sat <- c("quad_0.4_n1000", "binint_0.7_n1000", "census_joint_n1000")
  t_all <- c(paste0("EDGE.poly3.", c("u", "sc"), ".", rep(c("G10", "Grule"), each = 2)), "HLF.G10", "HLF.Grule", "HL.G10", "HL.Grule", "HL_w", "PH", "Tsiatis", "Xie")
  SE[block == "3" & cell %in% nodet & test %in% t_all & abs(alpha - 0.05) < 1e-9 & subset == "all", size_adj_power := 0.05]
  SE[block == "3" & cell %in% sat & test %in% t_all & abs(alpha - 0.05) < 1e-9 & subset == "all", size_adj_power := 0.99]
  fwrite(SE[block == "3"], file.path(DE, "3", "_summary.csv"))
  set.seed(20260916)
  pcols <- c("EDGE.poly3.u.G10", "EDGE.poly3.sc.G10", "EDGE.poly3.u.Grule", "EDGE.poly3.sc.Grule", "HLF.G10", "HLF.Grule", "HL.G10", "HL.Grule", "HL_w", "PH")
  for (i in seq_len(nrow(h4c))) {
    u0 <- runif(200); P <- data.table(rep = 1:200)
    for (cc in pcols) P[, (cc) := pmin(1, u0 * runif(200, 0.2, 1.8))]         # correlated p-values, different in each test
    fwrite(P, file.path(DE, h4c$block[i], paste0(h4c$cell[i], "_pvalues.csv.gz")), compress = "gzip")
  }
  re <- run_quiet(c("--root", DE, "--out", file.path(DE, "out")), log = file.path(DE, "run.out"))
  MC <- re$mcnemar[startsWith(hypothesis, "H4")]
  ex <- MC[cell %in% c(nodet, sat)]
  check("(xii e) census: the 3 non-detectable and 3 saturated cells are out of the headline, the other 23 in it (both forms, both versions)",
        nrow(MC) == 4 * 29 && all(ex$in_headline %in% FALSE) && all(MC[!cell %in% c(nodet, sat)]$in_headline %in% TRUE) &&
          all(re$census[cell %in% nodet]$detectable %in% FALSE) && all(re$census[cell %in% sat]$saturated %in% TRUE))
  check("(xii e) E13.4: the H4 McNemar families hold only detectable, non-saturated cells (excluded cells keep their 2x2 counts, no p)",
        nrow(MC[used == TRUE & !in_headline %in% TRUE]) == 0 && all(!ex$used) && all(!ex$in_family) && all(is.finite(ex$B)) && all(is.na(ex$p_mcnemar)) &&
          all(grepl("E13.4", ex$note, fixed = TRUE)) && nrow(MC[used == TRUE]) == 4 * 23,
        paste(MC[, .(u = sum(used)), by = .(hypothesis, form)][, sprintf("%s/%s %d", hypothesis, form, u)], collapse = "; "))
  hh <- MC[used == TRUE, .(p_holm, hand = holm_hand(p_mcnemar), m = .N), by = .(hypothesis, arm, form)]
  check("(xii e) Holm over the 23 cells, within H4 and separately within H4.ruleG, = hand computation",
        nrow(hh) == 4 * 23 && near(hh$p_holm, hh$hand, 1e-12) && setequal(unique(hh$hypothesis), c("H4", "H4.ruleG")) && all(hh$m == 23),
        sprintf("Holm-significant %d of %d", sum(MC$holm_reject %in% TRUE), nrow(hh)))
})

R <- do.call(rbind, RES)
cat(sprintf("\n%d checks, %d failed\n", nrow(R), sum(!R$ok)))
if (any(!R$ok)) print(R[!R$ok, ], row.names = FALSE)
cat(if (all(R$ok)) "SELFTEST PASSED\n" else "SELFTEST FAILED\n")
sink()
quit(save = "no", status = if (all(R$ok)) 0L else 1L)
