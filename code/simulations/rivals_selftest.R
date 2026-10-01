## rivals_selftest.R -- checks of run_M_rivals.R (block 8, PREDECLARATION_restructure_battery.md E9) with the dry-run battery
## files as the stored battery (--root battery/dryrun: 10 replicates of block 2 loglog_base_n610 and cauchit_auc_n460 and of
## block 1b null_base_n610). Kept light while the main battery runs: at most 2 R processes compute at any time, BLAS on one
## thread, 6 projection data sets (n = 610) and 2 BAGofT data sets (n = 460) in all.
##   1  the cell list equals E9 (30 cells: names, n, generators, fitted formulas, matched nulls)
##   2  the identity gate passes on the dry-run replicates, fails on perturbed stored values and draws no random number
##   3  refusals before any computation: missing battery file, missing null, too few stored replicates, a test run into
##      battery/8, a finished file with the wrong row count, a .part from another run
##   4  the projection test through the driver (3 replicates, 2 workers) equals direct proj_pvalue calls on the same data and
##      random-number state
##   5  resume: a finished cell is skipped unchanged, a complete .part is renamed, a partial .part is continued
##   6  BAGofT through the driver (1 replicate, 1 worker) gives a p-value in [0, 1] equal to a direct call
##   7  the summary: raw and size-adjusted rejection, McNemar and Holm against a hand computation on constructed p-values, the
##      rivals decided at p < alpha and the comparators at p <= alpha (review F1), and the level of that rule on the two p-value grids
##
##   Rscript rivals_selftest.R            every check (outputs under battery/8_test; an older one goes to battery/_review/rivals)
##   Rscript rivals_selftest.R --quick    no rival is computed: 4 and 6 are left out, and the driver checks of 3 and 5 use a copy of
##                                        the dry-run battery whose stored EDGE.poly3.u.G10 is shifted, so every replicate fails the
##                                        identity gate (outputs under battery/8_test/quick)
## Log <out>/selftest.log, checks <out>/selftest_checks.csv, timings <out>/timings.csv.
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
QUICK <- "--quick" %in% commandArgs(trailingOnly = TRUE)
BAT <- edge_battery()
T8R <- if (QUICK) "battery/8_test/quick" else "battery/8_test"             # as the driver's --out and --root see it
T8  <- file.path(SIMDIR, T8R)
REV <- edge_battery("_review", "rivals")
if (dir.exists(T8)) {
  dir.create(REV, recursive = TRUE, showWarnings = FALSE)
  invisible(file.rename(T8, file.path(REV, paste0(if (QUICK) "8_test_quick_prev_" else "8_test_prev_", format(Sys.time(), "%Y%m%d_%H%M%S")))))
}
dir.create(T8, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }

options(rivals.source_only = TRUE)
source(file.path(SIMDIR, "run_M_rivals.R"))
SETUP <- c("--root", "battery/dryrun", "--out", paste0(T8R, "/inproc"), "--reps", "10", "--workers", "2")
rv_setup(SETUP)

LOG <- file.path(T8, "selftest.log")
RES <- list(); TIM <- list()
say <- function(...) {
  line <- sprintf("%s  %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
  cat(line, "\n", sep = ""); cat(line, "\n", file = LOG, append = TRUE, sep = "")
}
check <- function(id, what, ok, detail = "") {
  ok <- isTRUE(ok)
  RES[[length(RES) + 1]] <<- data.frame(check = id, what = what, ok = ok, detail = detail, stringsAsFactors = FALSE)
  say("%-5s %-4s %s%s", id, if (ok) "ok" else "FAIL", what, if (nzchar(detail)) paste0("  [", detail, "]") else "")
  invisible(ok)
}
timing <- function(what, n, sec, note)
  TIM[[length(TIM) + 1]] <<- data.frame(what = what, n = n, seconds = sec, note = note, stringsAsFactors = FALSE)
RS <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
drv <- function(args, tag, wait = TRUE) {
  lf <- file.path(T8, paste0("cli_", tag, ".log"))
  t0 <- Sys.time()
  st <- system2(RS, c(shQuote(file.path(SIMDIR, "run_M_rivals.R")), args), stdout = lf, stderr = lf, wait = wait)
  list(status = st, log = lf, sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
       text = if (wait && file.exists(lf)) readLines(lf, warn = FALSE) else character(0))
}
args_of <- function(test, cells, reps, workers, root, out, extra = character(0))
  c("--test", test, "--cells", cells, "--reps", as.character(reps), "--workers", as.character(workers), "--root", root,
    "--out", paste0(T8R, "/", out), extra)
has <- function(r, pattern) any(grepl(pattern, r$text, fixed = TRUE))
## the log of a driver started in the background; Windows keeps it locked while that process runs, so read nothing until then
read_log <- function(f) tryCatch(if (file.exists(f)) readLines(f, warn = FALSE) else character(0),
                                 warning = function(w) character(0), error = function(e) character(0))
no_outputs <- function(dir) !length(list.files(dir, pattern = "_pvalues\\.csv\\.gz"))
## a copy of the dry-run files of loglog_base_n610 and its null, with the stored EDGE.poly3.u.G10 shifted from replicate 'from'
fake_root <- function(name, from) {
  fr <- file.path(T8, name); dir.create(file.path(fr, "2"), recursive = TRUE); dir.create(file.path(fr, "1b"), recursive = TRUE)
  Bf <- fread(edge_battery("dryrun", "2", "loglog_base_n610_pvalues.csv.gz"))
  Bf[rep >= from, EDGE.poly3.u.G10 := EDGE.poly3.u.G10 + 0.5]
  fwrite(Bf, file.path(fr, "2", "loglog_base_n610_pvalues.csv.gz"))
  file.copy(edge_battery("dryrun", "1b", "null_base_n610_pvalues.csv.gz"), file.path(fr, "1b"))
  paste0(T8R, "/", name)
}
say("rivals_selftest%s: R %s, BAGofT %s, root battery/dryrun, out %s", if (QUICK) " --quick" else "", getRversion(),
    as.character(utils::packageVersion("BAGofT")), T8)

## ---- 1: the cell list equals E9 ----------------------------------------------------------------------------------------------
E9_B3 <- sprintf("%s_n%d", rep(c("stk_long", "stk_short", "stk_asym", "cauchit", "t4", "loglog"), each = 2), c(500L, 1000L))
E9_B2 <- c("loglog_base_n610", "loglog_base_n1000", "cauchit_auc_n460", "cauchit_auc_n760", "loglog_auc_n380", "loglog_auc_n640",
           "loglog_e12_n1000")
E9 <- rbind(
  data.frame(block = "1a", cell = c("null_link_n500", "null_link_n1000"), null_block = NA, null_cell = NA),
  data.frame(block = "3", cell = E9_B3, null_block = "1a", null_cell = paste0("null_link_n", sub("^.*_n", "", E9_B3))),
  data.frame(block = "1b", cell = paste0("null_", sub("^[a-z]+_", "", E9_B2)), null_block = NA, null_cell = NA),
  data.frame(block = "2", cell = E9_B2, null_block = "1b", null_cell = paste0("null_", sub("^[a-z]+_", "", E9_B2))),
  data.frame(block = "1b", cell = "null_crossover_n1000", null_block = NA, null_cell = NA),
  data.frame(block = "4", cell = "crossover_n1000", null_block = "1b", null_cell = "null_crossover_n1000"), stringsAsFactors = FALSE)
C8 <- Cells8; CA <- battery_cells()
k8 <- paste(C8$block, C8$cell); kE <- paste(E9$block, E9$cell)
check("1.1", "30 cells, no duplicate name", nrow(C8) == 30L && !anyDuplicated(C8$cell) && nrow(E9) == 30L)
check("1.2", "the block and cell names equal E9", setequal(k8, kE), paste(c(setdiff(k8, kE), setdiff(kE, k8)), collapse = ", "))
mE <- match(k8, kE)
check("1.3", "every matched null equals E9 (block 3 -> 1a null_link, block 2 -> 1b null_<design>, crossover -> 1b null_crossover)",
      identical(ifelse(is.na(C8$null_cell), NA, paste(C8$null_block, C8$null_cell)),
                ifelse(is.na(E9$null_cell[mE]), NA, paste(E9$null_block[mE], E9$null_cell[mE]))))
check("1.4", "n in each name equals the cell's n, all n <= 1000",
      all(as.integer(sub("^.*_n", "", C8$cell)) == C8$n) && all(C8$n <= 1000L), paste(sort(unique(C8$n)), collapse = ", "))
check("1.5", "10 nulls and 20 alternatives", sum(C8$role == "null") == 10L && sum(C8$role == "alternative") == 20L)
gen_exp <- ifelse(C8$block == "1a", "dgp_null", ifelse(C8$cell == "null_crossover_n1000", "design_null",
           ifelse(C8$cell == "crossover_n1000", "crossover", "design")))
check("1.6", "generators as the battery (design, dgp_null, design_null, crossover) and every fitted formula y ~ x + d",
      identical(C8$generator, gen_exp) && all(C8$formula == "y ~ x + d"))
check("1.7", "every row is the battery_cells() row (all columns, seed_base included)",
      isTRUE(all.equal(C8, CA[match(k8, paste(CA$block, CA$cell)), ], check.attributes = FALSE)))
check("1.8", "block 2 cells with n <= 1000 in battery_cells() (from nplan.csv) are exactly E9's seven",
      setequal(CA$cell[CA$block == "2" & CA$n <= 1000L], E9_B2))
check("1.9", "every cell has B >= 500 in the battery (replicates 1-500 exist)", all(C8$B >= 500L), paste(range(C8$B), collapse = "-"))
g1 <- lapply(seq_len(nrow(C8)), function(i) rv_regen(1L, as.list(C8[i, ])))
check("1.10", "replicate 1 of every cell: formula, columns x, d, y and n as the table",
      all(vapply(seq_len(nrow(C8)), function(i) identical(paste(deparse(g1[[i]]$dat$f), collapse = ""), C8$formula[i]) &&
                   setequal(names(g1[[i]]$dat$d), c("x", "d", "y")) && nrow(g1[[i]]$dat$d) == C8$n[i], logical(1))))

## ---- 2: the identity gate on the dry-run battery files -------------------------------------------------------------------------
DR <- c("loglog_base_n610", "cauchit_auc_n460", "null_base_n610")
for (cn in DR) {
  cd <- as.list(C8[C8$cell == cn, ]); S <- rv_stored(cd, 10L)
  ok <- logical(10); dif <- numeric(10)
  for (r in 1:10) { g <- rv_regen(r, cd); ok[r] <- rv_identity(g, as.list(S[r, ]), cd$seed_base + r); dif[r] <- abs(g$p - S$p[r]) }
  check(sprintf("2.%d", match(cn, DR)), sprintf("identity gate passes on the 10 dry-run replicates of %s", cn), all(ok),
        sprintf("max |diff| %.1e", max(dif)))
}
ce <- as.list(C8[C8$cell == "loglog_base_n610", ]); S <- rv_stored(ce, 10L); g <- rv_regen(1L, ce); s1 <- as.list(S[1, ])
pert <- function(field, v) { z <- s1; z[[field]] <- v; rv_identity(g, z, ce$seed_base + 1) }
check("2.4", "the gate fails on p + 2e-8, the row of another replicate, p missing, events + 1; passes on p + 5e-9",
      !pert("p", s1$p + 2e-8) && !rv_identity(g, as.list(S[2, ]), ce$seed_base + 1) && !pert("p", NA_real_) &&
        !pert("events", s1$events + 1) && pert("p", s1$p + 5e-9))
RNGkind("L'Ecuyer-CMRG"); set.seed(ce$seed_base + 1); invisible(bt_data(ce)); rs1 <- .Random.seed
invisible(rv_regen(1L, ce)); rs2 <- .Random.seed
check("2.5", "the working fit and EDGE-poly3 draw no random number after the data (the rival continues the data's stream)",
      identical(rs1, rs2))
t0 <- proc.time()[["elapsed"]]; z <- s1; z$p <- s1$p + 1e-3
o <- rv_one(1L, ce, "proj", z); tf <- proc.time()[["elapsed"]] - t0
check("2.6", "an identity failure returns id_ok 0 without running the rival",
      o[["id_ok"]] == 0 && is.na(o[["proj"]]) && is.na(o[["sec"]]) && tf < 5, sprintf("%.2f s", tf))

## ---- 3: refusals before any computation ------------------------------------------------------------------------------------------
r <- drv(args_of("proj", "loglog_base_n1000", 3, 2, "battery/dryrun", "refuse_missing"), "refuse_missing")
check("3.1", "refuses a cell whose battery file does not exist yet", r$status != 0 && has(r, "does not exist yet") &&
        has(r, "block 8 not started") && no_outputs(file.path(T8, "refuse_missing")))
r <- drv(args_of("bagoft", "cauchit_auc_n460", 1, 1, "battery/dryrun", "refuse_null"), "refuse_null")
check("3.2", "refuses an alternative whose matched null's battery file does not exist (null_auc_n460)",
      r$status != 0 && has(r, "matched null null_auc_n460") && no_outputs(file.path(T8, "refuse_null")))
r <- drv(args_of("proj", "loglog_base_n610", 11, 2, "battery/dryrun", "refuse_reps"), "refuse_reps")
check("3.3", "refuses when the battery file has fewer replicates than planned", r$status != 0 && has(r, "replicates 1-11 are needed"))
e <- tryCatch({ rv_setup(c("--reps", "1", "--root", "battery/dryrun", "--out", "battery/8")); "no error" },
              error = function(e) conditionMessage(e))
check("3.4", "a test run cannot write to battery/8 (checked in-process; battery/8 not created)",
      grepl("cannot write to", e) && !dir.exists(edge_battery("8")), e)
rv_setup(SETUP)

## ---- 4: the projection test through the driver (the reference file of 3.5, 3.6 and 5) -------------------------------------------
if (!QUICK) {
  a4 <- args_of("proj", "loglog_base_n610", 3, 2, "battery/dryrun", "proj")
  r4 <- drv(a4, "proj")
  f4 <- file.path(T8, "proj", "loglog_base_n610_proj_pvalues.csv.gz")
  check("4.1", "driver: projection test, 3 replicates of loglog_base_n610, 2 workers, exit 0, file renamed from .part",
        r4$status == 0 && file.exists(f4) && !file.exists(paste0(f4, ".part")) && has(r4, "run finished"), sprintf("%.0f s wall", r4$sec))
  P4 <- if (file.exists(f4)) as.data.frame(fread(f4)) else data.frame()
  check("4.2", "3 rows, identity kept (|diff| <= 1e-8), p-value in [0, 1], no error",
        nrow(P4) == 3 && all(P4$rep == 1:3) && all(P4$id_ok == 1) && all(P4$flag.error == 0) && all(P4$proj >= 0 & P4$proj <= 1) &&
          all(P4$id.absdiff <= 1e-8),
        if (nrow(P4)) sprintf("p %s; %s s", paste(P4$proj, collapse = ", "), paste(round(P4$sec, 1), collapse = ", ")) else "")
  if (nrow(P4)) for (k in 1:3) timing("proj, driver (2 workers)", 610L, P4$sec[k], sprintf("replicate %d", k))
} else {
  a4 <- args_of("proj", "loglog_base_n610", 3, 2, fake_root("fakeroot_all", 1L), "proj")
  r4 <- drv(a4, "proj")
  f4 <- file.path(T8, "proj", "loglog_base_n610_proj_pvalues.csv.gz")
  P4 <- if (file.exists(f4)) as.data.frame(fread(f4)) else data.frame()
  check("Q.1", "driver with 2 workers on a battery copy that fails the gate: 3 rows, id_ok 0, no rival run, file renamed from .part",
        r4$status == 0 && nrow(P4) == 3 && all(P4$id_ok == 0) && all(is.na(P4$proj)) && all(is.na(P4$sec)) &&
          !file.exists(paste0(f4, ".part")) && has(r4, "identity failures 3") && has(r4, "run finished"), sprintf("%.1f s", r4$sec))
  lb <- file.path(T8, "cli_background.log")
  invisible(drv(args_of("proj", "loglog_base_n1000", 3, 2, "battery/dryrun", "refuse_background"), "background", wait = FALSE))
  dl <- Sys.time() + 120
  while (Sys.time() < dl && !any(grepl("not started", read_log(lb)))) Sys.sleep(1)
  check("Q.2", "a driver call started in the background writes its log, readable once it ends (the mechanism of 6)",
        any(grepl("not started", read_log(lb))))
}

## ---- 5: resume, and the refusals that need a finished file (no rival computed) ---------------------------------------------------
md5 <- unname(tools::md5sum(f4))
r <- drv(a4, "proj_resume")
check("5.1", "a finished cell is skipped: exit 0, file unchanged, nothing run", r$status == 0 && unname(tools::md5sum(f4)) == md5 &&
        has(r, "skipped") && has(r, "nothing to run"), sprintf("%.1f s", r$sec))
dir.create(file.path(T8, "proj_part"), showWarnings = FALSE); fp <- file.path(T8, "proj_part", basename(f4))
invisible(file.copy(f4, paste0(fp, ".part")))
a5 <- a4; a5[which(a5 == "--out") + 1L] <- paste0(T8R, "/proj_part")
r <- drv(a5, "proj_part_complete")
check("5.2", "a complete .part is renamed without computing", r$status == 0 && file.exists(fp) && !file.exists(paste0(fp, ".part")) &&
        has(r, "complete .part found") && isTRUE(all.equal(as.data.frame(fread(fp)), P4)), sprintf("%.1f s", r$sec))
## a partial .part (replicate 1 of the reference file) on a battery copy whose stored EDGE.poly3.u.G10 is shifted from replicate 2:
## the driver keeps replicate 1 from the .part and runs 2 and 3, which fail the identity gate, so no rival runs
dir.create(file.path(T8, "proj_mid"), showWarnings = FALSE); fm <- file.path(T8, "proj_mid", basename(f4))
fwrite(P4[1, , drop = FALSE], paste0(fm, ".part"), compress = "gzip")
r <- drv(args_of("proj", "loglog_base_n610", 3, 2, fake_root("fakeroot_from2", 2L), "proj_mid"), "proj_part_partial")
Pm <- if (file.exists(fm)) as.data.frame(fread(fm)) else data.frame()
check("5.3", "a partial .part is continued: replicate 1 kept from it, 2-3 run (identity failures, rival not run)",
      r$status == 0 && nrow(Pm) == 3 && has(r, "2 of 3 replicates to run (1 rows resumed from .part)") &&
        isTRUE(all.equal(unlist(Pm[1, ]), unlist(P4[1, ]), check.attributes = FALSE)) && all(Pm$id_ok[2:3] == 0) &&
        all(is.na(Pm$proj[2:3])) && all(is.na(Pm$sec[2:3])), sprintf("%.1f s", r$sec))
dir.create(file.path(T8, "wrong"), showWarnings = FALSE); fw <- file.path(T8, "wrong", basename(f4))
fwrite(P4[1:2, ], fw); md5w <- unname(tools::md5sum(fw))
a35 <- a4; a35[which(a35 == "--out") + 1L] <- paste0(T8R, "/wrong")
r <- drv(a35, "refuse_rows")
check("3.5", "refuses a finished file with the wrong row count (2 rows, 3 planned) and leaves it", r$status != 0 &&
        has(r, "exists with 2 rows but 3 replicates are planned") && unname(tools::md5sum(fw)) == md5w)
dir.create(file.path(T8, "wrongpart"), showWarnings = FALSE); fq <- file.path(T8, "wrongpart", basename(f4))
Pw <- P4[1, , drop = FALSE]; Pw$rep <- 7; fwrite(Pw, paste0(fq, ".part"), compress = "gzip")
a36 <- a4; a36[which(a36 == "--out") + 1L] <- paste0(T8R, "/wrongpart")
r <- drv(a36, "refuse_part")
check("3.6", "refuses a .part with a replicate outside the planned range", r$status != 0 && has(r, "does not belong to this run") &&
        !file.exists(fq))

if (!QUICK) {
  ## ---- 6 (start): BAGofT through the driver in the background, 1 worker ------------------------------------------------------------
  a6 <- args_of("bagoft", "cauchit_auc_n460", 1, 1, "battery/dryrun", "bagoft", "--allow-missing-null")
  invisible(drv(a6, "bagoft", wait = FALSE))
  say("BAGofT driver started in the background (1 worker); the direct projection calls run meanwhile (2 processes computing)")

  ## ---- 4 (cont.): direct projection calls on the same data and random-number state -----------------------------------------------
  dp <- ds <- numeric(3)
  for (k in 1:3) {
    RNGkind("L'Ecuyer-CMRG"); set.seed(ce$seed_base + k)
    dat <- bt_data(ce)
    fit <- suppressWarnings(glm(dat$f, data = dat$d, family = binomial()))
    X <- model.matrix(fit); y <- fit$y
    if (k == 1) {                                                     # the n x n kernel alone (draws no random number)
      tA <- system.time(invisible(proj_build_A(X)))[["elapsed"]]
      timing("proj_build_A alone, direct", 610L, tA, "the fixed n x n kernel; the rest of a data set is the 250 refits")
    }
    tk <- system.time(z <- proj_pvalue(y, X, B = 250))[["elapsed"]]
    dp[k] <- z$p_value; ds[k] <- z$stat
    timing("proj, direct (with 1 BAGofT worker running)", 610L, tk, sprintf("replicate %d", k))
  }
  check("4.3", "direct proj_pvalue(y, X, B = 250) on the regenerated data with the same RNG state gives the driver's p-values",
        nrow(P4) == 3 && all(abs(dp - P4$proj) < 1e-12), sprintf("direct %s", paste(dp, collapse = ", ")))
  check("4.4", "... and the driver's statistic (relative 1e-12)",
        nrow(P4) == 3 && all(abs(ds - P4$aux.proj_stat) <= 1e-12 * pmax(1, abs(ds))),
        sprintf("max rel diff %.1e", if (nrow(P4)) max(abs(ds - P4$aux.proj_stat) / pmax(1, abs(ds))) else NA))

  ## ---- 6 (cont.): a direct BAGofT call, then the driver's result -----------------------------------------------------------------
  ceb <- as.list(C8[C8$cell == "cauchit_auc_n460", ])
  RNGkind("L'Ecuyer-CMRG"); set.seed(ceb$seed_base + 1)
  datb <- bt_data(ceb)
  tb <- system.time(zb <- suppressWarnings(suppressMessages(
    BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = datb$f, link = "logit"), data = datb$d))))[["elapsed"]]
  timing("BAGofT, direct (with the BAGofT driver worker running)", 460L, tb, "package defaults")
  say("direct BAGofT: p %s in %.0f s; waiting for the driver", format(zb$p.value), tb)
  fb <- file.path(T8, "bagoft", "cauchit_auc_n460_bagoft_pvalues.csv.gz"); lfb <- file.path(T8, "cli_bagoft.log")
  lines_b <- function() read_log(lfb)
  deadline <- Sys.time() + 5400
  while (Sys.time() < deadline && !any(grepl("run finished|Execution halted|not started", lines_b()))) Sys.sleep(10)
  Pb <- if (file.exists(fb)) as.data.frame(fread(fb)) else data.frame()
  check("6.1", "driver: BAGofT on replicate 1 of cauchit_auc_n460 (1 worker) finishes with id_ok 1 and a p-value in [0, 1]",
        any(grepl("run finished", lines_b())) && nrow(Pb) == 1 && Pb$id_ok == 1 && Pb$flag.error == 0 && Pb$BAGofT >= 0 && Pb$BAGofT <= 1,
        if (nrow(Pb)) sprintf("p %s, %.0f s", format(Pb$BAGofT), Pb$sec) else "no file")
  check("6.2", "the driver's BAGofT p-value (and its mean, median and minimum split p-values) equals the direct call",
        nrow(Pb) == 1 && abs(Pb$BAGofT - zb$p.value) < 1e-12 && abs(Pb$aux.bagoft_p_median - zb$p.value2) < 1e-12 &&
          abs(Pb$aux.bagoft_p_min - zb$p.value3) < 1e-12 && abs(Pb$aux.bagoft_pmean - zb$pmean) < 1e-12,
        sprintf("direct %s / %s / %s", format(zb$p.value), format(zb$p.value2), format(zb$p.value3)))
  check("6.3", "the log records the allowed missing null of this test run", any(grepl("test run, allowed:", lines_b(), fixed = TRUE)))
  if (nrow(Pb)) timing("BAGofT, driver (1 worker, with the direct call running)", 460L, Pb$sec, "package defaults")
}

## ---- 7: the summary on constructed p-values ----------------------------------------------------------------------------------------
fk <- file.path(T8, "summary_fake"); dir.create(fk, showWarnings = FALSE)
mk <- function(cc, p, idok) {
  M <- as.data.frame(matrix(NA_real_, 10, length(rv_names("proj")), dimnames = list(NULL, rv_names("proj"))))
  M$rep <- 1:10; M$seed <- cc$seed_base + 1:10; M$n <- cc$n; M$events <- 1; M$id_ok <- idok; M$flag.degenerate <- 0
  M$proj <- ifelse(idok == 1, p, NA); M$sec <- ifelse(idok == 1, 1:10, NA); M$flag.error <- as.numeric(idok == 1 & is.na(p))
  M
}
cn_alt <- as.list(C8[C8$cell == "loglog_base_n610", ]); cn_nul <- as.list(C8[C8$cell == "null_base_n610", ])
## p-values equal to 0.01, 0.05 and 0.10 (grid points of a Monte Carlo p-value): the rival does not reject them (review F1)
p_nul <- c(0.004, 0.05, 0.05, 0.35, 0.4, 0.5, 0.6, NA, 0.8, 0.9); id_nul <- c(rep(1, 9), 0)
p_alt <- c(0.001, 0.01, 0.04, 0.05, 0.2, NA, 0.03, 0.5, 0.10, 0.7); id_alt <- c(rep(1, 9), 0)
fwrite(mk(cn_nul, p_nul, id_nul), file.path(fk, "null_base_n610_proj_pvalues.csv.gz"))
fwrite(mk(cn_alt, p_alt, id_alt), file.path(fk, "loglog_base_n610_proj_pvalues.csv.gz"))
rv_setup(c("--root", "battery/dryrun", "--out", paste0(T8R, "/summary_fake")))
invisible(capture.output(rv_summary(list(poly3 = "u", sym = "u"))))
SS <- as.data.frame(fread(file.path(fk, "_summary.csv"))); PP <- as.data.frame(fread(file.path(fk, "_paired.csv")))
## the same numbers by hand
ka <- 1:9; kn <- 1:9; pa <- p_alt[ka]; pn <- p_nul[kn]
Ba <- as.data.frame(fread(edge_battery("dryrun", "2", "loglog_base_n610_pvalues.csv.gz")))
Bn <- as.data.frame(fread(edge_battery("dryrun", "1b", "null_base_n610_pvalues.csv.gz")))
rowS <- function(t, a) SS[SS$cell == "loglog_base_n610" & SS$rival == "proj" & SS$test == t & abs(SS$alpha - a) < 1e-12, ]
rj <- function(p, a) is.finite(p) & p < a                                 # the rival's rule (review F1): p = alpha does not reject
ok7 <- TRUE; det <- character(0)
for (a in c(0.01, 0.05, 0.10)) {
  crit <- sort(ifelse(is.finite(pn), pn, 1))[max(1, ceiling(a * length(pn)))]
  z <- rowS("proj", a)
  e <- c(mean(rj(pa, a)), mean(is.finite(pa) & pa <= crit), mean(rj(pn, a)), mean(!is.finite(pa)), sum(rj(pa, a)) / sum(is.finite(pa)),
         mean(is.finite(pn) & pn <= crit))
  got <- if (nrow(z) == 1) c(z$rejection, z$size_adj_power, z$null_size, z$declined, z$rejection_given_p, z$size_adj_null_level) else rep(NA, 6)
  if (nrow(z) != 1 || any(abs(got - e) > 1e-12) || z$decision_rule != "p < alpha")
    { ok7 <- FALSE; det <- c(det, sprintf("proj alpha %g: got %s, hand %s", a, paste(got, collapse = "/"), paste(e, collapse = "/"))) }
  if (mean(rj(pa, a)) == mean(is.finite(pa) & pa <= a))                   # the constructed values separate the two rules at every alpha
    { ok7 <- FALSE; det <- c(det, sprintf("alpha %g does not separate p < alpha from p <= alpha", a)) }
  zc <- rowS("EDGE.poly3.u.G10", a)
  pc <- Ba$EDGE.poly3.u.G10[match(ka, Ba$rep)]; pcn <- Bn$EDGE.poly3.u.G10[match(kn, Bn$rep)]
  critc <- sort(ifelse(is.finite(pcn), pcn, 1))[max(1, ceiling(a * length(pcn)))]
  ec <- c(mean(is.finite(pc) & pc <= a), mean(is.finite(pc) & pc <= critc))
  if (nrow(zc) != 1 || any(abs(c(zc$rejection, zc$size_adj_power) - ec) > 1e-12) || zc$decision_rule != "p <= alpha")
    { ok7 <- FALSE; det <- c(det, sprintf("EDGE.poly3.u.G10 alpha %g", a)) }
}
z <- rowS("proj", 0.05)
check("7.1", "summary: raw rejection (rival at p < alpha, comparator at p <= alpha), size-adjusted power and its null level, null size, declined rate equal a hand computation",
      ok7 && nrow(z) == 1 && z$reps_kept == 9 && z$identity_failures == 1 && z$null_reps_kept == 9 && abs(z$rival_median_sec - 5) < 1e-12,
      paste(det, collapse = "; "))
hand <- do.call(rbind, lapply(RV_COMP, function(t) {
  r1 <- rj(pa, 0.05); pcm <- Ba[[t]][match(ka, Ba$rep)]; r2 <- is.finite(pcm) & pcm <= 0.05
  b <- sum(r1 & !r2); c <- sum(!r1 & r2)
  hc <- mean(is.finite(Bn[[t]]) & Bn[[t]] <= 0.05) <= 0.05 + 3 * sqrt(0.0475 / nrow(Bn))
  data.frame(comparator = t, b = b, c = c, p = if (b + c == 0) 1 else 2 * min(0.5, pbinom(min(b, c), b + c, 0.5)),
             inh = hc && !grepl("\\.sc\\.", t), stringsAsFactors = FALSE)
}))
hand$holm <- NA_real_; hand$holm[hand$inh] <- p.adjust(hand$p[hand$inh], "holm")
pr8 <- PP[PP$cell == "loglog_base_n610" & PP$rival == "proj", ]; pr8 <- pr8[match(hand$comparator, pr8$comparator), ]
check("7.2", "paired: discordant counts, exact McNemar p, Holm family (size held, forms u) and Holm p equal a hand computation",
      nrow(pr8) == length(RV_COMP) && all(pr8$rival_only == hand$b) && all(pr8$comparator_only == hand$c) &&
        all(abs(pr8$mcnemar_p - hand$p) < 1e-12) && identical(as.logical(pr8$in_holm), hand$inh) &&
        isTRUE(all.equal(pr8$holm_p, hand$holm)) && all(pr8$rival_holds_size) && sum(hand$inh) >= 1 &&
        all(abs(pr8$rival_size - 1 / 9) < 1e-12) && mean(is.finite(pn) & pn <= 0.05) > 0.05 + 3 * sqrt(0.0475 / 9) &&
        all(pr8$rival_rule == "p < alpha") && all(pr8$comparator_rule == "p <= alpha"),
      sprintf("%d in the Holm family; b/c %s; rival size 1/9 at p < 0.05 holds (p <= 0.05 would give 3/9 and fail)", sum(hand$inh),
              paste(sprintf("%d/%d", hand$b, hand$c), collapse = " ")))
pn8 <- PP[PP$cell == "null_base_n610", ]
check("7.3", "paired rows of the null cell are reported and kept out of the Holm family", nrow(pn8) == length(RV_COMP) && !any(pn8$in_holm))
## the rule on the two grids, with the p-values formed as the packages form them: BAGofT mean(observed > simulated) over nsim = 100,
## the projection test mean(T_boot >= T) over B = 250; every rank of the observed statistic equally likely (no computation of a rival)
kb <- 0:100; pb <- vapply(kb, function(k) mean(c(rep(TRUE, k), rep(FALSE, 100 - k))), 0)
jp <- 0:250; pp <- vapply(jp, function(j) mean(c(rep(TRUE, j), rep(FALSE, 250 - j))), 0)
lv <- function(p, mc) vapply(RV_ALPHAS, function(a) mean(rv_rej(p, a, mc)), 0)
l75 <- list(bag = lv(pb, TRUE) * 101, bag_le = lv(pb, FALSE) * 101, proj = lv(pp, TRUE) * 251, proj_le = lv(pp, FALSE) * 251)
check("7.5", "the rivals' rule on the grids: BAGofT rejects 1, 5, 10 of 101 ranks at 0.01, 0.05, 0.10 (p <= alpha: 2, 6, 11), the decisions of (k + 1)/101 <= alpha; projection 3, 13, 25 of 251 (p <= alpha: 3, 13, 26)",
      isTRUE(all.equal(l75$bag, c(1, 5, 10))) && isTRUE(all.equal(l75$bag_le, c(2, 6, 11))) &&
        all(vapply(RV_ALPHAS, function(a) identical(rv_rej(pb, a, TRUE), (kb + 1) / 101 <= a), TRUE)) &&
        isTRUE(all.equal(l75$proj, c(3, 13, 25))) && isTRUE(all.equal(l75$proj_le, c(3, 13, 26))) &&
        all(RV_RULE == c("p < alpha", "p <= alpha")),
      paste(vapply(names(l75), function(k) sprintf("%s %s", k, paste(round(l75[[k]], 6), collapse = "/")), ""), collapse = "; "))
rv_setup(SETUP)
r <- drv(c("--summary", "--root", "battery/dryrun", "--out", paste0(T8R, "/proj")), "summary_cli")
S7 <- if (file.exists(file.path(T8, "proj", "_summary.csv"))) as.data.frame(fread(file.path(T8, "proj", "_summary.csv"))) else data.frame()
st7 <- if (QUICK) "no replicate kept" else "matched null not run yet"      # quick: every replicate of the reference file failed the gate
check("7.4", sprintf("--summary on the driver's file: exit 0; the alternative's rows are marked '%s'", st7),
      r$status == 0 && nrow(S7) > 0 && all(S7$status[S7$test == "proj"] == st7))

## ---- done ---------------------------------------------------------------------------------------------------------------------
R <- do.call(rbind, RES)
fwrite(R, file.path(T8, "selftest_checks.csv"))
if (length(TIM)) fwrite(do.call(rbind, TIM), file.path(T8, "timings.csv"))
say("rivals_selftest%s: %d checks, %d failed", if (QUICK) " --quick" else "", nrow(R), sum(!R$ok))
quit(save = "no", status = if (all(R$ok)) 0L else 1L)
