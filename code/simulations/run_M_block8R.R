## run_M_block8R.R -- block 8R: le Cessie, the Liu projection test and BAGofT on omitted terms and rough misfit.
## Contract: paper_EDGE/theory/PREDECLARATION_block8R_rivals_rough_omitted.md (sha256 586d8f51...), frozen before any
## 8R scenario was run.
##
##   Rscript run_M_block8R.R [--tests lecessie,proj,bagoft] [--cells a,b] [--workers 20] [--reps N] [--out battery/8R] [--force]
##   Rscript run_M_block8R.R --summary [--out battery/8R]
##
## Scenarios (section 1): le Cessie on the 56 omitted-term and rough alternatives with n <= 1000 and their 11 nulls,
## replicates 1-500; the Liu projection test on the 19 of them at n = 500 and their 4 nulls, 1-500; BAGofT on the three
## rough shapes and the median rung of each omitted-term ladder at n = 500 and the same 4 nulls, 1-200.
##
## One replicate: block 8's path (rv_regen: set.seed(seed_base + rep), bt_data, bt_fit, EDGE-poly3 unit at G = 10) and
## block 8's identity gate (rv_identity); the rival then continues the replicate's random-number stream.
##   lecessie  block 8L's l8_lecessie() -- run.all.gof(fit, tests = "le-Cessie"), checked against its definition
##   proj      proj_pvalue(y, X, B = 250), as block 8
##   bagoft    BAGofT 1.0.0 at its defaults with parFun = the package's parRF with its one-covariate line given
##             drop = FALSE and nothing else changed (section 2; checked before the run)
## Output (<out> = battery/8R): <cell>_<test>_pvalues.csv.gz built in .part, _progress.log, _lock; --summary writes
## _summary.csv and _paired.csv (rules of section 3, written above r8_summary).
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
options(block8L.source_only = TRUE)
suppressMessages(source(file.path(SIMDIR, "run_M_block8L.R")))    # run_M_rivals.R too: rv_regen, rv_identity, proj_pvalue, ...
BAT  <- edge_battery()
ROOT <- BAT

R8_TESTS   <- c("lecessie", "proj", "bagoft")
R8_REPS    <- c(lecessie = 500L, proj = 500L, bagoft = 200L)
R8_BATCH   <- c(lecessie = 25L, proj = 4L, bagoft = 1L)                    # replicates per worker per batch
R8_GRID    <- c(lecessie = FALSE, proj = TRUE, bagoft = TRUE)              # grid p-values reject at p < alpha
R8_PCOL    <- c(lecessie = "lecessie", proj = "proj", bagoft = "BAGofT")
R8_AUX     <- list(lecessie = c("lc.stat", "lc.df"), proj = "aux.proj_stat", bagoft = "BAGofT.p3")
R8_BAG_ALT <- c("rough_osc2_n500", "rough_osc4_n500", "rough_sawtooth_n500", "binint_0.3_n500", "contint_0.5_n500", "quad_0.05_n500")
R8_PAIR    <- "EDGE.poly3.u.Grule"

## ---- BAGofT's partition function with its one-covariate line fixed (the DeepGOF study's parRF_dropfix.R) ----------------
r8_make_parRF_dropfix <- function() {
  src <- deparse(BAGofT::parRF)
  bad <- "datRf <- Train.data[, -which(names(Train.data) == Rsp)]"
  i <- grep(bad, src, fixed = TRUE)
  if (length(i) != 1L) stop("parRF source changed: expected exactly one match for the line to patch")
  src[i] <- sub("== Rsp)]", "== Rsp), drop = FALSE]", src[i], fixed = TRUE)
  eval(parse(text = src), envir = asNamespace("BAGofT"))
}
## the patched source equals the original with that one subset given drop = FALSE and nothing else changed. Compared with
## whitespace collapsed, because deparse() re-wraps the longer line.
r8_patch_check <- function() {
  flat <- function(f) gsub("\\s+", " ", paste(deparse(f), collapse = " "))
  a <- flat(BAGofT::parRF); b <- flat(r8_make_parRF_dropfix())
  bad <- "datRf <- Train.data[, -which(names(Train.data) == Rsp)]"
  good <- "datRf <- Train.data[, -which(names(Train.data) == Rsp), drop = FALSE]"
  lengths(regmatches(a, gregexpr(bad, a, fixed = TRUE))) == 1L &&
    identical(sub(bad, good, a, fixed = TRUE), b)
}

## ---- the scenarios (section 1) -------------------------------------------------------------------------------------------
r8_plan <- function() {
  C <- battery_cells()
  M <- unique(data.table::fread(edge_battery("analysis", "rule_A_membership.csv"))[, .(cell, family)], by = "cell")
  fam <- M$family[match(C$cell, M$cell)]
  alt <- C[!is.na(fam) & fam %in% c(4, 5) & C$n <= 1000 & C$role == "alternative", ]
  nk  <- unique(paste(alt$null_block, alt$null_cell))
  nul <- C[paste(C$block, C$cell) %in% nk, ]
  if (nrow(alt) != 56L || nrow(nul) != 11L || any(nul$role != "null"))
    stop("block 8R: expected 56 alternatives and 11 nulls, found ", nrow(alt), " and ", nrow(nul))
  mk <- function(rows, test) if (nrow(rows)) data.frame(test = test, block = rows$block, cell = rows$cell, role = rows$role,
                                                         n = rows$n, null_block = rows$null_block, null_cell = rows$null_cell,
                                                         R = R8_REPS[[test]], stringsAsFactors = FALSE)
  P <- rbind(mk(rbind(nul, alt), "lecessie"),
             mk(rbind(nul[nul$n == 500L, ], alt[alt$n == 500L, ]), "proj"),
             mk(rbind(nul[nul$n == 500L, ], alt[alt$cell %in% R8_BAG_ALT, ]), "bagoft"))
  cnt <- table(P$test, P$role)
  if (cnt["lecessie", "alternative"] != 56 || cnt["lecessie", "null"] != 11 || cnt["proj", "alternative"] != 19 ||
      cnt["proj", "null"] != 4 || cnt["bagoft", "alternative"] != 6 || cnt["bagoft", "null"] != 4)
    stop("block 8R: the plan does not match section 1")
  rownames(P) <- NULL
  P
}
r8_cell <- function(ce) as.list(battery_cells()[battery_cells()$cell == ce$cell & battery_cells()$block == ce$block, ])

## ---- one replicate ---------------------------------------------------------------------------------------------------------
r8_names <- function(test) c("rep", "seed", "n", "events", "id.stored", "id.regen", "id.absdiff", "id_ok", "flag.degenerate",
                             R8_PCOL[[test]], R8_AUX[[test]], "sec", "flag.error")
r8_one <- function(rep, cell, test, st) {
  nm <- r8_names(test)
  out <- setNames(rep(NA_real_, length(nm)), nm)
  seed <- cell$seed_base + rep
  out["rep"] <- rep; out["seed"] <- seed; out["id_ok"] <- 0; out["id.stored"] <- suppressWarnings(as.numeric(st$p))
  g <- tryCatch(rv_regen(rep, cell), error = function(e) NULL)
  if (is.null(g)) return(out)
  out["n"] <- g$n; out["events"] <- g$events; out["id.regen"] <- g$p; out["flag.degenerate"] <- as.numeric(g$degenerate)
  if (is.finite(out[["id.stored"]]) && is.finite(g$p)) out["id.absdiff"] <- abs(out[["id.stored"]] - g$p)
  if (!rv_identity(g, st, seed)) return(out)                         # identity failure: the rival is not run
  out["id_ok"] <- 1
  if (is.null(g$fq)) { out["flag.error"] <- as.numeric(!g$degenerate); return(out) }
  t0 <- proc.time()[["elapsed"]]
  if (test == "lecessie") {
    lc <- l8_lecessie(g$fq$fit)
    p <- lc[["p"]]; out["lc.stat"] <- lc[["stat"]]; out["lc.df"] <- lc[["df"]]
  } else if (test == "proj") {
    r <- tryCatch(proj_pvalue(g$fq$fit$y, stats::model.matrix(g$fq$fit), B = RV_BPROJ), error = function(e) NULL)
    p <- rv_num1(r$p_value); out["aux.proj_stat"] <- rv_num1(r$stat)
  } else {
    r <- tryCatch(suppressWarnings(suppressMessages(BAGofT::BAGofT(
      testModel = BAGofT::testGlmBi(formula = g$dat$f, link = "logit"), data = g$dat$d, parFun = PARRF_FIX()))),
      error = function(e) NULL)
    p <- rv_num1(r$p.value); out["BAGofT.p3"] <- rv_num1(r$p.value3)
  }
  out["sec"] <- proc.time()[["elapsed"]] - t0
  out[R8_PCOL[[test]]] <- if (is.finite(p) && p >= 0 && p <= 1) p else NA_real_
  out["flag.error"] <- as.numeric(!is.finite(out[[R8_PCOL[[test]]]]))
  out
}
r8_task <- function(task, cell, test) r8_one(task$rep, cell, test, task$st)

## ---- summary --------------------------------------------------------------------------------------------------------------
## Rules (section 3): a grid test rejects at p < alpha, le Cessie and EDGE at p <= alpha; no p-value = no rejection. Size
## band 0.05 +/- 3 sqrt(0.05 x 0.95 / B). Size-adjusted power: critical value = 5% quantile (type 1, no p-value = 1) of the
## same test's matched-null p-values on its own kept replicates, EDGE likewise on the same replicates; power = share with
## p <= critical value (block 8's rule). Lead or tie: EDGE > rival - 0.01. Pairs: exact McNemar on the raw decisions at 0.05,
## Holm across every pair of block 8R.
r8_rej <- function(p, grid) if (grid) is.finite(p) & p < 0.05 - 1e-9 else is.finite(p) & p <= 0.05
r8_summary <- function(OUT) {
  P <- r8_plan()
  rd <- function(f) if (file.exists(f)) as.data.frame(data.table::fread(f)) else NULL
  ## the battery files are large and each null is read once per alternative: keep the two columns needed, and cache
  cache <- new.env(parent = emptyenv())
  rd_bat <- function(ce) {
    f <- rv_bat_path(ce)
    key <- basename(f)
    if (!exists(key, envir = cache, inherits = FALSE))
      assign(key, if (file.exists(f)) as.data.frame(data.table::fread(f, select = c("rep", R8_PAIR))) else NULL, envir = cache)
    get(key, envir = cache)
  }
  path <- function(cell, test) file.path(OUT, sprintf("%s_%s_pvalues.csv.gz", cell, test))
  S <- list(); PR <- list()
  for (i in seq_len(nrow(P))) {
    ce <- as.list(P[i, ]); test <- ce$test; pc <- R8_PCOL[[test]]; grid <- R8_GRID[[test]]
    X <- rd(path(ce$cell, test)); if (is.null(X)) next
    K <- X[X$id_ok %in% 1, , drop = FALSE]
    B <- rd_bat(ce); pe <- B[[R8_PAIR]][match(K$rep, B$rep)]
    alt <- ce$role == "alternative"
    crit_r <- crit_e <- NA_real_
    if (alt) {
      NX <- rd(path(ce$null_cell, test))
      if (!is.null(NX)) {
        NK <- NX[NX$id_ok %in% 1, , drop = FALSE]
        crit_r <- l8_crit(NK[[pc]])
        BN <- rd_bat(list(block = ce$null_block, cell = ce$null_cell))
        crit_e <- l8_crit(BN[[R8_PAIR]][match(NK$rep, BN$rep)])
      }
    }
    rate <- mean(r8_rej(K[[pc]], grid)); band <- 0.05 + c(-3, 3) * sqrt(0.05 * 0.95 / nrow(K))
    S[[length(S) + 1]] <- data.frame(test = test, cell = ce$cell, role = ce$role, n = ce$n, null_cell = ce$null_cell,
      reps_in_file = nrow(X), reps_kept = nrow(K), identity_failures = sum(X$id_ok %in% 0),
      max_identity_absdiff = suppressWarnings(max(K$id.absdiff, na.rm = TRUE)), no_pvalue = sum(!is.finite(K[[pc]])),
      rejection = rate, band_lo = band[1], band_hi = band[2],
      in_band = if (!alt) rate >= band[1] && rate <= band[2] else NA,
      rival_size_adj = if (alt) mean(is.finite(K[[pc]]) & K[[pc]] <= crit_r) else NA_real_,
      edge_rejection = mean(r8_rej(pe, FALSE)),
      edge_size_adj = if (alt) mean(is.finite(pe) & pe <= crit_e) else NA_real_,
      median_sec = stats::median(K$sec, na.rm = TRUE), stringsAsFactors = FALSE)
    r1 <- r8_rej(K[[pc]], grid); r2 <- r8_rej(pe, FALSE); nb <- sum(r1 & !r2); nc <- sum(!r1 & r2)
    PR[[length(PR) + 1]] <- data.frame(test = test, cell = ce$cell, role = ce$role, shared_reps = nrow(K),
      rival_rejection = mean(r1), edge_rejection = mean(r2), rival_only = nb, edge_only = nc,
      mcnemar_p = if (nb + nc == 0) 1 else stats::binom.test(nb, nb + nc, 0.5)$p.value, stringsAsFactors = FALSE)
  }
  S <- data.table::rbindlist(S); PR <- data.table::rbindlist(PR)
  if (nrow(PR)) { PR$holm_p <- stats::p.adjust(PR$mcnemar_p, method = "holm"); PR$holm_reject_05 <- PR$holm_p <= 0.05 }
  data.table::fwrite(S, file.path(OUT, "_summary.csv")); data.table::fwrite(PR, file.path(OUT, "_paired.csv"))
  cat(sprintf("summary: %d rows; paired: %d rows (Holm across all of block 8R)\n", nrow(S), nrow(PR)))
  invisible(list(summary = S, paired = PR))
}

## ---- options, guards, cluster -------------------------------------------------------------------------------------------
r8_opts <- function(a = commandArgs(TRUE)) {
  opt <- list(); i <- 1L
  while (i <= length(a)) {
    if (!startsWith(a[i], "--")) stop("unexpected argument ", a[i])
    key <- sub("^--", "", a[i])
    if (key %in% c("force", "summary")) { opt[[key]] <- TRUE; i <- i + 1L }
    else { if (i == length(a)) stop("--", key, " needs a value"); opt[[key]] <- a[i + 1L]; i <- i + 2L }
  }
  bad <- setdiff(names(opt), c("tests", "cells", "workers", "reps", "out", "force", "summary"))
  if (length(bad)) stop("unknown option --", bad[1])
  opt
}
r8_setup <- function(opt) {
  OPT8 <<- opt
  TEST8 <<- !is.null(opt$reps)
  REPS8 <<- if (TEST8) as.integer(opt$reps) else NA_integer_
  if (TEST8 && (is.na(REPS8) || REPS8 < 1L)) stop("--reps needs a positive integer")
  out <- l8_abs(if (!is.null(opt$out)) opt$out else edge_battery(if (TEST8) "8R_test" else "8R"))
  if (TEST8 && identical(out, l8_abs(edge_battery("8R")))) stop("a test run (--reps) cannot write to the real block 8R folder")
  busy <- l8_abs(edge_battery(c(as.character(0:9), "1a", "1b", "8L", "8L_test", "9b", "9bL", "9c", "9d", "dryrun", "analysis")))
  if (any(startsWith(paste0(out, "/"), paste0(busy, "/")))) stop("--out cannot be a folder of another block: ", out)
  if (!startsWith(paste0(out, "/"), paste0(l8_abs(BAT), "/"))) stop("--out must lie inside ", BAT)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  OUT8 <<- out
  invisible(out)
}
r8_R <- function(ce) if (TEST8) min(REPS8, ce$R) else ce$R
r8_out_path <- function(ce) file.path(OUT8, sprintf("%s_%s_pvalues.csv.gz", ce$cell, ce$test))
r8_log <- function(msg) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "", msg)
  cat(line, "\n", sep = ""); cat(line, "\n", sep = "", file = file.path(OUT8, "_progress.log"), append = TRUE)
}
r8_check <- function(sel) {
  bad <- character(0)
  if ("lecessie" %in% sel$test) {
    if (utils::packageVersion("ebrahim.gof") < L8_MINVER) bad <- c(bad, "ebrahim.gof is older than 2.6.0")
    tc <- tryCatch(l8_transpose_check(), error = function(e) c(package = NA, declared = NA, vendor = NA))
    if (!isTRUE(abs(tc[["package"]] - tc[["declared"]]) < 1e-10) || !isTRUE(abs(tc[["package"]] - tc[["vendor"]]) > 1e-6))
      bad <- c(bad, "the installed le Cessie is not the corrected (I-H)'R(I-H)")
  }
  if ("bagoft" %in% sel$test) {
    if (!requireNamespace("BAGofT", quietly = TRUE) || utils::packageVersion("BAGofT") != "1.0.0")
      bad <- c(bad, "BAGofT 1.0.0 is required")
    else if (!isTRUE(tryCatch(r8_patch_check(), error = function(e) FALSE)))
      bad <- c(bad, "the parRF patch does not change exactly the one declared line")
  }
  for (i in seq_len(nrow(sel))) {
    ce <- as.list(sel[i, ]); f <- rv_bat_path(ce); R <- r8_R(ce)
    if (!file.exists(f)) { bad <- c(bad, sprintf("%s: battery file %s missing", ce$cell, f)); next }
    Sx <- data.table::fread(f, select = c("rep", RV_ID))
    if (!all(seq_len(R) %in% Sx$rep)) bad <- c(bad, sprintf("%s: the battery file lacks some of replicates 1-%d", ce$cell, R))
    o <- r8_out_path(ce); pt <- paste0(o, ".part")
    if (file.exists(o)) {
      nr <- nrow(data.table::fread(o, select = "rep"))
      if (nr != R) bad <- c(bad, sprintf("%s exists with %d rows but %d are planned: move it away", o, nr, R))
    } else if (file.exists(pt) || file.exists(paste0(pt, ".tmp"))) {
      Pp <- rv_read_part(pt)
      if (is.null(Pp) || !identical(names(Pp), r8_names(ce$test)) || anyDuplicated(Pp$rep) || any(!Pp$rep %in% seq_len(R)))
        bad <- c(bad, sprintf("%s does not belong to this run: move it away", pt))
    }
  }
  bad
}
r8_cluster <- function(w, need_bagoft) {
  rv_pin_blas()
  cl <- parallel::makePSOCKcluster(w)
  parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
  parallel::clusterCall(cl, function(simdir, bag) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
    options(block8L.source_only = TRUE)
    suppressPackageStartupMessages(library(data.table))
    suppressMessages(source(file.path(simdir, "run_M_block8L.R")))
    RNGkind("L'Ecuyer-CMRG")
    NULL
  }, SIMDIR, need_bagoft)
  parallel::clusterExport(cl, c("r8_one", "r8_task", "r8_names", "R8_PCOL", "R8_AUX", "r8_make_parRF_dropfix"))
  if (need_bagoft) parallel::clusterEvalQ(cl, { PARRF_FIX <- r8_make_parRF_dropfix(); NULL })
  cl
}
r8_run_cell <- function(cl, ce, prog) {
  R <- r8_R(ce); out <- r8_out_path(ce); part <- paste0(out, ".part")
  cell <- r8_cell(ce)
  S <- rv_stored(ce, R)
  M <- rv_read_part(part)
  todo <- setdiff(seq_len(R), M$rep)
  r8_log(sprintf("%-8s %-22s: %d of %d replicates to run%s", ce$test, ce$cell, length(todo), R,
                 if (is.null(M)) "" else sprintf(" (%d resumed from .part)", nrow(M))))
  bs <- length(cl) * R8_BATCH[[ce$test]]
  for (b in split(todo, ceiling(seq_along(todo) / bs))) {
    tasks <- lapply(b, function(r) list(rep = r, st = as.list(S[r, ])))
    rows <- parallel::clusterApplyLB(cl, tasks, r8_task, cell = cell, test = ce$test)
    M <- rbind(M, as.data.frame(do.call(rbind, rows)))
    M <- M[order(M$rep), , drop = FALSE]; rownames(M) <- NULL
    tmp <- paste0(part, ".tmp")
    data.table::fwrite(M, tmp, compress = "gzip")
    if (file.exists(part)) file.remove(part)
    if (!file.rename(tmp, part)) stop("could not rename ", tmp)
    prog$add(length(b) * r8_cost(ce$test, ce$n))
    r8_log(sprintf("  %-8s %-22s %d/%d written  ETA %s", ce$test, ce$cell, nrow(M), R, prog$eta()))
  }
  if (is.null(M) || nrow(M) != R || any(M$rep != seq_len(R))) stop("internal: ", part, " is incomplete")
  if (file.exists(out)) file.remove(out)
  if (!file.rename(part, out)) stop("could not rename ", part)
  M
}
## seconds per data set on one core (block 8's measurements), for the ETA only
r8_cost <- function(test, n) switch(test, lecessie = 0.2 * (n / 1000)^2 + 0.05,
                                    proj = if (n <= 500) 20 else 290, bagoft = if (n <= 500) 377 else 546)

## ---- main ---------------------------------------------------------------------------------------------------------------
r8_main <- function() {
  if (isTRUE(OPT8$summary)) { r8_summary(OUT8); return(0L) }
  W <- if (is.null(OPT8$workers)) 20L else as.integer(OPT8$workers)
  if (is.na(W) || W < 1L) stop("--workers needs a positive integer")
  sel <- r8_plan()
  if (!is.null(OPT8$tests)) {
    want <- trimws(strsplit(OPT8$tests, ",")[[1]])
    if (any(!want %in% R8_TESTS)) stop("--tests takes lecessie, proj, bagoft")
    sel <- sel[sel$test %in% want, ]
  }
  if (!is.null(OPT8$cells)) sel <- sel[sel$cell %in% trimws(strsplit(OPT8$cells, ",")[[1]]), ]
  sel <- sel[order(match(sel$test, R8_TESTS), sel$role != "null"), ]            # le Cessie first; nulls first within a test
  r8_log(sprintf("block 8R: %d scenario-test pairs, %d workers%s, out %s", nrow(sel), W,
                 if (TEST8) sprintf(", replicates 1-%d (TEST RUN)", REPS8) else "", OUT8))
  bad <- r8_check(sel)
  if (length(bad)) { for (b in bad) r8_log(paste("refused:", b)); r8_log("block 8R not started"); return(1L) }
  lock <- file.path(OUT8, "_lock")
  if (file.exists(lock) && !isTRUE(OPT8$force)) {
    pid <- suppressWarnings(as.integer(readLines(lock, n = 1L)))
    if (rv_pid_alive(pid)) stop("another block 8R run (process ", pid, ") holds ", lock)
    r8_log(sprintf("a stale lock from process %s is being replaced", format(pid)))
  }
  writeLines(c(as.character(Sys.getpid()), format(Sys.time())), lock)
  on.exit(if (file.exists(lock)) file.remove(lock), add = TRUE)
  done <- vapply(seq_len(nrow(sel)), function(i) file.exists(r8_out_path(as.list(sel[i, ]))), logical(1))
  todo <- which(!done)
  if (!length(todo)) { r8_log("block 8R: nothing to run"); r8_log("block 8R: run finished"); return(0L) }
  prog <- rv_progress(sum(vapply(todo, function(i) { ce <- as.list(sel[i, ]); r8_R(ce) * r8_cost(ce$test, ce$n) }, 1)))
  cl <- r8_cluster(W, any(sel$test[todo] == "bagoft"))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  kk <- 0L
  for (i in todo) {
    ce <- as.list(sel[i, ]); kk <- kk + 1L; t1 <- Sys.time()
    M <- r8_run_cell(cl, ce, prog)
    pc <- R8_PCOL[[ce$test]]
    r8_log(sprintf("[%d/%d] %-8s %-22s n=%-4d R=%-3d kept %d, identity failures %d, no p-value %d, rejection %.3f, median %.1f s, %.1f min  ETA %s",
                   kk, length(todo), ce$test, ce$cell, ce$n, nrow(M), sum(M$id_ok %in% 1), sum(M$id_ok %in% 0),
                   sum(M$id_ok %in% 1 & !is.finite(M[[pc]])), mean(M$id_ok %in% 1 & r8_rej(M[[pc]], R8_GRID[[ce$test]])),
                   stats::median(M$sec, na.rm = TRUE), as.numeric(difftime(Sys.time(), t1, units = "mins")), prog$eta()))
  }
  r8_log("block 8R: run finished")
  0L
}

if (!isTRUE(getOption("block8R.source_only"))) {
  r8_setup(r8_opts())
  quit(save = "no", status = r8_main())
}
