## run_M_rivals.R -- block 8 of the EDGE restructure battery (paper_EDGE/theory/PREDECLARATION_restructure_battery.md, E9 and
## E10.7): the projection test of Liu, Li, Chen, Haerdle and Liang (2024, Statistics and Computing 34:175) and BAGofT (Zhang,
## Ding and Yang 2023, JASA) on the battery's own data sets, regenerated from their stored seeds. Block 8 is not part of
## run_M_all.R: it runs after block 7 and changes nothing in blocks 0-7.
##
## Usage
##   Rscript run_M_rivals.R --test proj|bagoft|both [--cells a,b,...] [--workers 20] [--root battery] [--out battery/8]
##   Rscript run_M_rivals.R --summary [--root battery] [--out battery/8] [--forms poly3=u,sym=u]
##   Rscript run_M_rivals.R --list
##   Rscript run_M_rivals.R --test proj --cells loglog_base_n610 --reps 3 --workers 2 --root battery/dryrun --out battery/8_test/x
## --reps marks a test run: every test runs replicates 1 to --reps, the output goes to battery/8_test unless --out is given (never
## to battery/8), and --allow-missing-null (test runs only) lets an alternative start before its null's battery file exists.
##
## Cells (E9, 30): block 3 stk_long, stk_short, stk_asym, cauchit, t4 and loglog at n = 500 and 1000, with 1a null_link_n500
## and null_link_n1000; every block 2 cell with n <= 1000 and its 1b null; block 4 crossover_n1000 with 1b null_crossover_n1000.
## Generator, seed base, fitted formula and matched null of every cell come from battery_cells() (_battery_cells.R).
##
## One replicate (projection: replicates 1-500 of the cell; BAGofT: 1-200), as battery_one() of _battery_tests.R does it:
## L'Ecuyer-CMRG, set.seed(seed_base + rep), the cell's generator (bt_data), the working fit (bt_fit) and EDGE-poly3 unit at
## G = 10 (bt_edge_arm); nothing after the data draws a random number. The replicate is kept only if seed, n and events equal
## the cell's battery file and EDGE.poly3.u.G10 agrees with it to 1e-8 (or is missing in both); otherwise the rival is not run
## and the replicate is counted as an identity failure. The rival then continues the same random-number stream:
##   proj    proj_pvalue(y, X, B = 250), with y and X of the working fit, as grid_proj_power.R (_proj_test.R);
##   bagoft  BAGofT(testModel = testGlmBi(formula = <fitted formula>, link = "logit"), data = <data frame of the fit>) at the
##           package defaults (BAGofT 1.0.0: nsplits 100, ne floor(5 sqrt(n)), nsim 100, parFun parRF()).
## An error or no p-value counts as no rejection (E0.5); a sample with no event or no non-event has no working fit and gets no
## p-value (E8.1). Seconds per data set are the elapsed time of the rival's own call.
##
## Output (<out> = battery/8): <cell>_<test>_pvalues.csv.gz, one row per replicate. The finished rows go to <file>.part after
## every batch and the file is renamed when the cell is complete, so a run resumes inside a cell; a finished file with the wrong
## number of rows stops the run. A run refuses to start while any selected cell's battery file, or its matched null's, does not
## exist. _progress.log (with ETA). --summary writes _summary.csv (raw and size-adjusted rejection, declined rate, identity
## failures, seconds) and _paired.csv (exact McNemar at 0.05 against EDGE-poly3, EDGE-sym, Stk.joint and GiViTI, Holm within
## block 8). Both rivals' p-values lie on a Monte Carlo grid (BAGofT k/100, projection j/250), so a rival rejects at p < alpha; the
## comparators keep the battery's p <= alpha (the rule and its levels are written above rv_summary).
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
suppressPackageStartupMessages({ library(parallel); library(data.table) })
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))
source(file.path(SIMDIR, "_proj_test.R"))

BAT       <- edge_battery()
RV_TESTS  <- c("proj", "bagoft")
RV_REPS   <- c(proj = 500L, bagoft = 200L)                          # E9: replicates 1-500 and 1-200 of each cell
RV_PCOL   <- c(proj = "proj", bagoft = "BAGofT")                     # the p-value column of each test
RV_AUX    <- list(proj = "aux.proj_stat", bagoft = c("aux.bagoft_pmean", "aux.bagoft_p_median", "aux.bagoft_p_min"))
RV_BPROJ  <- 250L                                                    # bootstrap replicates of the projection test
RV_TOL    <- 1e-8                                                    # identity with the stored EDGE.poly3.u.G10
RV_ID     <- "EDGE.poly3.u.G10"
RV_ALPHAS <- c(0.01, 0.05, 0.10)
RV_BATCH  <- c(proj = 8L, bagoft = 2L)                               # replicates per batch, per worker (one .part write per batch)
RV_E9_B3  <- c("stk_long", "stk_short", "stk_asym", "cauchit", "t4", "loglog")
RV_COMP   <- c(c(t(outer(c("EDGE.poly3", "EDGE.sym"), c("u.G10", "u.Grule", "sc.G10", "sc.Grule"), paste, sep = "."))),
               "Stk.joint", "GiViTI")                                # the paired comparators of E9
## serial seconds per data set, for the ETA only (bench_slow_timing.csv: BAGofT 454 s at n = 1000; proj_timing.txt: 41.7 s at
## n = 1000 with B = 1000, most of it the n x n kernel)
rv_cost_guess <- function(test, n) if (test == "proj") 35 * (n / 1000)^3 + 2 * n / 1000 else 460 * n / 1000

## ---- options ----------------------------------------------------------------------------------------------------------
rv_parse <- function(a) {
  opt <- list(test = NULL, cells = NULL, reps = NULL, workers = NULL, root = NULL, out = NULL, forms = "poly3=u,sym=u",
              summary = FALSE, list = FALSE, "allow-missing-null" = FALSE, force = FALSE)
  flags <- c("summary", "list", "allow-missing-null", "force")
  i <- 1L
  while (i <= length(a)) {
    key <- sub("^--", "", a[i]); val <- NULL
    if (grepl("=", key)) { val <- sub("^[^=]*=", "", key); key <- sub("=.*$", "", key) }
    if (!key %in% names(opt)) stop("unknown option --", key)
    if (key %in% flags) { opt[[key]] <- TRUE; i <- i + 1L; next }
    if (is.null(val)) { i <- i + 1L; val <- a[i] }
    opt[[key]] <- val; i <- i + 1L
  }
  opt
}
rv_inside <- function(p, what) {
  if (!grepl("^([A-Za-z]:|/)", p)) p <- file.path(SIMDIR, p)
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  b <- normalizePath(BAT, winslash = "/", mustWork = FALSE)
  if (!startsWith(paste0(p, "/"), paste0(b, "/"))) stop("--", what, " must lie inside ", b)
  p
}
## sets OPT, TEST_MODE, ROOT, OUT and W; creates no directory
rv_setup <- function(a = commandArgs(trailingOnly = TRUE)) {
  opt <- rv_parse(a)
  tm <- !is.null(opt$reps)
  if (tm && (is.na(suppressWarnings(as.integer(opt$reps))) || as.integer(opt$reps) < 1L)) stop("--reps needs a positive integer")
  if (isTRUE(opt[["allow-missing-null"]]) && !tm) stop("--allow-missing-null is for test runs (--reps) only")
  root <- rv_inside(if (is.null(opt$root)) BAT else opt$root, "root")
  out <- rv_inside(if (!is.null(opt$out)) opt$out else edge_battery(if (tm) "8_test" else "8"), "out")
  real8 <- normalizePath(edge_battery("8"), winslash = "/", mustWork = FALSE)
  if (tm && startsWith(paste0(out, "/"), paste0(real8, "/"))) stop("a test run (--reps) cannot write to ", real8)
  busy <- normalizePath(edge_battery(c("0", "1a", "1b", "2", "3", "4", "5", "6", "7", "dryrun")), winslash = "/", mustWork = FALSE)
  if (any(startsWith(paste0(out, "/"), paste0(busy, "/")))) stop("--out cannot be a folder of blocks 0-7 or the dry run: ", out)
  w <- if (!is.null(opt$workers)) as.integer(opt$workers) else {
    f <- edge_battery("0", "workers_choice.txt")                   # the worker count block 0 chose (E0.8)
    if (file.exists(f)) as.integer(readLines(f, warn = FALSE)[1]) else max(1L, min(20L, detectCores() - 2L))
  }
  if (is.na(w) || w < 1L) stop("--workers needs a positive integer")
  OPT <<- opt; TEST_MODE <<- tm; ROOT <<- root; OUT <<- out; W <<- w
  invisible(opt)
}
rv_forms <- function(s) {
  kv <- strsplit(strsplit(s, ",")[[1]], "=")
  f <- setNames(vapply(kv, function(z) z[2], ""), vapply(kv, function(z) z[1], ""))
  if (!setequal(names(f), c("poly3", "sym")) || any(!f %in% c("u", "sc"))) stop("--forms needs poly3=u|sc,sym=u|sc")
  as.list(f)
}

## ---- the cells of E9 -------------------------------------------------------------------------------------------------------
## E9's groups in its order, each matched null before its alternatives; every row is the battery_cells() row itself
rv_cells <- function(C = battery_cells()) {
  key <- paste(C$block, C$cell)
  pick <- function(k) {
    i <- match(k, key)
    if (anyNA(i)) stop("battery_cells() has no cell ", paste(k[is.na(i)], collapse = ", "))
    C[i, ]
  }
  groups <- list(pick(paste("3", sprintf("%s_n%d", rep(RV_E9_B3, each = 2), c(500L, 1000L)))),
                 C[C$block == "2" & C$n %in% 1:1000, ],
                 pick("4 crossover_n1000"))
  out <- list()
  for (A in groups) {
    if (any(A$role != "alternative") || anyNA(A$null_cell)) stop("block 8: every listed cell must be an alternative with a matched null")
    out[[length(out) + 1]] <- pick(unique(paste(A$null_block, A$null_cell)))
    out[[length(out) + 1]] <- A
  }
  R <- do.call(rbind, out)
  rownames(R) <- NULL
  nul <- R$cell %in% R$null_cell
  if (nrow(R) != 30L || anyDuplicated(R$cell) || any(R$n > 1000L) || anyNA(R$seed_base) || any(R$role[nul] != "null") ||
      any(!is.na(R$null_cell[nul])))
    stop("the block 8 cell list does not match E9 (30 cells with n <= 1000, nulls and alternatives)")
  R
}
rv_bat_path <- function(ce) file.path(ROOT, ce$block, paste0(ce$cell, "_pvalues.csv.gz"))
rv_out_path <- function(ce, test) file.path(OUT, sprintf("%s_%s_pvalues.csv.gz", ce$cell, test))
rv_null_of  <- function(ce) if (is.na(ce$null_cell)) NULL else as.list(Cells8[Cells8$block == ce$null_block & Cells8$cell == ce$null_cell, ])

## the stored rows of replicates 1..R: rep, seed, n, events and p (= EDGE.poly3.u.G10)
rv_stored <- function(ce, R) {
  S <- as.data.frame(fread(rv_bat_path(ce), select = c("rep", "seed", "n", "events", RV_ID)))
  S <- S[S$rep %in% seq_len(R), , drop = FALSE]
  S <- S[order(S$rep), , drop = FALSE]
  if (nrow(S) != R || any(S$rep != seq_len(R))) stop(sprintf("%s: the battery file lacks some of replicates 1-%d", ce$cell, R))
  names(S)[names(S) == RV_ID] <- "p"
  rownames(S) <- NULL
  S
}

## ---- one replicate ------------------------------------------------------------------------------------------------------------
rv_num1 <- function(z) if (length(z) == 1L && is.numeric(z) && is.finite(z)) as.numeric(z) else NA_real_
rv_names <- function(test) c("rep", "seed", "n", "events", "id.stored", "id.regen", "id.absdiff", "id_ok", "flag.degenerate",
                              RV_PCOL[[test]], "sec", "flag.error", RV_AUX[[test]])

## the data set of one replicate and its EDGE-poly3 unit p-value at G = 10, by the path of battery_one() and battery_rep()
rv_regen <- function(rep, cell) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  set.seed(cell$seed_base + rep)
  dat <- bt_data(cell)
  n <- nrow(dat$d); ev <- sum(dat$d$y)
  deg <- min(ev, n - ev) == 0                                        # no working fit (E8.1)
  fq <- if (deg) NULL else bt_fit(dat)
  p <- if (is.null(fq)) NA_real_ else tryCatch(bt_edge_arm(fq, 10L, cell)[["EDGE.poly3.u"]], error = function(e) NA_real_)
  list(dat = dat, fq = fq, n = n, events = ev, degenerate = deg, p = p)
}
rv_identity <- function(g, st, seed) {
  sp <- suppressWarnings(as.numeric(st$p)); rp <- g$p
  same <- if (is.finite(sp) && is.finite(rp)) abs(sp - rp) <= RV_TOL else !is.finite(sp) && !is.finite(rp)
  isTRUE(same) && isTRUE(st$seed == seed) && isTRUE(st$n == g$n) && isTRUE(st$events == g$events)
}

rv_one <- function(rep, cell, test, st) {
  nm <- rv_names(test)
  out <- setNames(rep(NA_real_, length(nm)), nm)
  seed <- cell$seed_base + rep
  out["rep"] <- rep; out["seed"] <- seed; out["id_ok"] <- 0; out["id.stored"] <- suppressWarnings(as.numeric(st$p))
  g <- tryCatch(rv_regen(rep, cell), error = function(e) NULL)
  if (is.null(g)) return(out)                                        # the generator failed: the stored data cannot be reproduced
  out["n"] <- g$n; out["events"] <- g$events; out["id.regen"] <- g$p; out["flag.degenerate"] <- as.numeric(g$degenerate)
  if (is.finite(out[["id.stored"]]) && is.finite(g$p)) out["id.absdiff"] <- abs(out[["id.stored"]] - g$p)
  if (!rv_identity(g, st, seed)) return(out)                         # identity failure: the rival is not run
  out["id_ok"] <- 1
  if (is.null(g$fq)) { out["flag.error"] <- as.numeric(!g$degenerate); return(out) }
  ## the rival continues the replicate's stream: no random number has been drawn since the data
  t0 <- proc.time()[["elapsed"]]
  if (test == "proj") {
    r <- tryCatch(proj_pvalue(g$fq$fit$y, stats::model.matrix(g$fq$fit), B = RV_BPROJ), error = function(e) NULL)
    p <- rv_num1(r$p_value); aux <- rv_num1(r$stat)
  } else {
    r <- tryCatch(suppressWarnings(suppressMessages(
      BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = g$dat$f, link = "logit"), data = g$dat$d))), error = function(e) NULL)
    p <- rv_num1(r$p.value); aux <- c(rv_num1(r$pmean), rv_num1(r$p.value2), rv_num1(r$p.value3))
  }
  out["sec"] <- proc.time()[["elapsed"]] - t0
  out[RV_PCOL[[test]]] <- if (is.finite(p) && p >= 0 && p <= 1) p else NA_real_
  out["flag.error"] <- as.numeric(!is.finite(out[[RV_PCOL[[test]]]]))
  out[RV_AUX[[test]]] <- aux
  out
}
rv_task <- function(task, cell, test) rv_one(task$rep, cell, test, task$st)

## ---- workers, log, checks ------------------------------------------------------------------------------------------------------
rv_pin_blas <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
}
rv_cluster <- function(w, bagoft) {
  rv_pin_blas()
  cl <- tryCatch(makeCluster(w), error = function(e) {               # one retry on another port, as bt_cluster()
    message("makeCluster failed (", conditionMessage(e), "); retrying once on another port")
    Sys.sleep(3); makeCluster(w, port = 11000L + (Sys.getpid() %% 997L))
  })
  parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
  clusterCall(cl, function(simdir, bag) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
    source(file.path(simdir, "_battery_tests.R")); source(file.path(simdir, "_proj_test.R"))
    if (bag && !suppressWarnings(requireNamespace("BAGofT", quietly = TRUE))) stop("BAGofT is not installed")
    RNGkind("L'Ecuyer-CMRG")
    NULL
  }, SIMDIR, bagoft)
  clusterExport(cl, c("rv_names", "rv_num1", "rv_regen", "rv_identity", "rv_one", "rv_task", "RV_PCOL", "RV_AUX", "RV_BPROJ", "RV_TOL"))
  cl
}
rv_log <- function(msg) {
  line <- sprintf("%s  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", sep = "")
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  cat(line, "\n", file = file.path(OUT, "_progress.log"), append = TRUE, sep = "")
}
rv_progress <- function(total) {
  t0 <- Sys.time(); done <- 0
  list(add = function(w) done <<- done + w,
       eta = function() {
         if (done <= 0) return("unknown")
         el <- as.numeric(difftime(Sys.time(), t0, units = "secs")); rem <- el / done * max(total - done, 0)
         sprintf("%s (%.2f h)", format(Sys.time() + rem, "%Y-%m-%d %H:%M"), rem / 3600)
       })
}

## a .part file, or the complete .part.tmp left if a run stopped between writing and renaming; NULL if neither reads
rv_read_part <- function(part) {
  for (f in c(part, paste0(part, ".tmp"))) {
    if (!file.exists(f)) next
    P <- tryCatch(as.data.frame(fread(f)), error = function(e) NULL)
    if (!is.null(P)) return(P)
  }
  NULL
}

## is a process still running (Windows: tasklist; elsewhere: signal 0)
rv_pid_alive <- function(pid) {
  if (length(pid) != 1L || is.na(pid)) return(FALSE)
  if (.Platform$OS.type == "windows") {
    tl <- tryCatch(system2("tasklist", c("/FI", shQuote(paste("PID eq", pid), type = "cmd"), "/NH"), stdout = TRUE, stderr = FALSE),
                   error = function(e) character(0))
    return(any(grepl(paste0(" ", pid, " "), paste0(" ", tl, " "), fixed = TRUE)))
  }
  isTRUE(tools::pskill(pid, 0L))
}

## every reason not to start, found before any computation
rv_check <- function(sel, tests, reps_of, allow_missing_null) {
  bad <- character(0)
  need <- max(vapply(tests, reps_of, integer(1)))
  for (i in seq_len(nrow(sel))) {
    ce <- as.list(sel[i, ])
    f <- rv_bat_path(ce)
    if (!file.exists(f)) {
      bad <- c(bad, sprintf("%s: its battery file %s does not exist yet", ce$cell, f))
    } else {
      S <- fread(f, select = c("rep", RV_ID))
      if (!all(seq_len(need) %in% S$rep))
        bad <- c(bad, sprintf("%s: the battery file has %d replicates; replicates 1-%d are needed", ce$cell, nrow(S), need))
    }
    if (!is.na(ce$null_cell)) {
      nf <- file.path(ROOT, ce$null_block, paste0(ce$null_cell, "_pvalues.csv.gz"))
      if (!file.exists(nf)) {
        msg <- sprintf("%s: the battery file of its matched null %s does not exist yet (%s)", ce$cell, ce$null_cell, nf)
        if (allow_missing_null) rv_log(paste("test run, allowed:", msg)) else bad <- c(bad, msg)
      }
    }
    g <- tryCatch(rv_regen(1L, ce), error = function(e) NULL)       # the fitted formula and the data frame the rivals receive
    if (is.null(g)) bad <- c(bad, sprintf("%s: the generator fails", ce$cell))
    else if (!identical(paste(deparse(g$dat$f), collapse = ""), ce$formula) || !setequal(names(g$dat$d), all.vars(g$dat$f)) ||
             nrow(g$dat$d) != ce$n)
      bad <- c(bad, sprintf("%s: generator formula or data frame differs from the cell table (%s)", ce$cell, ce$formula))
    for (t in tests) {
      o <- rv_out_path(ce, t); R <- reps_of(t); pt <- paste0(o, ".part")
      if (file.exists(o)) {
        nr <- nrow(fread(o, select = "rep"))
        if (nr != R) bad <- c(bad, sprintf("%s exists with %d rows but %d replicates are planned: move the file away and start again", o, nr, R))
      } else if (file.exists(pt) || file.exists(paste0(pt, ".tmp"))) {
        P <- rv_read_part(pt)
        if (is.null(P))
          bad <- c(bad, sprintf("%s (and its .tmp) cannot be read, probably an interrupted write: move them away", pt))
        else if (!identical(names(P), rv_names(t)) || anyDuplicated(P$rep) || any(!P$rep %in% seq_len(R)))
          bad <- c(bad, sprintf("%s does not belong to this run (columns, or replicates outside 1-%d): move it away", pt, R))
      }
    }
  }
  bad
}

## ---- one cell -------------------------------------------------------------------------------------------------------------------
rv_run_cell <- function(cl, ce, test, R, prog) {
  out <- rv_out_path(ce, test); part <- paste0(out, ".part")
  S <- rv_stored(ce, R)
  M <- rv_read_part(part)
  todo <- setdiff(seq_len(R), M$rep)
  rv_log(sprintf("%s %s: %d of %d replicates to run%s", test, ce$cell, length(todo), R,
                 if (is.null(M)) "" else sprintf(" (%d rows resumed from .part)", nrow(M))))
  bs <- length(cl) * RV_BATCH[[test]]
  for (b in split(todo, ceiling(seq_along(todo) / bs))) {
    tasks <- lapply(b, function(r) list(rep = r, st = as.list(S[r, ])))
    rows <- clusterApplyLB(cl, tasks, rv_task, cell = ce, test = test)
    M <- rbind(M, as.data.frame(do.call(rbind, rows)))
    M <- M[order(M$rep), , drop = FALSE]; rownames(M) <- NULL
    tmp <- paste0(part, ".tmp")                                      # write whole, then swap in, so a crash never truncates .part
    fwrite(M, tmp, compress = "gzip")
    if (file.exists(part)) file.remove(part)
    if (!file.rename(tmp, part)) stop("could not rename ", tmp)
    prog$add(length(b) * rv_cost_guess(test, ce$n))
    rv_log(sprintf("  %s %-22s %d/%d replicates written to .part  ETA %s", test, ce$cell, nrow(M), R, prog$eta()))
  }
  if (is.null(M) || nrow(M) != R || any(M$rep != seq_len(R))) stop("internal: ", part, " is incomplete")
  if (file.exists(out)) file.remove(out)
  if (!file.rename(part, out)) stop("could not rename ", part)
  if (file.exists(paste0(part, ".tmp"))) file.remove(paste0(part, ".tmp"))
  M
}

rv_run <- function() {
  if (is.null(OPT$test) || !OPT$test %in% c(RV_TESTS, "both")) stop("--test must be proj, bagoft or both")
  tests <- if (OPT$test == "both") RV_TESTS else OPT$test
  reps_of <- function(t) if (TEST_MODE) as.integer(OPT$reps) else RV_REPS[[t]]
  if (!TEST_MODE && !isTRUE(OPT$force)) {                            # E11.4: block 8 starts after blocks 0-7
    lg <- edge_battery("launch.log")
    L <- if (file.exists(lg)) readLines(lg, warn = FALSE) else character(0)
    last_start <- max(c(0L, grep("launch from step", L, fixed = TRUE)))
    if (!any(grepl("launch finished", L[seq_along(L) > last_start], fixed = TRUE)))
      stop("blocks 0-7 have not finished (", lg, "); block 8 starts after them (E11.4), or pass --force")
  }
  sel <- Cells8
  if (!is.null(OPT$cells)) {
    want <- strsplit(OPT$cells, ",")[[1]]
    if (any(!want %in% sel$cell)) stop("not a block 8 cell: ", paste(setdiff(want, sel$cell), collapse = ", "))
    sel <- sel[sel$cell %in% want, ]
  }
  rv_log(sprintf("block 8, %s: %d cells, %s, %d workers, root %s, out %s%s", paste(tests, collapse = " then "), nrow(sel),
                 paste(sprintf("%s replicates 1-%d", tests, vapply(tests, reps_of, integer(1))), collapse = ", "), W, ROOT, OUT,
                 if (TEST_MODE) " [test run]" else ""))
  bad <- rv_check(sel, tests, reps_of, isTRUE(OPT[["allow-missing-null"]]))
  if (length(bad)) {
    for (b in bad) rv_log(paste("refused:", b))
    rv_log("block 8 not started")
    quit(save = "no", status = 1L)
  }
  if (!TEST_MODE) {                                                  # one real block 8 run at a time, taken only once nothing refused
    lock <- file.path(OUT, "_lock")
    if (file.exists(lock) && !isTRUE(OPT$force)) {
      pid <- suppressWarnings(as.integer(readLines(lock, warn = FALSE)[1]))
      if (rv_pid_alive(pid)) stop("another block 8 run (process ", pid, ") holds ", lock)
      rv_log(sprintf("a lock left by process %s, which is no longer running, is replaced", pid))
    }
    writeLines(c(as.character(Sys.getpid()), format(Sys.time())), lock)
    on.exit(unlink(lock), add = TRUE)                                # an Rscript error skips this; the next run then finds the process gone
  }
  plan <- do.call(rbind, lapply(tests, function(t) data.frame(test = t, i = seq_len(nrow(sel)), stringsAsFactors = FALSE)))
  plan$R <- vapply(plan$test, reps_of, integer(1))
  plan$out <- vapply(seq_len(nrow(plan)), function(k) rv_out_path(as.list(sel[plan$i[k], ]), plan$test[k]), "")
  plan$done <- file.exists(plan$out)
  plan$cost <- 0
  for (k in seq_len(nrow(plan))) {
    ce <- as.list(sel[plan$i[k], ]); part <- paste0(plan$out[k], ".part")
    if (plan$done[k]) { rv_log(sprintf("%s %s: finished file present (%d rows), skipped", plan$test[k], ce$cell, plan$R[k])); next }
    nres <- if (file.exists(part)) nrow(fread(part, select = "rep")) else 0L
    if (nres == plan$R[k]) {                                         # a complete .part (checked by rv_check): rename, nothing to run
      file.rename(part, plan$out[k]); plan$done[k] <- TRUE
      rv_log(sprintf("%s %s: complete .part found (%d rows), renamed", plan$test[k], ce$cell, nres)); next
    }
    plan$cost[k] <- (plan$R[k] - nres) * rv_cost_guess(plan$test[k], ce$n)
  }
  if (all(plan$done)) { rv_log("block 8: nothing to run"); rv_log("block 8: run finished"); return(invisible(NULL)) }
  prog <- rv_progress(sum(plan$cost))
  cl <- rv_cluster(W, any(plan$test[!plan$done] == "bagoft"))
  on.exit(stopCluster(cl))
  kk <- 0L; K <- sum(!plan$done)
  for (k in which(!plan$done)) {
    kk <- kk + 1L; ce <- as.list(sel[plan$i[k], ]); t1 <- Sys.time()
    M <- rv_run_cell(cl, ce, plan$test[k], plan$R[k], prog)
    ran <- M$id_ok %in% 1 & is.finite(M$sec)
    rv_log(sprintf("[%d/%d] %-6s %-22s n=%-5d R=%-3d kept %d, identity failures %d, errors %d, degenerate %d, median %.1f s per data set, %.1f min  ETA %s",
                   kk, K, plan$test[k], ce$cell, ce$n, plan$R[k], sum(M$id_ok %in% 1), sum(M$id_ok %in% 0), sum(M$flag.error %in% 1),
                   sum(M$flag.degenerate %in% 1), if (any(ran)) stats::median(M$sec[ran]) else NA_real_,
                   as.numeric(difftime(Sys.time(), t1, units = "mins")), prog$eta()))
  }
  rv_log("block 8: run finished")
  invisible(NULL)
}

## ---- summary ------------------------------------------------------------------------------------------------------------------
## Decision rule (review F1). Both rivals' p-values are Monte Carlo proportions on a grid: BAGofT's is k/100 (k = simulated
## statistics below the observed one, nsim = 100), the projection test's j/250 (j = bootstrap statistics at or above the observed
## one). A rival rejects at p < alpha, so a p-value equal to alpha does not reject. With k uniform on 0..100 under the null, BAGofT's
## level is then 1/101, 5/101 and 10/101 at 0.01, 0.05 and 0.10 (the decisions of (k + 1)/(nsim + 1) <= alpha), where p <= alpha
## would give 2/101, 6/101 and 11/101; the projection test changes only at 0.10, its one grid point (25/251 instead of 26/251).
## The rule is used for raw rejection, rejection among samples with a p-value, null size, the size gate and the McNemar
## decisions, as the earlier benchmarks did (_bench_bagoft.R, grid_proj_power.R: p < 0.05). The comparators' p-values from the
## battery files are continuous and keep block_summary's p <= alpha. Size-adjusted power keeps block_summary's p <= critical value
## for every test; size_adj_null_level is the matched null's own rejection rate at that critical value (ties on a grid can put it
## above alpha).
RV_RULE <- c(rival = "p < alpha", comparator = "p <= alpha")
RV_EPS  <- 1e-9                                                      # far below either grid step (1/250); absorbs the representation of k/100
rv_rej  <- function(p, a, mc) if (mc) is.finite(p) & p < a - RV_EPS else is.finite(p) & p <= a
## size at 0.05 within three Monte Carlo standard errors of 0.05, or below (mc: a rival's p-value, decided at p < alpha)
rv_holds <- function(p, mc, a = 0.05) {
  if (!length(p)) return(NA)
  mean(rv_rej(p, a, mc)) <= a + 3 * sqrt(a * (1 - a) / length(p))
}
rv_summary <- function(forms) {
  cache <- new.env()
  rd <- function(f) {
    if (!exists(f, envir = cache, inherits = FALSE)) assign(f, if (file.exists(f)) as.data.frame(fread(f)) else NULL, envir = cache)
    get(f, envir = cache)
  }
  S <- list(); PR <- list()
  for (test in RV_TESTS) for (i in seq_len(nrow(Cells8))) {
    ce <- as.list(Cells8[i, ]); pc <- RV_PCOL[[test]]
    P <- rd(rv_out_path(ce, test))
    if (is.null(P)) next
    K <- P[P$id_ok %in% 1, , drop = FALSE]
    alt <- !is.na(ce$null_cell); nce <- rv_null_of(ce)
    NK <- if (alt) { z <- rd(rv_out_path(nce, test)); if (is.null(z)) NULL else z[z$id_ok %in% 1, , drop = FALSE] } else NULL
    B  <- rd(rv_bat_path(ce)); BN <- if (alt) rd(rv_bat_path(nce)) else NULL
    ran <- K$flag.degenerate %in% 0 & is.finite(K$sec)
    kd <- P$id.absdiff[P$id_ok %in% 1]
    base <- data.frame(block = ce$block, cell = ce$cell, role = ce$role, null_cell = if (alt) ce$null_cell else NA_character_,
                       n = ce$n, rival = test, stringsAsFactors = FALSE)
    info <- data.frame(reps_in_file = nrow(P), reps_kept = nrow(K), identity_failures = sum(P$id_ok %in% 0),
                       identity_both_missing = sum(P$id_ok %in% 1 & !is.finite(P$id.stored)),
                       max_identity_absdiff = if (any(is.finite(kd))) max(kd[is.finite(kd)]) else NA_real_,
                       null_reps_kept = if (is.null(NK)) NA_integer_ else nrow(NK), rival_errors = sum(K$flag.error %in% 1),
                       rival_degenerate = sum(K$flag.degenerate %in% 1),
                       rival_median_sec = if (any(ran)) stats::median(K$sec[ran]) else NA_real_, stringsAsFactors = FALSE)
    ## the rival on its kept replicates; each comparator from the battery file on the same replicates, size-adjusted on the same
    ## replicates of the matched null
    for (t in c(pc, RV_COMP)) {
      mc <- t == pc                                                  # the rival's Monte Carlo p-value: p < alpha (review F1)
      if (mc) {
        p <- K[[pc]]; pn <- if (is.null(NK)) NULL else NK[[pc]]; src <- "block 8"
      } else {
        if (is.null(B) || !t %in% names(B)) next
        p <- B[[t]][match(K$rep, B$rep)]
        pn <- if (is.null(NK) || is.null(BN) || !t %in% names(BN)) NULL else BN[[t]][match(NK$rep, BN$rep)]
        src <- "battery file, the same replicates"
      }
      status <- if (!nrow(K)) "no replicate kept" else if (alt && is.null(pn)) "matched null not run yet" else ""
      for (a in RV_ALPHAS) {
        rate <- if (nrow(K)) mean(rv_rej(p, a, mc)) else NA_real_
        crit <- if (length(pn)) as.numeric(stats::quantile(ifelse(is.finite(pn), pn, 1), a, type = 1)) else NA_real_
        S[[length(S) + 1]] <- cbind(base, data.frame(test = t, source = src, alpha = a,
          decision_rule = RV_RULE[[if (mc) "rival" else "comparator"]], rejection = rate,
          mcse = if (nrow(K)) sqrt(rate * (1 - rate) / nrow(K)) else NA_real_,
          size_adj_power = if (is.finite(crit) && nrow(K)) mean(is.finite(p) & p <= crit) else NA_real_, null_crit = crit,
          size_adj_null_level = if (is.finite(crit)) mean(is.finite(pn) & pn <= crit) else NA_real_,
          null_size = if (length(pn)) mean(rv_rej(pn, a, mc)) else NA_real_,
          declined = if (nrow(K)) mean(!is.finite(p)) else NA_real_,
          rejection_given_p = if (any(is.finite(p))) sum(rv_rej(p, a, mc)) / sum(is.finite(p)) else NA_real_,
          status = status, stringsAsFactors = FALSE), info)
      }
    }
    ## paired decisions at 0.05 on the shared replicates (E9, C3). Size of the rival: its block 8 replicates of the matched null
    ## (of the cell itself for a null cell); size of a comparator: every replicate of the battery's null cell.
    if (!nrow(K) || is.null(B)) next
    size_r <- if (alt) (if (is.null(NK)) NULL else NK[[pc]]) else K[[pc]]
    hold_r <- rv_holds(size_r, TRUE)
    bsz <- if (alt) BN else B
    for (t in RV_COMP) {
      if (!t %in% names(B)) next
      m <- match(K$rep, B$rep); sh <- !is.na(m)
      pr <- K[[pc]][sh]; pcmp <- B[[t]][m[sh]]
      r1 <- rv_rej(pr, 0.05, TRUE); r2 <- rv_rej(pcmp, 0.05, FALSE)
      nb <- sum(r1 & !r2); nc <- sum(!r1 & r2)
      size_c <- if (is.null(bsz) || !t %in% names(bsz)) NULL else bsz[[t]]
      hold_c <- rv_holds(size_c, FALSE)
      form_ok <- if (grepl("^EDGE\\.", t)) {
        bs <- sub("^EDGE\\.([a-z0-9]+)\\..*$", "\\1", t); startsWith(t, sprintf("EDGE.%s.%s.", bs, forms[[bs]]))
      } else TRUE
      inh <- alt && isTRUE(hold_r) && isTRUE(hold_c) && form_ok
      note <- c(if (!alt) "null cell: size comparison, not in the Holm family",
                if (alt && is.na(hold_r)) "the rival's matched null has not run yet",
                if (isFALSE(hold_r)) "the rival does not hold size at 0.05", if (isFALSE(hold_c)) "the comparator does not hold size at 0.05",
                if (!form_ok) "not the weighting form chosen by Section A (--forms)")
      PR[[length(PR) + 1]] <- data.frame(base, comparator = t, shared_reps = sum(sh), rival_rejection = mean(r1),
        comparator_rejection = mean(r2), difference = mean(r1) - mean(r2), both = sum(r1 & r2), rival_only = nb, comparator_only = nc,
        neither = sum(!r1 & !r2), mcnemar_p = if (nb + nc == 0) 1 else stats::binom.test(nb, nb + nc, 0.5)$p.value,
        rival_size = if (length(size_r)) mean(rv_rej(size_r, 0.05, TRUE)) else NA_real_, rival_size_reps = length(size_r),
        rival_holds_size = hold_r, comparator_size = if (length(size_c)) mean(rv_rej(size_c, 0.05, FALSE)) else NA_real_,
        comparator_size_reps = length(size_c), comparator_holds_size = hold_c,
        rival_rule = RV_RULE[["rival"]], comparator_rule = RV_RULE[["comparator"]], chosen_form = form_ok,
        forms = sprintf("poly3=%s,sym=%s", forms$poly3, forms$sym), in_holm = inh, note = paste(note, collapse = "; "),
        stringsAsFactors = FALSE)
    }
  }
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  if (length(S)) { S <- rbindlist(S); fwrite(S, file.path(OUT, "_summary.csv")) }
  if (length(PR)) {
    PR <- rbindlist(PR)
    PR$holm_p <- NA_real_
    h <- which(PR$in_holm)
    if (length(h)) PR$holm_p[h] <- stats::p.adjust(PR$mcnemar_p[h], method = "holm")
    PR$holm_reject_05 <- PR$holm_p <= 0.05
    fwrite(PR, file.path(OUT, "_paired.csv"))
  }
  cat(sprintf("summary: %s (%d rows); paired: %s (%d rows, %d in the Holm family)\n", file.path(OUT, "_summary.csv"),
              if (is.data.frame(S)) nrow(S) else 0L, file.path(OUT, "_paired.csv"), if (is.data.frame(PR)) nrow(PR) else 0L,
              if (is.data.frame(PR)) sum(PR$in_holm) else 0L))
  invisible(list(summary = S, paired = PR))
}

## ---- main ---------------------------------------------------------------------------------------------------------------------
Cells8 <- rv_cells()
if (!isTRUE(getOption("rivals.source_only"))) {
  rv_setup()
  if (isTRUE(OPT$list)) {
    print(Cells8[, c("block", "cell", "role", "n", "null_block", "null_cell", "generator", "formula", "seed_base", "B")], row.names = FALSE)
    quit(save = "no", status = 0L)
  }
  if (isTRUE(OPT$summary)) { rv_summary(rv_forms(OPT$forms)); quit(save = "no", status = 0L) }
  rv_run()
  quit(save = "no", status = 0L)
}
