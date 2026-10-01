## run_M_block8L.R -- block 8L of the EDGE restructure: the le Cessie-van Houwelingen test on data sets the
## battery has already seen. Contract: paper_EDGE/theory/PREDECLARATION_block8L_lecessie.md (sha256 0c1a1731...),
## frozen 2026-09-19 before any 8L cell was run.
##
##   Rscript run_M_block8L.R [--part A|B|both] [--cells a,b,...] [--workers 20] [--reps N] [--out battery/8L] [--force]
##   Rscript run_M_block8L.R --summary [--out battery/8L]
##
## Part A: the 30 cells of block 8 (Cells8 of run_M_rivals.R), replicates 1-500, the replicates the projection test
## used. Part B: block 9's logistic cells at n = 1000, clean and C1 at k = 1, 2, 5, 10, replicates 1-1000.
##
## Nothing is simulated afresh. Each data set is regenerated from the seed its block stored -- Part A along block 8's
## own path (rv_regen), Part B along block 9's (b9_gen, then the battery's working fit) -- and it is used only if
## EDGE-poly3 unit at G = 10 recomputed on it equals the stored value to 1e-8 and seed, n and events agree: block 8's
## identity gate (rv_identity). A data set that fails the gate is not tested, and the cell reports it.
##
## le Cessie is the package's own, through its public battery: run.all.gof(fit, tests = "le-Cessie") of
## ebrahim.gof >= 2.6.0, whose moment reference is the corrected (I - H)' R (I - H). Before anything runs, the installed
## function is checked against that formula computed densely from its definition, and against the untransposed vendor
## form it must NOT match. The test draws no random number. Seconds per data set are the elapsed time of that call.
##
## Output (<out> = battery/8L): <cell>_lecessie_pvalues.csv.gz, one row per replicate, built in <file>.part and renamed
## when the cell is complete, so a run resumes inside a cell. --summary writes _summary.csv and _paired.csv (rules of
## section 3 of the declaration, written above l8_summary).
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
options(rivals.source_only = TRUE)
source(file.path(SIMDIR, "run_M_rivals.R"))                         # Cells8, rv_regen, rv_identity, rv_stored, rv_bat_path
source(file.path(SIMDIR, "_block9_contam.R"))                       # b9_cells, b9_gen
BAT  <- edge_battery()
ROOT <- BAT                                                          # rv_bat_path() reads the battery files from here

L8_REPS    <- c(A = 500L, B = 1000L)                                 # section 2
L8_MINVER  <- "2.6.0"                                                # section 1: the corrected moment reference
L8_BATCH   <- 25L                                                    # replicates per worker per batch
L8_PAIR    <- "EDGE.poly3.u.Grule"                                   # rule 5: EDGE-poly3, unit form, rule G
L8_NAMES   <- c("rep", "seed", "n", "events", "id.stored", "id.regen", "id.absdiff", "id_ok", "flag.degenerate",
                "lecessie", "lc.stat", "lc.df", "sec", "flag.error")

## ---- the cells ------------------------------------------------------------------------------------------------------
l8_cells <- function() {
  A <- Cells8
  A$part <- "A"; A$R <- L8_REPS[["A"]]; A$k <- NA_integer_
  B <- b9_cells()
  B <- B[B$truth == "logit" & B$n == 1000L & (B$corruption == "clean" | (B$corruption == "C1" & B$k %in% c(1L, 2L, 5L, 10L))), ]
  B$part <- "B"; B$R <- L8_REPS[["B"]]
  B$role <- ifelse(B$corruption == "clean", "clean", "corrupted")
  B$null_block <- NA_character_; B$null_cell <- NA_character_
  keep <- c("part", "block", "cell", "role", "n", "k", "null_block", "null_cell", "seed_base", "R")
  C <- rbind(A[, keep], B[, keep])
  rownames(C) <- NULL
  if (sum(C$part == "A") != 30L || sum(C$part == "B") != 5L || anyDuplicated(C$cell) || any(C$n > 1000L))
    stop("block 8L: the cell list does not match section 2 (30 + 5 cells, n <= 1000)")
  C
}
## the full generator row of a cell (the battery's or block 9's), with part and R added
l8_cell <- function(ce) {
  full <- if (ce$part == "A") as.list(Cells8[Cells8$cell == ce$cell, ]) else { z <- b9_cells(); as.list(z[z$cell == ce$cell, ]) }
  full$part <- ce$part; full$R <- ce$R
  full
}

## ---- one replicate --------------------------------------------------------------------------------------------------
## Part B's data set by block 9's path: set.seed(seed_base + rep), b9_gen, the battery's working fit, EDGE-poly3 unit at G = 10
l8_regen_b9 <- function(rep, cell) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  set.seed(cell$seed_base + rep)
  g <- b9_gen(cell$truth, cell$n, cell$corruption, cell$rate)
  dat <- list(d = g$d, f = g$f)
  n <- nrow(dat$d); ev <- sum(dat$d$y)
  deg <- min(ev, n - ev) == 0
  fq <- if (deg) NULL else bt_fit(dat)
  p <- if (is.null(fq)) NA_real_ else tryCatch(bt_edge_arm(fq, 10L, cell)[["EDGE.poly3.u"]], error = function(e) NA_real_)
  list(dat = dat, fq = fq, n = n, events = ev, degenerate = deg, p = p, k_corrupt = length(g$corrupt))
}
l8_regen <- function(rep, cell) if (cell$part == "A") rv_regen(rep, cell) else l8_regen_b9(rep, cell)

## le Cessie through the package's public battery; a missing row, an error or no p-value gives NA (rule 2)
l8_lecessie <- function(fit) {
  r <- tryCatch(suppressWarnings(suppressMessages(ebrahim.gof::run.all.gof(fit, tests = "le-Cessie", install = "no"))),
                error = function(e) NULL)
  row <- if (is.null(r)) NULL else r[r$Test == "le-Cessie", , drop = FALSE]
  if (is.null(row) || nrow(row) != 1L) return(c(p = NA_real_, stat = NA_real_, df = NA_real_))
  p <- rv_num1(row$p_value)
  c(p = if (is.finite(p) && p >= 0 && p <= 1) p else NA_real_, stat = rv_num1(row$Statistic), df = rv_num1(row$df))
}

l8_one <- function(rep, cell, st) {
  out <- setNames(rep(NA_real_, length(L8_NAMES)), L8_NAMES)
  seed <- cell$seed_base + rep
  out["rep"] <- rep; out["seed"] <- seed; out["id_ok"] <- 0; out["id.stored"] <- suppressWarnings(as.numeric(st$p))
  g <- tryCatch(l8_regen(rep, cell), error = function(e) NULL)
  if (is.null(g)) return(out)                                        # the generator failed: the data cannot be reproduced
  out["n"] <- g$n; out["events"] <- g$events; out["id.regen"] <- g$p; out["flag.degenerate"] <- as.numeric(g$degenerate)
  if (is.finite(out[["id.stored"]]) && is.finite(g$p)) out["id.absdiff"] <- abs(out[["id.stored"]] - g$p)
  if (!rv_identity(g, st, seed)) return(out)                         # identity failure: le Cessie is not run
  out["id_ok"] <- 1
  if (is.null(g$fq)) { out["flag.error"] <- as.numeric(!g$degenerate); return(out) }
  t0 <- proc.time()[["elapsed"]]
  lc <- l8_lecessie(g$fq$fit)
  out["sec"] <- proc.time()[["elapsed"]] - t0
  out["lecessie"] <- lc[["p"]]; out["lc.stat"] <- lc[["stat"]]; out["lc.df"] <- lc[["df"]]
  out["flag.error"] <- as.numeric(!is.finite(lc[["p"]]))
  out
}
l8_task <- function(task, cell) l8_one(task$rep, cell, task$st)

## ---- the installed le Cessie is the corrected one ---------------------------------------------------------------------
## A fixed, RNG-free data set; the statistic from its definition with dense n x n algebra. H = V X (X'VX)^-1 X', the
## declared moment matrix M = (I - H)' R (I - H), the vendor's (I - H) R (I - H). Returns the three p-values.
l8_transpose_check <- function() {
  n <- 60L; x <- seq(-3, 3, length.out = n); d <- rep(0:1, n / 2)
  u <- ((seq_len(n) * 37L) %% 59L + 0.5) / 59                        # a deterministic stand-in for uniforms
  y <- as.integer(u < plogis(-0.2 + 0.8 * x + 0.6 * d))
  fit <- suppressWarnings(glm(y ~ x + d, family = binomial()))
  ph <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6); r <- y - ph; V <- ph * (1 - ph)
  X <- model.matrix(fit); I <- diag(n)
  H <- (V * X) %*% solve(crossprod(X, V * X)) %*% t(X)
  D <- sqrt(0.5) * as.matrix(dist(scale(cbind(x, d))))               # the package's scaled distance
  K <- pmax(1 - D / mean(D), 0)
  Q <- as.numeric(r %*% K %*% r)
  mom <- function(M) {
    M <- (M + t(M)) / 2; dM <- diag(M)
    E <- sum(dM * V); Var <- sum(dM^2 * (V * (1 - 3 * V) - 3 * V^2)) + 2 * sum(M * M * outer(V, V))
    Tt <- Q * 2 * E / Var; df <- 2 * E^2 / Var
    stats::pchisq(Tt, df, lower.tail = FALSE)
  }
  c(package = l8_lecessie(fit)[["p"]], declared = mom(t(I - H) %*% K %*% (I - H)), vendor = mom((I - H) %*% K %*% (I - H)))
}

## ---- summary --------------------------------------------------------------------------------------------------------
## Rules (section 3). Rejection at p <= 0.05; no p-value = no rejection. Size band: 0.05 +/- 3 sqrt(0.05 x 0.95 / B) on
## the kept replicates. Size-adjusted power: the critical value is the 5% quantile (type 1) of the matched null's p-values,
## no p-value counted as 1, for le Cessie and for EDGE alike, each on its own kept replicates -- block 8's rule. Pairing:
## le Cessie against EDGE-poly3 unit at rule G on the same replicates, exact McNemar at 0.05. The Holm family is every cell
## of block 8L, as the declaration says ("within block 8L"); the null and clean cells are size comparisons and are in it,
## which can only make the correction more conservative for the others.
L8_A <- 0.05
l8_rej  <- function(p) is.finite(p) & p <= L8_A
l8_band <- function(B) L8_A + c(-3, 3) * sqrt(L8_A * (1 - L8_A) / B)
l8_crit <- function(pn) as.numeric(stats::quantile(ifelse(is.finite(pn), pn, 1), L8_A, type = 1))

l8_summary <- function(OUT) {
  C <- l8_cells()
  rd <- function(f) if (file.exists(f)) as.data.frame(fread(f)) else NULL
  lc_path <- function(ce) file.path(OUT, sprintf("%s_lecessie_pvalues.csv.gz", ce$cell))
  S <- list(); PR <- list()
  kept <- function(P) if (is.null(P)) NULL else P[P$id_ok %in% 1, , drop = FALSE]
  for (i in seq_len(nrow(C))) {
    ce <- as.list(C[i, ])
    P <- rd(lc_path(ce)); if (is.null(P)) next
    K <- kept(P); Bf <- rd(rv_bat_path(ce))
    pe <- Bf[[L8_PAIR]][match(K$rep, Bf$rep)]                        # EDGE on the same replicates
    alt <- ce$part == "A" && !is.na(ce$null_cell)
    crit_lc <- crit_ed <- NA_real_; lvl_lc <- lvl_ed <- NA_real_
    if (alt) {
      nce <- as.list(C[C$cell == ce$null_cell, ])
      NK <- kept(rd(lc_path(nce))); BN <- rd(rv_bat_path(nce))
      if (!is.null(NK)) {
        crit_lc <- l8_crit(NK$lecessie); lvl_lc <- mean(is.finite(NK$lecessie) & NK$lecessie <= crit_lc)
        pen <- BN[[L8_PAIR]][match(NK$rep, BN$rep)]
        crit_ed <- l8_crit(pen); lvl_ed <- mean(is.finite(pen) & pen <= crit_ed)
      }
    }
    band <- l8_band(nrow(K)); rate <- mean(l8_rej(K$lecessie)); ran <- is.finite(K$sec)
    S[[length(S) + 1]] <- data.frame(part = ce$part, cell = ce$cell, role = ce$role, n = ce$n, k = ce$k,
      null_cell = ce$null_cell, reps_in_file = nrow(P), reps_kept = nrow(K), identity_failures = sum(P$id_ok %in% 0),
      max_identity_absdiff = suppressWarnings(max(K$id.absdiff, na.rm = TRUE)), no_pvalue = sum(!is.finite(K$lecessie)),
      degenerate = sum(K$flag.degenerate %in% 1),
      lecessie_rejection = rate, mcse = sqrt(rate * (1 - rate) / nrow(K)), band_lo = band[1], band_hi = band[2],
      in_band = if (ce$role %in% c("null", "clean")) rate >= band[1] && rate <= band[2] else NA,
      lecessie_size_adj = if (alt) mean(is.finite(K$lecessie) & K$lecessie <= crit_lc) else NA_real_,
      lecessie_null_crit = crit_lc, lecessie_null_level = lvl_lc,
      edge_rejection = mean(l8_rej(pe)),
      edge_size_adj = if (alt) mean(is.finite(pe) & pe <= crit_ed) else NA_real_,
      edge_null_crit = crit_ed, edge_null_level = lvl_ed,
      median_sec = if (any(ran)) stats::median(K$sec[ran]) else NA_real_, stringsAsFactors = FALSE)
    r1 <- l8_rej(K$lecessie); r2 <- l8_rej(pe); nb <- sum(r1 & !r2); nc <- sum(!r1 & r2)
    PR[[length(PR) + 1]] <- data.frame(part = ce$part, cell = ce$cell, role = ce$role, shared_reps = nrow(K),
      lecessie_rejection = mean(r1), edge_rejection = mean(r2), difference = mean(r1) - mean(r2),
      both = sum(r1 & r2), lecessie_only = nb, edge_only = nc, neither = sum(!r1 & !r2),
      mcnemar_p = if (nb + nc == 0) 1 else stats::binom.test(nb, nb + nc, 0.5)$p.value, stringsAsFactors = FALSE)
  }
  S <- rbindlist(S); PR <- rbindlist(PR)
  if (nrow(PR)) { PR$holm_p <- stats::p.adjust(PR$mcnemar_p, method = "holm"); PR$holm_reject_05 <- PR$holm_p <= L8_A }
  fwrite(S, file.path(OUT, "_summary.csv")); fwrite(PR, file.path(OUT, "_paired.csv"))
  cat(sprintf("summary: %d cells -> %s; paired: %d rows (Holm over all of them)\n", nrow(S), file.path(OUT, "_summary.csv"), nrow(PR)))
  invisible(list(summary = S, paired = PR))
}

## ---- options, guards ----------------------------------------------------------------------------------------------------
l8_opts <- function(a = commandArgs(TRUE)) {
  opt <- list(); i <- 1L
  while (i <= length(a)) {
    if (!startsWith(a[i], "--")) stop("unexpected argument ", a[i])
    key <- sub("^--", "", a[i])
    if (key %in% c("force", "summary")) { opt[[key]] <- TRUE; i <- i + 1L }
    else { if (i == length(a)) stop("--", key, " needs a value"); opt[[key]] <- a[i + 1L]; i <- i + 2L }
  }
  bad <- setdiff(names(opt), c("part", "cells", "workers", "reps", "out", "force", "summary"))
  if (length(bad)) stop("unknown option --", bad[1])
  opt
}
l8_abs <- function(p) normalizePath(ifelse(grepl("^([A-Za-z]:|/)", p), p, file.path(SIMDIR, p)), winslash = "/", mustWork = FALSE)

l8_log <- function(msg) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "", msg)
  cat(line, "\n", sep = ""); cat(line, "\n", sep = "", file = file.path(OUT, "_progress.log"), append = TRUE)
}

## every reason not to start, found before any computation
l8_check <- function(sel) {
  bad <- character(0)
  if (utils::packageVersion("ebrahim.gof") < L8_MINVER)
    bad <- c(bad, sprintf("ebrahim.gof %s is installed; section 1 needs >= %s", utils::packageVersion("ebrahim.gof"), L8_MINVER))
  tc <- tryCatch(l8_transpose_check(), error = function(e) c(package = NA, declared = NA, vendor = NA))
  if (!isTRUE(abs(tc[["package"]] - tc[["declared"]]) < 1e-10) || !isTRUE(abs(tc[["package"]] - tc[["vendor"]]) > 1e-6))
    bad <- c(bad, sprintf("the installed le Cessie is not the corrected (I-H)'R(I-H): package %s, declared %s, vendor %s",
                          format(tc[["package"]]), format(tc[["declared"]]), format(tc[["vendor"]])))
  for (i in seq_len(nrow(sel))) {
    ce <- as.list(sel[i, ]); f <- rv_bat_path(ce); R <- l8_R(ce)
    if (!file.exists(f)) { bad <- c(bad, sprintf("%s: its battery file %s does not exist", ce$cell, f)); next }
    S <- fread(f, select = c("rep", RV_ID))
    if (!all(seq_len(R) %in% S$rep)) bad <- c(bad, sprintf("%s: the battery file lacks some of replicates 1-%d", ce$cell, R))
    o <- l8_out_path(ce); pt <- paste0(o, ".part")
    if (file.exists(o)) {
      nr <- nrow(fread(o, select = "rep"))
      if (nr != R) bad <- c(bad, sprintf("%s exists with %d rows but %d replicates are planned: move it away", o, nr, R))
    } else if (file.exists(pt) || file.exists(paste0(pt, ".tmp"))) {
      P <- rv_read_part(pt)
      if (is.null(P) || !identical(names(P), L8_NAMES) || anyDuplicated(P$rep) || any(!P$rep %in% seq_len(R)))
        bad <- c(bad, sprintf("%s does not belong to this run (or cannot be read): move it away", pt))
    }
  }
  bad
}

l8_cluster <- function(w) {
  rv_pin_blas()
  cl <- makePSOCKcluster(w)
  parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
  clusterCall(cl, function(simdir) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
    options(rivals.source_only = TRUE)
    suppressPackageStartupMessages({ library(data.table) })
    source(file.path(simdir, "run_M_rivals.R")); source(file.path(simdir, "_block9_contam.R"))
    if (!requireNamespace("ebrahim.gof", quietly = TRUE)) stop("ebrahim.gof is not installed")
    RNGkind("L'Ecuyer-CMRG")
    NULL
  }, SIMDIR)
  clusterExport(cl, c("l8_regen", "l8_regen_b9", "l8_lecessie", "l8_one", "l8_task", "L8_NAMES", "ROOT"))
  cl
}

l8_run_cell <- function(cl, ce, prog) {
  R <- l8_R(ce); out <- l8_out_path(ce); part <- paste0(out, ".part")
  cell <- l8_cell(ce)
  S <- rv_stored(ce, R)
  M <- rv_read_part(part)
  todo <- setdiff(seq_len(R), M$rep)
  l8_log(sprintf("%s (part %s): %d of %d replicates to run%s", ce$cell, ce$part, length(todo), R,
                 if (is.null(M)) "" else sprintf(" (%d resumed from .part)", nrow(M))))
  bs <- length(cl) * L8_BATCH
  for (b in split(todo, ceiling(seq_along(todo) / bs))) {
    tasks <- lapply(b, function(r) list(rep = r, st = as.list(S[r, ])))
    rows <- clusterApplyLB(cl, tasks, l8_task, cell = cell)
    M <- rbind(M, as.data.frame(do.call(rbind, rows)))
    M <- M[order(M$rep), , drop = FALSE]; rownames(M) <- NULL
    tmp <- paste0(part, ".tmp")                                      # write whole, then swap in: a crash never truncates .part
    fwrite(M, tmp, compress = "gzip")
    if (file.exists(part)) file.remove(part)
    if (!file.rename(tmp, part)) stop("could not rename ", tmp)
    prog$add(length(b) * ce$n / 1000)
    l8_log(sprintf("  %-24s %d/%d replicates written to .part  ETA %s", ce$cell, nrow(M), R, prog$eta()))
  }
  if (is.null(M) || nrow(M) != R || any(M$rep != seq_len(R))) stop("internal: ", part, " is incomplete")
  if (file.exists(out)) file.remove(out)
  if (!file.rename(part, out)) stop("could not rename ", part)
  M
}

## ---- main -----------------------------------------------------------------------------------------------------------
## inside a function so that the lock and the cluster are released by on.exit on every return; the exit status is 1 if the
## run refused to start
l8_main <- function() {
  if (!isTRUE(OPT$summary) && !is.null(OPT$part) && !OPT$part %in% c("A", "B", "both")) stop("--part must be A, B or both")
  if (isTRUE(OPT$summary)) { l8_summary(OUT); return(0L) }
  part <- if (is.null(OPT$part)) "both" else OPT$part
  W <- if (is.null(OPT$workers)) 20L else as.integer(OPT$workers)
  if (is.na(W) || W < 1L) stop("--workers needs a positive integer")
  sel <- l8_cells()
  if (part != "both") sel <- sel[sel$part == part, ]
  if (!is.null(OPT$cells)) {
    want <- trimws(strsplit(OPT$cells, ",")[[1]])
    if (any(!want %in% sel$cell)) stop("not a block 8L cell of part ", part, ": ", paste(setdiff(want, sel$cell), collapse = ", "))
    sel <- sel[sel$cell %in% want, ]
  }
  l8_log(sprintf("block 8L: %d cells (part %s), %s, %d workers, ebrahim.gof %s, out %s", nrow(sel), part,
                 if (TEST_RUN) sprintf("replicates 1-%d (TEST RUN)", REPS_OVERRIDE) else "A 1-500, B 1-1000", W,
                 as.character(utils::packageVersion("ebrahim.gof")), OUT))
  bad <- l8_check(sel)
  if (length(bad)) { for (b in bad) l8_log(paste("refused:", b)); l8_log("block 8L not started"); return(1L) }
  tc <- l8_transpose_check()
  l8_log(sprintf("le Cessie check: package p %.12f = declared (I-H)'R(I-H) %.12f; vendor (I-H)R(I-H) would give %.12f",
                 tc[["package"]], tc[["declared"]], tc[["vendor"]]))

  lock <- file.path(OUT, "_lock")
  if (file.exists(lock) && !isTRUE(OPT$force)) {
    pid <- suppressWarnings(as.integer(readLines(lock, n = 1L)))
    if (rv_pid_alive(pid)) stop("another block 8L run (process ", pid, ") holds ", lock)
    l8_log(sprintf("a stale lock from process %s is being replaced", format(pid)))
  }
  writeLines(c(as.character(Sys.getpid()), format(Sys.time())), lock)
  on.exit(if (file.exists(lock)) file.remove(lock), add = TRUE)

  done <- vapply(seq_len(nrow(sel)), function(i) file.exists(l8_out_path(as.list(sel[i, ]))), logical(1))
  for (i in which(done)) l8_log(sprintf("%s: finished file present, skipped", sel$cell[i]))
  todo <- which(!done)
  if (!length(todo)) { l8_log("block 8L: nothing to run"); l8_log("block 8L: run finished"); return(0L) }
  prog <- rv_progress(sum(vapply(todo, function(i) l8_R(as.list(sel[i, ])) * sel$n[i] / 1000, numeric(1))))
  cl <- l8_cluster(W)
  on.exit(stopCluster(cl), add = TRUE)
  kk <- 0L
  for (i in todo) {
    ce <- as.list(sel[i, ]); kk <- kk + 1L; t1 <- Sys.time()
    M <- l8_run_cell(cl, ce, prog)
    ran <- M$id_ok %in% 1 & is.finite(M$sec)
    l8_log(sprintf("[%d/%d] %-24s n=%-5d R=%-4d kept %d, identity failures %d, no p-value %d, degenerate %d, rejection %.3f, median %.2f s, %.1f min  ETA %s",
                   kk, length(todo), ce$cell, ce$n, nrow(M), sum(M$id_ok %in% 1), sum(M$id_ok %in% 0),
                   sum(M$id_ok %in% 1 & !is.finite(M$lecessie)), sum(M$flag.degenerate %in% 1),
                   mean(M$id_ok %in% 1 & is.finite(M$lecessie) & M$lecessie <= 0.05),
                   if (any(ran)) stats::median(M$sec[ran]) else NA_real_,
                   as.numeric(difftime(Sys.time(), t1, units = "mins")), prog$eta()))
  }
  l8_log("block 8L: run finished")
  0L
}

## the output folder and the replicate count; set here so that a sourcing script (the self-test) can set its own
l8_setup <- function(opt) {
  OPT <<- opt
  TEST_RUN <<- !is.null(opt$reps)
  REPS_OVERRIDE <<- if (TEST_RUN) as.integer(opt$reps) else NA_integer_
  if (TEST_RUN && (is.na(REPS_OVERRIDE) || REPS_OVERRIDE < 1L)) stop("--reps needs a positive integer")
  out <- l8_abs(if (!is.null(opt$out)) opt$out else edge_battery(if (TEST_RUN) "8L_test" else "8L"))
  if (TEST_RUN && identical(out, l8_abs(edge_battery("8L")))) stop("a test run (--reps) cannot write to the real block 8L folder")
  busy <- l8_abs(edge_battery(c(as.character(0:9), "1a", "1b", "9b", "9c", "9d", "dryrun", "analysis")))
  if (any(startsWith(paste0(out, "/"), paste0(busy, "/")))) stop("--out cannot be a folder of another block: ", out)
  if (!startsWith(paste0(out, "/"), paste0(l8_abs(BAT), "/"))) stop("--out must lie inside ", BAT)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  OUT <<- out
  invisible(out)
}
l8_R <- function(ce) if (TEST_RUN) min(REPS_OVERRIDE, ce$R) else ce$R
l8_out_path <- function(ce) file.path(OUT, sprintf("%s_lecessie_pvalues.csv.gz", ce$cell))

if (!isTRUE(getOption("block8L.source_only"))) {
  l8_setup(l8_opts())
  quit(save = "no", status = l8_main())
}
