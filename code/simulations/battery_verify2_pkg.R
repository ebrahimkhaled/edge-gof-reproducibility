## battery_verify2_pkg.R -- independent pre-launch checks (1) and (2).
## (1) For replicates 1..R of five real cells (seeds of the cell table), every EDGE form and G arm and Stukel's joint score
##     from battery_rep() against the installed ebrahim.gof (def.gof, gof_stukel form "joint"): equal to 1e-8 and the same
##     NA pattern, degenerate and guarded replicates included. battery_one() (the driver's path) must give the same row.
## (2) The same replicates through the harness snapshot of 2026-09-14 00:09 (before E8: no no-fit rule, no information
##     guard, unscaled 2.7.0 filter, clamped risks in the joint Stukel score, the A-solve shortcut for Z'Omega Z), and the
##     scaled against the unscaled kept columns inside the current harness (E8.3 scale invariance).
## Run: Rscript battery_verify2_pkg.R [R = 300] [workers = 16]
## Output: battery/_review/v2/pkg_compare_reps.csv.gz, pkg_compare.log
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

args <- commandArgs(trailingOnly = TRUE)
SIMDIR <- edge_path("code/simulations")
OUT <- edge_battery("_review", "v2")
OLD <- file.path(OUT, "old_harness_20260914_0009.R")
R  <- if (length(args) >= 1) as.integer(args[1]) else 300L
NW <- if (length(args) >= 2) as.integer(args[2]) else 16L
suppressPackageStartupMessages({ library(parallel); library(data.table) })
source(file.path(SIMDIR, "_battery_cells.R"))
C <- battery_cells()
WANT <- data.frame(block = c("1b", "1b", "2", "3", "4"),
                   cell = c("sparse49_n100", "sparse49_n200", "probit_auc_n3500", "stk_short_n1000", "sparse49_cloglog_n300"),
                   stringsAsFactors = FALSE)

worker_init <- function(simdir, old) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  source(file.path(simdir, "_battery_tests.R"))
  assign("OLDENV", new.env(parent = globalenv()), envir = globalenv())
  sys.source(old, envir = get("OLDENV", envir = globalenv()), keep.source = FALSE)
  RNGkind("L'Ecuyer-CMRG")
  invisible(NULL)
}

## bt_basis without the final unit-length scaling (the 2.7.0 filter alone)
basis_unscaled <- function(pbar, basis) {
  if (basis %in% c("poly2", "poly3")) {
    deg <- if (basis == "poly2") 2L else 3L
    if (length(unique(round(pbar, 8))) < deg + 1) return(NULL)
    Z <- as.matrix(stats::poly(pbar, deg))
  } else {
    e <- stats::qlogis(pbar)
    Z <- switch(basis, stk = cbind(e, e^2 * (e >= 0), -e^2 * (e < 0)), sym = cbind(e * abs(e)))
  }
  Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]
  if (!ncol(Z)) NULL else Z
}

verify_rep <- function(rep, cell) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(cell$seed_base + rep)
  dat <- bt_data(cell)
  n <- nrow(dat$d); ev <- sum(dat$d$y)
  h <- battery_rep(dat, cell)
  tests <- c(paste0("EDGE.", rep(BT_BASES, each = 4), ".", c("u", "sc"), ".", rep(c("G10", "Grule"), each = 2)), "Stk.joint")
  ## the driver's path: battery_one() draws the same data again
  b1 <- battery_one(rep, cell)[names(h)]
  one_diff <- if (identical(is.na(b1), is.na(h))) max(c(0, abs(b1 - h)), na.rm = TRUE) else Inf

  ## the pre-E8 harness snapshot on the same data
  o <- tryCatch(get("OLDENV", envir = globalenv())$battery_rep(dat, cell), error = function(e) NULL)
  otests <- setdiff(intersect(names(h), if (is.null(o)) character(0) else names(o)),
                    c("n", "events", "fit_ok", "glm_conv", grep("^flag\\.", names(h), value = TRUE)))

  ## the installed package on the same data
  fit <- suppressWarnings(stats::glm(dat$f, data = dat$d, family = stats::binomial()))
  G2 <- c(G10 = 10, Grule = bt_rule_G(n))
  pk <- c(); wn <- c()
  for (a in names(G2)) for (b in BT_BASES) for (w in c("unit", "score")) {
    key <- paste0("EDGE.", b, ".", if (w == "unit") "u" else "sc", ".", a)
    cls <- character(0)
    pk[key] <- tryCatch(withCallingHandlers(
      as.numeric(ebrahim.gof::def.gof(fit, G = G2[[a]], basis = if (b == "stk") "stukel" else b, weights = w)$p_value),
      warning = function(x) { cls <<- c(cls, class(x)[1]); invokeRestart("muffleWarning") }),
      error = function(e) { cls <<- c(cls, "ERROR"); NA_real_ })
    wn[key] <- paste(unique(cls), collapse = "+")
  }
  ns <- asNamespace("ebrahim.gof")
  ctx <- get(".gof_context", envir = ns)(fit, G = 10)
  sj <- tryCatch(get("gof_stukel", envir = ns)(ctx, list(form = "joint")),
                 error = function(e) list(p_value = NA_real_, df = NA_real_, Note = paste("ERROR:", conditionMessage(e))))
  pk["Stk.joint"] <- as.numeric(sj$p_value)

  ## diagnostics from the current harness's own pieces
  dg <- c(clamp = NA, g_edge_cols = NA, g_joint = NA, g_sym1 = NA, zero_info_cols = NA, joint_df = NA,
          scale_unit_maxdiff = NA, scale_score_maxdiff = NA, unscaled_unit_fail = NA, unscaled_score_fail = NA)
  if (ev > 0 && ev < n) {
    fq <- bt_fit(dat)
    if (!is.null(fq)) {
      dg["clamp"] <- as.numeric(any(fq$p_raw < 1e-6 | fq$p_raw > 1 - 1e-6))
      ge <- 0; zd <- 0; du <- 0; ds <- 0; fu <- 0; fs <- 0
      for (G in unique(G2)) {
        gs <- tryCatch(bt_groups(fq, G), error = function(e) NULL)
        if (is.null(gs)) next
        for (b in BT_BASES) {
          Z <- tryCatch(bt_basis(gs$pbar, b), error = function(e) NULL)
          if (is.null(Z)) next
          sc <- tryCatch(bt_edge_score(gs, fq$A, Z), error = function(e) NULL)
          if (!is.null(sc) && isTRUE(sc$guard)) ge <- ge + 1
          Zs <- Z * sqrt(gs$Vg); I <- crossprod(Zs, gs$Omega %*% Zs); zd <- zd + sum(diag(I) <= 0)
          Zu <- basis_unscaled(gs$pbar, b)
          pu_s <- tryCatch(bt_edge_unit(gs, fq$A, Z)$p, error = function(e) NA_real_)
          pu_u <- tryCatch(bt_edge_unit(gs, fq$A, Zu)$p, error = function(e) NA_real_)
          ps_s <- if (is.null(sc)) NA_real_ else sc$p
          ps_u <- tryCatch(bt_edge_score(gs, fq$A, Zu)$p, error = function(e) NA_real_)
          if (is.finite(pu_s) && is.finite(pu_u)) du <- max(du, abs(pu_s - pu_u)) else if (is.finite(pu_s)) fu <- fu + 1
          if (is.finite(ps_s) && is.finite(ps_u)) ds <- max(ds, abs(ps_s - ps_u)) else if (is.finite(ps_s)) fs <- fs + 1
        }
      }
      js <- bt_stukel_joint_stat(fq)
      W <- fq$ph * (1 - fq$ph); zs <- fq$eta * abs(fq$eta)
      I1 <- tryCatch({ zWX <- crossprod(fq$X, W * zs); sum(W * zs^2) - as.numeric(crossprod(zWX, solve(crossprod(fq$X, W * fq$X), zWX))) },
                     error = function(e) NA_real_)
      dg[c("g_edge_cols", "g_joint", "g_sym1", "zero_info_cols", "joint_df", "scale_unit_maxdiff", "scale_score_maxdiff",
           "unscaled_unit_fail", "unscaled_score_fail")] <-
        c(ge, as.numeric(isTRUE(js$guard)), as.numeric(is.finite(I1) && I1 <= 1e-10 * sum(W * zs^2)), zd, js$df, du, ds, fu, fs)
    }
  }
  num <- c(rep = rep, n = n, events = ev, one_diff = one_diff, old_ok = as.numeric(!is.null(o)),
           h[c("fit_ok", "glm_conv", "flag.degenerate", "flag.info_guard", "flag.stk_half0")],
           setNames(h[tests], paste0("h.", tests)), setNames(pk[tests], paste0("p.", tests)),
           if (!is.null(o)) setNames(o[otests], paste0("o.", otests)), setNames(h[otests], paste0("hx.", otests)), dg)
  list(num = num, warn = wn, stk_note = as.character(sj$Note), stk_df = as.numeric(sj$df))
}

sink(file.path(OUT, "pkg_compare.log"), split = TRUE)
cat("battery_verify2_pkg.R |", format(Sys.time()), "| ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "| R =", R, "| workers", NW, "\n")
cat("old harness snapshot:", OLD, "sha256", digest::digest(file = OLD, algo = "sha256"), "\n")
cl <- makeCluster(NW)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
invisible(clusterCall(cl, worker_init, SIMDIR, OLD))
clusterExport(cl, c("verify_rep", "basis_unscaled"))
ALL <- list()
for (i in seq_len(nrow(WANT))) {
  cell <- as.list(C[C$block == WANT$block[i] & C$cell == WANT$cell[i], ])
  t0 <- Sys.time()
  res <- parLapplyLB(cl, seq_len(R), verify_rep, cell = cell)
  M <- rbindlist(lapply(res, function(z) as.list(z$num)), fill = TRUE)
  W <- rbindlist(lapply(res, function(z) as.list(z$warn)))
  setnames(W, paste0("w.", names(W)))
  M <- cbind(data.table(block = cell$block, cell = cell$cell, seed = cell$seed_base + M$rep), M, W,
             stk_note = vapply(res, `[[`, "", "stk_note"), stk_df = vapply(res, `[[`, 0, "stk_df"))
  ALL[[i]] <- M
  cat(sprintf("%s %s: %d replicates in %.0f s\n", cell$block, cell$cell, R, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
stopCluster(cl)
A <- rbindlist(ALL, fill = TRUE)
fwrite(A, file.path(OUT, "pkg_compare_reps.csv.gz"))

tests <- sub("^h\\.", "", grep("^h\\.", names(A), value = TRUE))
cat("\n==== (1) harness battery_rep() against the installed package, tolerance 1e-8 ====\n")
S1 <- list()
for (cc in unique(A$cell)) {
  M <- A[cell == cc]
  cat(sprintf("\n-- %s %s: %d replicates, events %d-%d, degenerate %d, flag.info_guard %d, clamp active %d, Stukel half-column zero %d, battery_one() differs %d\n",
              M$block[1], cc, nrow(M), min(M$events), max(M$events), sum(M$flag.degenerate == 1), sum(M$flag.info_guard == 1, na.rm = TRUE),
              sum(M$clamp == 1, na.rm = TRUE), sum(M$flag.stk_half0 == 1, na.rm = TRUE), sum(M$one_diff != 0)))
  for (t in tests) {
    h <- M[[paste0("h.", t)]]; p <- M[[paste0("p.", t)]]
    both <- is.finite(h) & is.finite(p)
    md <- if (any(both)) max(abs(h[both] - p[both])) else NA_real_
    mism <- sum(is.finite(h) != is.finite(p))
    deg <- M$flag.degenerate == 1; grd <- M$flag.info_guard %in% 1
    S1[[length(S1) + 1]] <- data.frame(block = M$block[1], cell = cc, test = t, reps = nrow(M), both_finite = sum(both),
      both_na = sum(!is.finite(h) & !is.finite(p)), na_mismatch = mism, max_abs_diff = md,
      degenerate_reps = sum(deg), degenerate_both_na = sum(deg & !is.finite(h) & !is.finite(p)),
      guard_reps = sum(grd), guard_both_na = sum(grd & !is.finite(h) & !is.finite(p)), guard_both_finite_maxdiff =
        if (any(grd & both)) max(abs(h[grd & both] - p[grd & both])) else NA_real_,
      ok = mism == 0 && (is.na(md) || md <= 1e-8), stringsAsFactors = FALSE)
  }
}
S1 <- rbindlist(S1)
fwrite(S1, file.path(OUT, "pkg_compare_summary.csv"))
print(S1[, .(cells = uniqueN(cell), reps = sum(reps), both_finite = sum(both_finite), both_na = sum(both_na),
             na_mismatch = sum(na_mismatch), max_abs_diff = suppressWarnings(max(max_abs_diff, na.rm = TRUE)), not_ok = sum(!ok)), by = cell])
cat("\nper test (all five cells):\n")
print(S1[, .(both_finite = sum(both_finite), both_na = sum(both_na), na_mismatch = sum(na_mismatch),
             max_abs_diff = suppressWarnings(max(max_abs_diff, na.rm = TRUE)), degenerate_both_na = sum(degenerate_both_na),
             guard_reps = sum(guard_reps), guard_both_na = sum(guard_both_na)), by = test])
cat("rows not ok:", sum(!S1$ok), "\n")
if (any(!S1$ok)) print(S1[ok == FALSE])

cat("\n-- degenerate replicates (no event or no non-event)\n")
D <- A[flag.degenerate == 1]
wcols <- grep("^w\\.", names(A), value = TRUE)
cat(sprintf("replicates %d; harness all 17 NA: %s; package all 17 NA: %s; package def.gof calls with a def_degenerate warning: %d of %d; Stukel notes: %s\n",
            nrow(D), all(!is.finite(as.matrix(D[, paste0("h.", tests), with = FALSE]))),
            all(!is.finite(as.matrix(D[, paste0("p.", tests), with = FALSE]))),
            sum(grepl("def_degenerate", as.matrix(D[, wcols, with = FALSE]))), length(wcols) * nrow(D),
            paste(unique(D$stk_note), collapse = " | ")))

cat("\n-- guarded replicates (flag.info_guard = 1)\n")
G <- A[flag.info_guard %in% 1]
if (nrow(G)) {
  for (k in seq_len(nrow(G))) {
    z <- G[k]
    na_h <- tests[!is.finite(unlist(z[, paste0("h.", tests), with = FALSE]))]
    cat(sprintf("%s rep %d events %d: edge guard columns %s, joint guard %s, sym1 guard %s, joint df harness %s / package %s (%s); NA in harness: %s; package warnings with def_no_information: %s\n",
                z$cell, z$rep, z$events, z$g_edge_cols, z$g_joint, z$g_sym1, z$joint_df, z$stk_df, z$stk_note,
                if (length(na_h)) paste(na_h, collapse = ",") else "none",
                paste(sub("^w\\.", "", wcols[grepl("def_no_information", unlist(z[, wcols, with = FALSE]))]), collapse = ",")))
  }
}

cat("\n-- flag completeness: a score-form or joint p-value NA on a non-degenerate replicate without flag.info_guard\n")
sc_t <- grep("\\.sc\\.|Stk.joint", tests, value = TRUE)
F <- A[flag.degenerate == 0 & !(flag.info_guard %in% 1)]
nf <- F[, .(rep, cell, events, zero_info_cols, na = apply(!is.finite(as.matrix(.SD)), 1, function(r) paste(sc_t[r], collapse = ","))), .SDcols = paste0("h.", sc_t)]
nf <- nf[nzchar(na)]
cat("replicates:", nrow(nf), "\n"); if (nrow(nf)) print(nf)

cat("\n==== (2a) scaled against unscaled kept columns inside the current harness (E8.3) ====\n")
print(A[, .(reps = .N, unit_maxdiff = max(scale_unit_maxdiff, na.rm = TRUE), score_maxdiff = max(scale_score_maxdiff, na.rm = TRUE),
            unscaled_unit_fail = sum(unscaled_unit_fail, na.rm = TRUE), unscaled_score_fail = sum(unscaled_score_fail, na.rm = TRUE)), by = cell])

cat("\n==== (2b) current harness against the 00:09 snapshot on non-degenerate, non-guarded replicates ====\n")
otests <- sub("^o\\.", "", grep("^o\\.", names(A), value = TRUE))
S2 <- list()
for (cc in unique(A$cell)) for (t in otests) {
  M <- A[cell == cc & flag.degenerate == 0 & !(flag.info_guard %in% 1)]
  h <- M[[paste0("hx.", t)]]; o <- M[[paste0("o.", t)]]
  if (is.null(h) || is.null(o)) next
  both <- is.finite(h) & is.finite(o); d <- abs(h - o)
  cl <- M$clamp %in% 1
  S2[[length(S2) + 1]] <- data.frame(cell = cc, test = t, reps = nrow(M), both_finite = sum(both), na_mismatch = sum(is.finite(h) != is.finite(o)),
    max_abs_diff = if (any(both)) max(d[both]) else NA_real_, n_gt_1e8 = sum(both & d > 1e-8), n_gt_1e6 = sum(both & d > 1e-6),
    max_abs_diff_noclamp = if (any(both & !cl)) max(d[both & !cl]) else NA_real_, n_gt_1e6_noclamp = sum(both & !cl & d > 1e-6),
    clamp_reps = sum(cl), stringsAsFactors = FALSE)
}
S2 <- rbindlist(S2)
fwrite(S2, file.path(OUT, "old_harness_compare_summary.csv"))
print(S2[max_abs_diff > 1e-10 | na_mismatch > 0], nrow = 500)
cat("\nall (cell, test) rows with max |diff| <= 1e-10 and no NA mismatch:", sum(S2$max_abs_diff <= 1e-10 & S2$na_mismatch == 0, na.rm = TRUE),
    "of", nrow(S2), "\n")
cat("\n-- replicates with a difference above 1e-6 (test, cell, rep, events, clamp, current, old)\n")
for (cc in unique(A$cell)) for (t in otests) {
  M <- A[cell == cc & flag.degenerate == 0 & !(flag.info_guard %in% 1)]
  h <- M[[paste0("hx.", t)]]; o <- M[[paste0("o.", t)]]
  if (is.null(h)) next
  j <- which((is.finite(h) & is.finite(o) & abs(h - o) > 1e-6) | (is.finite(h) != is.finite(o)))
  for (k in j) cat(sprintf("  %-20s %-22s rep %3d events %3d clamp %s  current %s  old %s\n", t, cc, M$rep[k], M$events[k], M$clamp[k],
                           format(h[k], digits = 6), format(o[k], digits = 6)))
}
cat("\n-- the snapshot on degenerate replicates (tests with a finite p-value; share with p <= 0.05)\n")
if (nrow(D)) {
  fin <- sapply(otests, function(t) sum(is.finite(D[[paste0("o.", t)]])))
  rej <- sapply(otests, function(t) sum(is.finite(D[[paste0("o.", t)]]) & D[[paste0("o.", t)]] <= 0.05))
  print(data.frame(test = otests, finite = fin, p_le_05 = rej)[fin > 0, ], row.names = FALSE)
}
sink()
