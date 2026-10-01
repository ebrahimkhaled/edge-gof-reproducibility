## battery_verify2_guard.R -- check (1) on guarded replicates. The 1,500 replicates of battery_verify2_pkg.R held no replicate with
## flag.info_guard, so this script scans the full seed range of the sparse cells for replicates where the information guard
## fires (a light scan with the harness's own pieces), confirms flag.info_guard with battery_rep(), and compares those
## replicates with the installed package (def.gof 4 bases x 2 forms x 2 G arms, gof_stukel joint) at 1e-8 with the NA pattern.
## Run: Rscript battery_verify2_guard.R [workers = 16] [max compared per cell = 60]
## Output: battery/_review/v2/guard_scan.csv, guard_compare_reps.csv, guard_compare.log
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
NW <- if (length(args) >= 1) as.integer(args[1]) else 16L
KMAX <- if (length(args) >= 2) as.integer(args[2]) else 60L
suppressPackageStartupMessages({ library(parallel); library(data.table) })
source(file.path(SIMDIR, "_battery_cells.R"))
C <- battery_cells()
SCAN <- data.frame(block = c("1b", "1b", "1b", "1b", "1b", "4", "4", "4", "5", "5", "1b"),
                   cell = c("sparse49_n100", "sparse49_n150", "sparse49_n200", "sparse49_n300", "sparse49_n500",
                            "sparse49_cloglog_n200", "sparse49_cloglog_n300", "sparse49_cloglog_n500",
                            "runG_null_sparse_n500", "runG_cloglog_sparse_n500", "null_quad_n100"), stringsAsFactors = FALSE)

worker_init <- function(simdir) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  source(file.path(simdir, "_battery_tests.R")); RNGkind("L'Ecuyer-CMRG"); invisible(NULL)
}

## light scan: the guard conditions of bt_edge_score, bt_stukel_joint_stat and Stk.sym1, no rival tests
scan_rep <- function(rep, cell) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(cell$seed_base + rep)
  dat <- bt_data(cell); n <- nrow(dat$d); ev <- sum(dat$d$y)
  if (min(ev, n - ev) == 0) return(c(rep = rep, events = ev, degenerate = 1, g_edge = 0, g_joint = 0, g_sym1 = 0))
  fq <- bt_fit(dat)
  if (is.null(fq)) return(c(rep = rep, events = ev, degenerate = 0, g_edge = NA, g_joint = NA, g_sym1 = NA))
  ge <- 0
  for (G in unique(c(10, bt_rule_G(n)))) {
    gs <- tryCatch(bt_groups(fq, G), error = function(e) NULL); if (is.null(gs)) next
    for (b in BT_BASES) {
      Z <- tryCatch(bt_basis(gs$pbar, b), error = function(e) NULL); if (is.null(Z)) next
      sc <- tryCatch(bt_edge_score(gs, fq$A, Z), error = function(e) NULL)
      if (!is.null(sc) && isTRUE(sc$guard)) ge <- ge + 1
    }
  }
  js <- bt_stukel_joint_stat(fq)
  W <- fq$ph * (1 - fq$ph); zs <- fq$eta * abs(fq$eta)
  I1 <- tryCatch({ zWX <- crossprod(fq$X, W * zs); sum(W * zs^2) - as.numeric(crossprod(zWX, solve(crossprod(fq$X, W * fq$X), zWX))) },
                 error = function(e) NA_real_)
  c(rep = rep, events = ev, degenerate = 0, g_edge = ge, g_joint = as.numeric(isTRUE(js$guard)),
    g_sym1 = as.numeric(is.finite(I1) && I1 > 0 && I1 <= 1e-10 * sum(W * zs^2)))
}

compare_rep <- function(rep, cell) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(cell$seed_base + rep)
  dat <- bt_data(cell); n <- nrow(dat$d)
  h <- battery_rep(dat, cell)
  b1 <- battery_one(rep, cell)[names(h)]
  one_diff <- if (identical(is.na(b1), is.na(h))) max(c(0, abs(b1 - h)), na.rm = TRUE) else Inf
  tests <- c(paste0("EDGE.", rep(BT_BASES, each = 4), ".", c("u", "sc"), ".", rep(c("G10", "Grule"), each = 2)), "Stk.joint")
  fit <- suppressWarnings(stats::glm(dat$f, data = dat$d, family = stats::binomial()))
  G2 <- c(G10 = 10, Grule = bt_rule_G(n)); pk <- c(); wn <- c()
  for (a in names(G2)) for (b in BT_BASES) for (w in c("unit", "score")) {
    key <- paste0("EDGE.", b, ".", if (w == "unit") "u" else "sc", ".", a); cls <- character(0)
    pk[key] <- tryCatch(withCallingHandlers(
      as.numeric(ebrahim.gof::def.gof(fit, G = G2[[a]], basis = if (b == "stk") "stukel" else b, weights = w)$p_value),
      warning = function(x) { cls <<- c(cls, class(x)[1]); invokeRestart("muffleWarning") }),
      error = function(e) { cls <<- c(cls, "ERROR"); NA_real_ })
    wn[key] <- paste(unique(cls), collapse = "+")
  }
  ns <- asNamespace("ebrahim.gof")
  sj <- tryCatch(get("gof_stukel", envir = ns)(get(".gof_context", envir = ns)(fit, G = 10), list(form = "joint")),
                 error = function(e) list(p_value = NA_real_, df = NA_real_, Note = paste("ERROR:", conditionMessage(e))))
  pk["Stk.joint"] <- as.numeric(sj$p_value)
  list(num = c(rep = rep, n = n, events = sum(dat$d$y), one_diff = one_diff, h[c("flag.degenerate", "flag.info_guard", "flag.stk_half0")],
               h.Stk.sym1 = h[["Stk.sym1"]], setNames(h[tests], paste0("h.", tests)), setNames(pk[tests], paste0("p.", tests))),
       warn = wn, stk_note = as.character(sj$Note), stk_df = as.numeric(sj$df))
}

sink(file.path(OUT, "guard_compare.log"), split = TRUE)
cat("battery_verify2_guard.R |", format(Sys.time()), "| ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "| workers", NW, "\n")
cl <- makeCluster(NW)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
invisible(clusterCall(cl, worker_init, SIMDIR))
clusterExport(cl, c("scan_rep", "compare_rep"))
SC <- list(); CMP <- list()
for (i in seq_len(nrow(SCAN))) {
  cell <- as.list(C[C$block == SCAN$block[i] & C$cell == SCAN$cell[i], ])
  t0 <- Sys.time()
  S <- as.data.table(do.call(rbind, parLapplyLB(cl, splitIndices(cell$B, NW * 8), function(r, cell) do.call(rbind, lapply(r, scan_rep, cell = cell)), cell = cell)))
  S[, `:=`(block = cell$block, cell = cell$cell)]
  SC[[i]] <- S
  gr <- S[g_edge > 0 | g_joint > 0 | g_sym1 > 0]
  cat(sprintf("%-3s %-26s B %5d: degenerate %4d, guard any %3d (edge %d, joint %d, sym1 %d), scan %.0f s\n", cell$block, cell$cell, nrow(S),
              sum(S$degenerate == 1), nrow(gr), sum(gr$g_edge > 0), sum(gr$g_joint > 0), sum(gr$g_sym1 > 0),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (nrow(gr)) {
    reps <- head(gr$rep, KMAX)
    res <- parLapplyLB(cl, reps, compare_rep, cell = cell)
    M <- rbindlist(lapply(res, function(z) as.list(z$num)), fill = TRUE)
    Wn <- rbindlist(lapply(res, function(z) as.list(z$warn))); setnames(Wn, paste0("w.", names(Wn)))
    CMP[[length(CMP) + 1]] <- cbind(data.table(block = cell$block, cell = cell$cell), M, Wn,
                                    stk_note = vapply(res, `[[`, "", "stk_note"), stk_df = vapply(res, `[[`, 0, "stk_df"))
  }
}
stopCluster(cl)
SC <- rbindlist(SC); fwrite(SC, file.path(OUT, "guard_scan.csv"))
if (!length(CMP)) { cat("no guarded replicate in any scanned cell\n"); sink(); quit(save = "no") }
A <- rbindlist(CMP, fill = TRUE); fwrite(A, file.path(OUT, "guard_compare_reps.csv"))
tests <- sub("^h\\.", "", grep("^h\\.(EDGE|Stk\\.joint)", names(A), value = TRUE))
`%||%` <- function(a, b) if (is.null(a)) b else a
cat(sprintf("\ncompared replicates: %d; flag.info_guard = 1 in battery_rep(): %d; battery_one() differs: %d\n", nrow(A),
            sum(A$flag.info_guard == 1), sum(A$one_diff != 0)))
T <- rbindlist(lapply(tests, function(t) {
  h <- A[[paste0("h.", t)]]; p <- A[[paste0("p.", t)]]; both <- is.finite(h) & is.finite(p)
  data.table(test = t, both_finite = sum(both), both_na = sum(!is.finite(h) & !is.finite(p)), na_mismatch = sum(is.finite(h) != is.finite(p)),
             max_abs_diff = if (any(both)) max(abs(h[both] - p[both])) else NA_real_,
             pkg_no_information_warnings = sum(grepl("def_no_information", A[[paste0("w.", t)]] %||% "")))
}))
`%||%` <- function(a, b) if (is.null(a)) b else a
print(T)
cat("not ok (NA mismatch or diff > 1e-8):", sum(T$na_mismatch > 0 | (is.finite(T$max_abs_diff) & T$max_abs_diff > 1e-8)), "\n")
cat("\nper replicate:\n")
for (k in seq_len(nrow(A))) {
  z <- A[k]
  nh <- tests[!is.finite(unlist(z[, paste0("h.", tests), with = FALSE]))]
  cat(sprintf("  %-24s rep %4d events %2d guard %d half0 %s | harness NA: %s | Stk.joint %s (pkg %s, df %s, %s) | Stk.sym1 %s\n", z$cell, z$rep, z$events,
              z$flag.info_guard, z$flag.stk_half0, if (length(nh)) paste(nh, collapse = ",") else "none", format(z$h.Stk.joint, digits = 4),
              format(z$p.Stk.joint, digits = 4), z$stk_df, z$stk_note, format(z$h.Stk.sym1, digits = 4)))
}
sink()
