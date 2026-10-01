## battery_review_degenerate.R -- what the harness returns on degenerate sparse samples (0, 1 or 2 events), where the
## working fit is separated and the clamp ph >= 1e-6 is active. Uses the B = 50 review files in battery/_review/e and
## 50 serial replicates of 1b sparse49_n100 (about a fifth of them have no event). Writes battery/_review/review_degenerate.log.
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
OUT <- edge_battery("_review"); ROOT <- file.path(OUT, "e")
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
sink(file.path(OUT, "review_degenerate.log"), split = TRUE)
options(width = 250)
RNGkind("L'Ecuyer-CMRG")
Cells <- battery_cells()
KEY <- c("rep", "events", "glm_conv", "EDGE.poly3.u.G10", "EDGE.poly3.sc.G10", "EDGE.stk.u.G10", "EDGE.stk.sc.G10", "EDGE.sym.u.G10",
         "EDGE.sym.sc.G10", "HL.G10", "HLF.G10", "Stk.joint", "Stk.LR", "Stk.sym1", "Stk.marg", "GiViTI", "Cubic.LR", "HL_w", "PH",
         "Tsiatis", "Xie", "flag.stk_half0", "flag.stk_diag")
show <- function(P, lab) {
  cat(sprintf("\n== %s: events table\n", lab)); print(table(P$events))
  z <- P[events <= 2, intersect(KEY, names(P)), with = FALSE]
  print(format(as.data.frame(z), digits = 3), row.names = FALSE)
}
show(fread(file.path(ROOT, "4", "sparse49_cloglog_n200_pvalues.csv.gz")), "4 sparse49_cloglog_n200 (B = 50 file)")
show(fread(file.path(ROOT, "1b", "sparse49_n200_pvalues.csv.gz")), "1b sparse49_n200 (B = 50 file)")

ce <- as.list(Cells[Cells$block == "1b" & Cells$cell == "sparse49_n100", ])
H <- as.data.table(do.call(rbind, lapply(1:50, battery_one, cell = ce)))
show(H, "1b sparse49_n100 (50 serial replicates)")
cat("\nrejection at 0.05 among samples with 0 events, per test (NA = no p-value):\n")
z0 <- H[events == 0]
if (nrow(z0)) {
  tc <- setdiff(names(H), c("rep", "seed", "n", "events", "fit_ok", "glm_conv")); tc <- tc[!grepl("^flag", tc)]
  r <- sapply(tc, function(t) c(returned = sum(is.finite(z0[[t]])), reject05 = sum(is.finite(z0[[t]]) & z0[[t]] <= 0.05)))
  print(t(r)[t(r)[, "returned"] > 0, , drop = FALSE])
}

## inside one zero-event and one one-event replicate
for (target in c(0, 1)) {
  r <- H$rep[H$events == target][1]
  if (is.na(r)) next
  set.seed(ce$seed_base + r); dat <- bt_data(ce); fq <- bt_fit(dat)
  cat(sprintf("\n-- replicate %d (%d events): glm converged %s in %d iterations; coef %s; raw fitted p range %.2e-%.2e; clamped ph range %.2e-%.2e\n",
              r, target, fq$fit$converged, fq$fit$iter, paste(signif(coef(fq$fit), 4), collapse = " "), min(fq$p_raw), max(fq$p_raw), min(fq$ph), max(fq$ph)))
  gs <- bt_groups(fq, 10)
  cat("   group pbar:", signif(gs$pbar, 3), "\n   group V:", signif(gs$Vg, 3), "\n   r:", signif(gs$r, 3), "\n")
  for (b in BT_BASES) {
    Z <- bt_basis(gs$pbar, b)
    if (is.null(Z)) { cat(sprintf("   %-5s basis NULL\n", b)); next }
    un <- tryCatch(bt_edge_unit(gs, fq$A, Z), error = function(e) list(p = NA, S = NA, lam = conditionMessage(e)))
    Zs <- Z * sqrt(gs$Vg); I <- bt_zoz(gs, fq$A, Zs); d <- sqrt(pmax(diag(I), 0))
    sc <- tryCatch(bt_edge_score(gs, fq$A, Z), error = function(e) list(p = NA, S = NA, k = conditionMessage(e)))
    lamtxt <- if (is.numeric(un$lam)) paste(signif(un$lam, 3), collapse = " ") else paste("error:", un$lam)
    cat(sprintf("   %-5s unit S %.3g p %.3g lam %s | score: diag(I) %s, k %s, S %.3g p %.3g\n", b, as.numeric(un$S), as.numeric(un$p), lamtxt,
                paste(signif(diag(I), 3), collapse = " "), sc$k, sc$S, sc$p))
  }
  js <- bt_stukel_joint_stat(fq)
  cat(sprintf("   Stukel joint chi %.3g on %d df; harness p %.3g; sym1 p %.3g; GiViTI p %.3g; HL p %.3g\n", js$chi, js$df,
              H[rep == r]$Stk.joint, H[rep == r]$Stk.sym1, H[rep == r]$GiViTI, H[rep == r]$HL.G10))
}
sink()
