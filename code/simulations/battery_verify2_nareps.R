## battery_verify2_nareps.R -- follow-up of battery_verify2_pkg.R: why score-form or joint p-values are NA on replicates without
## flag.info_guard, and whether the raw-risk change of the joint Stukel score moves any rejection.
## Output: battery/_review/v2/nareps.log
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
OUT <- edge_battery("_review", "v2")
suppressPackageStartupMessages(library(data.table))
A <- fread(file.path(OUT, "pkg_compare_reps.csv.gz"))
sink(file.path(OUT, "nareps.log"), split = TRUE)
sc <- grep("^h\\.EDGE\\..*\\.sc\\.|^h\\.Stk\\.joint$", names(A), value = TRUE)
F <- A[flag.degenerate == 0 & !(flag.info_guard %in% 1)]
na_any <- apply(!is.finite(as.matrix(F[, sc, with = FALSE])), 1, any)
F <- F[na_any]
cat("non-degenerate replicates without the guard flag and with a score-form or joint NA:", nrow(F), "\n")
cat("  of these, poly3 unit form also NA at G10:", sum(!is.finite(F$h.EDGE.poly3.u.G10)), "; poly2 unit NA:", sum(!is.finite(F$h.EDGE.poly2.u.G10)),
    "; stk/sym score finite:", sum(is.finite(F$h.EDGE.stk.sc.G10) & is.finite(F$h.EDGE.sym.sc.G10)), "\n")
cat("  package warning text on EDGE.poly3.sc.G10:", paste(unique(F$w.EDGE.poly3.sc.G10), collapse = " | "), "\n")
cat("  clamp active in:", sum(F$clamp %in% 1), "\n")
J <- A[flag.degenerate == 0 & !is.finite(h.Stk.joint)]
cat("\nnon-degenerate replicates with Stk.joint NA:\n"); print(J[, .(cell, rep, events, flag.info_guard, flag.stk_half0, g_joint, joint_df, stk_df, stk_note, p.Stk.joint)])

cat("\njoint Stukel score, current (raw risks) against the 00:09 snapshot (clamped risks), non-degenerate replicates:\n")
Z <- A[flag.degenerate == 0 & is.finite(h.Stk.joint) & is.finite(o.Stk.joint)]
for (a in c(0.01, 0.05, 0.10))
  cat(sprintf("  alpha %.2f: rejections current %d, snapshot %d, replicates whose decision differs %d\n", a,
              sum(Z$h.Stk.joint <= a), sum(Z$o.Stk.joint <= a), sum((Z$h.Stk.joint <= a) != (Z$o.Stk.joint <= a))))
cat("  replicates that differ by more than 1e-6:", sum(abs(Z$h.Stk.joint - Z$o.Stk.joint) > 1e-6), "; all with the clamp active:",
    all(Z$clamp[abs(Z$h.Stk.joint - Z$o.Stk.joint) > 1e-6] == 1), "; smallest p among them (current / snapshot):",
    min(Z$h.Stk.joint[abs(Z$h.Stk.joint - Z$o.Stk.joint) > 1e-6]), "/", min(Z$o.Stk.joint[abs(Z$h.Stk.joint - Z$o.Stk.joint) > 1e-6]), "\n")
cat("  clamp-active replicates by cell:\n"); print(A[flag.degenerate == 0, .(reps = .N, clamp = sum(clamp %in% 1)), by = cell])
sink()
