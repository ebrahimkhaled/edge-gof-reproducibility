## battery_verify2_oldstored.R -- check (2), second route: the stored outputs of the pre-E8 harness (battery/_review/e, written
## 2026-09-14 05:46, B = 50) against battery_one() of the current harness on the same seeds.
## Output: battery/_review/v2/oldstored_compare.log, oldstored_compare.csv
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
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))
RNGkind("L'Ecuyer-CMRG")
C <- battery_cells()
E <- edge_battery("_review", "e")
files <- list.files(E, pattern = "_pvalues\\.csv\\.gz$", recursive = TRUE)
sink(file.path(OUT, "oldstored_compare.log"), split = TRUE)
cat("battery_verify2_oldstored.R |", format(Sys.time()), "\n")
res <- list(); detail <- list()
for (f in files) {
  blk <- dirname(f); cn <- sub("_pvalues\\.csv\\.gz$", "", basename(f))
  cell <- as.list(C[C$block == blk & C$cell == cn, ])
  S <- as.data.frame(fread(file.path(E, f)))
  t0 <- Sys.time()
  N <- as.data.frame(do.call(rbind, lapply(S$rep, battery_one, cell = cell)))
  cat(sprintf("\n-- %s %s: %d stored replicates, seeds = seed_base + rep: %s, same events: %s (%.0f s)\n", blk, cn, nrow(S),
              all(S$seed == cell$seed_base + S$rep & N$seed == S$seed), all(N$events == S$events),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  tests <- setdiff(intersect(names(S), names(N)), c("rep", "seed", "n", "events", "fit_ok", "glm_conv", grep("^flag\\.", names(S), value = TRUE)))
  cl_new <- N$flag.degenerate %in% 1 | N$flag.info_guard %in% 1
  cat(sprintf("   current harness: degenerate %d, info guard %d\n", sum(N$flag.degenerate %in% 1), sum(N$flag.info_guard %in% 1)))
  for (t in tests) {
    s <- S[[t]]; h <- N[[t]]; keep <- !cl_new
    both <- keep & is.finite(s) & is.finite(h)
    d <- abs(s - h)
    res[[length(res) + 1]] <- data.frame(block = blk, cell = cn, test = t, reps_kept = sum(keep), both_finite = sum(both),
      na_mismatch = sum(keep & is.finite(s) != is.finite(h)), max_abs_diff = if (any(both)) max(d[both]) else NA_real_,
      n_gt_1e6 = sum(both & d > 1e-6), excluded_reps = sum(!keep),
      excluded_changed = sum(!keep & !(is.finite(s) == is.finite(h) & (!is.finite(s) | d <= 1e-8))), stringsAsFactors = FALSE)
    j <- which((both & d > 1e-6) | (keep & is.finite(s) != is.finite(h)))
    for (k in j) detail[[length(detail) + 1]] <- sprintf("   %s %s rep %d events %d: stored %s current %s", cn, t, S$rep[k], S$events[k],
                                                        format(s[k], digits = 6), format(h[k], digits = 6))
    jx <- which(!keep & !(is.finite(s) == is.finite(h) & (!is.finite(s) | d <= 1e-8)))
    for (k in jx) detail[[length(detail) + 1]] <- sprintf("   [degenerate/guard] %s %s rep %d events %d: stored %s current %s", cn, t, S$rep[k],
                                                         S$events[k], format(s[k], digits = 6), format(h[k], digits = 6))
  }
}
R <- rbindlist(res)
fwrite(R, file.path(OUT, "oldstored_compare.csv"))
cat("\nper cell (non-degenerate, non-guarded replicates):\n")
print(R[, .(tests = .N, both_finite = sum(both_finite), na_mismatch = sum(na_mismatch), max_abs_diff = suppressWarnings(max(max_abs_diff, na.rm = TRUE)),
            n_gt_1e6 = sum(n_gt_1e6), excluded_reps = max(excluded_reps), excluded_changed = sum(excluded_changed)), by = .(block, cell)])
cat("\nrows with max |diff| > 1e-10 or an NA mismatch:\n")
print(R[max_abs_diff > 1e-10 | na_mismatch > 0], nrow = 400)
cat("\ndetail:\n"); cat(unlist(detail), sep = "\n")
sink()
