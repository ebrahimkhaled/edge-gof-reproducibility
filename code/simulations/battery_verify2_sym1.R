## battery_verify2_sym1.R -- follow-up of battery_verify2_guard.R: guarded replicates whose Stk.sym1 is NA although the sym1
## guard condition (0 < I1 <= 1e-10 * sum(W zs^2)) did not hold, i.e. post-fit information I1 <= 0 (dropped without a flag).
## Output: battery/_review/v2/sym1_check.log
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
source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
RNGkind("L'Ecuyer-CMRG")
C <- battery_cells()
S <- fread(file.path(OUT, "guard_scan.csv"))
sink(file.path(OUT, "sym1_check.log"), split = TRUE)
z <- S[g_joint == 1 & g_sym1 == 0]
cat("guarded replicates with the joint guard but not the sym1 guard:", nrow(z), "\n")
z2 <- S[g_sym1 == 1 & g_joint == 0]
cat("guarded replicates with the sym1 guard but not the joint guard:", nrow(z2), "\n")
for (k in seq_len(nrow(z))) {
  cell <- as.list(C[C$block == z$block[k] & C$cell == z$cell[k], ])
  set.seed(cell$seed_base + z$rep[k]); dat <- bt_data(cell); fq <- bt_fit(dat)
  W <- fq$ph * (1 - fq$ph); zs <- fq$eta * abs(fq$eta)
  zWX <- crossprod(fq$X, W * zs); I0 <- sum(W * zs^2)
  I1 <- I0 - as.numeric(crossprod(zWX, solve(crossprod(fq$X, W * fq$X), zWX)))
  h <- battery_rep(dat, cell)
  cat(sprintf("%s rep %d events %d: I1 %.3e, unadjusted %.3e, ratio %.3e; Stk.sym1 %s; flag.info_guard %d (set by the joint guard)\n",
              z$cell[k], z$rep[k], z$events[k], I1, I0, I1 / I0, format(h[["Stk.sym1"]]), as.integer(h[["flag.info_guard"]])))
}
sink()
