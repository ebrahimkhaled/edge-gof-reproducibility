## analyse_paired_families.R -- mean difference in size-adjusted power, EDGE (cubic basis, unit form, default partition)
## minus each comparator, by departure family, with a Monte Carlo 95% interval (referee M6.1, 2026-10-01). Replaces the
## "leads or ties within 0.01" counts as the main reading of the census.
## The interval ignores the positive correlation between two tests read on the same replicates, so it is conservative:
## var(diff in a cell) <= [p1(1 - p1) + p2(1 - p2)] / B, and a family mean averages its cells.
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
suppressMessages(library(data.table))
A <- edge_battery("analysis")
CEN <- fread(file.path(A, "census_extended.csv"))
MEM <- unique(fread(file.path(A, "rule_A_membership.csv"))[, .(cell, family_name)], by = "cell")
BB  <- unique(fread(edge_battery("_estimate_cells.csv"))[, .(cell, B)], by = "cell")
EDGE <- "EDGE.poly3.u.Grule"
COMP <- c("EDGE.poly3.u.G10" = "EDGE at ten groups", HL.G10 = "Hosmer-Lemeshow, ten groups", Stk.joint = "Stukel joint score",
          Stk.sym1 = "Stukel one-parameter", Cubic.LR = "cubic calibration LR", GiViTI = "GiViTI belt")
W <- dcast(CEN[test %in% c(EDGE, names(COMP))], cell ~ test, value.var = "power")
W <- merge(merge(W, MEM, by = "cell"), BB, by = "cell", all.x = TRUE)
stopifnot(!anyNA(W$B))
R <- rbindlist(lapply(names(COMP), function(cmp) {
  d <- W[is.finite(get(EDGE)) & is.finite(get(cmp))]
  d[, diff := get(EDGE) - get(cmp)]
  d[, v := (get(EDGE) * (1 - get(EDGE)) + get(cmp) * (1 - get(cmp))) / B]
  f <- function(x) x[, .(cells = .N, mean_diff = mean(diff), se = sqrt(sum(v)) / .N)]
  rbind(d[, f(.SD), by = family_name], cbind(family_name = "all families", f(d)))[, comparator := COMP[[cmp]]]
}))
R[, `:=`(lo = mean_diff - 1.96 * se, hi = mean_diff + 1.96 * se)]
setcolorder(R, c("comparator", "family_name", "cells", "mean_diff", "lo", "hi", "se"))
print(R, digits = 3)
fwrite(R, file.path(A, "paired_families.csv"))
