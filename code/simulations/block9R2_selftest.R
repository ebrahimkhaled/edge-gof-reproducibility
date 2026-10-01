## block9R2_selftest.R -- checks the hybrid partition before block 9R2 runs.
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
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
source(file.path(SIMDIR, "_block9_contam.R")); source(file.path(SIMDIR, "_block9R2_hybrid.R"))

ok <- 0L; bad <- 0L
chk <- function(label, cond) { if (isTRUE(cond)) { ok <<- ok + 1L; cat(sprintf("  ok   %s\n", label)) }
                               else { bad <<- bad + 1L; cat(sprintf("  FAIL %s\n", label)) } }

set.seed(11); g <- b9_gen("logit", 5000, "none", 0); fq <- bt_fit(list(d = g$d, f = g$f))
n <- fq$n

s10 <- fr2_sizes(fq, "HYB10"); s05 <- fr2_sizes(fq, "HYB05"); s20 <- fr2_sizes(fq, "HYB20")
chk("H1  HYB10 tail groups hold n/10", s10[1] == 500L && s10[length(s10)] == 500L)
chk("H2  HYB05 tail groups hold n/20", s05[1] == 250L && s05[length(s05)] == 250L)
chk("H3  HYB20 tail groups hold n/5",  s20[1] == 1000L && s20[length(s20)] == 1000L)
chk("H4  HYB10 middle groups near 25", { m <- s10[-c(1, length(s10))]; all(m >= 20L & m <= 30L) })
chk("H5  every record used once",      sum(s10) == n && sum(s05) == n && sum(s20) == n)
chk("H6  HYB10 group count is 2 + 4n/125", length(s10) == 2L + round(0.8 * n / 25))
chk("H7  G10 is ten equal groups",     { s <- fr2_sizes(fq, "G10"); length(s) == 10L && all(s == n / 10) })
chk("H8  V0 is the rule G",            length(fr2_sizes(fq, "V0")) == max(10L, ceiling(n / 25)))

## the partition is monotone in fitted risk: no group straddles another
grp <- fr2_partition(fq, "HYB10")$grp
rng <- tapply(fq$ph, grp, range)
chk("H9  groups are intervals of risk", all(diff(vapply(rng, `[`, numeric(1), 1)) >= 0))

p <- fr2_pvalues(fq)
chk("H10 five p-values, all in [0,1]", length(p) == 5L && all(is.finite(p)) && all(p >= 0 & p <= 1))
chk("H11 V0 equals bt_edge_unit at the rule G", {
  gs <- bt_groups(fq, max(10L, ceiling(n / 25))); Z <- bt_basis(gs$pbar, "poly3")
  abs(bt_edge_unit(gs, fq$A, Z)$p - p[["V0"]]) <= 1e-12 })
chk("H12 G10 equals bt_edge_unit at G = 10", {
  gs <- bt_groups(fq, 10L); Z <- bt_basis(gs$pbar, "poly3")
  abs(bt_edge_unit(gs, fq$A, Z)$p - p[["G10"]]) <= 1e-12 })
chk("H13 the hybrids differ from both anchors", {
  abs(p[["HYB10"]] - p[["V0"]]) > 1e-10 && abs(p[["HYB10"]] - p[["G10"]]) > 1e-10 })

## the corrupted records land in the top tail group, and it dilutes them: k = 50 of n = 5000
set.seed(12); gc <- b9_gen("logit", 5000, "C1", 50 / 5000); fqc <- bt_fit(list(d = gc$d, f = gc$f))
top10 <- fr2_sizes(fqc, "HYB10"); topV0 <- fr2_sizes(fqc, "V0")
chk("H14 default top group holds about 25", topV0[length(topV0)] <= 30L)
chk("H15 hybrid top group holds about 500", top10[length(top10)] >= 450L)

## small n: the guard, not a crash
set.seed(13); gs2 <- b9_gen("logit", 200, "none", 0); fq2 <- bt_fit(list(d = gs2$d, f = gs2$f))
p2 <- fr2_pvalues(fq2)
chk("H16 n = 200 returns five p-values", length(p2) == 5L && all(is.finite(p2)))
chk("H17 n = 200 HYB05 tails hold 10",   { s <- fr2_sizes(fq2, "HYB05"); s[1] == 10L })

cat(sprintf("\nblock 9R2 self-test: %d ok, %d failed\n", ok, bad))
quit(status = if (bad) 1L else 0L)
