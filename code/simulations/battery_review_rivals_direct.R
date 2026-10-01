## battery_review_rivals_direct.R -- independent review of block 8 (E9): direct rival calls on regenerated battery data.
## One R process, BLAS on one thread; run only after the driver test run has finished (never more than 2 computing processes).
##   (1) projection test, 2 data sets: cauchit_auc_n460 replicate 1 (compared with the driver's test run on the real battery
##       file) and loglog_base_n610 replicate 1 (compared with the builder's driver file). Data from my own generator, RNG
##       continued straight after the data (no EDGE step), the n x n kernel timed alone.
##   (2) BAGofT cost scaling in n: one BAGofT_multi pass (nsplits = 100, the package default; nsim = 0) at n = 380 and at
##       n = 1000. A full data set at the defaults is 1 + nsim = 101 such passes.
## Writes battery/_review/rivals/review_direct.csv and review_direct.log.
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
OUTD <- edge_battery("_review", "rivals")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_proj_test.R"))
LOG <- file.path(OUTD, "review_direct.log")
say <- function(...) { line <- sprintf("%s  %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)); cat(line, "\n", sep = "")
  cat(line, "\n", file = LOG, append = TRUE, sep = "") }
seed_base_of <- function(block, cell) { S <- fread(edge_battery(block, paste0(cell, "_pvalues.csv.gz")), select = c("rep", "seed"))
  unique(S$seed - S$rep) }
gen_design <- function(n, s, c0, pfun) { x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); data.frame(x = x, d = d, y = rbinom(n, 1, pfun(c0 + s * (0.6 * x + 0.5 * d)))) }
rows <- list()

## ---- (1) projection test --------------------------------------------------------------------------------------------------------
pj <- list(list(block = "2", cell = "cauchit_auc_n460", n = 460L, s = 2, pfun = pcauchy,
                ref = file.path(OUTD, "drv_real", "cauchit_auc_n460_proj_pvalues.csv.gz")),
           list(block = "2", cell = "loglog_base_n610", n = 610L, s = 1, pfun = function(e) exp(-exp(-e)),
                ref = edge_battery("8_test", "proj", "loglog_base_n610_proj_pvalues.csv.gz")))
for (z in pj) {
  sb <- seed_base_of(z$block, z$cell)
  RNGkind("L'Ecuyer-CMRG"); set.seed(sb + 1)
  D <- gen_design(z$n, z$s, 0, z$pfun)
  fit <- suppressWarnings(glm(y ~ x + d, data = D, family = binomial()))
  X <- model.matrix(fit); y <- fit$y
  s0 <- .Random.seed
  tA <- system.time(invisible(proj_build_A(X)))[["elapsed"]]
  kernel_no_rng <- identical(s0, .Random.seed)
  tP <- system.time(zz <- proj_pvalue(y, X, B = 250))[["elapsed"]]
  R <- if (file.exists(z$ref)) as.data.frame(fread(z$ref)) else NULL
  pd <- if (is.null(R)) NA_real_ else R$proj[R$rep == 1]
  sd <- if (is.null(R)) NA_real_ else R$aux.proj_stat[R$rep == 1]
  say("proj %s rep 1: direct p %.3f (stat %.6g), driver p %s; equal %s; kernel draws no RNG %s; kernel %.1f s, full call %.1f s",
      z$cell, zz$p_value, zz$stat, format(pd), isTRUE(abs(zz$p_value - pd) < 1e-12 && abs(zz$stat - sd) <= 1e-12 * abs(zz$stat)),
      kernel_no_rng, tA, tP)
  rows[[length(rows) + 1]] <- data.frame(what = "proj direct", cell = z$cell, n = z$n, p_direct = zz$p_value, p_driver = pd,
    equal = isTRUE(abs(zz$p_value - pd) < 1e-12), kernel_sec = tA, total_sec = tP, stringsAsFactors = FALSE)
}

## ---- (2) BAGofT: one pass of 100 splits at n = 380 and n = 1000 -------------------------------------------------------------------
bg <- list(list(block = "2", cell = "loglog_auc_n380", n = 380L, s = 2), list(block = "2", cell = "loglog_base_n1000", n = 1000L, s = 1))
for (z in bg) {
  sb <- seed_base_of(z$block, z$cell)
  RNGkind("L'Ecuyer-CMRG"); set.seed(sb + 1)
  D <- gen_design(z$n, z$s, 0, function(e) exp(-exp(-e)))
  t1 <- system.time(b <- suppressWarnings(suppressMessages(
    BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = y ~ x + d, link = "logit"), data = D, nsim = 0))))[["elapsed"]]
  say("BAGofT one pass (nsplits 100, nsim 0) %s: %.1f s; mean split p %.3f; x 101 = %.0f s per data set at the defaults",
      z$cell, t1, b$pmean, 101 * t1)
  rows[[length(rows) + 1]] <- data.frame(what = "BAGofT one pass (nsim 0)", cell = z$cell, n = z$n, p_direct = NA, p_driver = NA,
    equal = NA, kernel_sec = NA, total_sec = t1, stringsAsFactors = FALSE)
}
fwrite(do.call(rbind, rows), file.path(OUTD, "review_direct.csv"))
say("done")
