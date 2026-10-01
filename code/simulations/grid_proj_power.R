## grid_proj_power.R -- P0-3 (Statistics-and-Computing review, 2026-07-20):
## power head-to-head of the projection / RMEP test of Liu et al. (2024)
## against the pre-specified default EDGE-poly3 (+ EF, HL on the same draws).
##
## PRE-DECLARED grid (declared before running; no post-hoc selection):
##   * the 8 Table-2 scenarios ....... link/{cloglog,probit,stukel_heavy,stukel_light,stukel_asym},
##                                     quad/0.02, binint/0.3, contint/0.5
##   * the off-index scenario ........ omit_2cov (x~U(-2.5,2.5); eta=0.2+0.6x+0.8*z1+0.8*z2; fit y~x)
##                                     -- included deliberately: proj SHOULD win here (honest reporting).
##   * per-design NULLS .............. dgp_null for quad-, link-, and contint-design geometries
##   * n in {500, 1000};  REPS = 500 per cell (mcse <= 0.022);  MBB B = 250.
## proj implementation = the validated _proj_test.R (A built once per dataset, reused
## across all MBB replicates; reproduces Liu et al.'s UIS Example-2 rejection).
## Seed regime matches the paper harness: per-cell seed_base = GRID_SEED + cell_id*1e5,
## per-rep set.seed(seed_base + rep) under L'Ecuyer-CMRG cluster streams.
##
## Outputs (append-mode, RESUME-SAFE at the cell level -- rerun to continue):
##   proj_power_grid_pvalues.csv : scenario,param,n,rep,seed,p_proj,p_edge3,p_ef,p_hl
##   proj_power_grid.csv         : scenario,param,n,reps,B_mbb,reject_proj,reject_edge3,
##                                 reject_ef,reject_hl,mcse_proj
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

suppressMessages(library(parallel))
SIMDIR <- edge_path("code/simulations")
source(file.path(SIMDIR, "_harness.R"))
source(DGP)
source(file.path(SIMDIR, "_proj_test.R"))

GRID_SEED <- 20260720L
REPS      <- 500L
B_MBB     <- 250L
ALPHA     <- 0.05
PV_PATH   <- file.path(SIMDIR, "proj_power_grid_pvalues.csv")
SUM_PATH  <- file.path(SIMDIR, "proj_power_grid.csv")

## the pre-declared cells --------------------------------------------------------------
cells <- list(
  list(scn = "link",      param = "cloglog"),
  list(scn = "link",      param = "probit"),
  list(scn = "link",      param = "stukel_heavy"),
  list(scn = "link",      param = "stukel_light"),
  list(scn = "link",      param = "stukel_asym"),
  list(scn = "quad",      param = "0.02"),
  list(scn = "binint",    param = "0.3"),
  list(scn = "contint",   param = "0.5"),
  list(scn = "omit_2cov", param = "off-index"),
  list(scn = "null_quad",    param = "null"),
  list(scn = "null_link",    param = "null"),
  list(scn = "null_contint", param = "null")
)
NS <- c(500L, 1000L)

## dataset generator for one cell ------------------------------------------------------
gen_cell <- function(scn, param, n) {
  if (scn == "omit_2cov") {
    x <- runif(n, -2.5, 2.5); z1 <- rnorm(n); z2 <- rnorm(n)
    list(d = data.frame(x = x, y = rbinom(n, 1, plogis(0.2 + 0.6 * x + 0.8 * z1 + 0.8 * z2))),
         f = y ~ x)
  } else if (startsWith(scn, "null_")) {
    dgp_null(sub("^null_", "", scn), n)
  } else if (scn == "quad")   dgp_alt("quad",   as.numeric(param), n)
  else if (scn == "binint")   dgp_alt("binint", as.numeric(param), n)
  else if (scn == "contint")  dgp_alt("contint", as.numeric(param), n)
  else                        dgp_alt("link",   param, n)   # link family: param = link name
}

## one replicate: p-values of proj, EDGE-poly3, EF, HL on the SAME draw ----------------
one_rep <- function(rep, cell) {
  d   <- gen_cell(cell$scn, cell$param, cell$n)
  fit <- suppressWarnings(glm(d$f, data = d$d, family = binomial()))
  X   <- model.matrix(fit)
  y   <- fit$y
  p_proj <- tryCatch(proj_pvalue(y, X, B = B_MBB)$p_value, error = function(e) NA_real_)
  p_e3   <- tryCatch(suppressWarnings(edge.gof(fit, basis = "poly3"))$p_value, error = function(e) NA_real_)
  p_ef   <- tryCatch(suppressWarnings(ef.gof(fit))$p_value, error = function(e) NA_real_)
  p_hl   <- tryCatch({
    ph <- suppressWarnings(run.all.gof(fit, include_slow = FALSE))
    as.numeric(ph$p_value[ph$Test == "HL"][1])
  }, error = function(e) NA_real_)
  c(p_proj = p_proj, p_edge3 = p_e3, p_ef = p_ef, p_hl = p_hl)
}

## resume support ----------------------------------------------------------------------
done_keys <- character(0)
if (file.exists(SUM_PATH)) {
  ex <- read.csv(SUM_PATH)
  done_keys <- paste(ex$scenario, ex$param, ex$n)
}

cl <- ek_cluster(GRID_SEED)
clusterExport(cl, c("proj_build_A", "proj_pvalue", "gen_cell", "B_MBB"))
invisible(clusterEvalQ(cl, { source(file.path(
  edge_path("code/simulations"), "_proj_test.R")); NULL }))

cell_id <- 0L
for (cell0 in cells) for (n in NS) {
  cell_id <- cell_id + 1L
  key <- paste(cell0$scn, cell0$param, n)
  if (key %in% done_keys) { cat(sprintf("[skip] %-28s (done)\n", key)); next }
  cell <- c(cell0, list(n = n))
  seed_base <- GRID_SEED + cell_id * 100000L
  t0 <- Sys.time()
  cat(sprintf("[run ] %-28s reps=%d B=%d ...", key, REPS, B_MBB))
  M <- ek_run_cell(cl, REPS, seed_base, one_rep, cell)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  pv <- data.frame(scenario = cell0$scn, param = cell0$param, n = n,
                   rep = seq_len(REPS), seed = seed_base + seq_len(REPS),
                   p_proj = M$p_proj, p_edge3 = M$p_edge3, p_ef = M$p_ef, p_hl = M$p_hl)
  ek_append(pv, PV_PATH)
  ok <- !is.na(M$p_proj)
  sm <- data.frame(scenario = cell0$scn, param = cell0$param, n = n,
                   reps = REPS, B_mbb = B_MBB,
                   reject_proj  = mean(M$p_proj[ok] < ALPHA),
                   reject_edge3 = mean(M$p_edge3 < ALPHA, na.rm = TRUE),
                   reject_ef    = mean(M$p_ef   < ALPHA, na.rm = TRUE),
                   reject_hl    = mean(M$p_hl   < ALPHA, na.rm = TRUE),
                   n_proj_ok    = sum(ok),
                   mcse_proj    = ek_mcse(mean(M$p_proj[ok] < ALPHA), sum(ok)))
  ek_append(sm, SUM_PATH)
  cat(sprintf(" done in %.1f min | proj=%.3f edge3=%.3f ef=%.3f hl=%.3f\n",
              el, sm$reject_proj, sm$reject_edge3, sm$reject_ef, sm$reject_hl))
}
stopCluster(cl)
cat("ALL CELLS COMPLETE\n")
