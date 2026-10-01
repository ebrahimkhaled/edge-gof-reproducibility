## run_all.R — EDGE paper, Phase-A orchestrator.
##
## Runs the Phase-A simulation scripts IN DEPENDENCY ORDER. grid_null.R runs FIRST because it
## produces the size gate (sim_null.csv + sim_null_pvalues.csv) that grid_power_broad.R and
## null_calibration_checks.R consume.
##
## Each Phase-A script is a self-contained top-level program: it sources _harness.R + the DGP
## library itself and owns its own parallel cluster via ek_cluster() / on.exit(stopCluster).
## That on.exit teardown is ONLY correct when the file is the top-level program — sourcing or
## eval'ing it mid-session makes the on.exit fire against the wrong frame and tears the cluster
## down early ("invalid connection"), and bench_compute.R additionally self-guards its work
## behind `if (identical(environment(), globalenv()) && !interactive()) main()`, which a
## child-env source() would silently skip (a bogus PASS). To reproduce EXACT standalone
## semantics for every script we therefore run each one in its own Rscript subprocess,
## sequentially (one at a time, so each may use the full EK_NCORES box). The subprocess exit
## status is the ground-truth PASS/FAIL, and each script's own on.exit cleans up its cluster.
##
## This orchestrator session prints the sessionInfo()/packageVersion() header once, forwards
## REPS_SCALE (and any EK_* env vars already set) to every child, records PASS/FAIL + elapsed
## per script, then prints the CSV manifest (row counts) and cats _STOP_TRIGGERS.log.
##
## Honors Sys.getenv("REPS_SCALE") — set it BEFORE launching this script (e.g. "0.02" smoke,
## unset/"1" full). We do NOT override it here; children inherit the environment.
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
source(file.path(SIMDIR, "_harness.R"))   # SIMDIR, DGP, EK_NCORES, ek_reps, ek_cluster, ek_* helpers

STOP_LOG <- file.path(SIMDIR, "_STOP_TRIGGERS.log")
RSCRIPT  <- file.path(R.home("bin"), "Rscript")

## ---------------------------------------------------------------------------
## Header: sessionInfo() + ebrahim.gof version + run configuration
## ---------------------------------------------------------------------------
cat(strrep("=", 78), "\n", sep = "")
cat("run_all.R — Phase-A orchestrator\n")
cat("Started:", format(Sys.time()), "\n")
cat(strrep("=", 78), "\n", sep = "")

cat("\n---- sessionInfo() ----\n")
print(sessionInfo())

pkg_ver <- tryCatch(as.character(packageVersion("ebrahim.gof")),
                    error = function(e) paste0("<NOT INSTALLED: ", conditionMessage(e), ">"))
cat("\nebrahim.gof version:", pkg_ver, "\n")

reps_scale <- Sys.getenv("REPS_SCALE", "1")
cat("REPS_SCALE =", reps_scale, "\n")
cat("EK_NCORES  =", EK_NCORES, "\n")
cat("SIMDIR     =", SIMDIR, "\n")
cat("Rscript    =", RSCRIPT, "\n")

## ---------------------------------------------------------------------------
## Phase-A scripts IN ORDER (dependencies honored):
##   1. grid_null.R                -> size gate (MUST be first)
##   2. bench_compute.R            -> compute/scalability benchmark
##   3. grid_power_broad.R         -> headline power surfaces (consumes the size gate)
##   4. grid_asymmetry.R           -> EF-vs-HL asymmetry sweep
##   5. bench_slow_timing.R        -> slow-rival timing-only
##   6. grid_edge_loses.R          -> honest "where EDGE loses" grid
##   7. null_calibration_checks.R  -> null goodness / bootstrap / imhof (consumes the gate)
## ---------------------------------------------------------------------------
SCRIPTS <- c(
  "grid_null.R",
  "bench_compute.R",
  "grid_power_broad.R",
  "grid_asymmetry.R",
  "bench_slow_timing.R",
  "grid_edge_loses.R",
  "null_calibration_checks.R"
)

## CSVs each script is expected to (re)generate — used for the closing manifest.
EXPECTED_CSVS <- c(
  # grid_null.R
  "sim_null.csv", "sim_null_pvalues.csv",
  # bench_compute.R
  "bench_compute_time.csv", "bench_compute_time_summary.csv",
  "bench_compute_scaling_fits.csv", "bench_stukel_failure.csv",
  # grid_power_broad.R
  "sim_power_broad.csv", "sim_power_broad_pvalues.csv", "sim_g_sensitivity.csv",
  # grid_asymmetry.R
  "sim_asym_sweep.csv", "sim_asym_sweep_pvalues.csv", "sim_asym_delta.csv",
  # bench_slow_timing.R
  "bench_slow_timing.csv", "bench_slow_timing_summary.csv",
  # grid_edge_loses.R
  "sim_edge_loses.csv", "sim_edge_loses_pvalues.csv", "sim_edge_loses_Asearch.csv",
  # null_calibration_checks.R
  "null_ref_lambda.csv", "null_ks_table.csv", "null_bootstrap_agreement.csv",
  "imhof_satterthwaite.csv", "imhof_satterthwaite_rule.csv"
)

## ---------------------------------------------------------------------------
## Run each script in its own Rscript subprocess (exact standalone semantics),
## wrapped so a non-zero exit / crash is recorded as FAIL rather than aborting.
## ---------------------------------------------------------------------------
results <- data.frame(script = character(), status = character(),
                      elapsed_s = numeric(), detail = character(),
                      stringsAsFactors = FALSE)

for (scr in SCRIPTS) {
  path <- file.path(SIMDIR, scr)
  cat("\n", strrep("-", 78), "\n", sep = "")
  cat(">>> RUN ", scr, "   (", format(Sys.time()), ")\n", sep = "")
  cat(strrep("-", 78), "\n", sep = "")

  t0 <- proc.time()[["elapsed"]]
  status <- "FAIL"; detail <- NA_character_
  if (!file.exists(path)) {
    detail <- paste("script not found:", path)
    cat("    ERROR:", detail, "\n")
  } else {
    ## --vanilla child inherits our env (REPS_SCALE + any EK_* already set). Its stdout/stderr
    ## stream straight through to our console (and hence to the run log).
    rc <- tryCatch(
      system2(RSCRIPT, args = c("--vanilla", shQuote(path)),
              stdout = "", stderr = "", wait = TRUE),
      error = function(e) { detail <<- conditionMessage(e); 999L }
    )
    if (identical(rc, 0L)) {
      status <- "PASS"
    } else {
      status <- "FAIL"
      if (is.na(detail)) detail <- paste("Rscript exit status", rc)
    }
  }
  el <- proc.time()[["elapsed"]] - t0

  cat(sprintf("\n<<< %s: %s  (%.1fs)\n", scr, status, el))
  if (status == "FAIL") cat("    ", detail, "\n", sep = "")

  results <- rbind(results, data.frame(script = scr, status = status,
                                       elapsed_s = round(el, 1), detail = detail,
                                       stringsAsFactors = FALSE))
}

## ---------------------------------------------------------------------------
## Per-script PASS/FAIL + elapsed summary
## ---------------------------------------------------------------------------
cat("\n", strrep("=", 78), "\n", sep = "")
cat("PER-SCRIPT RESULTS\n")
cat(strrep("=", 78), "\n", sep = "")
for (i in seq_len(nrow(results))) {
  cat(sprintf("  %-28s %-4s  %8.1fs%s\n",
              results$script[i], results$status[i], results$elapsed_s[i],
              if (results$status[i] == "FAIL")
                paste0("   ", ifelse(is.na(results$detail[i]), "", results$detail[i])) else ""))
}

## ---------------------------------------------------------------------------
## CSV manifest with row counts
## ---------------------------------------------------------------------------
cat("\n", strrep("=", 78), "\n", sep = "")
cat("CSV MANIFEST (expected Phase-A outputs)\n")
cat(strrep("=", 78), "\n", sep = "")
row_count <- function(p) {
  tryCatch({
    d <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
    nrow(d)
  }, error = function(e) NA_integer_)
}
for (csv in EXPECTED_CSVS) {
  p <- file.path(SIMDIR, csv)
  if (file.exists(p)) {
    n <- row_count(p)
    cat(sprintf("  %-34s %s rows\n", csv,
                if (is.na(n)) "??" else format(n, big.mark = ",")))
  } else {
    cat(sprintf("  %-34s MISSING\n", csv))
  }
}

## ---------------------------------------------------------------------------
## STOP triggers
## ---------------------------------------------------------------------------
cat("\n", strrep("=", 78), "\n", sep = "")
cat("STOP TRIGGERS\n")
cat(strrep("=", 78), "\n", sep = "")
if (file.exists(STOP_LOG)) {
  cat("_STOP_TRIGGERS.log contents:\n")
  cat(readLines(STOP_LOG, warn = FALSE), sep = "\n")
  cat("\n")
} else {
  cat("No _STOP_TRIGGERS.log present — no STOP triggers fired.\n")
}

cat("\nrun_all.R DONE:", format(Sys.time()), "\n")

## Non-zero exit if any script FAILed (useful for CI / the orchestrator).
if (any(results$status == "FAIL")) quit(save = "no", status = 1L)
