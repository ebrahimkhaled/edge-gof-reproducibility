## battery_dryrun.R -- a dry run of every block of run_M_battery.R before the launch: block 0 in full at B = 20, then a few
## cells of each other block at B = 10, all with 8 workers, into battery/dryrun (never the result folders). A previous dry run
## is first moved to battery/_review/dryrun_prev_<time>, because the driver skips finished cells. Then the launcher's plan
## (run_M_all.R --plan) and the timing sample for the CPU-hour projection. Needs ebrahim.gof >= 2.8.0 installed.
## Writes battery/dryrun/dryrun.log and battery/dryrun/dryrun_check.csv; exits non-zero if any step fails.
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
ROOT_REL <- "battery/dryrun"
ROOT <- file.path(SIMDIR, ROOT_REL)
v <- tryCatch(suppressWarnings(utils::packageDescription("ebrahim.gof")$Version), error = function(e) NA_character_)
if (is.na(v) || utils::compareVersion(v, "2.8.0") < 0) {
  cat(sprintf("pending: package (ebrahim.gof %s installed, 2.8.0 or later needed)\n", v))
  quit(save = "no", status = 2L)
}
if (dir.exists(ROOT)) {
  prev <- edge_battery("_review", paste0("dryrun_prev_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(dirname(prev), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(ROOT, prev)) stop("could not move the previous dry run to ", prev)
}
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(ROOT, "dryrun.log")
RS <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
say <- function(...) { line <- sprintf("%s  %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)); cat(line, "\n"); cat(line, "\n", file = LOG, append = TRUE) }
## system2() truncates a file given as stdout, so every step writes its own log and dryrun.log only points to it
run <- function(args, tag, script = "run_M_battery.R") {
  t0 <- Sys.time()
  steplog <- file.path(ROOT, sprintf("step_%s.log", tag))
  st <- system2(RS, c(shQuote(file.path(SIMDIR, script)), args), stdout = steplog, stderr = steplog)
  if (st != 0) { tl <- utils::tail(readLines(steplog, warn = FALSE), 15); cat(tl, sep = "\n", file = LOG, append = TRUE) }
  list(status = st, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
say("dry run, ebrahim.gof %s, root %s", v, ROOT)

## cells per block, chosen to reach every generator type (external with a G list, EDGE-only, real data with a logit and a
## cloglog fit, the projection pairing, EDGE-ao, the new design nulls) and, where possible, an alternative with its matched
## null so that the size-adjusted summary is exercised. Block 1b runs in two calls, to check that the summary keeps both parts.
PICK <- list("1a" = list(c("null_quad_n200", "null_link_n1000")),
             "1b" = list("null_base_n610", c("null_quad_n100", "null_crossover_n1000", "null_census_int_binbin_n1000")),
             "2"  = list(c("loglog_base_n610", "cauchit_auc_n460")),
             "3"  = list(c("link_plateau_upper_n200", "stk_long_n1000", "census_int_binbin_n1000")),
             "4"  = list(c("omit_2cov_n500", "sparse49_cloglog_n200", "proj_link_plateau_lower_n500", "proj_null_link_n500", "crossover_n1000")),
             "5"  = list(c("runK2_cloglog_n500", "runH_temp_s0.85_n500", "runG_osc4_n500")),
             "6"  = list(c("glow_null_n500", "real_vaso_log", "real_beetle_logit", "real_beetle_cloglog", "real_diabetes_val")),
             "7"  = list(c("runJ_null_n2000", "runJ_cloglog_n2000")))

suppressPackageStartupMessages(library(data.table))
CHK <- list()
addchk <- function(block, cell, exit, rows, rep_errors, ok, secs, what = "")
  CHK[[length(CHK) + 1]] <<- data.frame(block = block, cell = cell, check = what, exit = exit, rows = rows, rep_errors = rep_errors, ok = ok, secs = secs)

r0 <- run(c("--block", "0", "--B", "20", "--workers", "8", "--worker-grid", "4,8", "--root", ROOT_REL), tag = "block_0")
say("block 0 (B = 20): exit %d in %.0f s", r0$status, r0$secs)
pk <- file.path(ROOT, "0", "package_check.csv")
PK <- if (file.exists(pk)) fread(pk) else NULL
addchk("0", "all", r0$status, NA, NA, r0$status == 0 && file.exists(file.path(ROOT, "0", "_PASSED")), r0$secs, "block 0 passed")
addchk("0", "constructed", r0$status, NA, NA, !is.null(PK) && sum(grepl("^constructed, all y", PK$cell)) == 6 &&
         all(PK$ok[grepl("^constructed, all y", PK$cell)]), r0$secs, "zero-event and all-event samples: every harness test NA")

for (b in names(PICK)) {
  for (k in seq_along(PICK[[b]])) {
    cells <- PICK[[b]][[k]]
    r <- run(c("--block", b, "--cells", paste(cells, collapse = ","), "--B", "10", "--workers", "8", "--root", ROOT_REL),
             tag = sprintf("block_%s_%d", b, k))
    say("block %s (%s): exit %d in %.0f s", b, paste(cells, collapse = ", "), r$status, r$secs)
    for (cn in cells) {
      f <- file.path(ROOT, b, paste0(cn, "_pvalues.csv.gz"))
      P <- if (file.exists(f)) fread(f) else NULL
      want <- if (grepl("^real_", cn)) 11L else 10L                  # a real fit: the stored order and 10 random orders
      addchk(b, cn, r$status, if (is.null(P)) 0L else nrow(P), if (is.null(P)) NA else sum(P$flag.rep_error %in% 1),
             r$status == 0 && !is.null(P) && nrow(P) == want && all(P$flag.rep_error %in% 0), r$secs, "output rows, no replicate error")
    }
  }
  sf <- file.path(ROOT, b, "_summary.csv")
  S <- if (file.exists(sf)) fread(sf) else NULL
  addchk(b, "summary", NA, NA, NA, !is.null(S) && all(unlist(PICK[[b]]) %in% S$cell), NA, "the summary holds every cell run in the block")
}

## specific checks of the fixes
S4 <- fread(file.path(ROOT, "4", "_summary.csv"))
addchk("4", "crossover_n1000", NA, NA, NA, any(is.finite(S4[cell == "crossover_n1000" & test == "EDGE.poly3.u.G10", size_adj_power])),
       NA, "crossover has a matched null: size-adjusted power computed")
pl <- readLines(file.path(ROOT, "4", "_progress.log"), warn = FALSE)
pd <- suppressWarnings(as.numeric(sub(".*stored p_edge3 max\\|diff\\| ([^ ]+).*", "\\1", grep("stored p_edge3", pl, value = TRUE))))
addchk("4", "projection pairing", NA, NA, NA, length(pd) == 3 && all(is.finite(pd) & pd <= 1e-6), NA,
       sprintf("stored p_edge3 reproduced in the 3 paired cells (max %s)", paste(format(pd, digits = 2), collapse = ", ")))
bc <- fread(file.path(ROOT, "6", "real_beetle_cloglog_pvalues.csv.gz"))
addchk("6", "real_beetle_cloglog", NA, NA, NA, all(is.na(unlist(bc[, grep("^(EDGE|Stk|Cubic)", names(bc)), with = FALSE]))) &&
         all(is.finite(bc$HL.G10)) && all(is.finite(bc$HLF.G10)), NA, "cloglog refit: EDGE and Stukel not run, HL and HL_F returned")
dv <- fread(file.path(ROOT, "6", "real_diabetes_val_pvalues.csv.gz"))
addchk("6", "real_diabetes_val", NA, NA, NA, all(c("EDGE.poly4.G2000", "stat.HL.ext.G50") %in% names(dv)) &&
         all(is.finite(dv$EDGE.poly4.G2000)) && all(is.finite(dv$stat.HL.ext.G50)), NA, "validation half: grouped tests at G = 10 ... 2000")
rof <- file.path(ROOT, "6", "_row_orders.csv")
RO <- if (file.exists(rof)) fread(rof) else NULL
bl <- fread(file.path(ROOT, "6", "real_beetle_logit_pvalues.csv.gz"))
addchk("6", "row orders", NA, NA, NA, !is.null(RO) && setequal(unique(RO$cell), grep("^real_", PICK[["6"]][[1]], value = TRUE)) &&
         all(RO$orders == 10L) && identical(as.integer(bl$rep), 0:10) && is.na(bl$seed[1]) && all(bl$seed[-1] == 20260914 + 1:10) &&
         bl$info.distinct_risks[1] == 8 && isTRUE(abs(RO[cell == "real_beetle_logit" & test == "HL.G10", p_stored] - bl$HL.G10[1]) < 1e-12),
       NA, "real fits: stored order + 10 random orders (seed 20260914 + k); _row_orders.csv per test (E8.5)")

rp <- run(c("--plan"), tag = "launch_plan", script = "run_M_all.R")
say("launcher plan: exit %d", rp$status)
addchk("launcher", "plan", rp$status, NA, NA, rp$status == 0 && file.exists(edge_battery("launch_plan.csv")), rp$secs,
       "run_M_all.R --plan")

re <- run(c("--estimate", "--block", "all", "--workers", "20"), tag = "estimate")
say("timing sample and projection: exit %d in %.0f s (battery/_estimate_blocks.csv)", re$status, re$secs)
addchk("estimate", "all", re$status, NA, NA, re$status == 0, re$secs, "timing sample")
K <- do.call(rbind, CHK)
fwrite(K, file.path(ROOT, "dryrun_check.csv"))
print(K, row.names = FALSE)
say("dry run %s", if (all(K$ok)) "PASSED" else "FAILED")
quit(save = "no", status = if (all(K$ok)) 0L else 1L)
