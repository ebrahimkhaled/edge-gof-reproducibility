## run_M_all.R -- runs the whole EDGE restructure battery in the order of PREDECLARATION_restructure_battery.md E4, one
## Rscript call of run_M_battery.R per step, and stops at the first step that fails:
##   0; 1b (the matched nulls of block 2); 2; 1a; the rest of 1b; 3 (the n = 1000 slice); the rest of 3; 4; 5; 6; 7;
##   then the summary of every block (all finished cells of the block).
## The driver skips finished cells, so the launcher can be started again after an interruption (with --from to skip steps).
##
##   Rscript run_M_all.R              the real run (block 0 times 12, 16 and 20 workers; later steps use its choice)
##   Rscript run_M_all.R --plan       print the steps and write battery/launch_plan.csv; nothing is run
##   Rscript run_M_all.R --from 6     start at step 6
## Log: battery/launch.log (each driver call also writes battery/<block>/_progress.log).
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
source(file.path(SIMDIR, "_battery_cells.R"))
a <- commandArgs(trailingOnly = TRUE)
plan_only <- "--plan" %in% a
from <- if ("--from" %in% a) suppressWarnings(as.integer(a[which(a == "--from") + 1L])) else 1L
if (is.na(from) || from < 1L) stop("--from needs a step number")
BAT <- edge_battery()

P <- launch_plan(battery_cells())
fwrite(P, edge_battery("launch_plan.csv"))
for (i in seq_len(nrow(P)))
  cat(sprintf("%2d  %-7s block %-2s  %s%s\n", P$step[i], P$action[i], P$block[i], P$note[i],
              if (nzchar(P$cells[i])) sprintf(" (%d cells)", length(strsplit(P$cells[i], ",")[[1]])) else ""))
if (plan_only) quit(save = "no", status = 0L)

RS  <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
LOG <- edge_battery("launch.log")
say <- function(...) {
  line <- sprintf("%s  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...))
  cat(line, "\n", sep = ""); cat(line, "\n", file = LOG, append = TRUE, sep = "")
}
say("launch from step %d of %d", from, nrow(P))
for (i in which(P$step >= from)) {
  args <- c(shQuote(file.path(SIMDIR, "run_M_battery.R")), "--block", P$block[i])
  if (P$action[i] == "summary") args <- c(args, "--summary")
  if (nzchar(P$cells[i])) args <- c(args, "--cells", P$cells[i])
  t0 <- Sys.time()
  say("step %d: %s block %s%s", P$step[i], P$action[i], P$block[i], if (nzchar(P$cells[i])) " (listed cells only)" else "")
  st <- system2(RS, args)
  say("step %d: exit %d after %.1f min", P$step[i], st, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  if (st != 0) {
    say("stopped at step %d; find the cause, then restart with: Rscript run_M_all.R --from %d", P$step[i], P$step[i])
    quit(save = "no", status = 1L)
  }
  if (P$block[i] == "0") {                                            # E0.8: the worker choice must come from 12, 16 and 20
    wf <- edge_battery("0", "workers.csv"); wc <- edge_battery("0", "workers_choice.txt")
    if (!file.exists(wc) || !file.exists(wf) || !all(c(12L, 16L, 20L) %in% fread(wf)$workers)) {
      say("block 0 did not record the worker timing at 12, 16 and 20 workers (E0.8)"); quit(save = "no", status = 1L)
    }
    say("workers chosen by block 0: %s", readLines(wc, warn = FALSE)[1])
  }
}
say("launch finished")
