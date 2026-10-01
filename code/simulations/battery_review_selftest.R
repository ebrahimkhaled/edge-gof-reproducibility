## battery_review_selftest.R -- re-runs battery_selftest.R unchanged except that its three output files are redirected
## to battery/_review/ (selftest_rerun.log, seed_overlaps_rerun.csv, old_id_check_rerun.csv), so the builder's outputs
## in battery/ are not overwritten. The substitutions are counted and must each hit exactly once.
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
dir.create(edge_battery("_review"), recursive = TRUE, showWarnings = FALSE)
txt <- paste(readLines(file.path(SIMDIR, "battery_selftest.R"), warn = FALSE), collapse = "\n")
subs <- list(
  c('edge_battery("selftest.log")', 'edge_battery("_review", "selftest_rerun.log")'),
  c('edge_battery("seed_overlaps.csv")', 'edge_battery("_review", "seed_overlaps_rerun.csv")'),
  c('edge_battery("old_id_check.csv")', 'edge_battery("_review", "old_id_check_rerun.csv")'))
for (s in subs) {
  hits <- lengths(regmatches(txt, gregexpr(s[1], txt, fixed = TRUE)))
  if (hits != 1L) stop("substitution did not hit exactly once: ", s[1], " (", hits, ")")
  txt <- sub(s[1], s[2], txt, fixed = TRUE)
}
stopifnot(!grepl('"battery", "selftest.log"', txt, fixed = TRUE))
eval(parse(text = txt), envir = globalenv())
