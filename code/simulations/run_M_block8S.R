## run_M_block8S.R -- block 8S: the three rivals on the link and tail scenarios at n = 200.
## Contract: paper_EDGE/theory/PREDECLARATION_block8S_smalln_rivals.md (sha256 712e944d...), frozen before any scenario ran.
##
##   Rscript run_M_block8S.R [--workers 20] [--tests lecessie,proj,bagoft] [--reps N] [--out battery/8S]
##   Rscript run_M_block8S.R --summary
##
## Block 8R's machinery is reused unchanged -- the same replicate function, the same identity gate against the battery's
## stored p-values, the same BAGofT repair, the same summary and pairing rules. Only the plan differs: the six link and
## tail alternatives at n = 200 and their matched null, which no rival has been run on.
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
options(block8R.source_only = TRUE)
suppressMessages(source(file.path(SIMDIR, "run_M_block8R.R")))

S8_ALT  <- c("cauchit_n200", "t4_n200", "loglog_n200", "stk_short_n200", "stk_long_n200", "stk_asym_n200")
S8_NULL <- "null_link_n200"

## the plan of section 1: every test on all six alternatives and on the matched null
r8_plan <- function() {
  C <- battery_cells()
  alt <- C[C$cell %in% S8_ALT & C$role == "alternative", ]
  nul <- C[C$cell == S8_NULL & C$role == "null", ]
  if (nrow(alt) != 6L || nrow(nul) != 1L)
    stop("block 8S: expected 6 alternatives and 1 null, found ", nrow(alt), " and ", nrow(nul))
  if (any(alt$n != 200L) || nul$n != 200L) stop("block 8S: a cell is not at n = 200")
  if (any(alt$null_cell != S8_NULL)) stop("block 8S: an alternative is matched to a different null")
  mk <- function(rows, test) data.frame(test = test, block = rows$block, cell = rows$cell, role = rows$role,
                                        n = rows$n, null_block = rows$null_block, null_cell = rows$null_cell,
                                        R = R8_REPS[[test]], stringsAsFactors = FALSE)
  P <- do.call(rbind, lapply(R8_TESTS, function(t) mk(rbind(nul, alt), t)))
  if (nrow(P) != 21L) stop("block 8S: the plan does not match section 1 (expected 21 rows, got ", nrow(P), ")")
  rownames(P) <- NULL
  P
}

## the block writes to battery/8S unless told otherwise; a test run (--reps) may not touch it
s8_opts <- function(a = commandArgs(TRUE)) {
  opt <- r8_opts(a)
  if (is.null(opt$out) && is.null(opt$reps)) opt$out <- file.path(BAT, "8S")
  opt
}

if (!isTRUE(getOption("block8S.source_only"))) {
  r8_setup(s8_opts())
  quit(save = "no", status = r8_main())
}
