## battery_review_zeroscan.R -- looks for p-values that are exactly 0 in null cells of every smoke-test output that exists
## (battery/dryrun, battery/_test, battery/_review/e, f4, f8), and reports the number of events in those replicates.
## A p-value of exactly 0 under a correctly specified model points to a numerical artefact such as finding D1.
## Writes battery/_review/review_zeroscan.log.
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
OUT <- edge_battery("_review")
suppressPackageStartupMessages(library(data.table))
sink(file.path(OUT, "review_zeroscan.log"), split = TRUE)
roots <- edge_battery(c("dryrun", "_test", "_review/e", "_review/f4", "_review/f8"))
files <- unlist(lapply(roots, function(r) list.files(r, pattern = "_pvalues\\.csv\\.gz$", recursive = TRUE, full.names = TRUE)))
META <- c("rep", "seed", "n", "events", "fit_ok", "glm_conv")
tot <- 0
for (f in files) {
  P <- fread(f)
  cn <- basename(f)
  null_like <- grepl("null|^id_null|sparse49_n|moderate1_n", cn)
  tests <- setdiff(names(P), META); tests <- tests[!grepl("^(flag\\.|lam[0-9]|chk\\.|proj$)", tests)]
  z <- sapply(tests, function(t) sum(P[[t]] %in% 0))
  if (any(z > 0)) {
    ev <- if ("events" %in% names(P)) paste(sort(unique(P$events[Reduce(`|`, lapply(tests[z > 0], function(t) P[[t]] %in% 0))])), collapse = ",") else "-"
    cat(sprintf("%-6s %-55s %s rows %d | p == 0: %s | events in those rows: %s\n", if (null_like) "NULL" else "alt", sub(paste0(SIMDIR, "/battery/"), "", f),
                "", nrow(P), paste(sprintf("%s %d", tests[z > 0], z[z > 0]), collapse = ", "), ev))
    if (null_like) tot <- tot + sum(z)
  }
}
cat(sprintf("\n%d files scanned; exact zeros in null-like files: %d\n", length(files), tot))
sink()
