## analysis_review_cells.R -- cell-table facts for the review of analyse_M_battery.R (names and roles only, no results):
## cell names used in two blocks, and the alternatives of blocks 2-4 that no family of E12.4 takes.
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
suppressPackageStartupMessages(library(data.table))
SIMDIR <- edge_path("code/simulations")
C <- fread(edge_battery("cells.csv"),
           colClasses = list(character = c("block", "cell", "role", "null_block", "null_cell", "generator", "family", "param",
                                           "link", "design", "seed_base")))
options(width = 200)
d <- C$cell[duplicated(C$cell)]
cat("cell names in more than one block:\n")
print(C[cell %in% d, .(block, cell, role, generator, n, B, seed_base, null_block, null_cell)])
cat("alternatives whose matched null carries one of those names:\n")
print(C[null_cell %in% d, .(block, cell, null_block, null_cell)])
src <- new.env()
sys.source(file.path(SIMDIR, "analyse_M_battery.R"), envir = src)
M <- src$an_membership(src$an_cell_table(edge_battery()), character(0))
A <- C[block %in% c("2", "3", "4") & role == "alternative"]
out <- A[!paste(block, cell) %in% paste(M$block, M$cell)]
cat(sprintf("\nblocks 2-4 alternatives: %d, in a family: %d, in none: %d\n", nrow(A), nrow(M), nrow(out)))
print(out[, .(block, cell, generator, family, param, link, design, n)], nrows = 200)
