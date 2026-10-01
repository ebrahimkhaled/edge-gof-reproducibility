## analysis_verify2_probe.R -- fread on this machine: a gzipped file against a plain file that carries the ".csv.gz"
## name (fwrite compresses by extension, so a plain file has to be written as .csv and renamed). Only synthetic files.
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
setDTthreads(1L)
V2 <- edge_battery("_review/analysis/verify2")
dir.create(file.path(V2, "probe"), recursive = TRUE, showWarnings = FALSE)
tm <- function(label, expr) { t0 <- proc.time(); v <- force(expr); el <- proc.time() - t0
  cat(sprintf("%-46s %6.2f s\n", label, el[["elapsed"]])); invisible(v) }
set.seed(1)
X <- data.table(rep = 1:120, seed = 1:120 + 1e6)
for (k in 1:20) set(X, j = paste0("t", k), value = runif(120))
f_gz <- file.path(V2, "probe", "real_pvalues.csv.gz")
f_tmp <- file.path(V2, "probe", "flat_pvalues.csv"); f_flat <- file.path(V2, "probe", "flat_pvalues.csv.gz")
fwrite(X, f_gz)                                        # fwrite compresses because the name ends in .gz
fwrite(X, f_tmp); invisible(file.rename(f_tmp, f_flat))   # plain text under the same kind of name
b <- function(f) paste(sprintf("0x%02X", as.integer(readBin(f, "raw", 2))), collapse = " ")
cat("gz file", file.size(f_gz), "bytes, magic", b(f_gz), "| flat file", file.size(f_flat), "bytes, magic", b(f_flat), "\n")
r1 <- tm("fread the gzipped file", tryCatch(fread(f_gz), error = function(e) conditionMessage(e)))
r2 <- tm("fread the plain file named .csv.gz", tryCatch(fread(f_flat), error = function(e) conditionMessage(e)))
r3 <- tm("fread the plain file, select 3 columns", tryCatch(fread(f_flat, select = c("rep", "t1", "t2")), error = function(e) conditionMessage(e)))
for (nm in c("gz", "flat", "flat select")) {
  z <- list(r1, r2, r3)[[match(nm, c("gz", "flat", "flat select"))]]
  cat(" ", nm, ":", if (is.data.table(z)) sprintf("%d rows x %d cols", nrow(z), ncol(z)) else z, "\n")
}
cat("values identical:", identical(all.equal(as.data.frame(r1), as.data.frame(r2)), TRUE), "\n")
