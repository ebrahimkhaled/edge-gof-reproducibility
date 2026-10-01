## battery_review_rivals_src.R -- review of block 8: reads the BAGofT source and the ebrahim.gof EDGE signature (no computation)
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
OUTD <- edge_battery("_review/rivals")
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUTD, "review_src.txt"))
ns <- asNamespace("BAGofT")
cat("==== BAGofT_multi\n"); print(ns$BAGofT_multi)
cat("==== BAGofT_sin\n"); print(ns$BAGofT_sin)
cat("==== dcPre\n"); print(ns$dcPre)
cat("==== functions of BAGofT that mention seeds or parallel back ends\n")
pat <- c("set.seed", "mclapply", "parLapply", "foreach", "makeCluster", "RNGkind", "detectCores")
for (f in ls(ns)) {
  b <- paste(deparse(get(f, ns)), collapse = " ")
  hit <- pat[vapply(pat, function(p) grepl(p, b, fixed = TRUE), logical(1))]
  if (length(hit)) cat(f, ":", paste(hit, collapse = ", "), "\n")
}
cat("==== imports of BAGofT\n"); print(packageDescription("BAGofT")[c("Version", "Imports", "Depends")])
cat("==== randomForest version\n"); print(packageVersion("randomForest"))
cat("==== ebrahim.gof exports matching gof\n"); print(grep("def|edge|stukel", getNamespaceExports("ebrahim.gof"), value = TRUE))
cat("==== def.gof args\n"); print(args(ebrahim.gof::def.gof))
sink()
