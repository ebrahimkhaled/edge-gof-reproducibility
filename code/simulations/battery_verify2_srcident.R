## battery_verify2_srcident.R -- is the installed ebrahim.gof the dev tree's R code? Compares the deparsed functions of the
## installed namespace with the same files sourced from the dev tree (branch fix-stukel-joint-sym).
## Output: battery/_review/v2/srcident.log
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
PKGDIR <- Sys.getenv("EBRAHIM_GOF_SRC")
OUT <- edge_battery("_review", "v2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "srcident.log"), split = TRUE)
cat("installed ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "|", packageDescription("ebrahim.gof")$Packaged, "\n")
cat("dev tree HEAD", system2("git", c("-C", shQuote(PKGDIR), "rev-parse", "--short", "HEAD"), stdout = TRUE), "\n")
cat("dev tree R/ changes not committed:", length(system2("git", c("-C", shQuote(PKGDIR), "status", "--porcelain", "R"), stdout = TRUE)), "\n")
ns <- asNamespace("ebrahim.gof")
env <- new.env()
for (f in c("def_gof.R", "run_all_gof.R", "ebrahim_farrington_test.R", "edge_gof.R", "def_ensemble_gof.R"))
  sys.source(file.path(PKGDIR, "R", f), envir = env, keep.source = FALSE)
fns <- c("def.gof", ".def_basis", ".def_auto_G", ".def_warn_few_events", ".def_warn_degenerate", ".def_warn_no_information",
         ".def_pvalue", "gof_stukel", ".gof_context", ".gof_groups_ef", ".gof_hl_stat", "gof_hl", "gof_hlw", "gof_ph_test", "gof_ef",
         "gof_def", "gof_tsiatis", "gof_xie", "gof_pr", ".gof_kmeans", ".gof_ginv", "ef.gof", "edge.gof")
for (f in fns) {
  a <- if (exists(f, envir = ns, inherits = FALSE)) get(f, envir = ns) else NULL
  b <- if (exists(f, envir = env, inherits = FALSE)) get(f, envir = env) else NULL
  same <- !is.null(a) && !is.null(b) && identical(deparse(a), deparse(b))
  cat(sprintf("%-26s installed %-5s dev %-5s identical %s\n", f, !is.null(a), !is.null(b), same))
}
sink()
