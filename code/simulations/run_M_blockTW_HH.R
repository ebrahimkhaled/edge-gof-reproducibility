## run_M_blockTW_HH.R -- block TW, Addendum 2: Hosmer and Hjort's (2002) weighted grouped tests on TW's own data sets.
## Contract: PREDECLARATION_blockTW_mechanism.md (sha256 80b5a83f) and ADDENDUM2 (sha256 de5ebb60), hashed before any
## Hosmer-Hjort replicate ran.
##
## The cells, the generator and the seeds are taken from run_M_blockTW.R itself (its definitions are evaluated, not
## copied), and every replicate must reproduce TW's EDGE.G10 p-value exactly, which shows the data set is the same.
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
PKG <- Sys.getenv("EBRAHIM_GOF_SRC")
suppressPackageStartupMessages({ library(parallel); library(data.table) })
src <- readLines(file.path(SIMDIR, "run_M_blockTW.R"))
stop_at <- grep("^## ---- arms", src)[1] - 1L                  # cells, BETA_B and tw_gen, nothing that runs
DEF <- src[grep("^A <- function", src)[1]:stop_at]
eval(parse(text = DEF))
tw_edge <- function(fit, G) tryCatch(suppressWarnings(edge.gof(fit, G = G, basis = "poly3")$p_value), error = function(e) NA_real_)
OUT <- edge_battery("TW_HH"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
REPS <- 1000L

hh_one <- function(rep, cell) {
  set.seed(980000000 + 10000 * cell$index + rep)
  g <- tw_gen(cell)
  fq <- bt_fit(g)
  if (is.null(fq)) return(c(rep = rep))
  h <- tryCatch(hh_test(fq$fit, g = 10), error = function(e) c(HL1 = NA, X2_1 = NA, HLnp = NA, X2np = NA))
  c(rep = rep, EDGE.G10.check = tw_edge(fq$fit, 10), HH.HL1 = h[["HL1"]], HH.X2_1 = h[["X2_1"]],
    HH.HLnp = h[["HLnp"]], HH.X2np = h[["X2np"]])
}

source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
source(file.path(SIMDIR, "_blockC1b_contam.R"))
cl <- makePSOCKcluster(16)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl))
clusterExport(cl, c("hh_one", "tw_gen", "tw_edge", "BETA_B", "PKG", "SIMDIR"))
invisible(clusterEvalQ(cl, {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  suppressMessages(library(ebrahim.gof))  # archive: was load_all() of the author's development tree; install ebrahim.gof (>= 2.9.0) from CRAN instead
  source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
  source(file.path(SIMDIR, "_blockC1b_contam.R")); source(file.path(SIMDIR, "_hosmer_hjort.R"))
  NULL
}))
for (j in seq_len(nrow(CELLS))) {
  cell <- as.list(CELLS[j])
  f <- file.path(OUT, sprintf("%s_pvalues.csv.gz", cell$cell))
  if (file.exists(f)) next
  M <- rbindlist(lapply(clusterApplyLB(cl, seq_len(REPS), hh_one, cell = cell), function(v) as.data.table(as.list(v))), fill = TRUE)
  ref <- fread(edge_battery("TW", sprintf("%s_pvalues.csv.gz", cell$cell)))
  same <- isTRUE(all.equal(M$EDGE.G10.check, ref$EDGE.G10[match(M$rep, ref$rep)]))
  if (!same) stop("block TW_HH: the data sets of ", cell$cell, " differ from block TW")
  fwrite(M, f, compress = "gzip")
  cat(sprintf("%-34s same data: %s | %s\n", cell$cell, same,
              paste(sprintf("%s %.3f", c("HH.HL1", "HH.X2_1", "HH.HLnp", "HH.X2np"),
                            vapply(c("HH.HL1", "HH.X2_1", "HH.HLnp", "HH.X2np"), function(t) mean(!is.na(M[[t]]) & M[[t]] <= 0.05), 0)),
                    collapse = "  ")))
}
cat("block TW_HH: run finished\n")
