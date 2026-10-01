## run_M_block9R2n.R -- the matched null of block 9R2's six battery scenarios, run through the same variants.
##
## Section 4.3 of the paper reads power size-adjusted, using each test's critical value from its own matched null
## scenario. Block 9R2 ran the six link and tail alternatives but not the null they are matched to, so the power rows of
## the EDGE-FR table could only be read raw. This adds `null_link_n1000`, the null of all six, for the same five
## variants and the same 1000 replicates, so that the paper's own reading rule can be applied to its own new results.
##
##   Rscript run_M_block9R2n.R [--workers 20] [--reps N] [--out battery/9R2n]
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
suppressPackageStartupMessages({ library(parallel); library(data.table) })
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))
source(file.path(SIMDIR, "_block9_contam.R"))
source(file.path(SIMDIR, "_block9R2_hybrid.R"))
BAT <- edge_battery()

N9_CELLS <- c("null_link_n1000")

n9_one <- function(rep, ce, cell) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  set.seed(ce$seed_base + rep)
  dat <- bt_data(cell)
  n <- nrow(dat$d); ev <- sum(dat$d$y)
  fq <- if (min(ev, n - ev) == 0) NULL else bt_fit(dat)
  c(rep = rep, n = n, events = ev, fr2_pvalues(fq))
}
n9_task <- function(rep, ce, cell) n9_one(rep, ce, cell)

OPT <- list(); a <- commandArgs(TRUE); i <- 1L
while (i <= length(a)) { OPT[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
W <- if (is.null(OPT$workers)) 20L else as.integer(OPT$workers)
REPS <- if (is.null(OPT$reps)) 1000L else as.integer(OPT$reps)
OUT <- normalizePath(if (!is.null(OPT$out)) OPT$out else edge_battery("9R2n"), winslash = "/", mustWork = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lg <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
                      cat(s, "\n"); cat(s, "\n", file = file.path(OUT, "_progress.log"), append = TRUE) }

C <- battery_cells()
cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl))
invisible(clusterCall(cl, function(simdir) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  SIMDIR <<- simdir
  source(file.path(simdir, "_battery_tests.R")); source(file.path(simdir, "_battery_cells.R"))
  source(file.path(simdir, "_block9_contam.R")); source(file.path(simdir, "_block9R2_hybrid.R"))
  RNGkind("L'Ecuyer-CMRG"); NULL
}, SIMDIR))
clusterExport(cl, c("n9_one", "n9_task"))

for (nm in N9_CELLS) {
  row <- C[C$cell == nm, ]
  if (nrow(row) != 1L) stop("expected one row for ", nm, ", found ", nrow(row))
  cell <- as.list(row); ce <- list(seed_base = row$seed_base)
  f <- file.path(OUT, sprintf("%s_variants.csv.gz", nm))
  if (file.exists(f)) { lg("%s: present, skipped", nm); next }
  t0 <- Sys.time()
  M <- as.data.table(do.call(rbind, clusterApplyLB(cl, seq_len(REPS), n9_task, ce = ce, cell = cell)))
  ## identity gate: V0 must be the stored EDGE-poly3 unit p-value at the rule G on the same data sets
  st <- fread(edge_battery(row$block, sprintf("%s_pvalues.csv.gz", nm)), select = c("rep", "EDGE.poly3.u.Grule"))
  m <- st$EDGE.poly3.u.Grule[match(M$rep, st$rep)]
  bad <- sum(!(abs(M$V0 - m) <= 1e-8 | (!is.finite(M$V0) & !is.finite(m))), na.rm = TRUE)
  M[, id_ok := as.integer(abs(V0 - m) <= 1e-8 | (!is.finite(V0) & !is.finite(m)))]
  fwrite(M, f, compress = "gzip")
  lg("%-18s n=%-5d R=%d  identity failures %d  size V0 %.3f HYB05 %.3f HYB10 %.3f HYB20 %.3f G10 %.3f  %.1f min",
     nm, row$n, REPS, bad,
     mean(M$V0 <= 0.05, na.rm = TRUE), mean(M$HYB05 <= 0.05, na.rm = TRUE), mean(M$HYB10 <= 0.05, na.rm = TRUE),
     mean(M$HYB20 <= 0.05, na.rm = TRUE), mean(M$G10 <= 0.05, na.rm = TRUE),
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
lg("block 9R2n: run finished")
