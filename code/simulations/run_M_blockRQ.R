## run_M_blockRQ.R -- block RQ: robust quasi-deviance tests of Stukel's alternative under corrupted covariates.
## Contract: paper_EDGE/theory/PREDECLARATION_blockRQ_robust_tests.md (sha256 f7354147...), frozen before any
## replicate ran.
##
##   Rscript run_M_blockRQ.R [--workers 20] [--reps 1000]
##
## Each cell's replicates are written to disk as soon as the cell finishes, and a cell already on disk is
## skipped, so an interrupted run resumes where it stopped.
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
OUT <- edge_battery("RQ")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
a <- commandArgs(TRUE); OPT <- list(); i <- 1L
while (i <= length(a)) { OPT[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
W <- if (is.null(OPT$workers)) 20L else as.integer(OPT$workers)
REPS <- if (is.null(OPT$reps)) 1000L else as.integer(OPT$reps)
rq_lg <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
                      cat(s, "\n"); cat(s, "\n", file = file.path(OUT, "_progress.log"), append = TRUE) }

CELLS <- rbind(
  data.table(cell = "RQ_logit_clean_n1000", truth = "logit", corruption = "clean", k = 0L),
  data.table(cell = sprintf("RQ_logit_C1_k%d_n1000", c(1L, 2L, 5L, 10L)), truth = "logit", corruption = "C1",
             k = c(1L, 2L, 5L, 10L)),
  data.table(cell = sprintf("RQ_logit_C1b_k%d_n1000", c(1L, 2L, 5L, 10L)), truth = "logit", corruption = "C1b",
             k = c(1L, 2L, 5L, 10L)),
  data.table(cell = "RQ_cloglog_clean_n1000", truth = "cloglog", corruption = "clean", k = 0L))
CELLS[, index := .I]

rq_one <- function(rep, cell) {
  suppressPackageStartupMessages({ library(robustbase); library(ebrahim.gof) })
  set.seed(970000000 + 100000 * cell$index + rep)
  g <- c1b_gen(cell$truth, 1000L, cell$corruption, cell$k / 1000)
  dat <- g$d
  qd <- function(wx) tryCatch({
    ctl <- glmrobMqle.control(tcc = 1.345, maxit = 100)
    f0 <- glmrob(y ~ x + d, family = binomial, data = dat, method = "Mqle", weights.on.x = wx, control = ctl)
    eta <- f0$linear.predictors
    d1 <- dat; d1$z1 <- 0.5 * eta^2 * (eta >= 0); d1$z2 <- -0.5 * eta^2 * (eta < 0)
    f1 <- glmrob(y ~ x + d + z1 + z2, family = binomial, data = d1, method = "Mqle", weights.on.x = wx,
                 control = ctl)
    tab <- anova(f0, f1, test = "QD")
    as.numeric(tab[2, ncol(tab)])
  }, error = function(e) NA_real_)
  fit <- stats::glm(y ~ x + d, family = stats::binomial(), data = dat)
  e_def <- tryCatch(edge.gof(fit, G = "auto", basis = "poly3")$p_value, error = function(e) NA_real_)
  e_g10 <- tryCatch(edge.gof(fit, G = 10, basis = "poly3")$p_value, error = function(e) NA_real_)
  c(rep = rep, k_corrupt = length(g$corrupt), RQD.Huber = qd("none"), RQD.Mallows = qd("hat"),
    EDGE.default = e_def, EDGE.G10 = e_g10)
}

source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
source(file.path(SIMDIR, "_blockC1b_contam.R"))
cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl))
invisible(clusterEvalQ(cl, {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  SIMDIR <- edge_path("code/simulations")
  source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
  source(file.path(SIMDIR, "_blockC1b_contam.R"))
  NULL
}))
clusterExport(cl, "rq_one")
rq_lg("block RQ: %d cells, %d replicates, %d workers (declaration sha256 f7354147)", nrow(CELLS), REPS, W)

for (j in seq_len(nrow(CELLS))) {
  cell <- as.list(CELLS[j])
  f <- file.path(OUT, sprintf("%s_pvalues.csv.gz", cell$cell))
  if (file.exists(f)) { rq_lg("%s: present, skipped", cell$cell); next }
  t0 <- Sys.time()
  M <- as.data.table(do.call(rbind, clusterApplyLB(cl, seq_len(REPS), rq_one, cell = cell)))
  fwrite(M, f, compress = "gzip")
  rt <- vapply(c("RQD.Huber", "RQD.Mallows", "EDGE.default", "EDGE.G10"),
               function(t) mean(!is.na(M[[t]]) & M[[t]] <= 0.05), 0)
  rq_lg("%-24s k=%2d  %s  missing %s  %.1f min", cell$cell, cell$k,
     paste(sprintf("%s %.3f", names(rt), rt), collapse = "  "),
     paste(vapply(c("RQD.Huber", "RQD.Mallows"), function(t) sum(is.na(M[[t]])), 0), collapse = "/"),
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
rq_lg("block RQ: run finished")
