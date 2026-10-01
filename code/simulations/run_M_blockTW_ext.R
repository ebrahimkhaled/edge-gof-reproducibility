## run_M_blockTW_ext.R -- block TW, Addendum 3: the ungrouped twin on frozen predictions, on block EXT's data sets.
## Contract: PREDECLARATION_blockTW_mechanism_ADDENDUM3.md (hashed before any replicate ran).
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
OUT <- edge_battery("TW_ext"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CELLS <- data.table(id = c(3L, 14:21), n = 1000L, fam = c("null", rep(c("exag", "sign"), each = 4)),
                    C = c(0, 4, 4, 4, 4, -4, -4, -4, -4), k = c(0L, 1L, 2L, 5L, 10L, 1L, 2L, 5L, 10L))
one <- function(rep, cell) {
  n <- cell$n; b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
  set.seed(if (cell$id <= 13) 7000000L + cell$id * 10000L + rep else 7600000L + cell$id * 10000L + rep)
  X <- matrix(stats::rnorm(n * 5), n, 5); eta_true <- as.numeric(cbind(1, X) %*% b0)
  if (cell$id <= 13) { y <- stats::rbinom(n, 1L, stats::plogis(eta_true)); eta0 <- eta_true } else {
    y <- stats::rbinom(n, 1L, stats::plogis(eta_true)); bad <- sample.int(n, cell$k)
    Xr <- X; Xr[bad, 2] <- cell$C * X[bad, 2]; eta0 <- as.numeric(cbind(1, Xr) %*% b0) }
  p <- stats::plogis(eta0); v <- p * (1 - p)
  Z <- stats::poly(p, 3); u <- colSums(Z * (y - p)); I <- crossprod(Z, v * Z)
  twin <- stats::pchisq(as.numeric(crossprod(u, solve(I, u))), 3, lower.tail = FALSE)
  zs <- sum((y - p) * (1 - 2 * p)) / sqrt(sum((1 - 2 * p)^2 * v))
  c(rep = rep, TWIN.ext = twin, SPZ.ext = 2 * stats::pnorm(-abs(zs)),
    EDGE.G10.ext = edge_external(y, p, G = 10)$p_value, EDGE.auto.ext = edge_external(y, p, G = "auto")$p_value)
}
cl <- makePSOCKcluster(16); on.exit(stopCluster(cl))
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl, c("one", "SIMDIR"))
invisible(clusterEvalQ(cl, { source(file.path(SIMDIR, "edge_external.R")); NULL }))
res <- list()
for (j in seq_len(nrow(CELLS))) {
  cell <- as.list(CELLS[j])
  M <- rbindlist(lapply(clusterApplyLB(cl, 1:1000, one, cell = cell), function(v) as.data.table(as.list(v))))
  ref <- fread(list.files(edge_battery("EXT"), sprintf("^c%02d_.*gz$", cell$id), full.names = TRUE))
  same <- isTRUE(all.equal(M$SPZ.ext, ref$spiegelhalter[match(M$rep, ref$rep)]))
  if (!same) stop("block TW_ext: the data sets of cell ", cell$id, " differ from block EXT")
  fwrite(M, file.path(OUT, sprintf("c%02d_%s_k%02d_pvalues.csv.gz", cell$id, cell$fam, cell$k)))
  res[[j]] <- data.table(cell = cell$id, fam = cell$fam, k = cell$k,
                         TWIN.ext = mean(M$TWIN.ext < .05), SPZ.ext = mean(M$SPZ.ext < .05),
                         EDGE.G10.ext = mean(M$EDGE.G10.ext < .05), EDGE.auto.ext = mean(M$EDGE.auto.ext < .05))
  print(res[[j]])
}
fwrite(rbindlist(res), file.path(OUT, "_summary.csv")); cat("done\n")
