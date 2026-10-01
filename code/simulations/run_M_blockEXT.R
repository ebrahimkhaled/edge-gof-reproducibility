## run_M_blockEXT.R -- block EXT: the directed test in external mode against the usual external-validation tests.
## Contract: paper_EDGE/theory/PREDECLARATION_blockEXT_external_validation.md (sha256 754cf309...), frozen before
## any replicate ran.
##
## Cells 1-13 are the design of paper_deepgof/theory/external/run_external.R, with its seeds, so replicate r is the
## same data set in both studies; its Cox, Spiegelhalter, GiViTI, Hosmer-Lemeshow and Stukel code is copied here
## unchanged. Cells 14-21 add what that study does not have: corrupted records under a correct frozen model.
##
##   Rscript run_M_blockEXT.R [--workers 20] [--reps 1000] [--smoke]
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
args <- commandArgs(TRUE)
opt <- function(k, d) { i <- match(k, args); if (is.na(i)) d else args[i + 1] }
SMOKE <- "--smoke" %in% args
W <- as.integer(opt("--workers", 20)); REPS <- if (SMOKE) 4L else as.integer(opt("--reps", 1000))
OUT <- edge_battery(if (SMOKE) "EXT_smoke" else "EXT")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "_progress.log")
lg <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
  cat(s, "\n"); cat(s, "\n", file = LOG, append = TRUE) }

CELLS <- rbind(
  data.frame(fam = "null", C = 0, n = c(250L, 500L, 1000L), k = 0L),
  data.frame(fam = rep(c("large", "slope", "ushape", "thresh", "inter"), each = 2),
             C = c(.2, .4, .8, .6, .4, .8, 2, 4, .6, 1), n = 500L, k = 0L),
  data.frame(fam = rep(c("exag", "sign"), each = 4), C = rep(c(4, -4), each = 4), n = 1000L,
             k = rep(c(1L, 2L, 5L, 10L), 2)))
CELLS$id <- seq_len(nrow(CELLS))

one_rep <- function(cc, r) {
  n <- cc$n; b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
  set.seed(if (cc$id <= 13) 7000000L + cc$id * 10000L + r else 7600000L + cc$id * 10000L + r)
  X <- matrix(stats::rnorm(n * 5), n, 5)
  eta_true <- as.numeric(cbind(1, X) %*% b0)
  if (cc$id <= 13) {
    eta0 <- eta_true
    d <- switch(cc$fam, null = 0, large = cc$C, slope = (cc$C - 1) * eta0,
                ushape = cc$C * (X[, 1]^2 - 1), thresh = cc$C * pmax(X[, 2] - 1, 0),
                inter = cc$C * X[, 1] * X[, 3])
    y <- stats::rbinom(n, 1L, stats::plogis(eta0 + d))
  } else {
    ## the model is correct at the true covariates; k records carry x2 recorded as C * x2, and their
    ## prediction is the frozen model's at the recorded value
    y <- stats::rbinom(n, 1L, stats::plogis(eta_true))
    bad <- sample.int(n, cc$k)
    Xr <- X; Xr[bad, 2] <- cc$C * X[bad, 2]
    eta0 <- as.numeric(cbind(1, Xr) %*% b0)
  }
  p <- stats::plogis(eta0); v <- p * (1 - p)
  ## --- the directed test, external mode (edge_external.R = edge.gof(external = TRUE) of ebrahim.gof 2.9.0)
  e_auto <- edge_external(y, p, G = "auto")$p_value
  e_10   <- edge_external(y, p, G = 10)$p_value
  ## --- run_external.R's standard external-validation tests, unchanged
  g <- pmin(ceiling(rank(p, ties.method = "first") / (n / 10)), 10)
  O <- tapply(y, g, sum); E <- tapply(p, g, sum); ng <- tabulate(g, 10)
  phl <- stats::pchisq(sum((O - E)^2 / (E * (1 - E / ng))), 10, lower.tail = FALSE)
  dev0 <- -2 * sum(y * log(p) + (1 - y) * log(1 - p))
  za <- 0.5 * eta0^2 * (p >= .5); zb <- -0.5 * eta0^2 * (p < .5)
  gs <- suppressWarnings(stats::glm(y ~ 0 + za + zb, offset = eta0, family = stats::binomial()))
  pstuk <- stats::pchisq(dev0 - gs$deviance, 2, lower.tail = FALSE)
  gc <- suppressWarnings(stats::glm(y ~ eta0, family = stats::binomial()))
  pcox <- stats::pchisq(dev0 - gc$deviance, 2, lower.tail = FALSE)
  zsp <- sum((y - p) * (1 - 2 * p)) / sqrt(sum((1 - 2 * p)^2 * v))
  pspz <- 2 * stats::pnorm(-abs(zsp))
  pgiv <- tryCatch(suppressWarnings(givitiR::givitiCalibrationBelt(o = y, e = p, devel = "external")$p.value),
                   error = function(e) NA_real_)
  c(id = cc$id, rep = r, events = sum(y), EDGE_auto = e_auto, EDGE_G10 = e_10, HL = phl, stukel = pstuk,
    cox = pcox, spiegelhalter = pspz, giviti = pgiv)
}

lg("block EXT: %d cells x %d replicates, %d workers%s", nrow(CELLS), REPS, W, if (SMOKE) " [SMOKE]" else "")
cl <- parallel::makeCluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
parallel::clusterExport(cl, c("one_rep", "SIMDIR"))
invisible(parallel::clusterEvalQ(cl, { source(file.path(SIMDIR, "edge_external.R")); NULL }))
for (i in seq_len(nrow(CELLS))) {
  cc <- CELLS[i, ]
  f <- file.path(OUT, sprintf("c%02d_%s_C%s_n%04d_k%02d_pvalues.csv.gz", cc$id, cc$fam, cc$C, cc$n, cc$k))
  if (file.exists(f)) { lg("cell %d present, skipped", cc$id); next }
  t0 <- proc.time()[["elapsed"]]
  R <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(REPS), function(r, cc) one_rep(cc, r), cc = cc))
  R <- as.data.frame(R)
  tmp <- paste0(f, ".part.gz"); con <- gzfile(tmp, "w"); utils::write.csv(R, con, row.names = FALSE); close(con)
  file.rename(tmp, f)
  rej <- vapply(R[, 4:10], function(x) mean(!is.na(x) & x < 0.05), 0)
  lg("cell %2d %-6s C=%-4s n=%4d k=%2d | %s | %.1f min", cc$id, cc$fam, cc$C, cc$n, cc$k,
     paste(sprintf("%s %.3f", names(rej), rej), collapse = "  "), (proc.time()[["elapsed"]] - t0) / 60)
}
parallel::stopCluster(cl)
lg("block EXT: run finished")
