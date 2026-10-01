## run_M_blockTW.R -- block TW: which property protects a test from corrupted records?
## Contract: paper_EDGE/theory/PREDECLARATION_blockTW_mechanism.md (sha256 80b5a83f) and its Addendum 1 (PRG2 arms),
## both hashed before any replicate of the full run.
##
##   Rscript run_M_blockTW.R [--workers 16] [--reps 1000] [--smoke]
##
## Every arm runs on the same data set in each replicate, so all comparisons are paired. Each cell is written to
## disk as soon as it finishes and a cell already on disk is skipped.
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
a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% a; a <- setdiff(a, "--smoke"); OPT <- list(); i <- 1L
while (i <= length(a)) { OPT[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
W <- if (is.null(OPT$workers)) 16L else as.integer(OPT$workers)
REPS <- if (SMOKE) 3L else if (is.null(OPT$reps)) 1000L else as.integer(OPT$reps)
OUT <- edge_battery(if (SMOKE) "TW_smoke" else "TW")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
tw_lg <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
                         cat(s, "\n"); cat(s, "\n", file = file.path(OUT, "_progress.log"), append = TRUE) }

A <- function(n, truth, corr, k) data.table(part = "A", n = n, truth = truth, corruption = corr, k = k)
B <- function(truth, corr, k) data.table(part = "B", n = 2000L, truth = truth, corruption = corr, k = k)
CELLS <- rbind(
  A(1000L, "logit", "clean", 0L), A(1000L, "logit", "C1", 1L), A(1000L, "logit", "C1", 2L), A(1000L, "logit", "C1", 5L),
  A(1000L, "logit", "C1", 10L), A(1000L, "logit", "C1b", 1L), A(1000L, "logit", "C1b", 2L), A(1000L, "logit", "C1b", 5L),
  A(1000L, "logit", "C1b", 10L), A(1000L, "cloglog", "clean", 0L), A(1000L, "logit", "C2", 10L),
  A(1000L, "tail0.6", "clean", 0L),
  A(5000L, "logit", "clean", 0L), A(5000L, "logit", "C1", 25L), A(5000L, "logit", "C1", 50L),
  A(5000L, "logit", "C1b", 10L), A(5000L, "logit", "C1b", 25L), A(5000L, "probit", "clean", 0L),
  A(5000L, "tail0.8", "clean", 0L),
  B("logit", "clean", 0L), B("cloglog", "clean", 0L), B("inter", "clean", 0L), B("logit", "C1", 10L),
  B("logit", "C1b", 10L), B("logit", "C1", 25L))
CELLS[, index := .I]
CELLS[, cell := sprintf("TW%s_%s_%s_k%02d_n%d", part, truth, corruption, k, n)]

## ---- data --------------------------------------------------------------------------------------------------------
BETA_B <- c(0.5, -0.4, 0.3, 0.3, -0.2, 0.2, 0.1, -0.1, 0.4, -0.3, 0.2, 0.2)
tw_gen <- function(cell) {
  n <- cell$n
  if (cell$part == "A") {
    if (cell$truth %in% c("logit", "cloglog", "probit")) {
      if (cell$corruption == "C2") return(b9_gen("logit", n, "C2", cell$k / n))
      return(c1b_gen(cell$truth, n, cell$corruption, cell$k / n))
    }
    ## localised tail misfit: the logistic model is right below the 98th percentile of true risk, and above it the
    ## true risk is a fraction of what the model says
    fac <- as.numeric(sub("tail", "", cell$truth))
    x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); p <- plogis(0.6 * x + 0.5 * d)
    top <- p > stats::quantile(p, 0.98)
    pt <- ifelse(top, fac * p, p)
    return(list(d = data.frame(y = rbinom(n, 1, pt), x = x, d = d), f = y ~ x + d, corrupt = integer(0), p_true = pt))
  }
  X <- cbind(matrix(rnorm(n * 8), n, 8), matrix(rbinom(n * 4, 1, 0.3), n, 4)); colnames(X) <- paste0("x", 1:12)
  eta <- -1 + as.numeric(X %*% BETA_B)
  p <- switch(cell$truth, logit = plogis(eta), cloglog = 1 - exp(-exp(eta)), inter = plogis(eta + 0.5 * X[, 1] * X[, 2]))
  y <- rbinom(n, 1, p)
  idx <- integer(0)
  if (cell$corruption %in% c("C1", "C1b")) {
    idx <- sample.int(n, cell$k); X[idx, 1] <- (if (cell$corruption == "C1") 4 else -4) * X[idx, 1]
  }
  list(d = data.frame(y = y, X), f = stats::reformulate(colnames(X), "y"), corrupt = idx, p_true = p)
}

## ---- arms --------------------------------------------------------------------------------------------------------
tw_score_twin <- function(fq) tryCatch({
  y <- fq$y; pr <- fq$p_raw; W <- pr * (1 - pr); X <- fq$X[, !is.na(stats::coef(fq$fit)), drop = FALSE]
  Z <- stats::poly(fq$ph, 3)
  u <- colSums(Z * (y - pr))
  ZWX <- crossprod(Z, W * X)
  I <- crossprod(Z, W * Z) - ZWX %*% solve(crossprod(X, W * X), t(ZWX))
  stats::pchisq(as.numeric(crossprod(u, solve(I, u))), 3, lower.tail = FALSE)
}, error = function(e) NA_real_)
tw_spz <- function(y, p) { v <- p * (1 - p)
  z <- sum((y - p) * (1 - 2 * p)) / sqrt(sum((1 - 2 * p)^2 * v)); 2 * stats::pnorm(-abs(z)) }
tw_hl10 <- function(fq) bt_hl_stat(fq$y, fq$ph, pmin(ceiling(rank(fq$ph, ties.method = "first") / (fq$n / 10)), 10))
tw_edge <- function(fit, G, method = "satterthwaite")
  tryCatch(suppressWarnings(edge.gof(fit, G = G, basis = "poly3", method = method)$p_value), error = function(e) NA_real_)

tw_one <- function(rep, cell) {
  set.seed(980000000 + 10000 * cell$index + rep)
  g <- tw_gen(cell)
  fq <- bt_fit(g)
  if (is.null(fq)) return(c(rep = rep))
  fit <- fq$fit
  out <- c(rep = rep, k_corrupt = length(g$corrupt),
           EDGE.default = tw_edge(fit, "auto"), EDGE.imhof = tw_edge(fit, "auto", "imhof"), EDGE.G10 = tw_edge(fit, 10),
           HL.G10 = tw_hl10(fq), TWIN.score = tw_score_twin(fq),
           EDGE.Gn = if (cell$n <= 1000) tw_edge(fit, cell$n) else NA_real_,
           SPZ = tw_spz(fq$y, fq$p_raw), Stk.joint = bt_stukel(fq)[["Stk.joint"]], Cubic.LR = bt_cubic(fq$y, fq$eta),
           GiViTI = bt_giviti(fq$y, fq$p_raw, 0.95))
  ## Pregibon's diagnostics, then the test on the refitted model without the flagged records
  h <- stats::hatvalues(fit); rp <- stats::residuals(fit, type = "pearson")
  drop <- which(rp^2 / (1 - h) > 4 | rp^2 * h / (1 - h)^2 > 1)
  g2 <- g; if (length(drop)) g2$d <- g$d[-drop, , drop = FALSE]
  fq2 <- bt_fit(g2)
  prg <- if (is.null(fq2)) rep(NA_real_, 4) else
    c(bt_stukel(fq2)[["Stk.joint"]], bt_cubic(fq2$y, fq2$eta), tw_edge(fq2$fit, "auto"), tw_edge(fq2$fit, 10))
  ## Addendum 1: the influence criterion alone, DeltaBeta > 1
  drop2 <- which(rp^2 * h / (1 - h)^2 > 1)
  g3 <- g; if (length(drop2)) g3$d <- g$d[-drop2, , drop = FALSE]
  fq3 <- if (length(drop2)) bt_fit(g3) else fq
  prg2 <- if (is.null(fq3)) rep(NA_real_, 4) else
    c(bt_stukel(fq3)[["Stk.joint"]], bt_cubic(fq3$y, fq3$eta), tw_edge(fq3$fit, "auto"), tw_edge(fq3$fit, 10))
  c(out, PRG.Stk = prg[1], PRG.Cubic = prg[2], PRG.EDGE.default = prg[3], PRG.EDGE.G10 = prg[4],
    PRG.dropped = length(drop), PRG.dropped_corrupt = sum(drop %in% g$corrupt),
    PRG2.Stk = prg2[1], PRG2.Cubic = prg2[2], PRG2.EDGE.default = prg2[3], PRG2.EDGE.G10 = prg2[4],
    PRG2.dropped = length(drop2), PRG2.dropped_corrupt = sum(drop2 %in% g$corrupt))
}

source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
source(file.path(SIMDIR, "_blockC1b_contam.R"))
stopifnot(c1b_identity_check())
cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl))
clusterExport(cl, c("tw_one", "tw_gen", "tw_score_twin", "tw_spz", "tw_hl10", "tw_edge", "BETA_B", "PKG", "SIMDIR"))
invisible(clusterEvalQ(cl, {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  suppressMessages(library(ebrahim.gof))  # archive: was load_all() of the author's development tree; install ebrahim.gof (>= 2.9.0) from CRAN instead
  source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
  source(file.path(SIMDIR, "_blockC1b_contam.R"))
  NULL
}))
tw_lg("block TW: %d cells, %d replicates, %d workers%s", nrow(CELLS), REPS, W, if (SMOKE) " [SMOKE]" else "")
ARMS <- c("EDGE.default", "EDGE.G10", "HL.G10", "TWIN.score", "EDGE.Gn", "SPZ", "Stk.joint", "Cubic.LR", "GiViTI",
          "PRG.Stk", "PRG.EDGE.G10", "PRG2.Stk", "PRG2.EDGE.G10")
for (j in seq_len(nrow(CELLS))) {
  cell <- as.list(CELLS[j])
  f <- file.path(OUT, sprintf("%s_pvalues.csv.gz", cell$cell))
  if (file.exists(f)) { tw_lg("%s: present, skipped", cell$cell); next }
  t0 <- Sys.time()
  M <- rbindlist(lapply(clusterApplyLB(cl, seq_len(REPS), tw_one, cell = cell), function(v) as.data.table(as.list(v))), fill = TRUE)
  fwrite(M, f, compress = "gzip")
  rt <- vapply(ARMS, function(t) if (t %in% names(M)) mean(!is.na(M[[t]]) & M[[t]] <= 0.05) else NA_real_, 0)
  tw_lg("%-34s %s | dropped %.1f  %.1f min", cell$cell, paste(sprintf("%s %.3f", ARMS, rt), collapse = " "),
        mean(M$PRG.dropped, na.rm = TRUE), as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
tw_lg("block TW: run finished")
