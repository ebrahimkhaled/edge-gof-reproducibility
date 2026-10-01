## bench_time_T_lecessie_n3.R -- block T addendum: le Cessie's test in the textbook O(n^3) form.
##
## ebrahim.gof computes le Cessie-van Houwelingen's moments through the rank-p factors of the hat matrix,
## O(n^2 p). The public code it was adapted from (smwrStats::leCessie.test) forms the n x n hat matrix and
## multiplies n x n matrices, O(n^3). Both give the same statistic (the transpose fix of 2026-07-29 is kept
## here, so the numbers agree to machine precision); only the time differs. Same data, same single core,
## same timer as bench_time_T.R.
##
##   Rscript bench_time_T_lecessie_n3.R   -> battery/T/timing_lecessie_n3.csv
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
PKG    <- Sys.getenv("EBRAHIM_GOF_SRC")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
suppressPackageStartupMessages(library(ebrahim.gof))  # archive: was load_all() of the author's development tree; install ebrahim.gof (>= 2.9.0) from CRAN instead

make_data <- function(n) {                               # bench_time_T.R's data, unchanged
  set.seed(20260920 + n)
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(0.6 * x + 0.5 * d))
  dat <- data.frame(y = y, x = x, d = d)
  list(dat = dat, fit = glm(y ~ x + d, family = binomial(), data = dat))
}

## the O(n^3) form: H = V X (X'VX)^{-1} X', M = (I-H)' R (I-H), Var = sum dM^2 (mu4 - 3 mu2^2) + 2 tr(M V M V)
lecessie_n3 <- function(fit) {
  y <- fit$y; p <- fitted(fit); r <- y - p; N <- length(y)
  covs <- model.frame(fit)[, -1, drop = FALSE]
  dl <- lapply(covs, function(x) if (is.numeric(x)) as.numeric(0.5 * stats::dist(scale(x))^2)
               else { xx <- as.numeric(as.factor(x)); nc <- length(unique(xx))
                      as.numeric((stats::dist(xx, method = "manhattan") != 0) * nc / (nc - 1)) })
  D <- matrix(0, N, N); D[lower.tri(D)] <- sqrt(rowSums(as.data.frame(dl))); D <- D + t(D)
  R <- pmax(1 - D / mean(D), 0)
  Q <- sum(as.numeric(r %*% R) * r)
  X <- model.matrix(fit); mu2 <- p * (1 - p)
  H <- (mu2 * X) %*% solve(crossprod(X, mu2 * X)) %*% t(X)
  IH <- diag(N) - H
  M <- t(IH) %*% R %*% IH; M <- (M + t(M)) / 2
  V <- diag(mu2)
  EQ <- sum(diag(M) * mu2)
  mu4 <- mu2 * (1 - 3 * mu2)
  VarQ <- sum(diag(M)^2 * (mu4 - 3 * mu2^2)) + 2 * sum(diag(M %*% V %*% M %*% V))
  stat <- Q * 2 * EQ / VarQ; df <- 2 * EQ^2 / VarQ
  c(stat = stat, df = df, p = stats::pchisq(stat, df, lower.tail = FALSE))
}

## agreement with the package at n = 500 before any timing
D5 <- make_data(500L)
a <- lecessie_n3(D5$fit)
b <- run.all.gof(D5$fit, tests = "le-Cessie")
bp <- as.numeric(b[grepl("Cessie", b[[1]]), grep("p", names(b), ignore.case = TRUE)[1]])
cat(sprintf("n=500: O(n^3) p = %.10f | ebrahim.gof p = %.10f\n", a[["p"]], bp))
stopifnot(abs(a[["p"]] - bp) < 1e-8)

OUT <- edge_battery("T", "timing_lecessie_n3.csv")
rows <- list()
for (n in c(500L, 1000L, 2000L, 5000L)) {
  D <- make_data(n); if (n <= 1000) invisible(lecessie_n3(D$fit))   # untimed warm-up at n <= 1000 only
  reps <- if (n <= 1000) 3L else 1L
  s <- vapply(seq_len(reps), function(i) system.time(lecessie_n3(D$fit))[["elapsed"]], 0)
  rows[[length(rows) + 1]] <- data.frame(test = "le Cessie, O(n^3) form", n = n, median_sec = median(s),
                                         min_sec = min(s), reps = reps)
  cat(sprintf("%s  n=%5d  median %.2f s\n", format(Sys.time(), "%H:%M:%S"), n, median(s)))
  write.csv(do.call(rbind, rows), OUT, row.names = FALSE)
}
cat("written:", OUT, "\n")
