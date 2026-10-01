## map_pkg_rankcheck.R -- the score form's "column rank" rule (MAP_package.md). When a Stukel half-column is non-zero
## only in a group whose mean risk sits just above 0.5, that column is tiny and the post-fit information I has an
## eigenvalue near 1e-15 relative to its largest (seen in map_pkg_nullcheck.R at c0 = -2). Two rules on null data:
##   relmax : keep eigenvalues of I above 1e-9 * max eigenvalue   (depends on column scale)
##   cor    : rescale I to a correlation matrix, keep eigenvalues above 1e-8 (depends on collinearity only)
## 5 workers, B = 2000 per cell. Not a paper result.
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
PROTO <- edge_path("code/simulations/map_pkg_proto.R")
source(PROTO)
suppressPackageStartupMessages(library(parallel))

score_two_rules <- function(fit, G, basis) {
  y <- as.numeric(fit$y); ph <- pmin(pmax(as.numeric(stats::fitted(fit)), 1e-6), 1 - 1e-6)
  X <- stats::model.matrix(fit); n <- length(y); V <- ph * (1 - ph)
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  og <- as.numeric(rowsum(y, grp)); eg <- as.numeric(rowsum(ph, grp)); Vg <- as.numeric(rowsum(V, grp))
  pbar <- eg / as.numeric(table(grp)); r <- (og - eg) / sqrt(Vg)
  U <- rowsum(V * X, grp) / sqrt(Vg)
  Omega <- diag(length(Vg)) - U %*% solve(crossprod(X, V * X)) %*% t(U)
  Z <- proto_basis(pbar, basis); Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]
  Zs <- Z * sqrt(Vg); u <- drop(crossprod(Zs, r)); Iz <- crossprod(Zs, Omega %*% Zs); Iz <- (Iz + t(Iz)) / 2
  e1 <- eigen(Iz, symmetric = TRUE); ok1 <- e1$values > 1e-9 * max(e1$values)
  T1 <- sum(drop(crossprod(e1$vectors[, ok1, drop = FALSE], u))^2 / e1$values[ok1]); k1 <- sum(ok1)
  d <- sqrt(pmax(diag(Iz), 0)); keep <- d > 0
  R <- Iz[keep, keep, drop = FALSE] / outer(d[keep], d[keep]); us <- u[keep] / d[keep]
  e2 <- eigen(R, symmetric = TRUE); ok2 <- e2$values > 1e-8
  T2 <- sum(drop(crossprod(e2$vectors[, ok2, drop = FALSE], us))^2 / e2$values[ok2]); k2 <- sum(ok2)
  c(p1 = stats::pchisq(T1, k1, lower.tail = FALSE), k1 = k1, p2 = stats::pchisq(T2, k2, lower.tail = FALSE), k2 = k2,
    ncol = ncol(Z), mind = min(d) / max(d))
}

BASES <- c("stukel", "stukel2", "poly3", "sym")
one_rep <- function(i, s, c0, n, G) {
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); y <- rbinom(n, 1, plogis(c0 + s * (0.6 * x + 0.5 * d)))
  fit <- tryCatch(suppressWarnings(glm(y ~ x + d, family = binomial())), error = function(e) NULL)
  out <- matrix(NA_real_, length(BASES), 6, dimnames = list(BASES, c("p1", "k1", "p2", "k2", "ncol", "mind")))
  if (is.null(fit)) return(as.vector(t(out)))
  for (b in BASES) out[b, ] <- tryCatch(score_two_rules(fit, G, b), error = function(e) rep(NA_real_, 6))
  as.vector(t(out))
}
CELLS <- data.frame(label = c("c0=-2 n500 G20", "c0=-2 n1000 G40", "s=2 c0=-2 n1000 G40", "base n1000 G40"),
                    s = c(1, 1, 2, 1), c0 = c(-2, -2, -2, 0), n = c(500, 1000, 1000, 1000), G = c(20, 40, 40, 40))
B <- 2000L
cl <- makeCluster(5L)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl, c("PROTO", "BASES", "one_rep", "score_two_rules"))
invisible(clusterEvalQ(cl, { source(PROTO); NULL }))
t0 <- Sys.time()
for (k in seq_len(nrow(CELLS))) {
  ce <- CELLS[k, ]; clusterSetRNGStream(cl, 20260920L + k)
  M <- do.call(rbind, parLapply(cl, seq_len(B), one_rep, s = ce$s, c0 = ce$c0, n = ce$n, G = ce$G))
  cat(sprintf("\n== %s (B = %d) ==\n", ce$label, B))
  for (j in seq_along(BASES)) {
    cols <- (j - 1) * 6 + 1:6; S <- M[, cols, drop = FALSE]
    cat(sprintf("  %-8s relmax: size05 %.4f  df table %s | cor: size05 %.4f  df table %s | ncol table %s | min col-scale ratio 1%% %.1e\n",
                BASES[j], mean(S[, 1] <= .05, na.rm = TRUE), paste(names(table(S[, 2])), table(S[, 2]), sep = ":", collapse = " "),
                mean(S[, 3] <= .05, na.rm = TRUE), paste(names(table(S[, 4])), table(S[, 4]), sep = ":", collapse = " "),
                paste(names(table(S[, 5])), table(S[, 5]), sep = ":", collapse = " "), quantile(S[, 6], .01, na.rm = TRUE)))
    dis <- which(S[, 2] != S[, 4])
    if (length(dis)) cat(sprintf("           rules disagree on df in %d reps; size05 in those reps relmax %.3f cor %.3f\n",
                                 length(dis), mean(S[dis, 1] <= .05), mean(S[dis, 3] <= .05)))
  }
}
stopCluster(cl)
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 2)), "\n")
