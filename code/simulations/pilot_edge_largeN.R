## pilot_edge_largeN.R -- PILOT (exploratory, not declared): a large-sample EDGE in the manner of Nattino, Pennell and
## Lemeshow (2020).
##
## In large samples every consistent test rejects a misfit too small to matter, because its noncentrality grows with n.
## Nattino et al. test H0: eps <= eps0 for the Hosmer-Lemeshow statistic, eps = sqrt(lambda / n). The EDGE statistic
## S = r' P_Z r has the same structure: under misfit its mean is sum(lam) + lambda with lambda proportional to n. With
## the Satterthwaite reference S ~ c chi2_nu (c = sum lam^2 / sum lam, nu = (sum lam)^2 / sum lam^2), a noncentral
## version S ~ c chi2_nu(ncp = lambda / c) keeps the mean right. The large-sample EDGE refers S to that law with
## lambda = eps0^2 n, eps0 chosen by Nattino's convention: the misfit whose expected statistic sits at the 5% critical
## value at n0 = 10^6, eps0^2 = (c qchisq(.95, nu) - sum lam) / n0.
##
## Design: x ~ N(0, 1), truth logit p = -1.5 + x + delta (x^2 - 1), model y ~ x; delta in {0, .02, .05, .10};
## n in {5k, 20k, 100k, 500k}; 400 replicates; ten groups. Tests: HL, HL-largeN, EDGE (G = 10), EDGE-largeN.
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
OUT <- edge_battery("PILOT_largeN"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
REPS <- 400L; N0 <- 1e6

## EDGE at G groups, unit form, internal mode: the statistic and its weights, as def.gof() computes them (rowsum for
## speed at n = 500,000; the check below compares the p-value with def.gof()).
edge_core <- function(y, ph, X, G = 10) {
  n <- length(y); V <- ph * (1 - ph)
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  og <- rowsum(y, grp)[, 1]; eg <- rowsum(ph, grp)[, 1]; Vg <- rowsum(V, grp)[, 1]
  pbar <- eg / rowsum(rep(1, n), grp)[, 1]
  r <- (og - eg) / sqrt(Vg)
  U <- rowsum(V * X, grp) / sqrt(Vg)                      # logit link: dmu = V
  Omega <- diag(G) - U %*% solve(crossprod(X, V * X)) %*% t(U)
  Z <- stats::poly(pbar, 3); Z <- Z / rep(sqrt(colSums(Z^2)), each = nrow(Z))
  ZtZ <- crossprod(Z); Zr <- crossprod(Z, r)
  S <- as.numeric(t(Zr) %*% solve(ZtZ) %*% Zr)
  lam <- Re(eigen(solve(ZtZ) %*% (t(Z) %*% Omega %*% Z), only.values = TRUE)$values); lam <- lam[lam > 1e-9]
  list(S = S, lam = lam)
}
edge_p <- function(S, lam, n, largeN = FALSE) {
  cc <- sum(lam^2) / sum(lam); nu <- sum(lam)^2 / sum(lam^2)
  if (!largeN) return(stats::pchisq(S / cc, nu, lower.tail = FALSE))
  eps0_sq <- (cc * stats::qchisq(0.95, nu) - sum(lam)) / N0
  stats::pchisq(S / cc, nu, ncp = eps0_sq * n / cc, lower.tail = FALSE)
}
hl_core <- function(y, ph, G = 10) {
  n <- length(y); grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  O <- rowsum(y, grp)[, 1]; E <- rowsum(ph, grp)[, 1]; m <- tabulate(grp, G)
  sum((O - E)^2 / (E * (1 - E / m)))
}

one <- function(rep, cell) {
  set.seed(9500000 + 10000 * cell$id + rep)
  n <- cell$n; x <- rnorm(n)
  y <- rbinom(n, 1, plogis(-1.5 + x + cell$delta * (x^2 - 1)))
  fit <- suppressWarnings(glm.fit(cbind(1, x), y, family = binomial()))
  ph <- pmin(pmax(fit$fitted.values, 1e-6), 1 - 1e-6); X <- cbind(1, x)
  e <- edge_core(y, ph, X, 10); C <- hl_core(y, ph, 10)
  c(rep = rep, HL = pchisq(C, 8, lower.tail = FALSE),
    HL.largeN = ebrahim.gof:::.hl_largeN_p(C, n, 8, N0)$p,
    EDGE = edge_p(e$S, e$lam, n), EDGE.largeN = edge_p(e$S, e$lam, n, TRUE),
    eps_EDGE = sqrt(max(e$S - sum(e$lam), 0) / n), eps_HL = sqrt(max(C - 8, 0) / n))
}

## the check: edge_core reproduces def.gof() at G = 10
suppressMessages(library(ebrahim.gof))  # archive: was load_all() of the author's development tree; install ebrahim.gof (>= 2.9.0) from CRAN instead
set.seed(1); xx <- rnorm(3000); yy <- rbinom(3000, 1, plogis(-1.5 + xx + 0.1 * (xx^2 - 1)))
ff <- glm(yy ~ xx, family = binomial())
ref <- def.gof(ff, G = 10)$p_value
mine <- with(edge_core(yy, pmin(pmax(fitted(ff), 1e-6), 1 - 1e-6), cbind(1, xx)), edge_p(S, lam, 3000))
cat(sprintf("check: def.gof p = %.10f, pilot p = %.10f\n", ref, mine)); stopifnot(abs(ref - mine) < 1e-8)

## the size of each misfit in absolute terms: ICI over a large sample at the fitted model's limit
ici <- vapply(c(0, .02, .05, .10), function(d) {
  set.seed(77); x <- rnorm(2e6); pt <- plogis(-1.5 + x + d * (x^2 - 1)); y <- rbinom(2e6, 1, pt)
  pm <- suppressWarnings(glm.fit(cbind(1, x), y, family = binomial()))$fitted.values; mean(abs(pt - pm)) }, 0)
cat("ICI by delta (0, .02, .05, .10):", sprintf("%.4f", ici), "\n")

CELLS <- CJ(delta = c(0, 0.02, 0.05, 0.10), n = c(5000L, 20000L, 100000L, 500000L)); CELLS[, id := .I]
cl <- makePSOCKcluster(16); on.exit(stopCluster(cl))
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl, c("one", "edge_core", "edge_p", "hl_core", "N0", "PKG"))
invisible(clusterEvalQ(cl, { Sys.setenv(OMP_NUM_THREADS = "1"); suppressMessages(library(ebrahim.gof)); NULL }))  # archive: was load_all() of the author's development tree; install ebrahim.gof (>= 2.9.0) from CRAN instead
res <- list()
for (j in seq_len(nrow(CELLS))) {
  cell <- as.list(CELLS[j])
  f <- file.path(OUT, sprintf("d%.2f_n%d.csv", cell$delta, cell$n))
  M <- if (file.exists(f)) fread(f) else {
    M <- rbindlist(lapply(clusterApplyLB(cl, seq_len(REPS), one, cell = cell), function(v) as.data.table(as.list(v))))
    fwrite(M, f); M }
  res[[j]] <- data.table(delta = cell$delta, n = cell$n, ICI = ici[match(cell$delta, c(0, .02, .05, .10))],
                         HL = mean(M$HL < .05), HL.largeN = mean(M$HL.largeN < .05),
                         EDGE = mean(M$EDGE < .05), EDGE.largeN = mean(M$EDGE.largeN < .05),
                         eps_HL = median(M$eps_HL), eps_EDGE = median(M$eps_EDGE))
  print(res[[j]], digits = 3)
}
R <- rbindlist(res); fwrite(R, file.path(OUT, "_summary.csv"))
cat("\n"); print(R, digits = 3)
