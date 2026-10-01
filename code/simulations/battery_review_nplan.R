## battery_review_nplan.R -- independent recomputation of the E0.7 sample-size rule, to check battery/nplan.csv.
## Own link functions (closed forms, not pnorm/pcauchy/pt where a closed form exists), own Newton-Raphson for the
## fractional logistic fit, own AUC. Three versions of D:
##   (1) the same population draw as battery_nplan.R (seed 20260913, N = 400,000): must reproduce nplan.csv;
##   (2) exact population values on a fine midpoint grid over x ~ U(-3, 3), d in {0, 1} (no Monte Carlo error);
##   (3) five other seeds at N = 400,000: Monte Carlo spread of n_env and how often the rounding would change.
## Writes battery/_review/review_nplan.log and review_nplan.csv.
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
OUT <- edge_battery("_review")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "review_nplan.log"), split = TRUE)
cat("battery_review_nplan.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")

## closed-form inverse links
t4cdf <- function(t) { u <- t^2 / 4; 0.5 + 0.375 * (t / sqrt(1 + u)) * (1 - (1 / 12) * t^2 / (1 + u)) }
FL <- list(logit = function(e) 1 / (1 + exp(-e)),
           probit = function(e) pnorm(e),
           cauchit = function(e) 0.5 + atan(e) / pi,
           t4 = t4cdf,
           loglog = function(e) exp(-exp(-e)),
           cloglog = function(e) -expm1(-exp(e)))
eg <- seq(-30, 30, by = 0.001)
cat(sprintf("closed-form checks: |t4cdf - pt(,4)| max %.1e; |cauchy - pcauchy| max %.1e\n\n",
            max(abs(t4cdf(eg) - pt(eg, 4))), max(abs(FL$cauchit(eg) - pcauchy(eg)))))

## fractional logistic fit by Newton-Raphson on weighted rows (weights = probability mass of each row)
frac_fit <- function(X, p, w) {
  b <- rep(0, ncol(X))
  for (it in 1:200) {
    pi <- 1 / (1 + exp(-drop(X %*% b)))
    g <- crossprod(X, w * (p - pi))
    H <- crossprod(X, (w * pi * (1 - pi)) * X)
    st <- drop(solve(H, g)); b <- b + st
    if (max(abs(st)) < 1e-14) break
  }
  pi <- 1 / (1 + exp(-drop(X %*% b)))
  list(b = b, pi = pi, it = it, score = max(abs(crossprod(X, w * (p - pi)))))
}
## AUC of the true probability as a score for y ~ Bernoulli(p), rows weighted by w; ties of p count one half
auc_w <- function(p, w) {
  o <- order(p); p <- p[o]; w <- w[o]
  cas <- w * p; con <- w * (1 - p)
  ## group exact ties
  key <- match(p, unique(p))
  cs <- tapply(cas, key, sum); cn <- tapply(con, key, sum)
  below <- cumsum(cn) - cn
  sum(cs * (below + cn / 2)) / (sum(cas) * sum(con))
}
rnd <- function(v) pmin(pmax(signif(v, 2), 300), 20000)

C0 <- c(logit = -2.64, probit = -2.01, t4 = -2.18, cauchit = -3.15, loglog = -1.60, cloglog = -2.76)
GRID <- list(base = c(1, 0), auc = c(2, 0), e12 = c(1, NA))
ONE <- function(x, d, w, dz, lk) {
  s <- GRID[[dz]][1]; c0 <- if (dz == "e12") C0[[lk]] else 0
  p <- FL[[lk]](c0 + s * (0.6 * x + 0.5 * d))
  ev <- sum(w * p) / sum(w); au <- auc_w(p, w)
  if (lk == "logit") return(c(event = ev, auc = au, D = NA, n_env = NA, it = NA, score = NA))
  f <- frac_fit(cbind(1, x, d), p, w)
  D <- sum(w * (p - f$pi)^2 / (f$pi * (1 - f$pi))) / sum(w)
  c(event = ev, auc = au, D = D, n_env = 7.85 / D, it = f$it, score = f$score)
}

NP <- read.csv(edge_battery("nplan.csv"), stringsAsFactors = FALSE)
todo <- NP[, c("design", "link")]

## (1) the builder's population draw, own code
set.seed(20260913L)
N <- 400000L
x1 <- runif(N, -3, 3); d1 <- rbinom(N, 1, 0.5); w1 <- rep(1, N)
## (2) exact: midpoint grid, M points per d
M <- 200000L
xg <- -3 + 6 * (seq_len(M) - 0.5) / M
x2 <- c(xg, xg); d2 <- rep(0:1, each = M); w2 <- rep(1, 2 * M)

rows <- list()
for (i in seq_len(nrow(todo))) {
  dz <- todo$design[i]; lk <- todo$link[i]
  a <- ONE(x1, d1, w1, dz, lk)
  e <- ONE(x2, d2, w2, dz, lk)
  rows[[i]] <- data.frame(design = dz, link = lk,
    event_nplan = NP$event_rate[i], event_draw = a[["event"]], event_exact = e[["event"]],
    auc_nplan = NP$auc[i], auc_draw = a[["auc"]], auc_exact = e[["auc"]],
    D_nplan = NP$D[i], D_draw = a[["D"]], D_exact = e[["D"]],
    nenv_nplan = NP$n_env[i], nenv_draw = a[["n_env"]], nenv_exact = e[["n_env"]],
    nlo_nplan = NP$n_lo[i], nhi_nplan = NP$n_hi[i],
    nlo_draw = if (is.na(a[["n_env"]])) NA else rnd(0.6 * a[["n_env"]]), nhi_draw = if (is.na(a[["n_env"]])) NA else rnd(a[["n_env"]]),
    nlo_exact = if (is.na(e[["n_env"]])) NA else rnd(0.6 * e[["n_env"]]), nhi_exact = if (is.na(e[["n_env"]])) NA else rnd(e[["n_env"]]),
    control_exact = !is.na(e[["n_env"]]) && e[["n_env"]] > 40000, newton_it = a[["it"]], score_draw = a[["score"]],
    stringsAsFactors = FALSE)
}
R <- do.call(rbind, rows)

## Gauss-Legendre cross-check of the exact D (statmod nodes), for the four cells named in the review
if (requireNamespace("statmod", quietly = TRUE)) {
  gq <- statmod::gauss.quad(2000, "legendre")
  xq <- 3 * gq$nodes; wq <- 3 * gq$weights / 6
  R$D_gauss <- NA_real_
  for (i in seq_len(nrow(R))) if (R$link[i] != "logit") {
    R$D_gauss[i] <- ONE(c(xq, xq), rep(0:1, each = length(xq)), c(wq, wq), R$design[i], R$link[i])[["D"]]
  }
}

## (3) Monte Carlo spread over other seeds
SEEDS <- c(11L, 22L, 33L, 44L, 55L)
mc <- matrix(NA_real_, nrow(R), length(SEEDS))
for (k in seq_along(SEEDS)) {
  set.seed(SEEDS[k]); xs <- runif(N, -3, 3); ds <- rbinom(N, 1, 0.5)
  for (i in seq_len(nrow(R))) if (R$link[i] != "logit") mc[i, k] <- ONE(xs, ds, w1, R$design[i], R$link[i])[["n_env"]]
}
R$nenv_mc_sd <- apply(mc, 1, sd)
R$nlo_mc_values <- apply(mc, 1, function(v) if (all(is.na(v))) "" else paste(unique(rnd(0.6 * v)), collapse = "/"))
R$nhi_mc_values <- apply(mc, 1, function(v) if (all(is.na(v))) "" else paste(unique(rnd(v)), collapse = "/"))

write.csv(R, file.path(OUT, "review_nplan.csv"), row.names = FALSE)

cat("(1) same draw, own code vs nplan.csv: max relative difference\n")
nz <- R$link != "logit"
cat(sprintf("  event %.1e  AUC %.1e  D %.1e  n_env %.1e; n_lo identical %s, n_hi identical %s\n",
            max(abs(R$event_draw / R$event_nplan - 1)), max(abs(R$auc_draw / R$auc_nplan - 1)),
            max(abs(R$D_draw[nz] / R$D_nplan[nz] - 1)), max(abs(R$nenv_draw[nz] / R$nenv_nplan[nz] - 1)),
            identical(R$nlo_draw, R$nlo_nplan), identical(R$nhi_draw, R$nhi_nplan)))
cat("\n(2) exact population values (midpoint grid, 400,000 nodes) and (3) Monte Carlo spread over 5 seeds\n")
cat("design  link     event(pl) event(ex)  AUC(pl) AUC(ex)   n_env(pl) n_env(ex)  D ex/gauss-1  sd(MC)  n_lo pl/ex/MC      n_hi pl/ex/MC\n")
for (i in seq_len(nrow(R))) with(R[i, ], cat(sprintf("%-6s  %-7s  %.4f    %.4f     %.3f   %.3f   %9.0f %9.0f  %9.1e  %6.0f  %5s/%5s/%-11s %5s/%5s/%s\n",
  design, link, event_nplan, event_exact, auc_nplan, auc_exact, ifelse(is.na(nenv_nplan), NA, nenv_nplan), ifelse(is.na(nenv_exact), NA, nenv_exact),
  ifelse(is.na(D_exact), NA, D_exact / D_gauss - 1), ifelse(is.na(nenv_mc_sd), NA, nenv_mc_sd),
  format(nlo_nplan), format(nlo_exact), nlo_mc_values, format(nhi_nplan), format(nhi_exact), nhi_mc_values)))

cat("\nE1 intercepts, exact event rate (target 0.120 +/- 0.005) and the intercept that gives exactly 12%\n")
for (lk in names(C0)) {
  evf <- function(c0) mean(FL[[lk]](c0 + 0.6 * x2 + 0.5 * d2))
  c12 <- uniroot(function(c0) evf(c0) - 0.12, c(-8, 2), tol = 1e-10)$root
  cat(sprintf("  %-7s c0 = %5.2f  exact event rate %.4f  exact 12%% intercept %.4f\n", lk, C0[[lk]], evf(C0[[lk]]), c12))
}
cat("\nReview targets (probit base, cauchit auc, t4 e12, loglog e12):\n")
for (k in list(c("base", "probit"), c("auc", "cauchit"), c("e12", "t4"), c("e12", "loglog"))) {
  z <- R[R$design == k[1] & R$link == k[2], ]
  cat(sprintf("  %-5s %-7s D nplan %.5e | own same draw %.5e | exact %.5e ; n_env %.0f | %.0f | %.0f ; n_lo/n_hi nplan %s/%s, exact %s/%s\n",
              k[1], k[2], z$D_nplan, z$D_draw, z$D_exact, z$nenv_nplan, z$nenv_draw, z$nenv_exact, z$nlo_nplan, z$nhi_nplan, z$nlo_exact, z$nhi_exact))
}
sink()
