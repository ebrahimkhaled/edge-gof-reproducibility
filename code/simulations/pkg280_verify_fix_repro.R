## pkg280_verify_fix_repro.R -- the findings of PKG280_verify_numerics.md and PKG280_verify_regression.md that the
## review fixes (5c029d7..906bd0c) address, re-run on the fixed branch loaded with pkgload::load_all(). The pre-fix
## code (0d63fba) and 2.7.0 (83c7e11) are parsed from git into their own environments for the "before" numbers.
## Also probes the new column rule of def.gof between "stopped before" and "unchanged".
## Run: Rscript pkg280_verify_fix_repro.R [package dir] > ../paper_EDGE/theory/pkg280_verify_fix_repro.log 2>&1
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

args <- commandArgs(trailingOnly = TRUE)
VF  <- file.path(tempdir(), "vfix")
PKG <- if (length(args)) args[1] else file.path(VF, "head")
SIM <- edge_path("code/simulations")
OUT <- edge_path("declarations")
options(width = 220, warn = 1)
suppressPackageStartupMessages(pkgload::load_all(PKG, export_all = FALSE, quiet = TRUE))
ns  <- asNamespace("ebrahim.gof")
sha <- system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE)
cat("ebrahim.gof", as.character(utils::packageVersion("ebrahim.gof")), "via load_all | HEAD", sha, "|",
    R.version.string, "|", format(Sys.time()), "\n")

old_env <- function(rev) {
  e <- new.env(parent = ns)
  KEEP <- c("def.gof", ".def_basis", ".def_pvalue", "def.ensemble.gof", ".combine_pvalues", "gof_stukel",
            ".gof_score_z", ".gof_context")
  for (f in c("R/def_gof.R", "R/def_ensemble_gof.R", "R/run_all_gof.R")) {
    txt <- system2("git", c("-C", shQuote(PKG), "show", paste0(rev, ":", f)), stdout = TRUE)
    for (ex in parse(text = txt, keep.source = FALSE))
      if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) && as.character(ex[[2]]) %in% KEEP)
        eval(ex, e)
  }
  e
}
pre  <- old_env("0d63fba")
v270 <- old_env("83c7e11")
stopifnot("weights" %in% names(formals(pre$def.gof)), !"weights" %in% names(formals(v270$def.gof)))

## helpers from pkg280_verify_num_exact.R (only the named assignments are evaluated)
SRC  <- file.path(SIM, "pkg280_verify_num_exact.R")
TAKE <- c("mk", "find_fit", "augment", "score_qr", "stukel_Z", "stk", "Gof", "step_cov", "count_warn", "quiet")
for (e in parse(SRC, keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% TAKE)
    eval(e, globalenv())
stopifnot(all(TAKE %in% ls(globalenv())))
warm <- function(fit) {                              # tight fit, then a fit started at its coefficients
  env <- environment(formula(fit))
  cl <- getCall(fit); cl$control <- glm.control(epsilon = 1e-13, maxit = 100)
  ft <- suppressWarnings(eval(cl, env))
  cl$start <- unname(coef(ft))
  suppressWarnings(eval(cl, env))
}
row <- function(r) sprintf("stat %-16s df %-3s p %-14s Note '%s'", format(r$Statistic, digits = 10), format(r$df),
                           format(r$p_value, digits = 8), r$Note)
old_stk <- function(env, fit, form = "joint") {
  ctx <- env$.gof_context(fit, NULL, NULL, G = 10)
  r <- tryCatch(env$gof_stukel(ctx, list(form = form)),
                error = function(e) list(Statistic = NA, df = NA, p_value = NA, Note = paste("ERROR", conditionMessage(e))))
  as.data.frame(r, stringsAsFactors = FALSE)
}
dg <- function(expr) {
  r <- tryCatch(suppressWarnings(expr), error = function(e) conditionMessage(e))
  if (is.character(r)) sprintf("ERROR: %s", r)
  else sprintf("S %-16s df %-10s p %-14s (%s)", format(r$Test_Statistic, digits = 10), format(r$df, digits = 7),
               format(r$p_value, digits = 8), r$Method)
}

## ---- 1a. aliased column, numerics R1 ----------------------------------------------------------------------------
cat("\n==== 1a. aliased column (numerics R1): xa2 = 2 * xa, n = 300 ====\n")
set.seed(3)
gd <- data.frame(xa = runif(300, -3, 3)); gd$out <- rbinom(300, 1, plogis(0.5 * gd$xa)); gd$xa2 <- 2 * gd$xa
fit_rd <- glm(out ~ xa + xa2, family = binomial(), data = gd)
fit_fr <- glm(out ~ xa, family = binomial(), data = gd)
cat("  coef(fit_rd):", format(coef(fit_rd), digits = 6), "\n")
for (f in c("joint", "lr", "marginal")) {
  cat(sprintf("  %-8s aliased fit, fixed  : %s\n", f, row(stk(fit_rd, f))))
  cat(sprintf("  %-8s full-rank fit, fixed: %s\n", f, row(stk(fit_fr, f))))
  cat(sprintf("  %-8s aliased fit, 0d63fba: %s\n", f, row(old_stk(pre, fit_rd, f))))
}
a_rd <- anova(fit_rd, augment(fit_rd, stukel_Z(fit_rd)), test = "Rao")
q_rd <- score_qr(fit_rd, stukel_Z(fit_rd))
cat(sprintf("  anova Rao on the aliased fit %.10f (Df %g), deviance drop %.10f | u'I^-1 u by pivoted weighted QR %.10f (df %d)\n",
            a_rd$Rao[2], a_rd$Df[2], a_rd$Deviance[2], q_rd$stat, q_rd$df))
j_rd <- stk(fit_rd, "joint")$Statistic; j_fr <- stk(fit_fr, "joint")$Statistic
cat(sprintf("  |joint aliased - joint full rank| %.2e | |joint aliased - QR| %.2e | |lr aliased - deviance drop| %.2e\n",
            abs(j_rd - j_fr), abs(j_rd - q_rd$stat), abs(stk(fit_rd, "lr")$Statistic - a_rd$Deviance[2])))
fw <- warm(fit_fr); a_w <- anova(fw, augment(fw, stukel_Z(fw)), test = "Rao")
cat(sprintf("  full-rank fit, warm restart: joint %.12f | anova Rao %.12f | |diff| %.2e\n",
            stk(fw, "joint")$Statistic, a_w$Rao[2], abs(stk(fw, "joint")$Statistic - a_w$Rao[2])))

## ---- 1b. aliased column, regression F1 design D5 -------------------------------------------------------------------
cat("\n==== 1b. aliased column (regression F1, design D5): x3 = 2 * x1, n = 800 ====\n")
set.seed(11); n <- 800
d <- data.frame(x1 = rnorm(n), x2 = runif(n, -2, 2), f = factor(sample(letters[1:3], n, TRUE)))
d$y <- rbinom(n, 1, plogis(-0.3 + 0.9 * d$x1 - 0.6 * d$x2 + c(0, 0.4, -0.5)[d$f]))
d$o <- runif(n, -1, 1); d$y2 <- rbinom(n, 1, plogis(0.2 + 0.7 * d$x1 + d$o))
d$x3 <- 2 * d$x1
fa <- glm(y ~ x1 + x2 + x3, data = d, family = binomial())
fb <- glm(y ~ x1 + x2, data = d, family = binomial())
for (f in c("joint", "lr", "marginal")) {
  cat(sprintf("  %-8s D5 aliased, fixed  : %s\n", f, row(stk(fa, f))))
  cat(sprintf("  %-8s D5 without x3      : %s\n", f, row(stk(fb, f))))
  cat(sprintf("  %-8s D5 aliased, 0d63fba: %s\n", f, row(old_stk(pre, fa, f))))
}
cat("  battery rows on D5 (DEF rows still stop on an aliased fit, as in 2.7.0):\n")
print(as.data.frame(run.all.gof(fa, tests = c("DEF.poly3", "DEF.sym", "Stukel"), install = "no"))[, c("Test", "Statistic", "df", "p_value", "Note")],
      row.names = FALSE)

## ---- 1c. clamp under near separation, regression F2 ---------------------------------------------------------------
cat("\n==== 1c. near separation (regression F2): slope 6, n = 300 ====\n")
set.seed(14); x <- runif(300, -3, 3); y <- rbinom(300, 1, plogis(6 * x))
dsep <- data.frame(x = x, y = y)
fs <- suppressWarnings(glm(y ~ x, family = binomial(), data = dsep))
p <- fitted(fs)
cat("  fitted risks outside [1e-6, 1 - 1e-6]:", sum(p < 1e-6 | p > 1 - 1e-6), "of 300 | converged", fs$converged, "\n")
cat("  joint, fixed  :", row(stk(fs, "joint")), "\n")
cat("  joint, 0d63fba:", row(old_stk(pre, fs, "joint")), "\n")
cat("  lr, fixed     :", row(stk(fs, "lr")), "\n")
q_s <- score_qr(fs, stukel_Z(fs)); a_s <- anova(fs, augment(fs, stukel_Z(fs)), test = "Rao")
cat(sprintf("  u'I^-1 u by QR on the raw fitted risks %.10f | anova Rao (default control) %.10f | |joint - QR| %.2e | rel |joint - anova| %.2e\n",
            q_s$stat, a_s$Rao[2], abs(stk(fs, "joint")$Statistic - q_s$stat),
            abs(stk(fs, "joint")$Statistic / a_s$Rao[2] - 1)))
fsw <- warm(fs); a_sw <- anova(fsw, augment(fsw, stukel_Z(fsw)), test = "Rao")
cat(sprintf("  warm restart: joint %.10f | anova Rao %.10f | rel diff %.2e | QR %.10f | converged %s\n",
            stk(fsw, "joint")$Statistic, a_sw$Rao[2], abs(stk(fsw, "joint")$Statistic / a_sw$Rao[2] - 1),
            score_qr(fsw, stukel_Z(fsw))$stat, fsw$converged))
cat("  complete separation (y = 1 when x > 0, n = 200), informational:\n")
set.seed(15); xc <- runif(200, -3, 3); dcs <- data.frame(x = xc, y = as.numeric(xc > 0))
fcs <- suppressWarnings(glm(y ~ x, family = binomial(), data = dcs))
cat("    fitted range", format(range(fitted(fcs)), digits = 3), "\n")
cat("    joint, fixed  :", row(stk(fcs, "joint")), "\n")
cat("    joint, 0d63fba:", row(old_stk(pre, fcs, "joint")), "\n")
cat("    anova Rao     :", tryCatch(format(anova(fcs, augment(fcs, stukel_Z(fcs)), test = "Rao")$Rao[2], digits = 10),
                                    error = function(e) conditionMessage(e)), "\n")

## ---- 1d. near-empty Stukel half-column, numerics R2 ----------------------------------------------------------------
cat("\n==== 1d. near-empty Stukel half-column (numerics R2): top group mean risk 0.5001, n = 500, G = 20 ====\n")
set.seed(4); n <- 500
ph2 <- c(sort(runif(475, 0.02, 0.45)), rep(0.5001, 25)); X2 <- cbind(1, qlogis(ph2)); y2 <- rbinom(n, 1, ph2)
for (m in c("satterthwaite", "imhof")) {
  cat(sprintf("  unit %-13s fixed  : %s\n", m, dg(def.gof(y2, ph2, X = X2, G = 20, basis = "stukel", method = m))))
  cat(sprintf("  unit %-13s 0d63fba: %s\n", m, dg(pre$def.gof(y2, ph2, X = X2, G = 20, basis = "stukel", method = m))))
  cat(sprintf("  unit %-13s 2.7.0  : %s\n", m, dg(v270$def.gof(y2, ph2, X = X2, G = 20, basis = "stukel", method = m))))
}
cat(sprintf("  score              fixed  : %s\n", dg(def.gof(y2, ph2, X = X2, G = 20, basis = "stukel", weights = "score"))))
cat(sprintf("  score              0d63fba: %s\n", dg(pre$def.gof(y2, ph2, X = X2, G = 20, basis = "stukel", weights = "score"))))
cat(sprintf("  edge.gof stukel    fixed  : %s\n", dg(edge.gof(y2, ph2, X = X2, G = 20, basis = "stukel"))))
## the unit form by hand, without the negligible column
g  <- pmin(ceiling(rank(ph2, ties.method = "first") / 25), 20)
V  <- ph2 * (1 - ph2); Vg <- as.numeric(tapply(V, g, sum))
r  <- as.numeric(tapply(y2 - ph2, g, sum)) / sqrt(Vg)
e  <- qlogis(as.numeric(tapply(ph2, g, mean)))
Zh <- cbind(e, -e^2 * (e < 0))
U  <- rowsum(V * X2, g) / sqrt(Vg)
Om <- diag(20) - U %*% solve(crossprod(X2, V * X2), t(U))
Sh <- drop(crossprod(r, Zh %*% solve(crossprod(Zh), crossprod(Zh, r))))
lam <- Re(eigen(solve(crossprod(Zh), crossprod(Zh, Om %*% Zh)), only.values = TRUE)$values); lam <- lam[lam > 1e-9]
ph_hand <- pchisq(Sh / (sum(lam^2) / sum(lam)), sum(lam)^2 / sum(lam^2), lower.tail = FALSE)
fx <- def.gof(y2, ph2, X = X2, G = 20, basis = "stukel")
cat(sprintf("  by hand without the column: S %.12f p %.12f | |fixed - hand| S %.2e p %.2e\n", Sh, ph_hand,
            abs(fx$Test_Statistic - Sh), abs(fx$p_value - ph_hand)))
## the score statistic with all three columns kept, scaled to unit length (a scale-free rank rule)
Z3 <- cbind(e, e^2 * (e >= 0), -e^2 * (e < 0)); Z3 <- Z3 / rep(sqrt(colSums(Z3^2)), each = 20)
Zs <- Z3 * sqrt(Vg); u3 <- drop(crossprod(Zs, r)); I3 <- crossprod(Zs, Om %*% Zs)
q3 <- qr(I3); cat(sprintf("  score form with all 3 columns (unit-length, pivoted QR rank): S %.10f, rank %d, p %.8g\n",
                          drop(crossprod(u3, qr.solve(I3, u3))), q3$rank, pchisq(drop(crossprod(u3, qr.solve(I3, u3))), q3$rank, lower.tail = FALSE)))

## ---- 1e. the new column rule between "stopped before" and "unchanged" ----------------------------------------------
cat("\n==== 1e. top group mean risk v just above 0.5 (the R2 construction, n = 500, G = 20) ====\n")
cat("  ratio = norm of the e^2 (e >= 0) column / largest column norm; the fixed rule drops the column when ratio < 1e-6\n")
for (v in 0.5 + c(1e-6, 1e-5, 5e-5, 1e-4, 1.7e-4, 2.5e-4, 5e-4, 1e-3, 1.5e-3, 1.6e-3, 2e-3, 5e-3)) {
  set.seed(4); ph <- c(sort(runif(475, 0.02, 0.45)), rep(v, 25)); X <- cbind(1, qlogis(ph)); yy <- rbinom(500, 1, ph)
  gg <- pmin(ceiling(rank(ph, ties.method = "first") / 25), 20); pbar <- as.numeric(tapply(ph, gg, mean))
  Z <- ns$.def_basis(pbar, "stukel"); nz <- sqrt(colSums(Z^2))
  cat(sprintf("  v = %-10s ratio %.2e\n", format(v, digits = 8), nz[2] / max(nz)))
  cat(sprintf("      unit  0d63fba %s\n      unit  fixed   %s\n", dg(pre$def.gof(yy, ph, X = X, G = 20, basis = "stukel")),
              dg(def.gof(yy, ph, X = X, G = 20, basis = "stukel"))))
  cat(sprintf("      score 0d63fba %s\n      score fixed   %s\n",
              dg(pre$def.gof(yy, ph, X = X, G = 20, basis = "stukel", weights = "score")),
              dg(def.gof(yy, ph, X = X, G = 20, basis = "stukel", weights = "score"))))
}

## ---- 1f. glm fits in that window: score form against anova Rao on the step covariates -----------------------------
cat("\n==== 1f. c0 = -2, n = 500, G = 20 fits whose groups with mean risk >= 0.5 all lie below 0.5016 ====\n")
found <- list(); tried <- 0
for (sd in 300:40000) {
  tried <- tried + 1
  f <- mk(sd, 500, c0 = -2); php <- pmin(pmax(fitted(f), 1e-6), 1 - 1e-6)
  gg <- pmin(ceiling(rank(php, ties.method = "first") / 25), 20); pbar <- tapply(php, gg, mean)
  up <- pbar[pbar >= 0.5]
  if (length(up) >= 1 && max(up) < 0.5016) { attr(f, "seed") <- sd; found[[length(found) + 1]] <- f
                                             if (length(found) == 4) break }
}
cat("  seeds tried:", tried, "| found:", length(found), "\n")
for (f in found) {
  fwin <- warm(f)
  sc  <- step_cov(fwin, 20, "stukel")                 # keeps every column with sum |z| > 1e-8, as anova would see it
  a   <- anova(fwin, augment(fwin, sc), test = "Rao"); q <- score_qr(fwin, sc)
  pb  <- attr(sc, "pbar")
  cat(sprintf("  seed %d | group means >= 0.5: %s\n", attr(f, "seed"), paste(format(pb[pb >= 0.5], digits = 8), collapse = " ")))
  cat(sprintf("      score fixed   %s\n      score 0d63fba %s\n", dg(def.gof(fwin, G = 20, basis = "stukel", weights = "score")),
              dg(pre$def.gof(fwin, G = 20, basis = "stukel", weights = "score"))))
  cat(sprintf("      anova Rao on the step covariates: %.10f, Df %g, p %.8g | u'I^-1 u by QR %.10f, rank %d\n",
              a$Rao[2], a$Df[2], a$`Pr(>Chi)`[2], q$stat, q$df))
  cat(sprintf("      unit  fixed   %s\n      unit  0d63fba %s\n", dg(def.gof(fwin, G = 20, basis = "stukel")),
              dg(pre$def.gof(fwin, G = 20, basis = "stukel"))))
}

## ---- 1g. the one null replicate that stopped (numerics finding 2: c0 = -2, n = 500, G = auto, rep 524) --------------
cat("\n==== 1g. size-run replicate 524 (c0 = -2, n = 500, G = auto) rebuilt from its RNG stream ====\n")
Pcsv <- read.csv(file.path(OUT, "pkg280_verify_num_size_pvalues.csv"), check.names = FALSE)
lab <- "c0=-2 n500 Gauto"; cellk <- 3L
bad <- Pcsv[Pcsv$cell == lab & !is.finite(Pcsv$stukel.unit), ]
cat("  replicates without a stukel unit p-value in the 2.8.0 size run:", paste(bad$rep, collapse = ", "), "\n")
RNGkind("L'Ecuyer-CMRG"); set.seed(20260914L + 100L * cellk)
seeds <- list(.Random.seed); for (i in 2:6) seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
chunks <- parallel::splitIndices(2000, 6)
regen <- function(r) {
  j <- which(vapply(chunks, function(ch) r %in% ch, logical(1)))
  assign(".Random.seed", seeds[[j]], envir = globalenv())
  for (qq in chunks[[j]]) {
    xa <- runif(500, -3, 3); db <- rbinom(500, 1, 0.5); resp <- rbinom(500, 1, plogis(-2 + (0.6 * xa + 0.5 * db)))
    if (qq == r) break
  }
  dd <- data.frame(resp = resp, xa = xa, db = db)
  suppressWarnings(glm(resp ~ xa + db, family = binomial(), data = dd))
}
for (r in bad$rep) {
  fr <- regen(r); rw <- Pcsv[Pcsv$cell == lab & Pcsv$rep == r, ]
  cat(sprintf("  rep %d rebuilt: sym unit p %.10g (csv %.10g), poly3 unit p %.10g (csv %.10g)\n", r,
              quiet(def.gof(fr, G = "auto", basis = "sym"))$p_value, rw$sym.unit,
              quiet(def.gof(fr, G = "auto", basis = "poly3"))$p_value, rw$poly3.unit))
  php <- pmin(pmax(fitted(fr), 1e-6), 1 - 1e-6); gg <- pmin(ceiling(rank(php, ties.method = "first") / 25), 20)
  cat("  top three group means:", format(tail(sort(tapply(php, gg, mean)), 3), digits = 7), "\n")
  for (m in c("satterthwaite", "imhof")) {
    cat(sprintf("  unit %-13s fixed  : %s\n", m, dg(def.gof(fr, G = "auto", basis = "stukel", method = m))))
    cat(sprintf("  unit %-13s 0d63fba: %s\n", m, dg(pre$def.gof(fr, G = "auto", basis = "stukel", method = m))))
  }
  cat(sprintf("  score              fixed  : %s\n", dg(def.gof(fr, G = "auto", basis = "stukel", weights = "score"))))
  cat(sprintf("  score              0d63fba: %s\n", dg(pre$def.gof(fr, G = "auto", basis = "stukel", weights = "score"))))
  cat("  battery DEF.stukel at G = auto:\n")
  print(as.data.frame(run.all.gof(fr, G = "auto", tests = c("DEF.stukel", "Ensemble.Vote(3DEF)"), install = "no")), row.names = FALSE)
}
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
cat("\ndone", format(Sys.time()), "\n")
