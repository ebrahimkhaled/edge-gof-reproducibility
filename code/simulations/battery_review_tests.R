## battery_review_tests.R -- battery_rep() against independent code on five data sets (review item c).
## For each data set and replicate: the harness vector from battery_one(), the same data re-drawn with own generator
## code, then
##   Stukel joint  vs own u'I^-1 u (QR residual form) on a tightly converged fit, and anova Rao on a warm-restarted fit
##   Stukel LR     vs the deviance drop of glm(. + za + zb)
##   Stukel sym1   vs anova Rao for the one column eta|eta|
##   Stukel marg   vs statmod::glm.scoretest on (za, zb)
##   EDGE score    vs anova Rao for the grouped step covariates (raw powers of pbar, own grouping), every basis and G arm
##   EDGE unit     vs edge_eig() of run_K_basis_score.R, every basis and G arm
##   GiViTI        vs a direct givitiR::givitiCalibrationTest call at thres 0.95 and 0.50 (called twice: determinism)
##   Cubic.LR      vs anova LRT of glm(y ~ eta + eta^2 + eta^3) against glm(y ~ eta)
##   flags         events < G (min of events and non-events), half-column zero, stuk_diag() of _dgp_library.R
## Writes battery/_review/review_tests.log and review_tests.csv.
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
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))
sink(file.path(OUT, "review_tests.log"), split = TRUE)
cat("battery_review_tests.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")

envK <- new.env()
for (e in parse(file.path(SIMDIR, "run_K_basis_score.R")))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% c("grp_stats", "edge_eig"))
    eval(e, envK)

Cells <- battery_cells()
PICK <- list(c("1b", "null_base_n1000"), c("2", "probit_auc_n3500"), c("1b", "sparse49_n300"),
             c("2", "cauchit_e12_n6400"), c("3", "stk_short_n1000"))
REPS <- 1:3
tight <- stats::glm.control(epsilon = 1e-14, maxit = 200)

## own inverse links
MYF <- list(logit = function(e) 1 / (1 + exp(-e)), probit = function(e) pnorm(e), cauchit = function(e) 0.5 + atan(e) / pi,
            stk_short = function(e) 1 / (1 + exp(-sign(e) * expm1(abs(e)))))
my_data <- function(ce, rep) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(ce$seed_base + rep); n <- ce$n
  if (ce$generator == "design") {
    x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
    y <- rbinom(n, 1, MYF[[ce$link]](ce$c0 + ce$s * (0.6 * x + 0.5 * d)))
    return(list(D = data.frame(x = x, d = d, y = y), f = y ~ x + d))
  }
  if (ce$generator == "sparse") {
    x <- rchisq(n, 4); x <- (x - mean(x)) / sd(x)
    y <- rbinom(n, 1, MYF[[ce$link]](ce$intercept + ce$slope * x))
    return(list(D = data.frame(x = x, y = y), f = y ~ x))
  }
  stop("generator")
}
rao <- function(fit, D, f, extra) {
  D2 <- D; nm <- paste0("zz", seq_len(ncol(extra))); for (j in seq_along(nm)) D2[[nm[j]]] <- extra[, j]
  f1 <- stats::update(f, stats::as.formula(paste(". ~ . +", paste(nm, collapse = " + "))))
  m1 <- suppressWarnings(stats::glm(f1, data = D2, family = stats::binomial(), control = tight))
  a <- stats::anova(fit, m1, test = "Rao")
  c(stat = a$Rao[2], df = a$Df[2], p = a[["Pr(>Chi)"]][2])
}
stk_exact <- function(fit, y) {                     # u'I^-1 u through weighted QR residuals
  X <- stats::model.matrix(fit); eta <- fit$linear.predictors; p <- fit$fitted.values; W <- p * (1 - p)
  Z <- cbind(0.5 * eta^2 * (eta >= 0), -0.5 * eta^2 * (eta < 0)); Z <- Z[, colSums(Z != 0) > 0, drop = FALSE]
  u <- crossprod(Z, y - p)
  Rz <- qr.resid(qr(sqrt(W) * X), sqrt(W) * Z)
  chi <- drop(crossprod(u, solve(crossprod(Rz), u)))
  c(stat = chi, df = ncol(Z), p = stats::pchisq(chi, ncol(Z), lower.tail = FALSE))
}
my_groups <- function(ph, G) { n <- length(ph); r <- integer(n); r[order(ph)] <- seq_len(n); pmin(ceiling(r * G / n), G) }
raw_basis <- function(pbar, b) {
  e <- log(pbar / (1 - pbar))
  Z <- switch(b, poly2 = cbind(pbar, pbar^2), poly3 = cbind(pbar, pbar^2, pbar^3),
              stk = cbind(e, e^2 * (e >= 0), -e^2 * (e < 0)), sym = cbind(e * abs(e)))
  Z[, colSums(abs(Z)) > 0, drop = FALSE]
}

ROWS <- list()
add <- function(ds, rep, test, h, ref, kind, note = "") {
  ROWS[[length(ROWS) + 1]] <<- data.frame(dataset = ds, rep = rep, test = test, harness = as.numeric(h), reference = as.numeric(ref),
    abs_diff = abs(as.numeric(h) - as.numeric(ref)), kind = kind, note = note, stringsAsFactors = FALSE)
}

for (pk in PICK) {
  ce <- as.list(Cells[Cells$block == pk[1] & Cells$cell == pk[2], ])
  stopifnot(length(ce$cell) == 1)
  cat(sprintf("== %s %s  generator %s link %s n %d seed_base %.0f\n", pk[1], pk[2], ce$generator, ce$link, ce$n, ce$seed_base))
  for (rep in REPS) {
    h <- battery_one(rep, ce)
    RNGkind("L'Ecuyer-CMRG"); set.seed(ce$seed_base + rep); dh <- bt_data(ce)
    md <- my_data(ce, rep); D <- md$D; f <- md$f; y <- D$y; n <- nrow(D)
    same_y <- identical(as.numeric(dh$d$y), as.numeric(y)); same_x <- max(abs(dh$d$x - D$x))
    ## fits: default (as the harness), tight, and tight warm-restarted
    fit_def <- suppressWarnings(stats::glm(f, data = D, family = stats::binomial()))
    fit0 <- suppressWarnings(stats::glm(f, data = D, family = stats::binomial(), control = tight))
    fitT <- suppressWarnings(stats::glm(f, data = D, family = stats::binomial(), control = tight, start = stats::coef(fit0)))
    ev <- sum(y); etaD <- fit_def$linear.predictors; etaT <- fitT$linear.predictors
    cat(sprintf("  rep %d: events %d of %d; own data identical y %s, max|dx| %.1e; default fit converged %s in %d it\n",
                rep, ev, n, same_y, same_x, fit_def$converged, fit_def$iter))
    add(pk[2], rep, "data.y_identical", as.numeric(same_y), 1, "identity")

    ## flags
    arms <- c(G10 = 10, Grule = max(10, round(n / 25)))
    for (a in names(arms)) add(pk[2], rep, paste0("flag.evlt.", a), h[[paste0("flag.evlt.", a)]], as.numeric(min(ev, n - ev) < arms[[a]]), "identity")
    add(pk[2], rep, "flag.stk_half0", h[["flag.stk_half0"]], as.numeric(!(any(etaD >= 0) && any(etaD < 0))), "identity")
    sd0 <- stuk_diag(fit_def)
    add(pk[2], rep, "flag.stk_diag", h[["flag.stk_diag"]], sd0[["fail"]], "identity", "stuk_diag() of _dgp_library.R")

    ## Stukel joint
    ex <- stk_exact(fitT, y); exD <- stk_exact(fit_def, y)
    Zs <- cbind(0.5 * etaT^2 * (etaT >= 0), -0.5 * etaT^2 * (etaT < 0)); Zs <- Zs[, colSums(Zs != 0) > 0, drop = FALSE]
    ra <- rao(fitT, D, f, Zs)
    add(pk[2], rep, "Stk.joint vs own u'I^-1u (default fit)", h[["Stk.joint"]], exD[["p"]], "same fit", sprintf("df %d", exD[["df"]]))
    add(pk[2], rep, "Stk.joint vs own u'I^-1u (tight fit)", h[["Stk.joint"]], ex[["p"]], "tight fit", sprintf("df %d", ex[["df"]]))
    add(pk[2], rep, "Stk.joint vs anova Rao (warm restart)", h[["Stk.joint"]], ra[["p"]], "tight fit", sprintf("df %d, stat own %.8f anova %.8f", ra[["df"]], ex[["stat"]], ra[["stat"]]))
    add(pk[2], rep, "Stk.joint own vs anova Rao (statistic)", ex[["stat"]], ra[["stat"]], "tight fit")

    ## Stukel LR
    ZD <- cbind(za = 0.5 * etaD^2 * (etaD >= 0), zb = -0.5 * etaD^2 * (etaD < 0))
    D3 <- D; D3$za <- ZD[, 1]; D3$zb <- ZD[, 2]
    mlr <- suppressWarnings(stats::glm(stats::update(f, . ~ . + za + zb), data = D3, family = stats::binomial()))
    k_lr <- sum(!is.na(stats::coef(mlr))) - sum(!is.na(stats::coef(fit_def)))
    plr <- if (isTRUE(mlr$converged)) stats::pchisq(max(stats::deviance(fit_def) - stats::deviance(mlr), 0), k_lr, lower.tail = FALSE) else NA
    add(pk[2], rep, "Stk.LR vs deviance drop", h[["Stk.LR"]], plr, "same fit", sprintf("df %d, refit converged %s", k_lr, mlr$converged))

    ## Stukel sym1
    rs <- rao(fitT, D, f, cbind(etaT * abs(etaT)))
    add(pk[2], rep, "Stk.sym1 vs one-column anova Rao", h[["Stk.sym1"]], rs[["p"]], "tight fit")

    ## Stukel marginal sum
    zm <- tryCatch(statmod::glm.scoretest(fit_def, ZD), error = function(e) c(NA, NA))
    pm <- stats::pchisq(sum(zm^2), 2, lower.tail = FALSE)
    add(pk[2], rep, "Stk.marg vs statmod::glm.scoretest", h[["Stk.marg"]], if (is.finite(pm)) pm else NA, "same fit")

    ## EDGE: grouping on the default fit's clamped probabilities, as the harness
    phD <- pmin(pmax(fit_def$fitted.values, 1e-6), 1 - 1e-6)
    for (a in names(arms)) {
      G <- arms[[a]]
      g <- my_groups(phD, G)
      fqh <- bt_fit(list(d = D, f = f)); gsh <- bt_groups(fqh, G)
      add(pk[2], rep, sprintf("grouping %s (G=%d) identical", a, G), as.numeric(identical(as.integer(gsh$grp), as.integer(g))), 1, "identity")
      pbar <- as.numeric(tapply(phD, g, mean))
      gk <- envK$grp_stats(fit_def, G)
      for (b in BT_BASES) {
        ## score form vs Rao
        Zr <- raw_basis(pbar, b)
        rr <- rao(fitT, D, f, Zr[g, , drop = FALSE])
        add(pk[2], rep, sprintf("EDGE.%s.sc.%s vs Rao step covariates", b, a), h[[sprintf("EDGE.%s.sc.%s", b, a)]], rr[["p"]], "tight fit",
            sprintf("G=%d df Rao %d", G, rr[["df"]]))
        ## unit form vs edge_eig
        Zk <- switch(b, poly3 = as.matrix(poly(gk$pb, 3)), poly2 = as.matrix(poly(gk$pb, 2)),
                     stk = cbind(gk$eb, gk$eb^2 * (gk$eb >= 0), -gk$eb^2 * (gk$eb < 0)), sym = cbind(gk$eb * abs(gk$eb)))
        add(pk[2], rep, sprintf("EDGE.%s.u.%s vs edge_eig", b, a), h[[sprintf("EDGE.%s.u.%s", b, a)]], envK$edge_eig(gk, Zk), "same fit",
            sprintf("G=%d", G))
      }
    }

    ## GiViTI, called twice
    pr <- fit_def$fitted.values
    for (th in c(0.95, 0.50)) {
      g1 <- suppressWarnings(givitiR::givitiCalibrationTest(y, pr, devel = "internal", thres = th))$p.value
      g2 <- suppressWarnings(givitiR::givitiCalibrationTest(y, pr, devel = "internal", thres = th))$p.value
      nm <- if (th == 0.95) "GiViTI" else "GiViTI.t50"
      add(pk[2], rep, sprintf("%s vs direct call", nm), h[[nm]], g1, "same fit", sprintf("second call differs by %.1e", abs(g1 - g2)))
    }

    ## Cubic LR
    e1 <- etaD
    m1 <- suppressWarnings(stats::glm(y ~ e1, family = stats::binomial()))
    m3 <- suppressWarnings(stats::glm(y ~ e1 + I(e1^2) + I(e1^3), family = stats::binomial()))
    a13 <- stats::anova(m1, m3, test = "LRT")
    add(pk[2], rep, "Cubic.LR vs anova LRT", h[["Cubic.LR"]], a13[["Pr(>Chi)"]][2], "same fit", sprintf("df %d, cubic converged %s", a13$Df[2], m3$converged))
  }
}

R <- do.call(rbind, ROWS)
fwrite(R, file.path(OUT, "review_tests.csv"))
cat("\nlargest absolute differences by test family and comparison kind\n")
R$family <- sub(" .*$", "", R$test); R$family <- sub("\\.(G10|Grule)$", "", R$family)
S <- as.data.table(R)[, .(checks = .N, max_abs_diff = if (all(is.na(abs_diff))) NA_real_ else max(abs_diff, na.rm = TRUE),
                          na_mismatch = sum(is.na(harness) != is.na(reference))), by = .(test = sub("\\.(G10|Grule)", ".<arm>", test), kind)]
print(S[order(kind, -max_abs_diff)], row.names = FALSE)
cat("\nrows with |diff| > 1e-6 or an NA mismatch:\n")
bad <- R[(is.finite(R$abs_diff) & R$abs_diff > 1e-6) | (is.na(R$harness) != is.na(R$reference)), ]
if (nrow(bad)) print(bad[, c("dataset", "rep", "test", "harness", "reference", "abs_diff", "note")], row.names = FALSE) else cat("  none\n")
sink()
