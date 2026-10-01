## battery_selftest.R -- checks of the battery code; 1-15 need no ebrahim.gof, 16 compares with the installed package (the full
## comparison is block 0).
##   1. Stukel joint score = anova(..., test = "Rao") on 5 designs (one of them with a half-column identically zero)
##   2. EDGE score form = the Rao score test for the grouped shapes added as step covariates
##   3. harness EDGE unit form = edge_eig() of run_K_basis_score.R on 3 data sets
##   4. stukel_h: a = 0 is the logistic, a1 = a2 gives h(-eta) = -h(eta), h is increasing
##   5. the Hosmer et al. (1997) Table V shapes give the tails named in E0.1 (long, short, long-short)
##   6. seed ranges of the result cells are disjoint, apart from overlaps inherited from the July grids; the row-order seeds
##      of the real fits lie outside every range
##   7. every E3 row is in the cell table (cells per block printed)
##   8. rebuilt old cell ids match the seeds stored in the July p-value files
##   9. ten replicates of the three block 0 identity cells reproduce the stored July p-values
##  10. EDGE column rule on the package's R2 construction: the tiny Stukel half-column is kept and scaled (E8.3), the unit
##      form gives the three-column value and the score form has 3 df
##  11. samples with no event or no non-event: no test returns a p-value, and the three score statistics decline on the fit
##  12. coefficients of the design nulls added before launch = pi* of the working model; event rates of their generators
##  13. the 24 projection-grid cells reproduce their stored seeds and EDGE-poly3 p-values
##  14. the launcher (launch_plan) follows the E4 order
##  15. the new cell types run (cloglog refit, external G list, design nulls, row orders of a real fit)
##  16. harness = installed ebrahim.gof on the R2 construction and on 6 samples with no event or no non-event
## Writes battery/selftest.log; exits non-zero if any check fails.
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
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))
dir.create(edge_battery(), showWarnings = FALSE)
sink(edge_battery("selftest.log"), split = TRUE)
cat("battery_selftest.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
RNGkind("L'Ecuyer-CMRG")
RES <- list()
check <- function(label, ok, detail = "") {
  RES[[length(RES) + 1]] <<- data.frame(check = label, ok = isTRUE(ok), detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("  [%s] %s %s\n", if (isTRUE(ok)) "ok" else "FAIL", label, detail))
}
tight <- stats::glm.control(epsilon = 1e-14, maxit = 100)
fit_tight <- function(dat) {
  ## anova's Rao statistic uses the working weights of the fit's last iteration, which are those of the previous
  ## iterate; refitting from the MLE makes that iterate the MLE, so the reference is exact
  fit0 <- suppressWarnings(stats::glm(dat$f, data = dat$d, family = stats::binomial(), control = tight))
  fit <- suppressWarnings(stats::glm(dat$f, data = dat$d, family = stats::binomial(), control = tight, start = stats::coef(fit0)))
  eta <- as.numeric(fit$linear.predictors); ph <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6)
  dmu <- fit$family$mu.eta(eta); X <- stats::model.matrix(fit); w <- dmu^2 / (ph * (1 - ph))
  list(fit = fit, y = as.numeric(fit$y), eta = eta, ph = ph, p_raw = as.numeric(fitted(fit)), dmu = dmu, X = X,
       A = crossprod(X, w * X), n = length(fit$y))
}
rao <- function(fit, dat, extra) {                  # anova Rao score for adding the columns of 'extra'
  d2 <- dat$d; nm <- paste0(".z", seq_len(ncol(extra))); d2[nm] <- as.data.frame(extra)
  f1 <- suppressWarnings(stats::glm(stats::update(dat$f, stats::as.formula(paste(". ~ . +", paste(nm, collapse = " + ")))),
                                    data = d2, family = stats::binomial(), control = tight))
  a <- stats::anova(fit, f1, test = "Rao")
  c(stat = a$Rao[2], df = a$Df[2])
}
mk <- function(generator, n, ...) c(list(generator = generator, n = n, type = "full", ao = FALSE, G_extra = ""), list(...))
designs <- list(
  base      = mk("design", 1000L, link = "logit", s = 1, c0 = 0, xdist = "uniform"),
  sparse    = mk("sparse", 500L, link = "logit", intercept = -4.9, slope = 1),
  auc_probit = mk("design", 2000L, link = "probit", s = 2, c0 = 0, xdist = "uniform"),
  e12_cauchit = mk("design", 3000L, link = "cauchit", s = 1, c0 = -3.15, xdist = "uniform"),
  skew_cauchit = mk("design", 1500L, link = "cauchit", s = 1, c0 = 0, xdist = "skewed"))

## ---- 1 ----
cat("1. Stukel joint score = anova Rao\n")
k <- 0L
for (nm in names(designs)) {
  k <- k + 1L; set.seed(9100 + k); dat <- bt_data(designs[[nm]]); fq <- fit_tight(dat)
  js <- bt_stukel_joint_stat(fq)
  ref <- rao(fq$fit, dat, js$Z[, js$keep, drop = FALSE])
  rel <- abs(js$chi - ref[["stat"]]) / max(1, ref[["stat"]])
  check(sprintf("joint score, %-12s", nm), rel < 1e-7 && js$df == ref[["df"]],
        sprintf("harness %.8f (%d df) anova Rao %.8f (%d df)  half-column zero: %s", js$chi, js$df, ref[["stat"]], ref[["df"]], js$half0 == 1))
}
check("a design with one half-column identically zero was exercised", any(vapply(names(designs), function(nm) {
  set.seed(9100 + match(nm, names(designs))); bt_stukel_joint_stat(fit_tight(bt_data(designs[[nm]])))$half0 == 1 }, logical(1))))

## ---- 2 ----
cat("\n2. EDGE score form = Rao score for grouped step covariates\n")
for (nm in c("base", "auc_probit", "sparse")) {
  set.seed(9200 + match(nm, names(designs))); dat <- bt_data(designs[[nm]]); fq <- fit_tight(dat)
  for (G in unique(c(10, bt_rule_G(fq$n)))) {
    gs <- bt_groups(fq, G)
    for (b in BT_BASES) {
      Z <- bt_basis(gs$pbar, b); sc <- bt_edge_score(gs, fq$A, Z)
      ref <- rao(fq$fit, dat, Z[gs$grp, , drop = FALSE])
      rel <- abs(sc$S - ref[["stat"]]) / max(1, ref[["stat"]])
      check(sprintf("score form, %-10s G=%-3d %-5s", nm, G, b), rel < 1e-7 && sc$k == ref[["df"]],
            sprintf("harness %.6f (%d df) anova Rao %.6f (%d df)", sc$S, sc$k, ref[["stat"]], ref[["df"]]))
    }
  }
}

## ---- 3 ----
cat("\n3. harness EDGE unit form = edge_eig() of run_K_basis_score.R\n")
ex <- parse(file.path(SIMDIR, "run_K_basis_score.R"))                 # only the two definitions are evaluated
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
                  as.character(e[[2]]) %in% c("grp_stats", "edge_eig")) eval(e, globalenv())
for (j in 1:3) {
  ce <- list(designs$base, designs$e12_cauchit, designs$skew_cauchit)[[j]]
  set.seed(9300 + j); dat <- bt_data(ce); fq <- bt_fit(dat)
  for (G in unique(c(10, bt_rule_G(fq$n)))) {
    gk <- grp_stats(fq$fit, G); gs <- bt_groups(fq, G)
    ref3 <- edge_eig(gk, as.matrix(poly(gk$pb, 3))); refs <- edge_eig(gk, cbind(gk$eb, gk$eb^2 * (gk$eb >= 0), -gk$eb^2 * (gk$eb < 0)))
    h3 <- bt_edge_unit(gs, fq$A, bt_basis(gs$pbar, "poly3"))$p; hs <- bt_edge_unit(gs, fq$A, bt_basis(gs$pbar, "stk"))$p
    check(sprintf("unit form, data set %d G=%d", j, G), abs(h3 - ref3) < 1e-10 && abs(hs - refs) < 1e-10,
          sprintf("poly3 %.10f vs %.10f; stk %.10f vs %.10f", h3, ref3, hs, refs))
  }
}

## ---- 4 ----
cat("\n4. stukel_h\n")
eg <- seq(-12, 12, by = 0.01)
check("a1 = a2 = 0 gives the logistic (h = eta)", identical(stukel_h(eg, 0, 0), eg))
sym_ok <- all(vapply(c(-2, -1, -0.5, 0.5, 1, 2), function(a) max(abs(stukel_h(-eg, a, a) + stukel_h(eg, a, a))) < 1e-10, logical(1)))
check("a1 = a2 gives h(-eta) = -h(eta) for a in {-2, -1, -0.5, 0.5, 1, 2}", sym_ok)
mono_ok <- all(apply(expand.grid(a1 = c(-2, -1, 0, 1, 2), a2 = c(-2, -1, 0, 1, 2)), 1, function(a) all(diff(stukel_h(eg[abs(eg) < 6], a[1], a[2])) > 0)))
check("h increasing in eta for (a1, a2) in {-2..2}^2", mono_ok)
check("h continuous at 0 with slope 1 on both sides", all(abs(stukel_h(c(-1e-7, 1e-7), 1, -1) - c(-1e-7, 1e-7)) < 1e-12))

## ---- 5 ----
cat("\n5. Hosmer et al. (1997) Table V shapes\n")
tail_up <- function(lk, e) 1 - BT_LINKINV[[lk]](e); tail_lo <- function(lk, e) BT_LINKINV[[lk]](-e)
e1 <- c(3, 6); lg_up <- 1 - plogis(e1); lg_lo <- plogis(-e1)
long_ok  <- all(tail_up("stk_long", e1) > lg_up) && all(tail_lo("stk_long", e1) > lg_lo) &&
            diff(tail_up("stk_long", e1) / lg_up) > 0 && diff(tail_lo("stk_long", e1) / lg_lo) > 0
short_ok <- all(tail_up("stk_short", e1) < lg_up) && all(tail_lo("stk_short", e1) < lg_lo) &&
            diff(tail_up("stk_short", e1) / lg_up) < 0 && diff(tail_lo("stk_short", e1) / lg_lo) < 0
asym_ok  <- all(tail_up("stk_asym", e1) > lg_up) && all(tail_lo("stk_asym", e1) < lg_lo)
check("(-1, -1): both tails longer than logistic, ratio growing", long_ok,
      sprintf("P(y=0 | eta=6) %.2e vs logistic %.2e", tail_up("stk_long", 6), lg_up[2]))
check("(1, 1): both tails shorter than logistic, ratio shrinking", short_ok,
      sprintf("P(y=0 | eta=6) %.2e vs logistic %.2e", tail_up("stk_short", 6), lg_up[2]))
check("(-1, 1): long upper tail, short lower tail", asym_ok)
pg <- seq(0.6, 2.3, by = 0.1)
check("plateau_upper is flat above eta = 0.5 while stk_long keeps rising",
      max(abs(diff(BT_LINKINV$plateau_upper(pg)))) < 1e-12 && all(diff(BT_LINKINV$stk_long(pg)) > 0),
      sprintf("plateau p = %.3f", BT_LINKINV$plateau_upper(1)))
check("plateau names map to the old _dgp_library.R links (name mapping only; July identity is block 0 (i))", identical(BT_LINKINV$plateau_upper(eg), linkp("stukel_heavy", eg)) &&
      identical(BT_LINKINV$plateau_lower(eg), linkp("stukel_light", eg)) && identical(BT_LINKINV$plateau_both(eg), linkp("stukel_asym", eg)))

## ---- 6 ----
cat("\n6. seed ranges\n")
Cells <- battery_cells()
O <- seed_overlaps(Cells)
new <- Cells[Cells$seed_family == "new", ]
check("no overlap that involves a new cell", sum(O$kind == "conflict") == 0, sprintf("%d conflicting pairs", sum(O$kind == "conflict")))
check("new seed bases within integer range, B <= 10,000 and ids <= 999", all(new$seed_base + new$B < .Machine$integer.max) &&
      all(new$B <= 10000) && all(new$cell_id <= 999))
check("no two cells share a seed base", !anyDuplicated(Cells$seed_base[!is.na(Cells$seed_base) & Cells$block != "0"]))
ro <- BT_ORDER_SEED + c(1, max(Cells$orders, na.rm = TRUE)); sb <- Cells[!is.na(Cells$seed_base), ]
check("row-order seeds of the real fits (E8.5) lie outside every cell's seed range", !any(sb$seed_base + 1 <= ro[2] & sb$seed_base + sb$B >= ro[1]),
      sprintf("orders seeded %.0f-%.0f", ro[1], ro[2]))
cat(sprintf("  inherited from the July grids: %d pairs of old cells share %d seed values (same id in two grids, offset 2-15)\n",
            nrow(O), sum(O$shared_seeds)))
if (nrow(O)) print(utils::head(O[order(-O$shared_seeds), c("cell_a", "cell_b", "offset", "shared_seeds")], 8), row.names = FALSE)
fwrite(O, edge_battery("seed_overlaps.csv"))

## ---- 7 ----
cat("\n7. E3 rows in the cell table\n")
E <- cells_expected(Cells)
for (i in seq_len(nrow(E))) check(E$check[i], E$ok[i])
cnt <- as.data.frame(table(block = Cells$block)); cnt$replicates <- as.numeric(tapply(Cells$B, Cells$block, sum)[as.character(cnt$block)])
print(cnt, row.names = FALSE); cat(sprintf("  total: %d cells, %s replicates\n", nrow(Cells), format(sum(Cells$B), big.mark = ",")))

## ---- 8 ----
cat("\n8. old cell ids against the stored July seeds\n")
V <- verify_old_ids(Cells)
fwrite(V, edge_battery("old_id_check.csv"))
check("every reused cell: stored seed - rep is the table's seed base for all replicates, stored B = table B", all(V$ok),
      sprintf("%d of %d cells", sum(V$ok), nrow(V)))

## ---- 9 ----
cat("\n9. block 0 identity cells, 10 replicates, against the stored July p-values\n")
b0 <- Cells[Cells$block == "0" & !is.na(Cells$stored), ]
for (i in seq_len(nrow(b0))) {
  ce <- as.list(b0[i, ])
  H <- as.data.frame(do.call(rbind, lapply(1:10, battery_one, cell = ce)))
  if (ce$stored == "sim_null") {
    S <- fread(file.path(SIMDIR, "sim_null_pvalues.csv"))[family == ce$family & n == ce$n & G == 10L & rep <= 10L][order(rep)]
    setnames(S, sub("^p\\.", "", names(S)))
  } else {
    S <- if (ce$stored == "sim_power_broad")
      fread(file.path(SIMDIR, "sim_power_broad_pvalues.csv"), select = c("family", "param", "n", "rep", "seed", "test", "p_value"),
            colClasses = list(character = "param"))[family == ce$family & param == ce$param & n == ce$n & rep <= 10L]
    else fread(file.path(SIMDIR, "sim_edge_loses_pvalues.csv"))[scenario == "crossover" & n == ce$n & rep <= 10L]
    S <- dcast(S, rep + seed ~ test, value.var = "p_value")[order(rep)]
  }
  pairs <- rbind(c("EDGE.poly2.u.G10", "DEF.poly2"), c("EDGE.poly3.u.G10", "DEF.poly3"), c("EDGE.stk.u.G10", "DEF.stukel"),
                 c("HL.G10", "HL"), c("HLF.G10", "EF"), c("HL_w", "HL-equalwidth"), c("PH", "Pigeon-Heyse"),
                 c("Tsiatis", "Tsiatis"), c("Xie", "Xie"), c("Stk.marg", "Stukel"))
  md <- vapply(seq_len(nrow(pairs)), function(j) max(abs(H[[pairs[j, 1]]] - S[[pairs[j, 2]]])), numeric(1))
  check(sprintf("%-22s seeds and 10 p-value columns", ce$cell), all(H$seed == S$seed) && all(md < 1e-6),
        sprintf("max|diff| %.1e", max(md)))
}

## ---- 10 ----
cat("\n10. EDGE column rule on the package's R2 construction: a Stukel half-column that is tiny but not zero\n")
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(4)        # as tests/testthat/test-def-gof.R of ebrahim.gof 2.8.0
ph <- c(sort(runif(475, 0.02, 0.45)), rep(0.5001, 25)); Xk <- cbind(1, qlogis(ph)); yk <- rbinom(500, 1, ph)
RNGkind("L'Ecuyer-CMRG")
fk <- list(n = 500L, ph = ph, y = yk, dmu = ph * (1 - ph), X = Xk, A = crossprod(Xk, ph * (1 - ph) * Xk))
gk <- bt_groups(fk, 20L)
ek <- qlogis(gk$pbar)
Zfull <- cbind(ek, ek^2 * (ek >= 0), -ek^2 * (ek < 0))
Zh <- bt_basis(gk$pbar, "stk")
nzk <- sqrt(colSums(Zfull^2))
check("stk basis, G = 20: the tiny half-column is kept (sum |z| > 1e-8) and every column has unit length",
      ncol(Zh) == 3 && all(abs(sqrt(colSums(Zh^2)) - 1) < 1e-12),
      sprintf("half-column length %.2e, longest %.2f, sum |z| %.2e", nzk[2], max(nzk), sum(abs(Zfull[, 2]))))
## the three-column value, built as in the package test: tapply group sums, rowsum U, Omega from solve, columns scaled
g3 <- pmin(ceiling(rank(ph, ties.method = "first") / 25), 20)
V3 <- ph * (1 - ph); Vg3 <- as.numeric(tapply(V3, g3, sum)); r3 <- as.numeric(tapply(yk - ph, g3, sum)) / sqrt(Vg3)
e3 <- qlogis(as.numeric(tapply(ph, g3, mean)))
Z3 <- cbind(e3, e3^2 * (e3 >= 0), -e3^2 * (e3 < 0)); Z3 <- sweep(Z3, 2, sqrt(colSums(Z3^2)), "/")
U3 <- rowsum(V3 * Xk, g3) / sqrt(Vg3); Om3 <- diag(20) - U3 %*% solve(crossprod(Xk, V3 * Xk), t(U3))
S3 <- drop(crossprod(r3, Z3 %*% solve(crossprod(Z3), crossprod(Z3, r3))))
l3 <- Re(eigen(solve(crossprod(Z3), crossprod(Z3, Om3 %*% Z3)), only.values = TRUE)$values); l3 <- l3[l3 > 1e-9]
p3 <- pchisq(S3 / (sum(l3^2) / sum(l3)), sum(l3)^2 / sum(l3^2), lower.tail = FALSE)
## score form by an independent route: u'I^-1 u by a pivoted QR solve on the weighted, scaled columns
Zw <- (Zfull / rep(nzk, each = nrow(Zfull))) * sqrt(gk$Vg); uq <- drop(crossprod(Zw, gk$r)); Iq <- crossprod(Zw, gk$Omega %*% Zw)
qi <- qr(Iq, tol = 1e-12); Ssq <- sum(uq * qr.coef(qi, uq), na.rm = TRUE)
hu <- bt_edge_unit(gk, fk$A, Zh); hs <- bt_edge_score(gk, fk$A, Zh)
ou <- tryCatch(bt_edge_unit(gk, fk$A, Zfull)$p, error = function(e) NA_real_)
r2 <- list(ph = ph, X = Xk, y = yk, hu = hu, hs = hs)                  # for check 16
check("unit form = the three-column value of the package test (S and p; p 0.4169424779 as installed)",
      isTRUE(abs(hu$p - p3) < 1e-10) && isTRUE(abs(hu$S - S3) < 1e-10 * max(1, S3)) && isTRUE(abs(hu$p - 0.4169424779) < 1e-10),
      sprintf("p %.10f vs %.10f, S %.10f vs %.10f (unscaled columns: %s)", hu$p, p3, hu$S, S3,
              if (is.na(ou)) "the solve stops" else paste("p", format(ou, digits = 6))))
check("score form with all three columns = the QR route, on 3 df",
      isTRUE(abs(hs$S - Ssq) < 1e-8 * max(1, Ssq)) && hs$k == 3L && qi$rank == 3L,
      sprintf("S %.6f vs %.6f on %d df (information rank %d)", hs$S, Ssq, hs$k, qi$rank))

## ---- 11 ----
cat("\n11. samples with no event or no non-event\n")
ce_s <- mk("sparse", 200L, link = "logit", intercept = -4.9, slope = 1)
for (yv in c(0, 1)) {
  set.seed(9400 + yv); xz <- as.numeric(scale(rchisq(200, 4)))
  dz <- list(d = data.frame(x = xz, y = rep(yv, 200)), f = y ~ x)
  h <- battery_rep(dz, ce_s)
  tcz <- names(h)[!grepl("^(flag\\.|n$|events$|fit_ok$|glm_conv$)", names(h))]
  check(sprintf("all y = %d: no test returns a p-value, flag.degenerate = 1", yv), all(!is.finite(h[tcz])) && h[["flag.degenerate"]] == 1,
        sprintf("%d test columns", length(tcz)))
  fz <- bt_fit(dz); jz <- bt_stukel_joint_stat(fz); sz <- bt_stukel(fz); gz <- bt_groups(fz, 10)
  scz <- vapply(c("stk", "sym"), function(b) { Z <- bt_basis(gz$pbar, b); if (is.null(Z)) NA_real_ else bt_edge_score(gz, fz$A, Z)$p }, numeric(1))
  check(sprintf("all y = %d, separated fit: Stk.joint, Stk.sym1, EDGE-stk and EDGE-sym score forms decline (were p = 0)", yv),
        !is.finite(jz$chi) && !is.finite(sz[["Stk.sym1"]]) && all(!is.finite(scz)),
        sprintf("joint %s, sym1 %s, stk score %s, sym score %s", jz$chi, sz[["Stk.sym1"]], scz[1], scz[2]))
}
set.seed(9402); x1 <- as.numeric(scale(rchisq(200, 4))); y1 <- rep(0, 200); y1[which.max(x1)] <- 1
h1 <- battery_rep(list(d = data.frame(x = x1, y = y1), f = y ~ x), ce_s)
check("one event: not degenerate, HL returns a p-value", h1[["flag.degenerate"]] == 0 && is.finite(h1[["HL.G10"]]))

## ---- 12 ----
cat("\n12. design nulls added before launch: coefficients = pi* of the working model, event rates of the generators\n")
set.seed(20260913L); Np <- 400000L
pis <- function(X, p) suppressWarnings(glm.fit(X, p, family = quasibinomial()))$coefficients
REF <- list(); EV <- list()
x <- runif(Np, -2.5, 2.5); d <- rbinom(Np, 1, .5); p <- plogis(0.9 * x - 0.5 * x * d); REF$crossover <- pis(cbind(1, x, d), p); EV$crossover <- mean(p)
x <- runif(Np, 0.3, 6); p <- plogis(-1 + 1.5 * log(x)); REF$logx <- pis(cbind(1, x), p); EV$logx <- mean(p)
d1 <- rbinom(Np, 1, .5); d2 <- rbinom(Np, 1, .5); p <- plogis(-0.6 + 0.6 * d1 + 1.0 * d2 + 1.4 * d1 * d2)
REF$int_binbin <- pis(cbind(1, d1, d2), p); EV$int_binbin <- mean(p)
x <- rchisq(Np, 4); xs <- (x - 4) / sqrt(8); p <- plogis(-0.3 + 0.7 * xs + 0.3 * (xs^2 - 1)); REF$skew <- pis(cbind(1, x), p); EV$skew <- mean(p)
x <- runif(Np, -2.5, 2.5); d <- rbinom(Np, 1, .5); p <- plogis(-0.2 + 0.8 * x + 1.0 * d + 0.8 * x * d); REF$joint <- pis(cbind(1, x), p); EV$joint <- mean(p)
x1 <- rnorm(Np); x2 <- 0.5 * x1 + sqrt(0.75) * rnorm(Np); p <- plogis(-0.2 + 0.7 * x1 + 0.7 * x2 + 0.6 * x1 * x2)
REF$corr <- pis(cbind(1, x1, x2), p); EV$corr <- mean(p)
x <- runif(Np, -2.5, 2.5); p <- plogis(0.3 + 0.7 * x + 0.2 * (x^2 - 2.083) + 0.09 * x^3); REF$omit_x2x3 <- pis(cbind(1, x), p); EV$omit_x2x3 <- mean(p)
x <- runif(Np, -2.5, 2.5); z1 <- rnorm(Np); z2 <- rnorm(Np); p <- plogis(0.2 + 0.6 * x + 0.8 * z1 + 0.8 * z2)
REF$omit_2cov <- pis(cbind(1, x), p); EV$omit_2cov <- mean(p)
x <- runif(Np, -2.5, 2.5); d <- rbinom(Np, 1, .5); z <- rnorm(Np); p <- plogis(0.2 + 0.6 * x + 0.5 * d + 0.4 * z + 0.7 * x * d + 0.7 * x * z)
REF$omit_2int <- pis(cbind(1, x, d, z), p); EV$omit_2int <- mean(p)
for (sc in names(REF)) {
  yb <- bt_gen_design_null(sc, 200000L)$d$y
  check(sprintf("%-10s coefficients within 0.015 of pi*, generator event rate within 0.006 of the alternative's", sc),
        length(REF[[sc]]) == length(BT_NULL_COEF[[sc]]) && max(abs(REF[[sc]] - BT_NULL_COEF[[sc]])) <= 0.015 && abs(mean(yb) - EV[[sc]]) <= 0.006,
        sprintf("pi* %s, used %s; events %.4f vs %.4f", paste(sprintf("%.3f", REF[[sc]]), collapse = " "),
                paste(sprintf("%.2f", BT_NULL_COEF[[sc]]), collapse = " "), mean(yb), EV[[sc]]))
}
rm(x, d, p, d1, d2, xs, x1, x2, z1, z2, z)

## ---- 13 ----
cat("\n13. the 24 projection-grid cells, replicates 1-3, against proj_power_grid_pvalues.csv\n")
pj <- Cells[Cells$pair_proj, ]
SP <- fread(file.path(SIMDIR, "proj_power_grid_pvalues.csv"), colClasses = list(character = "param"))
mxp <- 0; seeds_ok <- TRUE; nam <- 0L
for (i in seq_len(nrow(pj))) {
  ce <- as.list(pj[i, ]); sk <- strsplit(ce$stored_key, "|", fixed = TRUE)[[1]]
  s <- SP[scenario == sk[1] & param == sk[2] & n == ce$n & rep <= 3][order(rep)]
  H <- do.call(rbind, lapply(1:3, battery_one, cell = ce))
  seeds_ok <- seeds_ok && nrow(s) == 3 && all(H[, "seed"] == s$seed)
  nam <- nam + sum(is.finite(H[, "EDGE.poly3.u.G10"]) != is.finite(s$p_edge3))
  dd <- abs(H[, "EDGE.poly3.u.G10"] - s$p_edge3); if (any(is.finite(dd))) mxp <- max(mxp, dd[is.finite(dd)])
}
check("24 cells: stored seeds and EDGE.poly3.u.G10 = stored p_edge3", nrow(pj) == 24 && seeds_ok && nam == 0 && mxp < 1e-6,
      sprintf("max|diff| %.1e, NA mismatches %d", mxp, nam))

## ---- 14 ----
cat("\n14. the launcher follows E4\n")
PL <- launch_plan(Cells); RL <- PL[PL$action == "run", ]
check("run order 0; 1b (block 2 nulls); 2; 1a; 1b; 3 (n = 1000); 3; 4; 5; 6; 7",
      identical(RL$block, c("0", "1b", "2", "1a", "1b", "3", "3", "4", "5", "6", "7")) && all(RL$cells[-c(2, 6)] == ""))
check("step 2 holds exactly the matched nulls of block 2", all(Cells$null_block[Cells$block == "2"] == "1b") &&
      setequal(strsplit(RL$cells[2], ",")[[1]], Cells$null_cell[Cells$block == "2"]))
check("step 6 holds exactly the n = 1000 cells of block 3", setequal(strsplit(RL$cells[6], ",")[[1]], Cells$cell[Cells$block == "3" & Cells$n %in% 1000L]))
check("then a summary of every block, without --cells", setequal(PL$block[PL$action == "summary"], setdiff(BC_BLOCKS, "0")) &&
      all(PL$cells[PL$action == "summary"] == "") && min(PL$step[PL$action == "summary"]) > max(RL$step))

## ---- 15 ----
cat("\n15. new cell types\n")
hb <- battery_one(0L, as.list(Cells[Cells$cell == "real_beetle_cloglog", ]))
check("real_beetle_cloglog: EDGE, Stukel and Cubic.LR not run; HL and HL_F returned",
      all(!is.finite(hb[grep("^(EDGE|Stk|Cubic)", names(hb))])) && is.finite(hb[["HL.G10"]]) && is.finite(hb[["HLF.G10"]]),
      sprintf("HL %.4f, HL_F %.4f", hb[["HL.G10"]], hb[["HLF.G10"]]))
ext <- as.list(Cells[Cells$cell == "runH_temp_s0.85_n500", ]); ext$G_extra <- "10,50"
he <- battery_one(1L, ext)
check("external cell with a G list: grouped tests and statistics at each G",
      all(c("EDGE.poly4.G50", "stat.HL.ext.G10", "HL.ext.G50") %in% names(he)) && is.finite(he[["EDGE.poly4.G50"]]) && is.finite(he[["stat.HL.ext.G10"]]))
for (cn in c("null_crossover_n1000", "null_census_int_binbin_n1000", "null_census_skew_n1000", "null_census_omit_2int_n1000")) {
  h1 <- battery_one(1L, as.list(Cells[Cells$block == "1b" & Cells$cell == cn, ]))
  check(sprintf("%s: one replicate without error, HL returned", cn), h1[["flag.rep_error"]] == 0 && is.finite(h1[["HL.G10"]]))
}
bl <- as.list(Cells[Cells$cell == "real_beetle_logit", ])
o0 <- battery_one(0L, bl); o1 <- battery_one(1L, bl); o1b <- battery_one(1L, bl)
oo <- do.call(rbind, lapply(2:6, battery_one, cell = bl))
hlo <- c(o1[["HL.G10"]], oo[, "HL.G10"])
check("real_beetle_logit row orders (E8.5): rep 0 unseeded, order k seeded 20260914 + k and reproducible, 8 distinct risks, ties split at G 10, HL moves",
      is.na(o0[["seed"]]) && o1[["seed"]] == 20260915 && identical(o1, o1b) && all(oo[, "seed"] == 20260914 + 2:6) &&
        o0[["info.distinct_risks"]] == 8 && o0[["info.tied_boundaries.G10"]] > 0 && o0[["flag.rep_error"]] == 0 &&
        length(unique(round(hlo, 12))) > 1,
      sprintf("HL G10 stored %.4f, orders 1-6 %s; tied boundaries at G 10 %d (HL_F blocks %d)", o0[["HL.G10"]],
              paste(sprintf("%.3f", hlo), collapse = " "), as.integer(o0[["info.tied_boundaries.G10"]]),
              as.integer(o0[["info.tied_boundaries_hlf.G10"]])))

## ---- 16 ----
cat("\n16. harness against the installed ebrahim.gof: the R2 construction and 6 samples with no event or no non-event\n")
pkv <- tryCatch(as.character(utils::packageVersion("ebrahim.gof")), error = function(e) NA_character_)
if (is.na(pkv) || utils::compareVersion(pkv, "2.8.0") < 0) {
  check("ebrahim.gof 2.8.0 or later is installed", FALSE, sprintf("installed: %s", pkv))
} else {
  pu <- ebrahim.gof::def.gof(r2$y, r2$ph, X = r2$X, G = 20, basis = "stukel")
  ps <- ebrahim.gof::def.gof(r2$y, r2$ph, X = r2$X, G = 20, basis = "stukel", weights = "score")
  check("R2: package unit S and p = harness; package score S = harness, both on 3 df (1e-10)",
        isTRUE(abs(pu$p_value - r2$hu$p) < 1e-10) && isTRUE(abs(pu$Test_Statistic - r2$hu$S) < 1e-10 * max(1, r2$hu$S)) &&
          isTRUE(abs(ps$Test_Statistic - r2$hs$S) < 1e-10 * max(1, r2$hs$S)) && isTRUE(ps$df == 3) && r2$hs$k == 3L,
        sprintf("ebrahim.gof %s: unit p %.10f (harness %.10f), score S %.8f on %s df (harness %.8f on %d)", pkv, pu$p_value, r2$hu$p,
                ps$Test_Statistic, format(ps$df), r2$hs$S, r2$hs$k))
  nsp <- asNamespace("ebrahim.gof"); ce_d <- mk("sparse", 200L, link = "logit", intercept = -4.9, slope = 1)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")                 # the six samples of block 0 (ii), drawn as the driver does
  for (k in 1:6) {
    set.seed(20260914L + k); xd <- as.numeric(scale(rchisq(200, 4))); yv <- if (k <= 3) 0 else 1
    dd <- data.frame(x = xd, y = rep(yv, 200))
    h <- battery_rep(list(d = dd, f = y ~ x), ce_d)
    tcd <- names(h)[!grepl("^(flag\\.|n$|events$|fit_ok$|glm_conv$)", names(h))]
    fd <- suppressWarnings(glm(y ~ x, data = dd, family = binomial()))
    cls <- character(0); pd <- c()
    for (b in c("poly3", "poly2", "stukel", "sym")) for (w in c("unit", "score"))
      pd[paste(b, w)] <- tryCatch(withCallingHandlers(as.numeric(ebrahim.gof::def.gof(fd, G = 10, basis = b, weights = w)$p_value),
        warning = function(x) { cls <<- c(cls, class(x)[1]); invokeRestart("muffleWarning") }), error = function(e) -1)
    ctx <- get(".gof_context", envir = nsp)(fd, G = 10)
    pj <- vapply(c("joint", "lr", "marginal"), function(f) tryCatch(as.numeric(suppressWarnings(
      get("gof_stukel", envir = nsp)(ctx, list(form = f))$p_value)), error = function(e) -1), numeric(1))
    check(sprintf("all y = %d, set %d: no harness p-value; package def.gof NA with def_degenerate (8), gof_stukel NA (3)", yv, k),
          all(!is.finite(h[tcd])) && h[["flag.degenerate"]] == 1 && all(is.na(pd)) && sum(cls == "def_degenerate") == 8 &&
            all(cls == "def_degenerate") && all(is.na(pj)),
          sprintf("%d harness columns; def.gof NA %d of 8, warnings %d; gof_stukel NA %d of 3", length(tcd), sum(is.na(pd)), length(cls), sum(is.na(pj))))
  }
  RNGkind("L'Ecuyer-CMRG")
}

R <- do.call(rbind, RES)
cat(sprintf("\n%d checks, %d failed\n", nrow(R), sum(!R$ok)))
if (any(!R$ok)) print(R[!R$ok, ], row.names = FALSE)
cat(if (all(R$ok)) "SELFTEST PASSED\n" else "SELFTEST FAILED\n")
sink()
quit(save = "no", status = if (all(R$ok)) 0L else 1L)
