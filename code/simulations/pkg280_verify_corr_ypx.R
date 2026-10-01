## pkg280_verify_corr_ypx.R -- the (y, predicted_probs, X) path of def.gof after the column-rule correction, and the R2
## construction of the correction report.
## 1. R2: set.seed(4); 475 risks U(0.02, 0.45) sorted, 25 at v; n = 500, G = 20, basis stukel, for v = 0.5 + 1e-6 ...
##    2e-3. Installed (new) against 0d63fba (pre): unit Satterthwaite, unit Imhof, score. The report gives S 1.2850859905,
##    df 1.811174, p 0.4169424779 at v = 0.5001, where pre stops.
## 2. 1000 new c0 = -2, n = 500 fits: def.gof(y, fitted, X = model.matrix) for every basis, unit (Satterthwaite, Imhof)
##    and score, G 10 and 20, installed against 0d63fba; equal = S, df and p within 1e-9 (relative, floor 1).
## Run: Rscript pkg280_verify_corr_ypx.R <clone dir> > ../paper_EDGE/theory/pkg280_verify_corr_ypx.log 2>&1

args <- commandArgs(trailingOnly = TRUE)
PKG <- args[1]
options(width = 200, warn = 1)
suppressPackageStartupMessages(library(ebrahim.gof))
ns <- asNamespace("ebrahim.gof")
PRE <- new.env(parent = ns)
for (ex in parse(text = system2("git", c("-C", shQuote(PKG), "show", "0d63fba:R/def_gof.R"), stdout = TRUE), keep.source = FALSE))
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]])) eval(ex, PRE)
cat("installed ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "| pre = 0d63fba |", format(Sys.time()), "\n")

val <- function(expr) {
  r <- tryCatch(suppressWarnings(expr), error = function(e) NULL)
  if (is.null(r)) c(S = NA_real_, df = NA_real_, p = NA_real_, err = 1) else c(S = r$Test_Statistic, df = r$df, p = r$p_value, err = 0)
}
show <- function(x) if (x["err"] == 1) "stops" else sprintf("S %.10f df %.6f p %.10f", x["S"], x["df"], x["p"])
eqv <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-9 * pmax(1, abs(a)))

cat("\n==== 1. R2 construction (n = 500, G = 20, stukel) ====\n")
for (v in 0.5 + c(1e-6, 5e-5, 1e-4, 1.7e-4, 5e-4, 1e-3, 2e-3)) {
  set.seed(4); ph <- c(sort(runif(475, 0.02, 0.45)), rep(v, 25)); X <- cbind(1, qlogis(ph)); yy <- rbinom(500, 1, ph)
  for (fm in c("sat", "imh", "sco")) {
    call_it <- function(f) switch(fm, sat = f(yy, ph, X = X, G = 20, basis = "stukel"),
                                      imh = f(yy, ph, X = X, G = 20, basis = "stukel", method = "imhof"),
                                      sco = f(yy, ph, X = X, G = 20, basis = "stukel", weights = "score"))
    a <- val(call_it(PRE$def.gof)); b <- val(call_it(def.gof))
    cat(sprintf("  v = %-9s %s  pre %-52s | new %s\n", format(v, digits = 7), fm, show(a), show(b)))
  }
}

cat("\n==== 2. (y, predicted_probs, X) path on 1000 c0 = -2 fits, installed against 0d63fba ====\n")
res <- list()
for (i in 1:1000) {
  set.seed(71500000 + i)
  dat <- data.frame(xa = runif(500, -3, 3), db = rbinom(500, 1, 0.5))
  dat$out <- rbinom(500, 1, plogis(-2 + 0.6 * dat$xa + 0.5 * dat$db))
  f <- suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat))
  y <- f$y; ph <- fitted(f); X <- model.matrix(f)
  for (G in c(10, 20)) for (b in c("poly2", "poly3", "stukel", "sym")) for (fm in c("sat", "imh", "sco")) {
    call_it <- function(g) switch(fm, sat = g(y, ph, X = X, G = G, basis = b),
                                      imh = g(y, ph, X = X, G = G, basis = b, method = "imhof"),
                                      sco = g(y, ph, X = X, G = G, basis = b, weights = "score"))
    a <- val(call_it(PRE$def.gof)); bb <- val(call_it(def.gof))
    res[[length(res) + 1]] <- data.frame(i = i, G = G, basis = b, form = fm, pre_S = a["S"], pre_df = a["df"], pre_p = a["p"],
                                         pre_err = a["err"], new_S = bb["S"], new_df = bb["df"], new_p = bb["p"], new_err = bb["err"])
  }
}
R <- do.call(rbind, res); rownames(R) <- NULL
R$eq <- eqv(R$pre_S, R$new_S) & eqv(R$pre_df, R$new_df) & eqv(R$pre_p, R$new_p)
R$imh_noise <- !R$eq & R$form == "imh" & eqv(R$pre_S, R$new_S) & eqv(R$pre_df, R$new_df) & abs(R$pre_p - R$new_p) <= 1e-6
tab <- aggregate(cbind(n = 1, pre_ran = pre_err == 0, equal = pre_err == 0 & eq, imhof_noise = pre_err == 0 & imh_noise,
                       differ = pre_err == 0 & !eq & !imh_noise, pre_stop = pre_err == 1,
                       stop_new_finite = pre_err == 1 & new_err == 0 & !is.na(new_p), new_err = new_err == 1) ~ G + basis + form,
                 data = R, FUN = sum)
print(tab, row.names = FALSE)
cat(sprintf("\nTOTAL: %d comparisons; pre ran %d; differ beyond Imhof noise %d; Imhof noise %d; pre stopped %d, new finite on %d; new errors %d\n",
            nrow(R), sum(R$pre_err == 0), sum(R$pre_err == 0 & !R$eq & !R$imh_noise), sum(R$pre_err == 0 & R$imh_noise),
            sum(R$pre_err == 1), sum(R$pre_err == 1 & R$new_err == 0 & !is.na(R$new_p)), sum(R$new_err == 1)))
if (any(R$pre_err == 1 | (!R$eq & !R$imh_noise))) print(R[R$pre_err == 1 | (!R$eq & !R$imh_noise), ], digits = 10, row.names = FALSE)
cat("\ndone", format(Sys.time()), "\n")
