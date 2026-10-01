## pkg280_verify_corr_colrule.R -- verification of the corrected column rule of def.gof (commit 18a3293).
## Three copies of def.gof run on the same fits:
##   new = the installed ebrahim.gof 2.8.0 (2.7.0 filter colSums(abs(Z)) > 1e-8, kept columns scaled, score guard)
##   pre = 0d63fba:R/def_gof.R, every top-level function sourced into its own environment (2.7.0 filter, no scaling)
##   rel = 906bd0c:R/def_gof.R the same way (the withdrawn relative 1e-6 rule), a positive control for the comparison
## Parts: A. 2000 fits at c0 = -2, n = 500 (G 10, 20) and 2000 base-design fits, n = 1000 (G 10, 20, 40), every basis
## (poly2, poly3, stukel, sym) in the unit form (Satterthwaite and Imhof) and the score form; B. the score form of new on
## warm-restarted fits against anova(..., test = "Rao"); C. a seed scan for fits whose top group mean risk lies just above
## one half (window 0.5002-0.5012) and the same comparisons there. Seeds are new (71e6, 72e6, 73e6 + i).
## Verification of the package, not a paper result.
## Run: Rscript pkg280_verify_corr_colrule.R <clone dir> <workers> > ../paper_EDGE/theory/pkg280_verify_corr_colrule.log 2>&1
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
PKG <- args[1]; NW <- as.integer(args[2])
OUT <- edge_path("declarations")
SCR <- file.path(tempdir(), "vcorr")
options(width = 220, warn = 1)
library(parallel)
t0 <- Sys.time()
lap <- function(label) cat(sprintf("  [%s done at %.1f min]\n", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

src_of <- function(sha) system2("git", c("-C", shQuote(PKG), "show", paste0(sha, ":R/def_gof.R")), stdout = TRUE)
SRC <- list(pre = src_of("0d63fba"), rel = src_of("906bd0c"))

init <- function(SRC) {
  suppressPackageStartupMessages(library(ebrahim.gof))
  ns <- asNamespace("ebrahim.gof")
  mk_env <- function(txt) {
    e <- new.env(parent = ns)
    for (ex in parse(text = txt, keep.source = FALSE))
      if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]])) eval(ex, e)
    e
  }
  assign("PRE", mk_env(SRC$pre), envir = globalenv())
  assign("REL", mk_env(SRC$rel), envir = globalenv())
  NULL
}

## ---- helpers (identical on master and workers) ----
gen_data <- function(n, c0, s) {
  dat <- data.frame(xa = runif(n, -3, 3), db = rbinom(n, 1, 0.5))
  dat$out <- rbinom(n, 1, plogis(c0 + s * (0.6 * dat$xa + 0.5 * dat$db)))
  dat
}
fit_of <- function(dat, ctl = glm.control()) suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat, control = ctl))
warm_of <- function(dat) {
  ctl <- glm.control(epsilon = 1e-13, maxit = 100)
  f0 <- fit_of(dat, ctl)
  suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat, start = unname(coef(f0)), control = ctl))
}
grp_means <- function(fit, G) {
  n <- length(fit$y); ph <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6)
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  list(grp = grp, pbar = as.numeric(tapply(ph, grp, mean)))
}
## group-level step covariates for anova; only exactly-zero columns are removed (independent of the package threshold)
step_cov <- function(fit, G, basis) {
  gm <- grp_means(fit, G); pbar <- gm$pbar; e <- qlogis(pbar)
  Zg <- switch(basis, poly2 = as.matrix(poly(pbar, 2)), poly3 = as.matrix(poly(pbar, 3)),
               stukel = cbind(e, e^2 * (e >= 0), -e^2 * (e < 0)), sym = cbind(e * abs(e)))
  Zg <- Zg[, colSums(Zg != 0) > 0, drop = FALSE]
  Zg[gm$grp, , drop = FALSE]
}
augment <- function(fit, M) {
  e <- new.env(parent = globalenv())
  e$out <- fit$y; e$Xm <- model.matrix(fit); e$M <- as.matrix(M)
  fo <- as.formula("out ~ Xm + M - 1", env = e)
  suppressWarnings(glm(fo, family = binomial(), control = fit$control))
}
run1 <- function(fun, f, G, basis, form) {
  cls <- character(0)
  r <- tryCatch(withCallingHandlers(
         switch(form, sat = fun(f, G = G, basis = basis),
                      imh = fun(f, G = G, basis = basis, method = "imhof"),
                      sco = fun(f, G = G, basis = basis, weights = "score")),
         warning = function(w) { cls <<- c(cls, class(w)[1]); invokeRestart("muffleWarning") }),
       error = function(e) NULL)
  if (is.null(r)) return(c(S = NA_real_, df = NA_real_, p = NA_real_, err = 1, nw = length(cls)))
  c(S = r$Test_Statistic, df = r$df, p = r$p_value, err = 0, nw = length(cls))
}
three <- function(f, G, basis, form) {
  a <- run1(PRE$def.gof, f, G, basis, form); b <- run1(def.gof, f, G, basis, form); c <- run1(REL$def.gof, f, G, basis, form)
  c(pre = a, new = b, rel = c)
}
BASES <- c("poly2", "poly3", "stukel", "sym")
FORMS <- c("sat", "imh", "sco")

one_fit <- function(i, design, anova_too = TRUE) {
  if (design == "c0m2") { set.seed(71000000 + i); dat <- gen_data(500, -2, 1); Gs <- c(10, 20) }
  else                  { set.seed(72000000 + i); dat <- gen_data(1000, 0, 1); Gs <- c(10, 20, 40) }
  f <- fit_of(dat)
  rows <- list()
  for (G in Gs) {
    top <- max(grp_means(f, G)$pbar)
    for (b in BASES) for (fm in FORMS) {
      v <- three(f, G, b, fm)
      rows[[length(rows) + 1]] <- data.frame(part = "A", design = design, i = i, G = G, basis = b, form = fm, top = top,
                                             t(v), rao = NA_real_, rao_df = NA_real_)
    }
  }
  if (anova_too) {
    fw <- warm_of(dat)
    for (G in Gs) {
      top <- max(grp_means(fw, G)$pbar)
      for (b in BASES) {
        v <- three(fw, G, b, "sco")
        a <- tryCatch(anova(fw, augment(fw, step_cov(fw, G, b)), test = "Rao"), error = function(e) NULL)
        rows[[length(rows) + 1]] <- data.frame(part = "B", design = design, i = i, G = G, basis = b, form = "sco", top = top,
                                               t(v), rao = if (is.null(a)) NA_real_ else a$Rao[2],
                                               rao_df = if (is.null(a)) NA_real_ else a$Df[2])
      }
    }
  }
  do.call(rbind, rows)
}

scan_seed <- function(sd) {
  set.seed(73000000 + sd); dat <- gen_data(500, -2, 1); f <- fit_of(dat)
  c(sd = sd, top10 = max(grp_means(f, 10)$pbar), top20 = max(grp_means(f, 20)$pbar))
}
window_fit <- function(sd) {
  set.seed(73000000 + sd); dat <- gen_data(500, -2, 1)
  fw <- warm_of(dat); f <- fit_of(dat)
  rows <- list()
  for (G in c(10, 20)) {
    topw <- max(grp_means(fw, G)$pbar); topf <- max(grp_means(f, G)$pbar)
    if (!(topw > 0.5 && topw < 0.5015) && !(topf > 0.5 && topf < 0.5015)) next
    for (b in c("stukel", "sym", "poly3")) {
      for (fm in FORMS) {
        v <- three(f, G, b, fm)
        rows[[length(rows) + 1]] <- data.frame(part = "C", design = "window", i = sd, G = G, basis = b, form = fm, top = topf,
                                               t(v), rao = NA_real_, rao_df = NA_real_)
      }
      v <- three(fw, G, b, "sco")
      a <- tryCatch(anova(fw, augment(fw, step_cov(fw, G, b)), test = "Rao"), error = function(e) NULL)
      rows[[length(rows) + 1]] <- data.frame(part = "Cw", design = "window", i = sd, G = G, basis = b, form = "sco", top = topw,
                                             t(v), rao = if (is.null(a)) NA_real_ else a$Rao[2],
                                             rao_df = if (is.null(a)) NA_real_ else a$Df[2])
    }
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}

## ---- run ----
cat("pkg280_verify_corr_colrule | clone", system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE),
    "| workers", NW, "|", format(t0), "\n")
init(SRC)
cat("installed ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "| pre = 0d63fba, rel = 906bd0c sourced\n")
cl <- makeCluster(NW)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
invisible(clusterCall(cl, init, SRC))
clusterExport(cl, c("gen_data", "fit_of", "warm_of", "grp_means", "step_cov", "augment", "run1", "three", "BASES", "FORMS",
                    "one_fit", "scan_seed", "window_fit"))
clusterSetRNGStream(cl, 1)

A1 <- do.call(rbind, parLapplyLB(cl, 1:2000, one_fit, design = "c0m2")); lap("c0 = -2 fits")
A2 <- do.call(rbind, parLapplyLB(cl, 1:2000, one_fit, design = "base")); lap("base-design fits")
SC <- as.data.frame(do.call(rbind, parLapply(cl, 1:40000, scan_seed))); lap("seed scan")
cand <- SC$sd[(SC$top10 > 0.5 & SC$top10 < 0.5015) | (SC$top20 > 0.5 & SC$top20 < 0.5015)]
cat(sprintf("  seed scan: %d of 40000 fits have a top group mean in (0.5, 0.5015) at G 10 or 20\n", length(cand)))
W  <- do.call(rbind, parLapplyLB(cl, cand, window_fit)); lap("window fits")
stopCluster(cl)
ALL <- rbind(A1, A2, W)
saveRDS(list(ALL = ALL, SC = SC), file.path(SCR, "colrule_all.rds"))

## ---- summaries ----
eqv <- function(a, b, tol = 1e-9) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol * pmax(1, abs(a)))
same3 <- function(d, x, y) eqv(d[[paste0(x, ".S")]], d[[paste0(y, ".S")]]) & eqv(d[[paste0(x, ".df")]], d[[paste0(y, ".df")]]) &
                            eqv(d[[paste0(x, ".p")]], d[[paste0(y, ".p")]])
ALL$pre_ran  <- ALL$pre.err == 0
ALL$new_fin  <- ALL$new.err == 0 & is.finite(ALL$new.p)
ALL$eq_new   <- same3(ALL, "pre", "new")
ALL$eq_rel   <- same3(ALL, "pre", "rel")
ALL$imh_noise <- !ALL$eq_new & ALL$form == "imh" & eqv(ALL$pre.S, ALL$new.S) & eqv(ALL$pre.df, ALL$new.df) &
                 !is.na(ALL$pre.p) & !is.na(ALL$new.p) & abs(ALL$pre.p - ALL$new.p) <= 1e-6

cat("\n==== A. every fit: new (installed) against pre (0d63fba); rel (906bd0c) as the positive control ====\n")
cat("equal = S, df and p within 1e-9 (relative, floor 1); imhof noise = S and df equal, |dp| <= 1e-6\n")
A <- ALL[ALL$part == "A", ]
key <- interaction(A$design, A$G, A$basis, A$form, drop = TRUE, lex.order = TRUE)
tab <- do.call(rbind, lapply(split(A, key), function(d) data.frame(
  design = d$design[1], G = d$G[1], basis = d$basis[1], form = d$form[1], fits = nrow(d),
  pre_ran = sum(d$pre_ran), new_eq_pre = sum(d$pre_ran & d$eq_new), imh_noise = sum(d$pre_ran & d$imh_noise),
  new_differs = sum(d$pre_ran & !d$eq_new & !d$imh_noise),
  max_dS = if (any(d$pre_ran & !is.na(d$new.S))) max(abs(d$pre.S - d$new.S), na.rm = TRUE) else NA,
  max_dp = if (any(d$pre_ran & !is.na(d$new.p))) max(abs(d$pre.p - d$new.p), na.rm = TRUE) else NA,
  pre_stop = sum(!d$pre_ran), stop_new_finite = sum(!d$pre_ran & d$new_fin), new_err = sum(d$new.err == 1),
  rel_differs = sum(d$pre_ran & !d$eq_rel), new_warn = sum(d$new.nw > 0))))
rownames(tab) <- NULL
print(tab, digits = 3)
cat(sprintf("\nTOTAL part A: %d comparisons; pre ran %d; new differs (beyond Imhof noise) %d; Imhof noise %d (max |dp| %.2e);\n",
            nrow(A), sum(A$pre_ran), sum(A$pre_ran & !A$eq_new & !A$imh_noise), sum(A$pre_ran & A$imh_noise),
            if (any(A$pre_ran & A$imh_noise)) max(abs(A$pre.p - A$new.p)[A$pre_ran & A$imh_noise]) else 0))
cat(sprintf("              pre stopped %d, new finite on %d of them; new errors %d; rel (906bd0c) moves %d results that ran under pre\n",
            sum(!A$pre_ran), sum(!A$pre_ran & A$new_fin), sum(A$new.err == 1), sum(A$pre_ran & !A$eq_rel)))
cat(sprintf("              max |dS| over non-Imhof equal cases %.2e\n",
            max(abs(A$pre.S - A$new.S)[A$pre_ran & A$eq_new & !is.na(A$new.S)])))
bad <- A[(A$pre_ran & !A$eq_new & !A$imh_noise) | !A$pre_ran | A$imh_noise, ]
write.csv(bad, file.path(OUT, "pkg280_verify_corr_colrule_cases.csv"), row.names = FALSE)
if (nrow(bad)) { cat("\ncases (pre stopped, differing, or Imhof noise):\n"); print(bad[, c("design", "i", "G", "basis", "form", "top", "pre.S", "pre.df",
  "pre.p", "pre.err", "new.S", "new.df", "new.p", "new.err", "rel.S", "rel.df", "rel.p")], digits = 10) }

cat("\n==== B. score form of new on warm-restarted fits against anova(fit, augmented, test = 'Rao') ====\n")
B <- ALL[ALL$part %in% c("B", "Cw"), ]
B$rel_err <- abs(B$new.S - B$rao) / pmax(1, abs(B$rao))
keyB <- interaction(B$part, B$design, B$G, B$basis, drop = TRUE, lex.order = TRUE)
tabB <- do.call(rbind, lapply(split(B, keyB), function(d) data.frame(
  part = d$part[1], design = d$design[1], G = d$G[1], basis = d$basis[1], fits = nrow(d),
  rao_na = sum(is.na(d$rao)), new_na = sum(is.na(d$new.S)), df_mismatch = sum(d$new.df != d$rao_df, na.rm = TRUE),
  max_rel_err = max(d$rel_err, na.rm = TRUE), over_1e8 = sum(d$rel_err > 1e-8, na.rm = TRUE),
  new_eq_pre = sum(eqv(d$pre.S, d$new.S) & eqv(d$pre.df, d$new.df)), rel_df_lt = sum(d$rel.df < d$new.df, na.rm = TRUE))))
rownames(tabB) <- NULL
print(tabB, digits = 3)
cat(sprintf("\nTOTAL part B + window: %d score forms; df mismatches %d; relative |S - Rao| > 1e-8: %d; max %.2e\n",
            nrow(B), sum(B$new.df != B$rao_df, na.rm = TRUE), sum(B$rel_err > 1e-8, na.rm = TRUE), max(B$rel_err, na.rm = TRUE)))

cat("\n==== C. window fits (top group mean in (0.5, 0.5015) at G 10 or 20), stukel basis ====\n")
Cw <- ALL[ALL$part %in% c("C", "Cw") & ALL$basis == "stukel", ]
Cw$inwin <- Cw$top >= 0.5002 & Cw$top <= 0.5012
print(Cw[order(Cw$i, Cw$G, Cw$part, Cw$form), c("part", "i", "G", "form", "top", "inwin", "pre.S", "pre.df", "pre.p", "pre.err", "new.S", "new.df",
  "new.p", "rel.df", "rel.p", "rao", "rao_df")], digits = 8, row.names = FALSE)
cw_w <- Cw[Cw$part == "Cw" & Cw$inwin, ]
cat(sprintf("\nwarm window fits with top mean in [0.5002, 0.5012]: %d (G-specific); score df = Rao df on %d; max rel |S - Rao| %.2e; rel (906bd0c) df smaller on %d\n",
            nrow(cw_w), sum(cw_w$new.df == cw_w$rao_df, na.rm = TRUE), max(abs(cw_w$new.S - cw_w$rao) / pmax(1, cw_w$rao), na.rm = TRUE),
            sum(cw_w$rel.df < cw_w$new.df, na.rm = TRUE)))
Cc <- ALL[ALL$part == "C", ]
cat(sprintf("window, unwarmed fits, all bases and forms: %d comparisons; pre ran %d; new differs beyond Imhof noise %d; pre stopped %d, new finite %d\n",
            nrow(Cc), sum(Cc$pre_ran), sum(Cc$pre_ran & !Cc$eq_new & !Cc$imh_noise), sum(!Cc$pre_ran), sum(!Cc$pre_ran & Cc$new_fin)))
cat("\ndone", format(Sys.time()), sprintf("(%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
