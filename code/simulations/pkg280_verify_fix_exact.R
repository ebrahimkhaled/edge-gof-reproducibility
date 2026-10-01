## pkg280_verify_fix_exact.R -- on the fixed branch (load_all): the joint Stukel statistic against the exact u'I^-1 u
## (pivoted weighted QR) and anova(..., test = "Rao") on warm-restarted fits; lr against the deviance drop; the
## marginal form against 2.7.0 (parsed from main) and statmod; and the def.gof score form against anova Rao on the
## grouped step covariates, on fixtures and on 100 random warm-restarted fits per design.
## Run: Rscript pkg280_verify_fix_exact.R [package dir] > ../paper_EDGE/theory/pkg280_verify_fix_exact.log 2>&1
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
    R.version.string, "| statmod", as.character(packageVersion("statmod")), "|", format(Sys.time()), "\n")

v270 <- new.env(parent = ns)
for (f in c("R/def_gof.R", "R/run_all_gof.R")) {
  txt <- system2("git", c("-C", shQuote(PKG), "show", paste0("83c7e11:", f)), stdout = TRUE)
  for (ex in parse(text = txt, keep.source = FALSE))
    if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) &&
        as.character(ex[[2]]) %in% c("gof_stukel", ".gof_score_z", ".gof_context"))
      eval(ex, v270)
}
SRC  <- file.path(SIM, "pkg280_verify_num_exact.R")
TAKE <- c("mk", "find_fit", "augment", "score_qr", "stukel_Z", "stk", "Gof", "step_cov", "count_warn", "quiet")
for (e in parse(SRC, keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% TAKE)
    eval(e, globalenv())
stopifnot(all(TAKE %in% ls(globalenv())))
warm <- function(fit) {
  env <- environment(formula(fit))
  cl <- getCall(fit); cl$control <- glm.control(epsilon = 1e-13, maxit = 100)
  ft <- suppressWarnings(eval(cl, env))
  cl$start <- unname(coef(ft))
  suppressWarnings(eval(cl, env))
}
ROWS <- list()
rec <- function(part, design, quantity, value) ROWS[[length(ROWS) + 1]] <<-
  data.frame(part = part, design = design, quantity = quantity, value = value, stringsAsFactors = FALSE)

## ---- 2a. fixtures ----
cat("\n==== 2a. Stukel on fixtures ====\n")
FX <- list(
  base      = mk(101, 1000),
  c0_m2     = mk(102, 500, c0 = -2),
  s_2       = mk(103, 1000, s = 2),
  cloglog   = mk(104, 1500, truth = "cloglog"),
  offset    = local({ set.seed(108); dat <- data.frame(xa = runif(800, -3, 3), db = rbinom(800, 1, 0.5))
                      dat$out <- rbinom(800, 1, plogis(-0.5 + 0.7 * dat$xa + 0.4 * dat$db))
                      glm(out ~ xa + offset(0.4 * db), family = binomial(), data = dat) }),
  all_below = find_fit(5:500, function(sd) mk(sd, 400, eta_fun = function(d) -4 + 0.5 * d$xa),
                       function(f) max(fitted(f)) < 0.5)
)
for (nm in names(FX)) {
  f  <- FX[[nm]]; fw <- warm(f); Zw <- stukel_Z(fw)
  j  <- stk(fw, "joint"); a <- anova(fw, augment(fw, Zw), test = "Rao"); q <- score_qr(fw, Zw)
  j0 <- stk(f, "joint"); q0 <- score_qr(f, stukel_Z(f))
  l0 <- stk(f, "lr"); A0 <- augment(f, stukel_Z(f))
  m  <- stk(f, "marginal"); o <- v270$gof_stukel(v270$.gof_context(f, NULL, NULL, G = 10))
  Zf <- stukel_Z(f)
  sm <- if (all(colSums(Zf != 0) > 0)) sum(statmod::glm.scoretest(f, Zf)^2) else NA_real_
  same270 <- identical(unname(m$Statistic), unname(o$Statistic)) && identical(unname(m$df), unname(o$df)) &&
             identical(unname(m$p_value), unname(o$p_value))
  cat(sprintf("  %-9s warm: joint %.12f df %g | anova Rao %.12f Df %g | |diff| %.1e | |joint - QR| %.1e\n",
              nm, j$Statistic, j$df, a$Rao[2], a$Df[2], abs(j$Statistic - a$Rao[2]), abs(j$Statistic - q$stat)))
  cat(sprintf("            default fit: joint %.10f (Note '%s') | |joint - QR| %.1e | lr %.10f df %g | deviance drop %.10f rank diff %d | |lr - drop| %.1e\n",
              j0$Statistic, j0$Note, abs(j0$Statistic - q0$stat), l0$Statistic, l0$df, f$deviance - A0$deviance,
              A0$rank - f$rank, abs(l0$Statistic - (f$deviance - A0$deviance))))
  cat(sprintf("            marginal %s | identical() to 2.7.0: %s | statmod sum %s\n",
              format(m$Statistic, digits = 10), same270, format(sm, digits = 10)))
  rec("fixture", nm, "joint vs anova (warm), abs", abs(j$Statistic - a$Rao[2]))
  rec("fixture", nm, "joint df == anova Df", as.numeric(j$df == a$Df[2]))
  rec("fixture", nm, "joint vs QR (warm), abs", abs(j$Statistic - q$stat))
  rec("fixture", nm, "joint vs QR (default fit), abs", abs(j0$Statistic - q0$stat))
  rec("fixture", nm, "lr vs deviance drop, abs", abs(l0$Statistic - (f$deviance - A0$deviance)))
  rec("fixture", nm, "marginal identical to 2.7.0", as.numeric(same270))
  rec("fixture", nm, "marginal vs statmod sum, abs", abs(m$Statistic - sm))
}

## ---- 2b. random fits ----
cat("\n==== 2b. 100 random warm-restarted fits per design ====\n")
gen <- function(n, c0, s) {
  dat <- data.frame(xa = runif(n, -3, 3), db = rbinom(n, 1, 0.5))
  dat$out <- rbinom(n, 1, plogis(c0 + s * (0.6 * dat$xa + 0.5 * dat$db)))
  suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat))
}
CELLS <- data.frame(label = c("base n1000", "c0=-2 n500", "s=2 n1000"), n = c(1000, 500, 1000), c0 = c(0, -2, 0),
                    s = c(1, 1, 2), stringsAsFactors = FALSE)
for (k in seq_len(nrow(CELLS))) {
  ce <- CELLS[k, ]; set.seed(990000 + k)
  D <- matrix(NA_real_, 100, 14, dimnames = list(NULL, c("joint", "joint_df", "joint_qr", "marg270", "onedf",
                                                         "sym", "sym_df", "poly3", "poly3_df", "poly2", "poly2_df",
                                                         "stukel", "stukel_df", "stukel_k2")))
  for (i in 1:100) {
    f <- gen(ce$n, ce$c0, ce$s); fw <- warm(f); Zw <- stukel_Z(fw)
    j <- stk(fw, "joint"); a <- anova(fw, augment(fw, Zw), test = "Rao"); q <- score_qr(fw, Zw)
    m <- stk(f, "marginal"); o <- v270$gof_stukel(v270$.gof_context(f, NULL, NULL, G = 10))
    D[i, c("joint", "joint_df", "joint_qr", "marg270", "onedf")] <-
      c(abs(j$Statistic - a$Rao[2]) / max(a$Rao[2], 1e-8), j$df != a$Df[2], abs(j$Statistic - q$stat),
        identical(unname(c(m$Statistic, m$df, m$p_value)), unname(c(o$Statistic, o$df, o$p_value))), j$df == 1)
    for (b in c("sym", "poly3", "poly2", "stukel")) {
      r  <- quiet(def.gof(fw, G = "auto", basis = b, weights = "score"))
      ab <- anova(fw, augment(fw, step_cov(fw, "auto", b)), test = "Rao")
      D[i, paste0(b, c("", "_df"))] <- c(abs(r$Test_Statistic - ab$Rao[2]) / max(ab$Rao[2], 1e-8), r$df != ab$Df[2])
      if (b == "stukel") D[i, "stukel_k2"] <- r$df == 2
    }
  }
  cat(sprintf("  %-11s joint vs anova: max rel %.2e, df mismatches %d, 1-df fits %d | joint vs QR max abs %.2e | marginal identical to 2.7.0: %d/100\n",
              ce$label, max(D[, "joint"]), sum(D[, "joint_df"]), sum(D[, "onedf"]), max(D[, "joint_qr"]), sum(D[, "marg270"])))
  cat(sprintf("              score vs anova, max rel (df mismatches): sym %.2e (%d) | poly3 %.2e (%d) | poly2 %.2e (%d) | stukel %.2e (%d), stukel df 2 in %d fits\n",
              max(D[, "sym"]), sum(D[, "sym_df"]), max(D[, "poly3"]), sum(D[, "poly3_df"]), max(D[, "poly2"]),
              sum(D[, "poly2_df"]), max(D[, "stukel"]), sum(D[, "stukel_df"]), sum(D[, "stukel_k2"])))
  for (v in c("joint", "joint_qr", "sym", "poly3", "poly2", "stukel")) rec("random", ce$label, paste(v, "max"), max(D[, v]))
  for (v in c("joint_df", "sym_df", "poly3_df", "poly2_df", "stukel_df")) rec("random", ce$label, paste(v, "mismatches"), sum(D[, v]))
  rec("random", ce$label, "marginal identical to 2.7.0 (of 100)", sum(D[, "marg270"]))
}
W <- do.call(rbind, ROWS); W$head <- sha
write.csv(W, file.path(OUT, "pkg280_verify_fix_exact_checks.csv"), row.names = FALSE)
cat("\ndone", format(Sys.time()), "\n")
