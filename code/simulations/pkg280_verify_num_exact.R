## pkg280_verify_num_exact.R -- exact numerical checks of the INSTALLED ebrahim.gof 2.8.0 against base R, statmod,
## LogisticDx and the 2.7.0 code on main. Verification of the package, not a paper result.
##
##  (a) Stukel: joint = anova(fit, augmented, test = "Rao"); lr = deviance drop with the right df;
##      marginal = 2.7.0 gof_stukel = LogisticDx SstBoth = sum of squared statmod::glm.scoretest z's.
##  (b) def.gof weights = "score" = Rao score test for adding the grouped step covariates (sym, poly3, poly2,
##      stukel); unit form identical to 2.7.0 def.gof (and def.ensemble.gof / edge.gof unchanged).
##  (d) G = "auto" = max(10, round(n/25)) end to end; the few-events warning fires exactly when
##      min(events, non-events) < G.
##
## Run: Rscript pkg280_verify_num_exact.R > ../paper_EDGE/theory/pkg280_verify_num_exact.log 2>&1
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

suppressPackageStartupMessages(library(ebrahim.gof))
PKGV <- as.character(packageVersion("ebrahim.gof"))
stopifnot(PKGV == "2.8.0")
PKG <- Sys.getenv("EBRAHIM_GOF_SRC")
OUT <- edge_path("declarations")
options(width = 220, warn = 1, digits = 10)
main_sha <- system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "main"), stdout = TRUE)
cat("ebrahim.gof", PKGV, "(installed) | 2.7.0 code from main", main_sha, "|", R.version.string,
    "| statmod", as.character(packageVersion("statmod")), "| LogisticDx", as.character(packageVersion("LogisticDx")),
    "|", format(Sys.time()), "\n")

## ---- the 2.7.0 functions, parsed from main into one environment (parent: the installed namespace) ----------------
old <- new.env(parent = asNamespace("ebrahim.gof"))
KEEP <- c("def.gof", ".def_basis", ".def_pvalue", "def.ensemble.gof", "edges.gof", ".combine_pvalues",
          "gof_stukel", ".gof_score_z", ".gof_context")
for (f in c("R/def_gof.R", "R/def_ensemble_gof.R", "R/run_all_gof.R")) {
  txt <- system2("git", c("-C", shQuote(PKG), "show", paste0("main:", f)), stdout = TRUE)
  for (e in parse(text = txt, keep.source = FALSE))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% KEEP)
      eval(e, old)
}
stopifnot(all(KEEP %in% ls(old, all.names = TRUE)))
stopifnot(!"weights" %in% names(formals(old$def.gof)))          # really the old signature

## ---- bookkeeping --------------------------------------------------------------------------------------------------
ROWS <- list()
add <- function(part, design, quantity, pkg, ref, tol, note = "") {
  pkg <- suppressWarnings(as.numeric(pkg)); ref <- suppressWarnings(as.numeric(ref))
  ad <- abs(pkg - ref); rd <- ad / max(abs(ref), 1e-12)
  pass <- if (is.na(pkg) && is.na(ref)) TRUE else if (is.na(ad)) FALSE else ad <= tol
  ROWS[[length(ROWS) + 1]] <<- data.frame(part = part, design = design, quantity = quantity, pkg = pkg, ref = ref,
                                          abs_diff = ad, rel_diff = rd, tol = tol, pass = pass, note = note,
                                          stringsAsFactors = FALSE)
}
addlog <- function(part, design, quantity, ok, note = "")
  ROWS[[length(ROWS) + 1]] <<- data.frame(part = part, design = design, quantity = quantity, pkg = NA_real_,
                                          ref = NA_real_, abs_diff = NA_real_, rel_diff = NA_real_, tol = NA_real_,
                                          pass = isTRUE(ok), note = note, stringsAsFactors = FALSE)

## run an expression, muffling and counting def_few_events warnings; other warnings muffled and counted apart
count_warn <- function(expr) {
  few <- 0L; other <- 0L; msgs <- character(0)
  val <- tryCatch(withCallingHandlers(expr,
           def_few_events = function(w) { few <<- few + 1L; msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
           warning = function(w) { other <<- other + 1L; invokeRestart("muffleWarning") }),
         error = function(e) structure(list(msg = conditionMessage(e)), class = "caught_error"))
  list(value = val, few = few, other = other, msgs = msgs)
}
quiet <- function(expr) count_warn(expr)$value

## ---- fixtures -----------------------------------------------------------------------------------------------------
## the run-L / MAP design: xa ~ U(-3,3), db ~ Bern(0.5), eta = c0 + s(0.6 xa + 0.5 db)
mk <- function(seed, n, c0 = 0, s = 1, truth = "logit", eta_fun = NULL, control = glm.control()) {
  set.seed(seed)
  dat <- data.frame(xa = runif(n, -3, 3), db = rbinom(n, 1, 0.5))
  eta <- if (is.null(eta_fun)) c0 + s * (0.6 * dat$xa + 0.5 * dat$db) else eta_fun(dat)
  dat$out <- rbinom(n, 1, if (truth == "logit") plogis(eta) else 1 - exp(-exp(eta)))
  suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat, control = control))
}
find_fit <- function(seeds, make, cond) {
  for (sd in seeds) { f <- make(sd); if (cond(f)) { attr(f, "seed") <- sd; return(f) } }
  stop("find_fit: no seed satisfied the condition")
}
tighten <- function(fit) {
  cl <- update(fit, control = glm.control(epsilon = 1e-13, maxit = 100), evaluate = FALSE)
  suppressWarnings(eval(cl, environment(formula(fit))))
}
## augmented fit with extra columns M, response named as in fit (so anova() accepts the pair)
augment <- function(fit, M) {
  resp <- deparse(formula(fit)[[2L]])
  e <- new.env(parent = globalenv())
  assign(resp, fit$y, envir = e); e$Xm <- model.matrix(fit); e$M <- as.matrix(M); e$off <- fit$offset
  e$fam <- fit$family; e$ctl <- fit$control
  e$fo <- as.formula(paste(resp, "~ Xm + M - 1"), env = e)
  suppressWarnings(eval(quote(glm(fo, family = fam, offset = off, control = ctl)), e))
}
## u'I^-1 u at the fitted values, information by weighted QR (an independent route to the score statistic)
score_qr <- function(fit, M) {
  M <- as.matrix(M); M <- M[, colSums(M != 0) > 0, drop = FALSE]
  p <- as.numeric(fitted(fit)); W <- p * (1 - p); X <- model.matrix(fit); y <- fit$y
  u <- colSums(M * (y - p))
  sw <- sqrt(W); Mr <- qr.resid(qr(sw * X), sw * M)
  q <- qr(Mr); k <- q$rank
  if (k < ncol(M)) { piv <- q$pivot[seq_len(k)]; u <- u[piv]; Mr <- Mr[, piv, drop = FALSE] }
  list(stat = drop(crossprod(u, solve(crossprod(Mr), u))), df = k)
}
## anova's Rao, rebuilt: weighted regression of working residuals on the big design; stale = glm's own weights
rao_regress <- function(fit, M, fresh = FALSE) {
  Xb <- cbind(model.matrix(fit), as.matrix(M))
  if (fresh) { p <- as.numeric(fitted(fit)); w <- p * (1 - p); r <- (fit$y - p) / w }
  else { w <- fit$weights; r <- fit$residuals }
  z <- lm.wfit(Xb, r, w)
  sum(w * (r - sum(w * r) / sum(w))^2) - sum(w * z$residuals^2)
}
stukel_Z <- function(fit) {
  eta <- as.numeric(predict(fit, type = "link"))
  ph  <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6)
  cbind(za = 0.5 * eta^2 * (ph >= 0.5), zb = -0.5 * eta^2 * (ph < 0.5))
}
stk <- function(fit, form = "joint")
  as.data.frame(suppressWarnings(run.all.gof(fit, tests = "Stukel", install = "no",
                                             control = list(Stukel = list(form = form)))))

data("gof_demo", package = "ebrahim.gof")
demo_wrong <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial())
demo_right <- glm(outcome ~ poly(age, 2) + bmi + sex + treatment, data = gof_demo, family = binomial())

A <- list(
  base      = mk(101, 1000),
  c0_m2     = mk(102, 500, c0 = -2),
  s_2       = mk(103, 1000, s = 2),
  cloglog   = mk(104, 1500, truth = "cloglog"),
  all_below = find_fit(5:500, function(sd) mk(sd, 400, eta_fun = function(d) -4 + 0.5 * d$xa),
                       function(f) max(fitted(f)) < 0.5),
  all_above = find_fit(6:500, function(sd) mk(sd, 400, eta_fun = function(d) 4 + 0.5 * d$xa),
                       function(f) min(fitted(f)) >= 0.5),
  two_above = find_fit(57:2000, function(sd) mk(sd, 400, eta_fun = function(d) -2.2 + 0.6 * d$xa),
                       function(f) sum(fitted(f) >= 0.5) %in% 1:3),
  demo_wrong = demo_wrong,
  offset    = local({ set.seed(108); dat <- data.frame(xa = runif(800, -3, 3), db = rbinom(800, 1, 0.5))
                      dat$out <- rbinom(800, 1, plogis(-0.5 + 0.7 * dat$xa + 0.4 * dat$db))
                      glm(out ~ xa + offset(0.4 * db), family = binomial(), data = dat) })
)
cat("\nfixtures:\n")
for (nm in names(A)) cat(sprintf("  %-10s n %5d  events %4d  fitted in [%.4f, %.4f]  #fitted>=0.5 %4d  seed %s\n", nm,
                                 length(A[[nm]]$y), as.integer(sum(A[[nm]]$y)), min(fitted(A[[nm]])),
                                 max(fitted(A[[nm]])), sum(fitted(A[[nm]]) >= 0.5),
                                 if (is.null(attr(A[[nm]], "seed"))) "-" else attr(A[[nm]], "seed")))

## =================================================================================================================
## (a) Stukel
## =================================================================================================================
cat("\n==== (a) Stukel ====\n")
for (nm in names(A)) {
  fit <- A[[nm]]; ft <- tighten(fit); Z <- stukel_Z(fit); Zt <- stukel_Z(ft)
  j <- stk(fit, "joint"); jt <- stk(ft, "joint"); l <- stk(fit, "lr"); m <- stk(fit, "marginal")
  f1 <- augment(fit, Z); f1t <- augment(ft, Zt)
  a  <- anova(fit, f1, test = "Rao"); at <- anova(ft, f1t, test = "Rao")
  sq <- score_qr(fit, Z)
  kz <- sum(colSums(Z != 0) > 0)
  cat(sprintf("\n[%s] joint %.8f df %g p %.6g Note '%s' | anova Rao %.8f df %g | tight: joint %.10f anova %.10f\n",
              nm, j$Statistic, j$df, j$p_value, j$Note, a$Rao[2], a$Df[2], jt$Statistic, at$Rao[2]))
  add("a", nm, "joint stat vs anova Rao (epsilon 1e-13 refit)", jt$Statistic, at$Rao[2], 1e-5)
  add("a", nm, "joint df vs anova Df (epsilon 1e-13 refit)", jt$df, at$Df[2], 0)
  add("a", nm, "joint stat vs anova Rao (default glm control)", j$Statistic, a$Rao[2], 1e-5,
      "informational: anova uses glm's one-iteration-old weights")
  add("a", nm, "joint df vs anova Df (default control)", j$df, a$Df[2], 0)
  add("a", nm, "joint stat vs u'I^-1u by weighted QR", j$Statistic, sq$stat, 1e-8)
  add("a", nm, "joint stat vs working-residual regression, fresh weights", j$Statistic, rao_regress(fit, Z, TRUE), 1e-6)
  add("a", nm, "anova Rao vs working-residual regression, glm weights", a$Rao[2], rao_regress(fit, Z, FALSE), 1e-8,
      "reproduces anova, so the default-control gap is anova's stale weights")
  add("a", nm, "joint p = pchisq(stat, df)", j$p_value, pchisq(j$Statistic, j$df, lower.tail = FALSE), 1e-14)
  add("a", nm, "joint df = number of non-zero Stukel columns", j$df, kz, 0)
  addlog("a", nm, "joint Note names the side iff 1 df",
         if (kz == 2) identical(j$Note, "") else grepl("^1 df: no fitted risk", j$Note), j$Note)
  ## lr
  add("a", nm, "lr stat vs deviance(fit) - deviance(augmented)", l$Statistic, fit$deviance - f1$deviance, 1e-8)
  add("a", nm, "lr stat vs anova Deviance column", l$Statistic, anova(fit, f1)$Deviance[2], 1e-8)
  add("a", nm, "lr df vs rank(augmented) - rank(fit)", l$df, f1$rank - fit$rank, 0)
  add("a", nm, "lr p = pchisq(stat, df)", l$p_value, pchisq(l$Statistic, l$df, lower.tail = FALSE), 1e-14)
  cat(sprintf("       lr %.8f df %g | dev drop %.8f rank diff %d || marginal %.8f p %.6g\n", l$Statistic, l$df,
              fit$deviance - f1$deviance, f1$rank - fit$rank, m$Statistic, m$p_value))
  ## marginal against 2.7.0, statmod and LogisticDx
  o <- old$gof_stukel(old$.gof_context(fit, NULL, NULL, G = 10))
  addlog("a", nm, "marginal Statistic/df/p_value identical() to 2.7.0 gof_stukel",
         identical(m$Statistic, o$Statistic) && identical(m$df, o$df) && identical(m$p_value, o$p_value),
         sprintf("2.8.0 %.10g / 2.7.0 %.10g", m$Statistic, o$Statistic))
  zst <- suppressWarnings(statmod::glm.scoretest(fit, Z))
  add("a", nm, "marginal stat vs sum(statmod::glm.scoretest z^2)", m$Statistic, sum(zst^2), 1e-10)
  sst <- tryCatch(suppressMessages(suppressWarnings({
    g <- LogisticDx::gof(fit, plotROC = FALSE); as.data.frame(g$gof) })), error = function(e) NULL)
  if (!is.null(sst)) {
    add("a", nm, "marginal stat vs LogisticDx::gof SstBoth val", m$Statistic, sst$val[sst$test == "SstBoth"], 1e-6)
    add("a", nm, "marginal p vs LogisticDx::gof SstBoth pVal", m$p_value, sst$pVal[sst$test == "SstBoth"], 1e-6)
  } else addlog("a", nm, "LogisticDx::gof ran", FALSE, "LogisticDx::gof failed on this fit")
}

## guards and odd fits (behaviour, recorded; not tolerance checks)
cat("\n---- (a) guards ----\n")
set.seed(3); gd <- data.frame(xa = runif(300, -3, 3)); gd$out <- rbinom(300, 1, pnorm(0.5 * gd$xa))
probit <- glm(out ~ xa, family = binomial("probit"), data = gd)
wtd    <- glm(out ~ xa, family = binomial(), data = gd, weights = rep(2, 300))
gd$tri <- rbinom(300, 3, plogis(0.5 * gd$xa))
agg    <- glm(cbind(tri, 3 - tri) ~ xa, family = binomial(), data = gd)
gd$xa2 <- 2 * gd$xa
rankdef <- glm(out ~ xa + xa2, family = binomial(), data = gd)
for (nm in c("probit", "wtd", "agg", "rankdef")) {
  fit <- get(nm)
  for (fm in c("joint", "lr", "marginal")) {
    r <- tryCatch(stk(fit, fm), error = function(e) data.frame(Statistic = NA, df = NA, p_value = NA,
                                                                Note = paste("ERROR:", conditionMessage(e))))
    cat(sprintf("  %-8s %-8s stat %-12s df %-3s p %-12s Note: %s\n", nm, fm, format(r$Statistic, digits = 8),
                format(r$df), format(r$p_value, digits = 6), r$Note))
  }
}
addlog("a", "probit", "joint and lr return NA with a logit Note",
       is.na(stk(probit, "joint")$p_value) && grepl("logit", stk(probit, "lr")$Note))
addlog("a", "weighted", "joint returns NA with a Note", is.na(stk(wtd, "joint")$p_value))
addlog("a", "aggregated", "joint returns NA with a Note", is.na(stk(agg, "joint")$p_value))
rd_rao <- anova(rankdef, augment(rankdef, stukel_Z(rankdef)), test = "Rao")
cat(sprintf("  rankdef: anova Rao %.6f df %g; score_qr %.6f df %d\n", rd_rao$Rao[2], rd_rao$Df[2],
            score_qr(rankdef, stukel_Z(rankdef))$stat, score_qr(rankdef, stukel_Z(rankdef))$df))

## =================================================================================================================
## (b) score form = Rao for the grouped step covariates; unit form identical to 2.7.0
## =================================================================================================================
cat("\n==== (b) def.gof weights = 'score' ====\n")
Gof <- function(n, G) if (identical(G, "auto")) max(10, round(n / 25)) else G
step_cov <- function(fit, G, basis) {
  y <- fit$y; n <- length(y); G <- Gof(n, G)
  ph  <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6)
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  lev <- sort(unique(grp))
  pbar <- vapply(lev, function(g) mean(ph[grp == g]), numeric(1))
  e <- qlogis(pbar)
  Zg <- switch(basis, poly2 = as.matrix(poly(pbar, 2)), poly3 = as.matrix(poly(pbar, 3)),
               stukel = cbind(e, e^2 * (e >= 0), -e^2 * (e < 0)), sym = cbind(e * abs(e)))
  Zg <- Zg[, colSums(abs(Zg)) > 1e-8, drop = FALSE]
  structure(Zg[match(grp, lev), , drop = FALSE], pbar = pbar)
}
one_group_half <- find_fit(300:3000, function(sd) mk(sd, 500, c0 = -2),
                           function(f) sum(attr(step_cov(f, 20, "sym"), "pbar") >= 0.5) == 1)
BD <- list(
  list(name = "base n1000 G10",     fit = mk(201, 1000),          G = 10),
  list(name = "c0=-2 n500 G20",     fit = mk(202, 500, c0 = -2),  G = 20),
  list(name = "s=2 n1000 G40",      fit = mk(203, 1000, s = 2),   G = 40),
  list(name = "base n1000 Gauto",   fit = mk(201, 1000),          G = "auto"),
  list(name = "cloglog n1500 Gauto", fit = A$cloglog,             G = "auto"),
  list(name = "demo_wrong G10",     fit = demo_wrong,             G = 10),
  list(name = sprintf("c0=-2 n500 G20, one group >= 0.5 (seed %d)", attr(one_group_half, "seed")),
       fit = one_group_half, G = 20)
)
for (d in BD) {
  fit <- d$fit; ft <- tighten(fit)
  for (b in c("sym", "poly3", "poly2", "stukel")) {
    r  <- quiet(def.gof(fit, G = d$G, basis = b, weights = "score"))
    rt <- quiet(def.gof(ft,  G = d$G, basis = b, weights = "score"))
    M  <- step_cov(fit, d$G, b); Mt <- step_cov(ft, d$G, b)
    sq <- score_qr(fit, M)
    at <- anova(ft, augment(ft, Mt), test = "Rao"); a <- anova(fit, augment(fit, M), test = "Rao")
    cat(sprintf("  %-44s %-6s S %.8f df %g p %.6g | QR %.8f df %d | anova(tight) %.8f df %g | anova(default) %.8f\n",
                d$name, b, r$Test_Statistic, r$df, r$p_value, sq$stat, sq$df, at$Rao[2], at$Df[2], a$Rao[2]))
    add("b", d$name, paste(b, "score S vs anova Rao, step covariates (epsilon 1e-13 refit)"), rt$Test_Statistic, at$Rao[2], 1e-5)
    add("b", d$name, paste(b, "score df vs anova Df (epsilon 1e-13 refit)"), rt$df, at$Df[2], 0)
    add("b", d$name, paste(b, "score S vs anova Rao (default control)"), r$Test_Statistic, a$Rao[2], 1e-5,
        "informational: anova's stale weights")
    add("b", d$name, paste(b, "score S vs u'I^-1u by weighted QR"), r$Test_Statistic, sq$stat, 1e-8)
    add("b", d$name, paste(b, "score df vs QR rank"), r$df, sq$df, 0)
    add("b", d$name, paste(b, "score p = pchisq(S, df)"), r$p_value, pchisq(r$Test_Statistic, r$df, lower.tail = FALSE), 1e-14)
    if (inherits(fit, "glm") && fit$family$link == "logit") {
      ry <- quiet(def.gof(as.numeric(fit$y), fitted(fit), X = model.matrix(fit), G = d$G, basis = b, weights = "score"))
      add("b", d$name, paste(b, "score: (y, ph, X) path vs glm path"), ry$Test_Statistic, r$Test_Statistic, 1e-9)
    }
  }
}
## the score form on a probit fit is not the probit Rao test (documentation wording check)
pr_def <- quiet(def.gof(probit, G = 10, basis = "sym", weights = "score"))
pr_rao <- anova(tighten(probit), augment(tighten(probit), step_cov(tighten(probit), 10, "sym")), test = "Rao")
cat(sprintf("  probit fit, sym score: def.gof S %.6f p %.4g | probit anova Rao %.6f\n",
            pr_def$Test_Statistic, pr_def$p_value, pr_rao$Rao[2]))

## unit form identical to 2.7.0
cat("\n---- (b) unit form against 2.7.0 def.gof ----\n")
UF <- c(A[setdiff(names(A), "offset")], list(demo_right = demo_right, offset = A$offset,
        mf1_500 = local({ set.seed(1); x <- runif(500, -3, 3); dd <- data.frame(x = x, yv = rbinom(500, 1, plogis(0.6 * x)))
                          glm(yv ~ x, family = binomial(), data = dd) })))
same_obj <- function(x, y) {
  if (inherits(x, "caught_error") || inherits(y, "caught_error")) return(identical(unclass(x), unclass(y)))
  identical(x, y)
}
n_cmp <- 0L; n_same <- 0L; bad <- character(0); n_few <- 0L
for (nm in names(UF)) {
  fit <- UF[[nm]]; n <- length(fit$y)
  for (b in c("poly2", "poly3", "stukel")) for (mth in c("satterthwaite", "imhof")) for (G in list(10, 5, 20, "auto")) {
    Gn <- Gof(n, G)
    new <- count_warn(def.gof(fit, G = G, basis = b, method = mth))
    ref <- count_warn(old$def.gof(fit, G = Gn, basis = b, method = mth))
    n_few <- n_few + new$few
    n_cmp <- n_cmp + 1L; ok <- same_obj(new$value, ref$value); n_same <- n_same + ok
    if (!ok) bad <- c(bad, paste(nm, b, mth, format(G)))
    newe <- count_warn(edge.gof(fit, G = G, basis = b, method = mth))
    if (!inherits(newe$value, "caught_error") && !inherits(ref$value, "caught_error")) {
      rr <- ref$value; rr$Test <- "EDGE"
      n_cmp <- n_cmp + 1L; ok <- identical(newe$value, rr); n_same <- n_same + ok
      if (!ok) bad <- c(bad, paste("edge", nm, b, mth, format(G)))
    }
  }
  X <- model.matrix(fit); yv <- as.numeric(fit$y); pv <- as.numeric(fitted(fit))
  for (b in c("poly2", "poly3", "stukel")) {
    n_cmp <- n_cmp + 2L
    ok1 <- same_obj(quiet(def.gof(yv, pv, X = X, basis = b)), quiet(old$def.gof(yv, pv, X = X, basis = b)))
    ok2 <- same_obj(quiet(def.gof(yv, pv, basis = b)), quiet(old$def.gof(yv, pv, basis = b)))
    n_same <- n_same + ok1 + ok2
    if (!ok1) bad <- c(bad, paste(nm, b, "yphX")); if (!ok2) bad <- c(bad, paste(nm, b, "yph naive"))
  }
  for (cmb in c("cct", "minp", "fisher")) for (ef in c(FALSE, TRUE)) {
    n_cmp <- n_cmp + 1L
    ok <- same_obj(quiet(def.ensemble.gof(fit, combine = cmb, add_ef = ef)),
                   quiet(old$def.ensemble.gof(fit, combine = cmb, add_ef = ef)))
    n_same <- n_same + ok; if (!ok) bad <- c(bad, paste("ensemble", nm, cmb, ef))
  }
}
cat(sprintf("  %d comparisons (def.gof glm x 3 bases x 2 methods x G in 10/5/20/auto, edge.gof, (y,ph,X), (y,ph), def.ensemble.gof x 3 combiners x add_ef) on %d fits: %d identical()\n",
            n_cmp, length(UF), n_same))
if (length(bad)) cat("  NOT identical:", paste(bad, collapse = " | "), "\n")
addlog("b", "unit form", sprintf("identical() to 2.7.0 in %d of %d comparisons", n_same, n_cmp), n_same == n_cmp,
       paste(bad, collapse = " | "))

## =================================================================================================================
## (d) G = "auto" and the few-events warning
## =================================================================================================================
cat("\n==== (d) G = 'auto' ====\n")
for (n in c(100, 249, 250, 262, 263, 287, 288, 312, 313, 500, 1000, 1012, 1013, 1500, 2000, 5000)) {
  fit <- mk(700 + n, n); Ge <- max(10, round(n / 25))
  ok <- TRUE
  for (b in c("poly3", "sym", "stukel")) for (wt in c("unit", "score"))
    ok <- ok && same_obj(quiet(def.gof(fit, G = "auto", basis = b, weights = wt)), quiet(def.gof(fit, G = Ge, basis = b, weights = wt)))
  ok_edge <- same_obj(quiet(edge.gof(fit, G = "auto", basis = "sym")), quiet(edge.gof(fit, G = Ge, basis = "sym")))
  ok_ens  <- same_obj(quiet(def.ensemble.gof(fit, G = "auto", add_ef = TRUE)), quiet(def.ensemble.gof(fit, G = Ge, add_ef = TRUE)))
  ok_yph  <- same_obj(quiet(def.gof(as.numeric(fit$y), fitted(fit), X = model.matrix(fit), G = "auto", weights = "score")),
                      quiet(def.gof(as.numeric(fit$y), fitted(fit), X = model.matrix(fit), G = Ge, weights = "score")))
  bat <- quiet(as.data.frame(run.all.gof(fit, tests = c("DEF.sym", "DEF.poly3"), install = "no",
                                         control = list(DEF.sym = list(G = "auto", weights = "score"),
                                                        DEF.poly3 = list(G = "auto")))))
  ref_sym <- quiet(def.gof(fit, G = Ge, basis = "sym", weights = "score"))
  ref_p3  <- quiet(def.gof(fit, G = Ge, basis = "poly3"))
  ok_bat <- identical(bat$Statistic[bat$Test == "DEF.sym"], ref_sym$Test_Statistic) &&
            identical(bat$p_value[bat$Test == "DEF.poly3"], ref_p3$p_value) &&
            all(grepl(sprintf("G = %d (auto)", as.integer(Ge)), bat$Note, fixed = TRUE)) &&
            grepl("score form", bat$Note[bat$Test == "DEF.sym"])
  cat(sprintf("  n %5d  expected G %3d  def.gof %s  edge %s  ensemble %s  (y,ph,X) %s  battery %s  Notes: %s\n", n, as.integer(Ge),
              ok, ok_edge, ok_ens, ok_yph, ok_bat, paste(bat$Note, collapse = " / ")))
  addlog("d", sprintf("n = %d", n), sprintf("G = 'auto' identical to G = %d in def.gof/edge/ensemble/(y,ph,X)/battery", as.integer(Ge)),
         ok && ok_edge && ok_ens && ok_yph && ok_bat)
}
e8 <- count_warn(def.gof(mk(9, 8), G = "auto"))
cat("  n = 8, G = 'auto':", if (inherits(e8$value, "caught_error")) paste("error:", e8$value$msg) else "no error", "\n")
top <- tryCatch(as.data.frame(suppressWarnings(run.all.gof(A$base, G = "auto", tests = c("HL", "DEF.sym"), install = "no"))),
                error = function(e) paste("error:", conditionMessage(e)))
cat("  run.all.gof(fit, G = 'auto') at top level:\n"); print(top)

cat("\n==== (d) few-events warning ====\n")
## sweep every event count with the (y, ph, X) interface
sweep <- function(n, G, seed) {
  set.seed(seed); xa <- runif(n, -3, 3); X <- cbind(1, xa); ph <- plogis(0.6 * xa)
  Gn <- Gof(n, G); mism <- 0L; nerr <- 0L; multi <- 0L
  for (k in 0:n) {
    yv <- numeric(n); if (k > 0) yv[sample(n, k, prob = ph)] <- 1
    w <- count_warn(def.gof(yv, ph, X = X, G = G, basis = "poly3"))
    if (inherits(w$value, "caught_error")) nerr <- nerr + 1L
    if ((w$few > 0) != (min(k, n - k) < Gn)) mism <- mism + 1L
    if (w$few > 1) multi <- multi + 1L
  }
  cat(sprintf("  sweep n %d G %s (= %d): %d event counts, mismatches %d, calls warning more than once %d, errors %d\n",
              n, format(G), as.integer(Gn), n + 1, mism, multi, nerr))
  addlog("d", sprintf("sweep n %d G %s", n, format(G)), "warning iff min(events, non-events) < G, every k in 0..n",
         mism == 0 && multi == 0, sprintf("errors (still counted): %d", nerr))
}
sweep(60, 10, 11); sweep(60, "auto", 12); sweep(300, 13, 13); sweep(1000, "auto", 14)

## glm path at the boundary, every entry point
fit_k <- function(n, k, seed) {
  set.seed(seed); dat <- data.frame(xa = runif(n, -3, 3)); yv <- numeric(n)
  yv[sample(n, k, prob = plogis(0.6 * dat$xa))] <- 1; dat$out <- yv
  suppressWarnings(glm(out ~ xa, family = binomial(), data = dat))
}
for (cs in list(list(n = 400, G = 10), list(n = 400, G = 20), list(n = 1000, G = "auto"))) {
  Gn <- Gof(cs$n, cs$G)
  for (k in c(Gn - 1, Gn, Gn + 1, cs$n - Gn - 1, cs$n - Gn, cs$n - Gn + 1)) {
    fit <- fit_k(cs$n, k, 900 + k); expect <- min(k, cs$n - k) < Gn
    cnt <- integer(0)
    for (b in c("poly2", "poly3", "stukel", "sym")) for (wt in c("unit", "score"))
      cnt <- c(cnt, count_warn(def.gof(fit, G = cs$G, basis = b, weights = wt))$few)
    w1  <- count_warn(def.gof(fit, G = cs$G))
    we  <- count_warn(edge.gof(fit, G = cs$G))$few
    wen <- count_warn(def.ensemble.gof(fit, G = cs$G, add_ef = TRUE))$few
    wen4 <- count_warn(def.ensemble.gof(fit, G = cs$G, components = c("poly2", "poly3", "stukel", "sym")))$few
    bt  <- count_warn(as.data.frame(run.all.gof(fit, tests = c("DEF.poly3", "DEF.sym"), G = Gn, install = "no",
                                                control = list(DEF.sym = list(G = cs$G)))))
    notes_ok <- all(grepl("than groups", bt$value$Note) == expect)
    what <- if (k <= cs$n - k) "events" else "non-events"
    msg_ok <- if (expect) grepl(sprintf("%d %s for G = %d groups", as.integer(min(k, cs$n - k)), what, as.integer(Gn)),
                                w1$msgs[1], fixed = TRUE) else length(w1$msgs) == 0
    ok <- all(cnt == as.integer(expect)) && we == as.integer(expect) && wen == as.integer(expect) &&
          wen4 == as.integer(expect) && bt$few == 0L && notes_ok && msg_ok
    cat(sprintf("  n %4d G %-4s k %4d  expect %-5s | def.gof counts %s | edge %d | ensemble %d / 4-basis %d | battery escaped %d, Note %s | msg %s\n",
                cs$n, format(cs$G), as.integer(k), expect, paste(cnt, collapse = ""), we, wen, wen4, bt$few, notes_ok, msg_ok))
    addlog("d", sprintf("glm n %d G %s k %d", cs$n, format(cs$G), as.integer(k)),
           "warning count 1/0 in def.gof (4 bases x 2 weights), edge.gof, def.ensemble.gof; battery puts it in Note", ok,
           if (length(w1$msgs)) w1$msgs[1] else "")
  }
}
## the full fast battery must not leak the warning (DEF rows to Note, ensemble rows silent)
for (k in c(9, 10)) {
  fit <- fit_k(400, k, 950 + k)
  bt <- count_warn(as.data.frame(run.all.gof(fit, include_slow = FALSE, install = "no")))
  nn <- bt$value$Note[grepl("^DEF", bt$value$Test)]
  cat(sprintf("  full fast battery n 400 G 10 k %d: def_few_events escaped %d, other warnings %d, DEF Notes: %s | ensemble p: %s\n",
              k, bt$few, bt$other, paste(unique(nn), collapse = " / "),
              paste(format(bt$value$p_value[bt$value$Family == "Ensemble"], digits = 4), collapse = ", ")))
  addlog("d", sprintf("full battery n 400 k %d", k), "no def_few_events warning escapes run.all.gof; DEF Note iff k < G",
         bt$few == 0L && all(grepl("than groups", nn) == (k < 10)))
}

## =================================================================================================================
CH <- do.call(rbind, ROWS)
CH$pkg_version <- PKGV
write.csv(CH, file.path(OUT, "pkg280_verify_num_exact_checks.csv"), row.names = FALSE)
cat("\n==== summary ====\n")
info <- grepl("informational", CH$note)
cat(sprintf("checks: %d | pass %d | fail %d (of which informational %d)\n", nrow(CH), sum(CH$pass), sum(!CH$pass),
            sum(!CH$pass & info)))
agg <- aggregate(cbind(n = 1, pass = CH$pass, max_abs = ifelse(is.na(CH$abs_diff), 0, CH$abs_diff),
                       max_rel = ifelse(is.na(CH$rel_diff), 0, CH$rel_diff)) ~ part + sub("^(sym|poly3|poly2|stukel) ", "", quantity),
                 data = CH, FUN = function(v) c(sum = sum(v), max = max(v)))
print(CH[!CH$pass, c("part", "design", "quantity", "pkg", "ref", "abs_diff", "rel_diff", "tol", "note")], row.names = FALSE)
for (q in unique(sub("^(sym|poly3|poly2|stukel) ", "", CH$quantity))) {
  s <- CH[sub("^(sym|poly3|poly2|stukel) ", "", CH$quantity) == q, ]
  cat(sprintf("  [%s] %-80s n %3d pass %3d  max|diff| %.2e  max rel %.2e\n", s$part[1], q, nrow(s), sum(s$pass),
              suppressWarnings(max(s$abs_diff, na.rm = TRUE)), suppressWarnings(max(s$rel_diff, na.rm = TRUE))))
}
cat("\ndone", format(Sys.time()), "\n")
