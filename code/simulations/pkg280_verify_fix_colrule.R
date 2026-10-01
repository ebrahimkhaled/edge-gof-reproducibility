## pkg280_verify_fix_colrule.R -- c44034b drops a def.gof basis column whose length is below 1e-6 times the longest
## one. That also moves results that did not stop before (a Stukel half-column reaching one group whose mean risk is
## within about 0.0016 of one half). This probe compares, on the same fits:
##   pre   = 0d63fba (2.7.0 filter colSums(abs(Z)) > 1e-8, no scaling)
##   fixed = 906bd0c (relative 1e-6 rule, scaled columns)
##   alt   = the 2.7.0 filter kept, and the kept columns scaled to unit length before any solve
## alt is built in memory from the fixed R/def_gof.R with those three lines replaced; the package is not edited.
## Expectation for alt: equal to pre wherever pre ran, finite wherever pre stopped.
## Run: Rscript pkg280_verify_fix_colrule.R [package dir] > ../paper_EDGE/theory/pkg280_verify_fix_colrule.log 2>&1
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
cat("ebrahim.gof", as.character(utils::packageVersion("ebrahim.gof")), "via load_all | HEAD", sha, "|", format(Sys.time()), "\n")

KEEP <- c("def.gof", ".def_basis", ".def_pvalue")
pre <- new.env(parent = ns)
for (ex in parse(text = system2("git", c("-C", shQuote(PKG), "show", "0d63fba:R/def_gof.R"), stdout = TRUE), keep.source = FALSE))
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) && as.character(ex[[2]]) %in% KEEP) eval(ex, pre)
txt <- readLines(file.path(PKG, "R", "def_gof.R"))
i <- grep("nz <- sqrt(colSums(Z^2))", txt, fixed = TRUE)
stopifnot(length(i) == 1, grepl("ok <- nz > 1e-6 * max(nz)", txt[i + 1], fixed = TRUE),
          grepl("Z  <- Z[, ok, drop = FALSE] / rep(nz[ok], each = nrow(Z))", txt[i + 2], fixed = TRUE))
txt[i]     <- "  Z  <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]"
txt[i + 1] <- "  Z  <- Z / rep(sqrt(colSums(Z^2)), each = nrow(Z))"
txt <- txt[-(i + 2)]
alt <- new.env(parent = ns)
for (ex in parse(text = txt, keep.source = FALSE))
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) && as.character(ex[[2]]) %in% KEEP) eval(ex, alt)
cat("alt column rule:\n", paste(txt[i:(i + 1)], collapse = "\n"), "\n")

SRC  <- file.path(SIM, "pkg280_verify_num_exact.R")
TAKE <- c("mk", "augment", "score_qr", "Gof", "step_cov", "count_warn", "quiet")
for (e in parse(SRC, keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% TAKE)
    eval(e, globalenv())
warm <- function(fit) {
  env <- environment(formula(fit))
  cl <- getCall(fit); cl$control <- glm.control(epsilon = 1e-13, maxit = 100)
  ft <- suppressWarnings(eval(cl, env))
  cl$start <- unname(coef(ft))
  suppressWarnings(eval(cl, env))
}
val <- function(expr) {
  r <- tryCatch(suppressWarnings(expr), error = function(e) NULL)
  if (is.null(r)) c(S = NA_real_, df = NA_real_, p = NA_real_) else c(S = r$Test_Statistic, df = r$df, p = r$p_value)
}
show <- function(x) if (is.na(x["S"])) "stops" else sprintf("S %.8f df %.6f p %.8f", x["S"], x["df"], x["p"])
three <- function(call_fun) list(pre = val(call_fun(pre$def.gof)), fixed = val(call_fun(def.gof)), alt = val(call_fun(alt$def.gof)))
same <- function(a, b) !anyNA(c(a, b)) && all(abs(a - b) <= 1e-10 * pmax(1, abs(a)))

## ---- A. the R2 construction, top group mean risk v ----
cat("\n==== A. R2 construction (n = 500, G = 20), top group mean risk v ====\n")
for (v in 0.5 + c(1e-6, 5e-5, 1e-4, 1.7e-4, 5e-4, 1e-3, 1.5e-3, 2e-3)) {
  set.seed(4); ph <- c(sort(runif(475, 0.02, 0.45)), rep(v, 25)); X <- cbind(1, qlogis(ph)); yy <- rbinom(500, 1, ph)
  u <- three(function(f) f(yy, ph, X = X, G = 20, basis = "stukel"))
  s <- three(function(f) f(yy, ph, X = X, G = 20, basis = "stukel", weights = "score"))
  cat(sprintf("  v = %-9s unit : pre %-44s | fixed %-44s | alt %s\n", format(v, digits = 7), show(u$pre), show(u$fixed), show(u$alt)))
  cat(sprintf("  %-13s score: pre %-44s | fixed %-44s | alt %s\n", "", show(s$pre), show(s$fixed), show(s$alt)))
}

## ---- B. glm fits in the window, with anova Rao for the score form ----
cat("\n==== B. c0 = -2, n = 500, G = 20 fits from pkg280_verify_fix_repro.R (warm-restarted) ====\n")
for (sd in c(392, 423, 588, 605)) {
  fw <- warm(mk(sd, 500, c0 = -2)); sc <- step_cov(fw, 20, "stukel")
  a <- anova(fw, augment(fw, sc), test = "Rao")
  u <- three(function(f) f(fw, G = 20, basis = "stukel"))
  s <- three(function(f) f(fw, G = 20, basis = "stukel", weights = "score"))
  cat(sprintf("  seed %d unit : pre %-44s | fixed %-44s | alt %s\n", sd, show(u$pre), show(u$fixed), show(u$alt)))
  cat(sprintf("  %-8s score: pre %-44s | fixed %-44s | alt %s\n", "", show(s$pre), show(s$fixed), show(s$alt)))
  cat(sprintf("  %-8s anova Rao %.8f Df %g p %.8f | |alt score - anova| %.2e\n", "", a$Rao[2], a$Df[2], a$`Pr(>Chi)`[2],
              abs(s$alt["S"] - a$Rao[2])))
}

## ---- C. sweeps ----
gen_data <- function(n, c0, s) {
  dat <- data.frame(xa = runif(n, -3, 3), db = rbinom(n, 1, 0.5))
  dat$out <- rbinom(n, 1, plogis(c0 + s * (0.6 * dat$xa + 0.5 * dat$db)))
  dat
}
fit_of <- function(dat) suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat))
sweep <- function(label, D, Gs) {
  for (G in Gs) {
    cnt <- list()
    add <- function(k) cnt[[k]] <<- (if (is.null(cnt[[k]])) 0 else cnt[[k]]) + 1
    changed <- character(0)
    for (i in seq_along(D)) {
      f <- fit_of(D[[i]])
      for (form in c("unit satterthwaite", "unit imhof", "score")) {
        cf <- switch(form,
          "unit satterthwaite" = function(g) g(f, G = G, basis = "stukel"),
          "unit imhof"         = function(g) g(f, G = G, basis = "stukel", method = "imhof"),
          "score"              = function(g) g(f, G = G, basis = "stukel", weights = "score"))
        r <- three(cf)
        if (is.na(r$pre["S"])) {
          add(paste(form, "| pre stops"))
          add(paste(form, "| pre stops, fixed finite:", !is.na(r$fixed["S"])))
          add(paste(form, "| pre stops, alt finite:", !is.na(r$alt["S"])))
        } else {
          add(paste(form, "| pre runs"))
          if (!same(r$pre, r$fixed)) { add(paste(form, "| pre runs, fixed differs")); if (form != "unit imhof")
            changed <- c(changed, sprintf("    fit %4d %-18s pre %-44s fixed %s", i, form, show(r$pre), show(r$fixed))) }
          if (!same(r$pre, r$alt)) add(paste(form, "| pre runs, alt differs"))
        }
      }
    }
    cat(sprintf("\n  %s, G = %d, %d fits:\n", label, G, length(D)))
    for (k in sort(names(cnt))) cat(sprintf("    %-58s %d\n", k, cnt[[k]]))
    if (length(changed)) { cat("  cases where the fixed build moves a result that ran before (first 12):\n")
                           cat(head(changed, 12), sep = "\n") }
  }
}
cat("\n==== C. sweeps: stukel basis, pre vs fixed vs alt ====\n")
set.seed(880004); D1 <- replicate(2000, gen_data(500, -2, 1), simplify = FALSE)   # the engine's 'extra' data sets
sweep("c0 = -2, n = 500 (engine 'extra' data)", D1, c(20, 10))
set.seed(880005); D2 <- replicate(1000, gen_data(500, 2, 1), simplify = FALSE)
sweep("c0 = +2, n = 500", D2, 20)
cat("\ndone", format(Sys.time()), "\n")
