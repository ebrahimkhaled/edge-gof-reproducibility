## pkg280_verify_fix_imhof.R -- size of the differences in Imhof p-values after the columns of Z are scaled to unit
## length (906bd0c, and the alternative rule of pkg280_verify_fix_colrule.R), against 0d63fba. The Satterthwaite p and
## the statistic agree to about 1e-14; this shows whether the Imhof differences are only integration noise, and
## compares them with CompQuadForm::imhof's own accuracy (epsabs = epsrel = 1e-6 by default).
## Run: Rscript pkg280_verify_fix_imhof.R [package dir] > ../paper_EDGE/theory/pkg280_verify_fix_imhof.log 2>&1
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
options(width = 200, warn = 1)
suppressPackageStartupMessages(pkgload::load_all(PKG, export_all = FALSE, quiet = TRUE))
ns <- asNamespace("ebrahim.gof")
cat("HEAD", system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE),
    "| CompQuadForm", as.character(packageVersion("CompQuadForm")), "\n")
print(formals(CompQuadForm::imhof)[c("epsabs", "epsrel", "limit")])
cat(".def_pvalue:\n"); print(ns$.def_pvalue)
KEEP <- c("def.gof", ".def_basis", ".def_pvalue")
pre <- new.env(parent = ns)
for (ex in parse(text = system2("git", c("-C", shQuote(PKG), "show", "0d63fba:R/def_gof.R"), stdout = TRUE), keep.source = FALSE))
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) && as.character(ex[[2]]) %in% KEEP) eval(ex, pre)
txt <- readLines(file.path(PKG, "R", "def_gof.R"))
i <- grep("nz <- sqrt(colSums(Z^2))", txt, fixed = TRUE)
txt[i] <- "  Z  <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]"
txt[i + 1] <- "  Z  <- Z / rep(sqrt(colSums(Z^2)), each = nrow(Z))"
txt <- txt[-(i + 2)]
alt <- new.env(parent = ns)
for (ex in parse(text = txt, keep.source = FALSE))
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) && as.character(ex[[2]]) %in% KEEP) eval(ex, alt)

gen_data <- function(n, c0, s) {
  dat <- data.frame(xa = runif(n, -3, 3), db = rbinom(n, 1, 0.5))
  dat$out <- rbinom(n, 1, plogis(c0 + s * (0.6 * dat$xa + 0.5 * dat$db)))
  dat
}
fit_of <- function(dat) suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat))
p_of <- function(f, fit, ...) tryCatch(suppressWarnings(f(fit, method = "imhof", ...))$p_value, error = function(e) NA_real_)
s_of <- function(f, fit, ...) tryCatch(suppressWarnings(f(fit, method = "imhof", ...))$Test_Statistic, error = function(e) NA_real_)
report <- function(label, D, G, basis) {
  dp_fixed <- dp_alt <- dS_alt <- rep(NA_real_, length(D))
  for (k in seq_along(D)) {
    f  <- fit_of(D[[k]])
    p0 <- p_of(pre$def.gof, f, G = G, basis = basis)
    if (is.na(p0)) next
    dp_fixed[k] <- abs(p_of(def.gof, f, G = G, basis = basis) - p0)
    dp_alt[k]   <- abs(p_of(alt$def.gof, f, G = G, basis = basis) - p0)
    dS_alt[k]   <- abs(s_of(alt$def.gof, f, G = G, basis = basis) - s_of(pre$def.gof, f, G = G, basis = basis))
  }
  q <- function(x) paste(format(quantile(x, c(0.5, 0.99, 1), na.rm = TRUE), digits = 3), collapse = " / ")
  cat(sprintf("  %-22s %-6s G %2d | Imhof |p - pre|, median / 99%% / max: fixed %s | alt %s | alt > 1e-10: %d | alt |S - pre| max %.2e\n",
              label, basis, G, q(dp_fixed), q(dp_alt), sum(dp_alt > 1e-10, na.rm = TRUE), max(dS_alt, na.rm = TRUE)))
}
set.seed(880005); D2 <- replicate(1000, gen_data(500, 2, 1), simplify = FALSE)
report("c0 = +2, n 500", D2, 20, "stukel")
set.seed(880001); D0 <- replicate(200, gen_data(1000, 0, 1), simplify = FALSE)
report("base, n 1000", D0, 10, "sym")
report("base, n 1000", D0, 40, "sym")
report("base, n 1000", D0, 40, "stukel")
cat("done\n")
