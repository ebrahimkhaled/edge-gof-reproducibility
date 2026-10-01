## pkg280_verify_num_followup.R -- follow-ups to pkg280_verify_num_exact.R and pkg280_verify_num_size.R, using the
## INSTALLED ebrahim.gof 2.8.0. Verification of the package, not a paper result.
##   A. marginal Stukel vs LogisticDx::gof SstBoth on every (a) fixture, refitted at global scope
##   B. the two_above fixture, where form = "lr" returned NA
##   C. a rank-deficient glm: joint vs lr vs anova Rao
##   D. the one replicate of the size run (c0 = -2, n = 500, G = auto) where DEF stukel unit gave no p-value
## Run: Rscript pkg280_verify_num_followup.R > ../paper_EDGE/theory/pkg280_verify_num_followup.log 2>&1
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
PKGV <- as.character(packageVersion("ebrahim.gof")); stopifnot(PKGV == "2.8.0")
PKG <- Sys.getenv("EBRAHIM_GOF_SRC")
OUT <- edge_path("declarations")
SRC <- edge_path("code/simulations/pkg280_verify_num_exact.R")
options(width = 200, warn = 1)
cat("ebrahim.gof", PKGV, "(installed) |", R.version.string, "|", format(Sys.time()), "\n")
data("gof_demo", package = "ebrahim.gof")
TAKE <- c("mk", "find_fit", "augment", "score_qr", "stukel_Z", "stk", "count_warn", "quiet", "demo_wrong", "A")
for (e in parse(SRC, keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% TAKE)
    eval(e, globalenv())
old <- new.env(parent = asNamespace("ebrahim.gof"))
for (f in c("R/def_gof.R", "R/def_ensemble_gof.R", "R/run_all_gof.R")) {
  txt <- system2("git", c("-C", shQuote(PKG), "show", paste0("main:", f)), stdout = TRUE)
  for (e in parse(text = txt, keep.source = FALSE))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% c("def.gof", ".def_basis", ".def_pvalue", "gof_stukel", ".gof_score_z", ".gof_context"))
      eval(e, old)
}

## ---- A. LogisticDx ------------------------------------------------------------------------------------------------
cat("\n==== A. marginal Stukel vs LogisticDx::gof SstBoth ====\n")
for (nm in names(A)) {
  fit <- A[[nm]]
  why <- tryCatch({ suppressMessages(suppressWarnings(LogisticDx::gof(fit, plotROC = FALSE))); "ran" },
                  error = function(e) conditionMessage(e))
  Dn <- paste0("D_", nm); assign(Dn, fit$data, envir = globalenv())
  fo <- formula(fit); environment(fo) <- globalenv()
  g <- suppressWarnings(eval(call("glm", formula = fo, family = quote(binomial()), data = as.name(Dn)), globalenv()))
  same_fit <- identical(unname(coef(g)), unname(coef(fit))) && identical(unname(fitted(g)), unname(fitted(fit)))
  lx <- tryCatch(as.data.frame(suppressMessages(suppressWarnings(LogisticDx::gof(g, plotROC = FALSE)))$gof),
                 error = function(e) conditionMessage(e))
  m <- stk(g, "marginal"); o <- old$gof_stukel(old$.gof_context(g, NULL, NULL, G = 10))
  if (is.data.frame(lx)) {
    v <- lx$val[lx$test == "SstBoth"]; p <- lx$pVal[lx$test == "SstBoth"]
    cat(sprintf("  %-10s original fit: %-45s | global refit identical coef/fitted %s | pkg marginal %.10f p %.8g | 2.7.0 %.10f | SstBoth %.10f p %.8g | |diff| %.2e\n",
                nm, substr(why, 1, 45), same_fit, m$Statistic, m$p_value, o$Statistic, v, p, abs(m$Statistic - v)))
  } else cat(sprintf("  %-10s LogisticDx failed on the global refit too: %s\n", nm, lx))
}

## ---- B. two_above, form = "lr" -----------------------------------------------------------------------------------
cat("\n==== B. two_above: form = 'lr' ====\n")
fit <- A$two_above
print(stk(fit, "lr")); print(stk(fit, "joint"))
Z <- stukel_Z(fit); X <- model.matrix(fit); keep <- colSums(Z != 0) > 0
cat("  observations with fitted >= 0.5:", sum(fitted(fit) >= 0.5), "| their y:", fit$y[fitted(fit) >= 0.5], "\n")
for (mx in c(25, 100, 1000)) {
  f1 <- suppressWarnings(glm.fit(cbind(X, Z[, keep, drop = FALSE]), fit$y, family = binomial(),
                                 control = glm.control(maxit = mx)))
  cat(sprintf("  glm.fit maxit %4d: converged %s, iter %d, deviance %.8f, drop %.8f, coef %s, max fitted of the 2 obs %.10f\n",
              mx, f1$converged, f1$iter, f1$deviance, fit$deviance - f1$deviance,
              paste(format(coef(f1), digits = 4), collapse = " "), max(f1$fitted.values[fitted(fit) >= 0.5])))
}
f1g <- augment(fit, Z)
cat("  the reference augment() fit in the exact script: converged", f1g$converged, "iter", f1g$iter, "\n")

## ---- C. rank-deficient glm ------------------------------------------------------------------------------------------
cat("\n==== C. rank-deficient glm (xa2 = 2 xa) ====\n")
set.seed(3); gd <- data.frame(xa = runif(300, -3, 3)); gd$out <- rbinom(300, 1, pnorm(0.5 * gd$xa)); gd$xa2 <- 2 * gd$xa
rankdef <- glm(out ~ xa + xa2, family = binomial(), data = gd)
full    <- glm(out ~ xa, family = binomial(), data = gd)
cat("  coef(rankdef):", format(coef(rankdef)), "| rank", rankdef$rank, "\n")
for (nm in c("rankdef", "full")) {
  f <- get(nm); j <- stk(f, "joint"); l <- stk(f, "lr"); m <- stk(f, "marginal")
  a <- anova(f, augment(f, stukel_Z(f)), test = "Rao")
  cat(sprintf("  %-8s joint %s (df %s, Note '%s') | lr %.8f df %g | marginal %.8f | anova Rao %.8f df %g | score_qr %.8f\n",
              nm, format(j$Statistic), format(j$df), j$Note, l$Statistic, l$df, m$Statistic, a$Rao[2], a$Df[2],
              score_qr(f, stukel_Z(f))$stat))
}
cat("  def.gof(rankdef):", tryCatch(format(def.gof(rankdef)$p_value), error = function(e) paste("error:", conditionMessage(e))), "\n")
cat("  2.7.0 def.gof(rankdef):", tryCatch(format(old$def.gof(rankdef)$p_value), error = function(e) paste("error:", conditionMessage(e))), "\n")
cat("  2.7.0 marginal Stukel(rankdef):", format(old$gof_stukel(old$.gof_context(rankdef, NULL, NULL, G = 10))$Statistic), "\n")

## ---- D. the missing DEF stukel unit p-value in the size run ----------------------------------------------------------
cat("\n==== D. size run, c0 = -2, n = 500, G = auto: replicates without a stukel unit p-value ====\n")
P <- read.csv(file.path(OUT, "pkg280_verify_num_size_pvalues.csv"), check.names = FALSE)
cellk <- 3L; lab <- "c0=-2 n500 Gauto"
bad <- P[P$cell == lab & !is.finite(P$stukel.unit), ]
print(bad[, c("cell", "rep", "ev_min", "stk.joint", "stk.joint.df", "sym.unit", "poly3.unit", "stukel.unit", "stukel.score", "stukel.score.df")])
## rebuild the worker streams exactly as clusterSetRNGStream(cl, 20260914 + 100 k) and parLapply over 6 nodes
RNGkind("L'Ecuyer-CMRG"); set.seed(20260914L + 100L * cellk)
seeds <- list(.Random.seed); for (i in 2:6) seeds[[i]] <- parallel::nextRNGStream(seeds[[i - 1]])
chunks <- parallel::splitIndices(2000, 6)
regen <- function(r) {
  j <- which(vapply(chunks, function(ch) r %in% ch, logical(1)))
  assign(".Random.seed", seeds[[j]], envir = globalenv())
  for (q in chunks[[j]]) {
    xa <- runif(500, -3, 3); db <- rbinom(500, 1, 0.5); resp <- rbinom(500, 1, plogis(-2 + (0.6 * xa + 0.5 * db)))
    if (q == r) break
  }
  dd <- data.frame(resp = resp, xa = xa, db = db)
  suppressWarnings(glm(resp ~ xa + db, family = binomial(), data = dd))
}
for (r in bad$rep) {
  fit <- regen(r); row <- P[P$cell == lab & P$rep == r, ]
  chk <- c(sym.unit = quiet(def.gof(fit, G = "auto", basis = "sym"))$p_value,
           poly3.unit = quiet(def.gof(fit, G = "auto", basis = "poly3"))$p_value)
  cat(sprintf("  rep %d rebuilt: sym unit p %.10g (csv %.10g), poly3 unit p %.10g (csv %.10g)\n", r,
              chk[1], row$sym.unit, chk[2], row$poly3.unit))
  n <- length(fit$y); G <- max(10, round(n / 25)); ph <- fitted(fit)
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  pbar <- tapply(ph, grp, mean)
  cat("  group mean risks >= 0.5:", sum(pbar >= 0.5), "| top three pbar:", format(tail(sort(pbar), 3), digits = 6), "\n")
  cat("  2.8.0 def.gof stukel unit:\n"); print(tryCatch(def.gof(fit, G = "auto", basis = "stukel"), error = function(e) conditionMessage(e)))
  cat("  2.7.0 def.gof stukel unit (G = 20):\n"); print(tryCatch(old$def.gof(fit, G = 20, basis = "stukel"), error = function(e) conditionMessage(e)))
  cat("  2.8.0 def.gof stukel score:\n"); print(tryCatch(def.gof(fit, G = "auto", basis = "stukel", weights = "score"), error = function(e) conditionMessage(e)))
  cat("  2.8.0 def.gof stukel unit, imhof:\n"); print(tryCatch(def.gof(fit, G = "auto", basis = "stukel", method = "imhof"), error = function(e) conditionMessage(e)))
}
cat("\ndone", format(Sys.time()), "\n")
