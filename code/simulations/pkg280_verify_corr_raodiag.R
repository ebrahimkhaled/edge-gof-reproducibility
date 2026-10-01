## pkg280_verify_corr_raodiag.R -- the score forms that differ from anova Rao in pkg280_verify_corr_colrule.R, and the
## fits on which the installed build raised a warning. For each Rao mismatch: the top half-column of the grouped stukel
## basis (sum |z| over groups with mean risk >= 0.5), whether it is below the 2.7.0 filter (1e-8), and whether the Rao
## statistic without that column equals the installed score form.
## Run: Rscript pkg280_verify_corr_raodiag.R > ../paper_EDGE/theory/pkg280_verify_corr_raodiag.log 2>&1
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

SCR <- file.path(tempdir(), "vcorr")
options(width = 220, warn = 1)
suppressPackageStartupMessages(library(ebrahim.gof))
A <- readRDS(file.path(SCR, "colrule_all.rds"))$ALL
## the colrule workers ran under clusterSetRNGStream(), i.e. RNGkind("L'Ecuyer-CMRG"); set.seed() must use the same kind
## to regenerate their data sets. Each regenerated fit is checked against the stored top group mean.
RNGkind("L'Ecuyer-CMRG")

gen_data <- function(n, c0, s) {
  dat <- data.frame(xa = runif(n, -3, 3), db = rbinom(n, 1, 0.5))
  dat$out <- rbinom(n, 1, plogis(c0 + s * (0.6 * dat$xa + 0.5 * dat$db)))
  dat
}
data_of <- function(design, i) {
  if (design == "c0m2") { set.seed(71000000 + i); gen_data(500, -2, 1) }
  else if (design == "base") { set.seed(72000000 + i); gen_data(1000, 0, 1) }
  else { set.seed(73000000 + i); gen_data(500, -2, 1) }
}
fit_of <- function(dat, ctl = glm.control()) suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat, control = ctl))
warm_of <- function(dat) {
  ctl <- glm.control(epsilon = 1e-13, maxit = 100); f0 <- fit_of(dat, ctl)
  suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat, start = unname(coef(f0)), control = ctl))
}
grp <- function(fit, G) {
  n <- length(fit$y); ph <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6)
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G); list(g = g, pbar = as.numeric(tapply(ph, g, mean)))
}
rao_with <- function(fw, M) {
  e <- new.env(parent = globalenv()); e$out <- fw$y; e$Xm <- model.matrix(fw); e$M <- as.matrix(M)
  big <- suppressWarnings(glm(as.formula("out ~ Xm + M - 1", env = e), family = binomial(), control = fw$control))
  a <- anova(fw, big, test = "Rao"); c(S = a$Rao[2], df = a$Df[2])
}

cat("==== score forms that differ from anova Rao (df or relative |S - Rao| > 1e-8) ====\n")
B <- A[A$part %in% c("B", "Cw"), ]
B$rel_err <- abs(B$new.S - B$rao) / pmax(1, abs(B$rao))
M <- B[(B$new.df != B$rao_df) | B$rel_err > 1e-8, ]
out <- list()
for (k in seq_len(nrow(M))) {
  fw <- warm_of(data_of(M$design[k], M$i[k])); G <- M$G[k]
  gm <- grp(fw, G); e <- qlogis(gm$pbar)
  if (abs(max(gm$pbar) - M$top[k]) > 1e-12) stop(sprintf("regenerated fit %s %d G %d: top %.12f, stored %.12f",
                                                         M$design[k], M$i[k], G, max(gm$pbar), M$top[k]))
  Zg <- cbind(e, e^2 * (e >= 0), -e^2 * (e < 0)); sz <- colSums(abs(Zg))
  keep_pkg <- sz > 1e-8; keep_nz <- colSums(Zg != 0) > 0
  r_pkg <- rao_with(fw, Zg[gm$g, keep_pkg, drop = FALSE])
  sc <- suppressWarnings(def.gof(fw, G = G, basis = "stukel", weights = "score"))
  out[[k]] <- data.frame(part = M$part[k], design = M$design[k], i = M$i[k], G = G, top = max(gm$pbar),
                         groups_ge_half = sum(e >= 0), sum_abs_top_half = sz[2], below_filter = sz[2] <= 1e-8,
                         cols_nonzero = sum(keep_nz), cols_pkg = sum(keep_pkg), new_S = sc$Test_Statistic, new_df = sc$df,
                         rao_all_S = M$rao[k], rao_all_df = M$rao_df[k], rao_pkgcols_S = r_pkg[["S"]], rao_pkgcols_df = r_pkg[["df"]],
                         rel_err_pkgcols = abs(sc$Test_Statistic - r_pkg[["S"]]) / max(1, abs(r_pkg[["S"]])),
                         new_eq_pre = isTRUE(all.equal(M$new.S[k], M$pre.S[k], tolerance = 1e-9)) && M$new.df[k] == M$pre.df[k])
}
R <- do.call(rbind, out)
print(R, digits = 6, row.names = FALSE)
cat(sprintf("\nall mismatches below the 1e-8 filter: %s; Rao on the package's columns equals the installed score form (rel <= 1e-8) on %d of %d; new = 0d63fba on %d of %d\n",
            all(R$below_filter), sum(R$rel_err_pkgcols <= 1e-8 & R$new_df == R$rao_pkgcols_df), nrow(R), sum(R$new_eq_pre), nrow(R)))

cat("\n==== fits on which the installed build warned (part A) ====\n")
Wn <- A[A$part == "A" & A$new.nw > 0, c("design", "i", "G", "basis", "form", "pre.nw", "new.nw", "pre.p", "new.p")]
print(Wn, digits = 10, row.names = FALSE)
for (k in seq_len(nrow(Wn))) {
  f <- fit_of(data_of(Wn$design[k], Wn$i[k]))
  msgs <- character(0)
  withCallingHandlers(def.gof(f, G = Wn$G[k], basis = Wn$basis[k], method = "imhof"),
                      warning = function(w) { msgs <<- c(msgs, paste0("[", class(w)[1], "] ", conditionMessage(w))); invokeRestart("muffleWarning") })
  cat(sprintf("  %s fit %d G %d %s imhof: %s\n", Wn$design[k], Wn$i[k], Wn$G[k], Wn$basis[k], paste(msgs, collapse = " / ")))
}
cat("\ndone", format(Sys.time()), "\n")
