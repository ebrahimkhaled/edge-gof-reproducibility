## pkg280_verify_corr_flatlimit.R -- when the information guard removes a joint Stukel statistic, was it rounding noise?
## For each sample: the pre-guard statistic (2e0f6d1 gof_stukel) on the default glm fit and on a warm-restarted fit
## (epsilon 1e-13); the same u'I^-1 u by a residualised QR route on the warm fit (u = Mr'(y - p)/sqrt(W), I = Mr'Mr, with
## Mr = sqrt(W) Z minus its projection on sqrt(W) X), which avoids both cancellations of the normal-equation route; and
## anova(warm, warm + Stukel columns, test = "Rao"). When the fitted slope is near zero and every risk is below one half,
## the Stukel column -eta^2/2 is, after removing (1, x), proportional to x^2, so the statistic tends to the score test for
## adding x^2: that limit is also given.
## Samples: the near-flat constructions of pkg280_verify_corr_extra.R, and every joint-Stukel firing of the guard sweep
## (scratchpad guard_partB.rds, regenerated from its seed) when that file exists.
## Run: Rscript pkg280_verify_corr_flatlimit.R <clone dir> > ../paper_EDGE/theory/pkg280_verify_corr_flatlimit.log 2>&1
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
PKG <- args[1]
SIM <- edge_path("code/simulations")
SCR <- file.path(tempdir(), "vcorr")
options(width = 230, warn = 1)
suppressPackageStartupMessages(library(ebrahim.gof))
source(file.path(SIM, "_dgp_library.R"))
ns <- asNamespace("ebrahim.gof")
PRES <- new.env(parent = ns)
for (ex in parse(text = system2("git", c("-C", shQuote(PKG), "show", "2e0f6d1:R/run_all_gof.R"), stdout = TRUE), keep.source = FALSE))
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) && as.character(ex[[2]]) == "gof_stukel") eval(ex, PRES)
cat("installed ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "|", format(Sys.time()), "\n")

warm <- function(fo, d) {
  ctl <- glm.control(epsilon = 1e-13, maxit = 200)
  f0 <- suppressWarnings(glm(fo, family = binomial(), data = d, control = ctl))
  suppressWarnings(glm(fo, family = binomial(), data = d, control = ctl, start = unname(coef(f0))))
}
stk_cols <- function(fit) {
  ph <- pmin(pmax(fitted(fit), 1e-6), 1 - 1e-6); eta <- fit$linear.predictors
  Z <- cbind(za = 0.5 * eta^2 * (ph >= 0.5), zb = -0.5 * eta^2 * (ph < 0.5)); Z[, colSums(Z != 0) > 0, drop = FALSE]
}
qr_route <- function(fit) {
  p <- as.numeric(fitted(fit)); W <- p * (1 - p); sw <- sqrt(W); X <- model.matrix(fit); Z <- stk_cols(fit)
  Mr <- qr.resid(qr(sw * X), sw * Z)
  u  <- crossprod(Mr, (fit$y - p) / sw)
  list(S = drop(crossprod(u, solve(crossprod(Mr), u))), df = ncol(Z),
       ratio = min(colSums(Mr^2) / colSums(W * Z^2)))
}
rao_stk <- function(fit, d, fo) {
  Z <- stk_cols(fit); d2 <- d; nm <- colnames(Z); for (k in nm) d2[[k]] <- Z[, k]
  fo2 <- update(fo, as.formula(paste(". ~ . +", paste(nm, collapse = " + "))))
  big <- suppressWarnings(glm(fo2, family = binomial(), data = d2, control = fit$control))
  a <- anova(fit, big, test = "Rao"); c(S = a$Rao[2], df = a$Df[2])
}
rao_x2 <- function(fit, d, fo) {
  d2 <- d; d2$x2_ <- d$x^2
  big <- suppressWarnings(glm(update(fo, . ~ . + x2_), family = binomial(), data = d2, control = fit$control))
  a <- anova(fit, big, test = "Rao"); c(S = a$Rao[2], df = a$Df[2])
}
one <- function(label, d, fo) {
  f  <- suppressWarnings(glm(fo, family = binomial(), data = d)); fw <- warm(fo, d)
  s0 <- PRES$gof_stukel(ns$.gof_context(f)); s0w <- PRES$gof_stukel(ns$.gof_context(fw))
  s1 <- ns$gof_stukel(ns$.gof_context(f))
  qr <- qr_route(fw); ra <- tryCatch(rao_stk(fw, d, fo), error = function(e) c(S = NA, df = NA))
  x2 <- tryCatch(rao_x2(fw, d, fo), error = function(e) c(S = NA, df = NA))
  data.frame(sample = label, events = sum(d$y), slope = unname(coef(f)["x"]), ratio_qr = qr$ratio,
             pre_S = as.numeric(s0$Statistic), pre_S_warm = as.numeric(s0w$Statistic), qr_S_warm = qr$S, rao_S_warm = ra[["S"]],
             x2_S = x2[["S"]], df = qr$df, pre_p = as.numeric(s0$p_value), qr_p = pchisq(qr$S, qr$df, lower.tail = FALSE),
             new_p = as.numeric(s1$p_value))
}

rows <- list()
x0 <- qnorm(ppoints(200))
for (ev in list(c(1, 200), c(1, 100, 200))) for (delta in c(0, 1e-6, 1e-3, 1e-2, 1e-1, 0.3)) {
  x <- x0; y <- rep(0, 200); y[ev] <- 1; x[ev[1]] <- x[ev[1]] + delta
  rows[[length(rows) + 1]] <- one(sprintf("rows %s delta %g", paste(ev, collapse = ","), delta), data.frame(x = x, y = y), y ~ x)
}
cat("\n==== constructed near-flat samples ====\n")
print(do.call(rbind, rows), digits = 6, row.names = FALSE)

rds <- file.path(SCR, "guard_partB.rds")
if (file.exists(rds)) {
  PB <- do.call(rbind, readRDS(rds))
  Fz <- PB[PB$test == "stukel.joint" & PB$fire_obs %in% TRUE, ]
  cat(sprintf("\n==== joint-Stukel firings of the guard sweep (%d) ====\n", nrow(Fz)))
  gen <- function(design, r) switch(design,
    sparse100 = { set.seed(79000000 + r); gen_sparse_link(100, -4.9) }, sparse300 = { set.seed(75000000 + r); gen_sparse_link(300, -4.9) },
    sparse500 = { set.seed(76000000 + r); gen_sparse_link(500, -4.9) }, base1000 = { set.seed(77000000 + r); dgp_null("link", 1000) })
  rr <- list()
  for (k in seq_len(nrow(Fz))) {
    g <- gen(Fz$design[k], Fz$r[k])
    rr[[k]] <- one(sprintf("%s r=%d", Fz$design[k], Fz$r[k]), g$d, g$f)
  }
  if (length(rr)) {
    R <- do.call(rbind, rr)
    print(R, digits = 5, row.names = FALSE)
    rel <- function(a, b) abs(a - b) / pmax(1e-12, abs(b))
    cat(sprintf("\nrelative |pre - qr (warm)|: median %.2e, max %.2e; relative |qr - anova Rao| (warm): median %.2e, max %.2e\n",
                median(rel(R$pre_S, R$qr_S_warm), na.rm = TRUE), max(rel(R$pre_S, R$qr_S_warm), na.rm = TRUE),
                median(rel(R$qr_S_warm, R$rao_S_warm), na.rm = TRUE), max(rel(R$qr_S_warm, R$rao_S_warm), na.rm = TRUE)))
    cat(sprintf("stable routes (qr and Rao within 1e-3): %d of %d; among them qr p < 0.2: %d, qr p < 0.05: %d\n",
                sum(rel(R$qr_S_warm, R$rao_S_warm) < 1e-3, na.rm = TRUE), nrow(R),
                sum(rel(R$qr_S_warm, R$rao_S_warm) < 1e-3 & R$qr_p < 0.2, na.rm = TRUE),
                sum(rel(R$qr_S_warm, R$rao_S_warm) < 1e-3 & R$qr_p < 0.05, na.rm = TRUE)))
  }
}
cat("\ndone", format(Sys.time()), "\n")
