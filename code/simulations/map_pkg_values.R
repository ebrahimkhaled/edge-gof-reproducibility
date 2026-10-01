## map_pkg_values.R -- current outputs of the ebrahim.gof dev tree at every fixture a test, an example or a vignette
## uses for Stukel / EDGE, and the same fixtures under the PROTOTYPES of map_pkg_proto.R. Read-only. Single core.
## Run:  Rscript map_pkg_values.R > <scratchpad>/map_pkg_values.log 2>&1
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
source(edge_path("code/simulations/map_pkg_proto.R"))
options(width = 220, warn = 1)

fmt <- function(x) if (is.null(x) || length(x) == 0) "NULL" else formatC(as.numeric(x)[1], digits = 7, format = "g")
row <- function(label, stat, df, p, note = "")
  cat(sprintf("  %-62s stat=%-13s df=%-10s p=%-13s %s\n", label, fmt(stat), fmt(df), fmt(p), note))
dres <- function(label, d) row(label, d$Test_Statistic, d$df, d$p_value, d$Method)
sres <- function(label, s) row(label, s$Statistic, s$df, s$p_value,
  paste0(if (!is.null(s$Note)) s$Note else "", if (!is.null(s$rho) && isTRUE(is.finite(s$rho))) sprintf(" rho=%.4f", s$rho) else ""))

all_def <- function(tag, fit, G = 10) {
  cur <- pk$def.gof(fit, G = G)
  dres(paste(tag, "CURRENT def.gof poly3"), cur)
  pu <- proto_def(fit, G = G, basis = "poly3")
  cat(sprintf("  %-62s %s\n", paste(tag, "proto unit poly3 identical to current?"),
              isTRUE(all.equal(c(pu$Test_Statistic, pu$p_value), c(cur$Test_Statistic, cur$p_value), tolerance = 1e-12))))
  for (b in c("poly2", "stukel")) dres(paste(tag, "CURRENT def.gof", b), pk$def.gof(fit, G = G, basis = b))
  for (b in c("poly3", "poly2", "stukel", "sym", "stukel2")) for (wt in c("unit", "score")) {
    if (b %in% c("poly3", "poly2", "stukel") && wt == "unit") next
    if (b == "stukel2" && wt == "unit") next
    r <- tryCatch(proto_def(fit, G = G, basis = b, weights = wt, diag = TRUE), error = function(e) e)
    if (inherits(r, "error")) { cat("  ", tag, b, wt, "ERROR:", conditionMessage(r), "\n"); next }
    dres(paste(tag, "PROTO", b, wt), r)
    ev <- attr(r, "eig_rel")
    if (!is.null(ev)) cat(sprintf("  %-62s %s\n", "      relative eigenvalues of I",
                                  paste(formatC(ev, digits = 3, format = "g"), collapse = " ")))
  }
  e3 <- tryCatch(pk$def.ensemble.gof(fit, G = G)$p_value, error = function(e) NA)
  row(paste(tag, "CURRENT def.ensemble.gof (poly2+poly3+stukel)"), NA, NA, e3)
}

stk <- function(tag, fit) {
  ctx <- pk$.gof_context(fit)
  cur <- withCallingHandlers(
    tryCatch(pk$gof_stukel(ctx), error = function(e) list(Statistic = NA, df = NA, p_value = NA,
                                                          Note = paste("ERROR", conditionMessage(e)))),
    warning = function(w) { cat("      warning:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  sres(paste(tag, "CURRENT gof_stukel (statmod branch)"), cur)
  cat(sprintf("  %-62s Statistic is.nan=%s p is.nan=%s\n", "", is.nan(cur$Statistic), is.nan(cur$p_value)))
  ph <- ctx$ph; eta <- as.numeric(stats::predict(fit, type = "link"))
  za <- 0.5 * eta^2 * (ph >= 0.5); zb <- -0.5 * eta^2 * (ph < 0.5)
  zz <- c(pk$.gof_score_z(za, ctx$y, ph, ctx$X), pk$.gof_score_z(zb, ctx$y, ph, ctx$X))
  row(paste(tag, "CURRENT no-statmod branch (.gof_score_z)"), sum(zz^2), 2,
      stats::pchisq(sum(zz^2), 2, lower.tail = FALSE), sprintf("z=(%s, %s)", fmt(zz[1]), fmt(zz[2])))
  zs <- tryCatch(statmod::glm.scoretest(fit, cbind(za, zb)), error = function(e) c(NA, NA))
  cat(sprintf("  %-62s statmod z=(%s, %s)  n(p>=.5)=%d of %d\n", "", fmt(zs[1]), fmt(zs[2]), sum(ph >= 0.5), length(ph)))
  sres(paste(tag, "PROTO joint (with one-df fallback)"), proto_stukel(fit, "joint"))
}

rao_joint <- function(fit) {
  y <- fit$y; Xm <- stats::model.matrix(fit); eta <- stats::predict(fit, type = "link")
  ph <- pmin(pmax(stats::fitted(fit), 1e-6), 1 - 1e-6)
  za <- 0.5 * eta^2 * (ph >= 0.5); zb <- -0.5 * eta^2 * (ph < 0.5)
  f0 <- stats::glm(y ~ Xm - 1, family = stats::binomial())
  f1 <- suppressWarnings(stats::glm(y ~ Xm + za + zb - 1, family = stats::binomial()))
  a <- stats::anova(f0, f1, test = "Rao"); c(Rao = a$Rao[2], df = a$Df[2], p = a[["Pr(>Chi)"]][2])
}
rao_step <- function(fit, basis, G = 10) {
  y <- fit$y; Xm <- stats::model.matrix(fit); n <- length(y)
  ph <- pmin(pmax(stats::fitted(fit), 1e-6), 1 - 1e-6)
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  pbar <- as.numeric(tapply(ph, grp, mean))
  S <- proto_basis(pbar, basis)[grp, , drop = FALSE]
  f0 <- stats::glm(y ~ Xm - 1, family = stats::binomial())
  f1 <- suppressWarnings(stats::glm(y ~ Xm + S - 1, family = stats::binomial()))
  a <- stats::anova(f0, f1, test = "Rao"); c(Rao = a$Rao[2], df = a$Df[2], p = a[["Pr(>Chi)"]][2])
}
harness_probes <- function(fit, G = 10) {   # the run_L3_closing.R formulas, clamp set to the package's 1e-6
  y <- fit$y; X <- stats::model.matrix(fit); n <- length(y)
  ph <- pmin(pmax(stats::fitted(fit), 1e-6), 1 - 1e-6); w <- ph * (1 - ph)
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  o <- as.numeric(tapply(y, g, sum)); e <- as.numeric(tapply(ph, g, sum)); V <- as.numeric(tapply(w, g, sum))
  pb <- as.numeric(tapply(ph, g, mean)); eb <- stats::qlogis(pb); r <- (o - e) / sqrt(V)
  U <- rowsum(w * X, g) / sqrt(V); XWX <- crossprod(X, w * X)
  ZOZ <- function(Z) crossprod(Z) - crossprod(Z, U) %*% solve(XWX, crossprod(U, Z))
  Zg <- cbind(sym = eb * abs(eb)); zg <- drop(crossprod(Zg, r)) / sqrt(diag(ZOZ(Zg)))
  Zt <- Zg * sqrt(V); zt <- drop(crossprod(Zt, r)) / sqrt(diag(ZOZ(Zt)))
  Zk <- cbind(eb^2 * (eb >= 0), -eb^2 * (eb < 0)) * sqrt(V)
  sc <- if (all(apply(Zk, 2, stats::sd) > 1e-10)) { sk <- drop(crossprod(Zk, r))
    stats::pchisq(drop(t(sk) %*% solve(ZOZ(Zk), sk)), 2, lower.tail = FALSE) } else NA
  c(g.sym = 2 * stats::pnorm(-abs(zg)), g.sym.sc = 2 * stats::pnorm(-abs(zt)), EDGE.stk.sc = sc)
}

cat("=== 0. environment ===\n")
cat(R.version.string, "\n")
cat("installed ebrahim.gof:", format(utils::packageVersion("ebrahim.gof")), "\n")
db <- tryCatch(tools::CRAN_package_db(), error = function(e) NULL)
if (!is.null(db)) print(db[db$Package == "ebrahim.gof", c("Package", "Version", "Published")], row.names = FALSE)
if (requireNamespace("LogisticDx", quietly = TRUE)) {
  src <- deparse(utils::getFromNamespace("gof.glm", "LogisticDx"))
  i <- grep("stuk|Sst|scoretest", src, ignore.case = TRUE)
  idx <- sort(unique(unlist(lapply(i, function(k) max(1, k - 3):min(length(src), k + 3)))))
  cat("--- LogisticDx:::gof.glm, Stukel lines ---\n"); cat(paste(idx, src[idx]), sep = "\n")
}

cat("\n=== 1. tests/testthat/test-def-gof.R fixtures ===\n")
mf_def <- function(seed = 1, n = 600, link = "logit") {
  set.seed(seed); x <- runif(n, -3, 3); eta <- 0.6 * x
  p <- if (link == "cloglog") 1 - exp(-exp(eta)) else 1 / (1 + exp(-eta))
  glm(rbinom(n, 1, p) ~ x, family = binomial())
}
all_def("[def make_fit(1,600)]", mf_def())
f2 <- mf_def(2)
a <- pk$def.gof(f2, basis = "poly2")
b <- pk$def.gof(as.numeric(f2$y), stats::fitted(f2), X = stats::model.matrix(f2), basis = "poly2")
row("[def make_fit(2)] CURRENT poly2 glm form", a$Test_Statistic, a$df, a$p_value)
row("[def make_fit(2)] CURRENT poly2 (y, ph, X) form", b$Test_Statistic, b$df, b$p_value)
for (wt in c("score")) for (bs in c("poly2", "sym")) {
  a2 <- proto_def(f2, basis = bs, weights = wt)
  b2 <- proto_def(as.numeric(f2$y), stats::fitted(f2), X = stats::model.matrix(f2), basis = bs, weights = wt)
  row(sprintf("[def make_fit(2)] PROTO %s %s glm vs (y,ph,X) |dp|", bs, wt), NA, NA, abs(a2$p_value - b2$p_value))
}
fc <- mf_def(7, 1500, "cloglog")
all_def("[def cloglog(7,1500)]", fc)

cat("\n=== 2. tests/testthat/test-def-ensemble-gof.R fixtures ===\n")
mf_ens <- function(seed = 1, n = 600) { set.seed(seed); x <- runif(n, -3, 3)
  glm(rbinom(n, 1, 1 / (1 + exp(-(0.6 * x)))) ~ x, family = binomial()) }
e1 <- mf_ens(); row("[ens make_fit(1)] CURRENT cct", NA, NA, pk$def.ensemble.gof(e1)$p_value)
row("[ens make_fit(1)] CURRENT cct + EF", NA, NA, pk$def.ensemble.gof(e1, add_ef = TRUE)$p_value)
row("[ens make_fit(3)] CURRENT cct", NA, NA, pk$def.ensemble.gof(mf_ens(3))$p_value)

cat("\n=== 3. tests/testthat/test-run-all-gof.R fixtures (Stukel) ===\n")
mf_ra <- function(seed = 1, n = 500) { set.seed(seed); x <- runif(n, -3, 3)
  glm(rbinom(n, 1, 1 / (1 + exp(-(0.6 * x)))) ~ x, family = binomial()) }
fr <- mf_ra()
stk("[ra make_fit(1,500)]", fr)
print(rao_joint(fr))
set.seed(7); xc <- runif(1500, -3, 3)
fcl <- glm(rbinom(1500, 1, 1 - exp(-exp(0.6 * xc))) ~ xc, family = binomial())
stk("[ra cloglog seed7 n1500]", fcl)
bat <- pk$run.all.gof(fr, include_slow = FALSE, install = "no")
cat(sprintf("  run.all.gof(make_fit(1,500), include_slow = FALSE): %d rows\n  %s\n", nrow(bat), paste(bat$Test, collapse = ", ")))
print(as.data.frame(bat)[bat$Test %in% c("Stukel", "DEF.poly2", "DEF.poly3", "DEF.stukel", "Ensemble.Vote(3DEF)", "Ensemble.Univ(3DEF+EF)"), ], row.names = FALSE)

cat("\n=== 4. gof_demo (examples of def.gof / def.ensemble.gof; test-exports-smoke.R) ===\n")
wrong <- glm(outcome ~ age + bmi + sex + treatment, data = pk$gof_demo, family = binomial())
right <- glm(outcome ~ poly(age, 2) + bmi + sex + treatment, data = pk$gof_demo, family = binomial())
all_def("[gof_demo wrong]", wrong); stk("[gof_demo wrong]", wrong)
all_def("[gof_demo right]", right); stk("[gof_demo right]", right)
row("[gof_demo wrong] CURRENT def.ensemble.gof add_ef", NA, NA, pk$def.ensemble.gof(wrong, add_ef = TRUE)$p_value)

cat("\n=== 5. roxygen examples (edge.gof, edges.gof, run.all.gof, package page) ===\n")
set.seed(1); x <- runif(500, -3, 3); y <- rbinom(500, 1, plogis(0.6 * x)); fex <- glm(y ~ x, family = binomial())
all_def("[edge.gof example]", fex); stk("[edge.gof / run.all.gof example fit]", fex)
set.seed(1); n <- 500; x <- runif(n, -3, 3); y <- rbinom(n, 1, 1 / (1 + exp(-(0.6 * x))))
fit <- glm(y ~ x, family = binomial())
y2 <- rbinom(n, 1, 1 / (1 + exp(-(0.6 * x + 0.5 * x^2)))); bad <- glm(y2 ~ x, family = binomial())
cat(sprintf("  run.all.gof example fit identical to edge.gof example fit? %s\n", isTRUE(all.equal(fitted(fit), fitted(fex)))))
stk("[run.all.gof example bad (omitted x^2)]", bad)
set.seed(1); n <- 200; x <- rnorm(n); y <- rbinom(n, 1, plogis(0.3 + 0.9 * x)); fpk <- glm(y ~ x, family = binomial())
all_def("[package-page example n200]", fpk)

cat("\n=== 6. degenerate Stukel columns ===\n")
set.seed(5); x <- runif(400, -3, 3)
ylo <- rbinom(400, 1, plogis(-4 + 0.5 * x)); flo <- glm(ylo ~ x, family = binomial())
cat(sprintf("  all-low fit: fitted range %.4f-%.4f, events %d\n", min(fitted(flo)), max(fitted(flo)), sum(ylo)))
stk("[all fitted < 0.5]", flo)
dres("[all fitted < 0.5] CURRENT def.gof stukel (zero column dropped)", pk$def.gof(flo, basis = "stukel"))
print(as.data.frame(pk$run.all.gof(flo, tests = c("Stukel", "DEF.stukel"), install = "no")), row.names = FALSE)
yhi <- rbinom(400, 1, plogis(4 + 0.5 * x)); fhi <- glm(yhi ~ x, family = binomial())
cat(sprintf("  all-high fit: fitted range %.4f-%.4f, non-events %d\n", min(fitted(fhi)), max(fitted(fhi)), 400 - sum(yhi)))
stk("[all fitted >= 0.5]", fhi)
ynr <- rbinom(400, 1, plogis(-1.6 + 0.5 * x)); fnr <- glm(ynr ~ x, family = binomial())
stk("[few fitted >= 0.5]", fnr)

cat("\n=== 7. independent golden references ===\n")
for (nm in c("fr", "wrong", "fex")) {
  f <- get(nm); pj <- proto_stukel(f)
  rj <- rao_joint(f)
  cat(sprintf("  %-6s joint proto T=%.8f p=%.8f | anova Rao T=%.8f df=%d p=%.8f\n", nm, pj$Statistic, pj$p_value, rj["Rao"], as.integer(rj["df"]), rj["p"]))
  for (bs in c("sym", "poly3", "stukel2", "stukel")) {
    ps <- proto_def(f, basis = bs, weights = "score"); rs <- rao_step(f, bs)
    cat(sprintf("  %-6s %-7s score proto T=%.8f df=%d p=%.8f | anova Rao(step covariate) T=%.8f df=%d p=%.8f\n",
                nm, bs, ps$Test_Statistic, as.integer(ps$df), ps$p_value, rs["Rao"], as.integer(rs["df"]), rs["p"]))
  }
  h <- harness_probes(f)
  cat(sprintf("  %-6s harness g.sym=%.8f proto sym unit=%.8f | harness g.sym.sc=%.8f proto sym score=%.8f | harness EDGE.stk.sc=%.8f proto stukel2 score=%.8f\n",
              nm, h["g.sym"], proto_def(f, basis = "sym")$p_value, h["g.sym.sc"], proto_def(f, basis = "sym", weights = "score")$p_value,
              h["EDGE.stk.sc"], proto_def(f, basis = "stukel2", weights = "score")$p_value))
}
if (requireNamespace("LogisticDx", quietly = TRUE)) {
  cap <- utils::capture.output(gg <- tryCatch(LogisticDx::gof(fex, g = 10, plotROC = FALSE), error = function(e) e))
  if (inherits(gg, "error")) cat("  LogisticDx::gof error:", conditionMessage(gg), "\n") else {
    tb <- as.data.frame(gg$gof); cat("  LogisticDx::gof(fex) Stukel rows:\n")
    print(tb[apply(tb, 1, function(r) any(grepl("Sst", r))), ], row.names = FALSE)
  }
}
cat("\ndone\n")
