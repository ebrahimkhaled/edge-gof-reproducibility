## map_pkg_values2.R -- follow-up to map_pkg_values.R: harness identity with names fixed, LogisticDx SstBoth identity,
## a near-degenerate Stukel sample, standalone fixtures for the proposed new tests, and timing. Read-only, single core.
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
options(width = 200, warn = 1)

harness_probes <- function(fit, G = 10) {
  y <- fit$y; X <- stats::model.matrix(fit); n <- length(y)
  ph <- pmin(pmax(stats::fitted(fit), 1e-6), 1 - 1e-6); w <- ph * (1 - ph)
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  o <- as.numeric(tapply(y, g, sum)); e <- as.numeric(tapply(ph, g, sum)); V <- as.numeric(tapply(w, g, sum))
  pb <- as.numeric(tapply(ph, g, mean)); eb <- stats::qlogis(pb); r <- (o - e) / sqrt(V)
  U <- rowsum(w * X, g) / sqrt(V); XWX <- crossprod(X, w * X)
  ZOZ <- function(Z) crossprod(Z) - crossprod(Z, U) %*% solve(XWX, crossprod(U, Z))
  Zg <- cbind(eb * abs(eb)); zg <- unname(drop(crossprod(Zg, r)) / sqrt(diag(ZOZ(Zg))))
  Zt <- Zg * sqrt(V); zt <- unname(drop(crossprod(Zt, r)) / sqrt(diag(ZOZ(Zt))))
  Zs <- cbind(eb, eb^2 * (eb >= 0), -eb^2 * (eb < 0)); Zs <- Zs[, apply(Zs, 2, stats::sd) > 1e-10, drop = FALSE]
  lam <- Re(eigen(solve(crossprod(Zs), ZOZ(Zs)), only.values = TRUE)$values); lam <- lam[lam > 1e-9]
  Sx <- sum(qr.fitted(qr(Zs), r)^2); cc <- sum(lam^2) / sum(lam); nu <- sum(lam)^2 / sum(lam^2)
  c(g.sym = 2 * stats::pnorm(-abs(zg)), g.sym.sc = 2 * stats::pnorm(-abs(zt)),
    EDGE.stk = stats::pchisq(Sx / cc, nu, lower.tail = FALSE))
}

cat("=== 1. harness formulas (run_L3_closing.R, clamp 1e-6) vs prototypes ===\n")
wrong <- glm(outcome ~ age + bmi + sex + treatment, data = pk$gof_demo, family = binomial())
set.seed(1); dd <- data.frame(x = runif(500, -3, 3)); dd$y <- rbinom(500, 1, plogis(0.6 * dd$x))
fex <- glm(y ~ x, family = binomial(), data = dd)
for (nm in c("wrong", "fex")) {
  f <- get(nm); h <- harness_probes(f)
  cat(sprintf("  %-5s g.sym %.10f vs proto sym unit %.10f | g.sym.sc %.10f vs proto sym score %.10f | EDGE.stk %.10f vs def.gof stukel %.10f\n",
              nm, h["g.sym"], proto_def(f, basis = "sym")$p_value, h["g.sym.sc"],
              proto_def(f, basis = "sym", weights = "score")$p_value, h["EDGE.stk"], pk$def.gof(f, basis = "stukel")$p_value))
  cat(sprintf("        sym unit, imhof method p = %.10f\n", proto_def(f, basis = "sym", method = "imhof")$p_value))
}

cat("\n=== 2. LogisticDx::gof SstBoth vs the package's marginal sum ===\n")
cap <- utils::capture.output(gg <- tryCatch(LogisticDx::gof(fex, g = 10, plotROC = FALSE), error = function(e) e))
if (inherits(gg, "error")) cat("  LogisticDx::gof error:", conditionMessage(gg), "\n") else {
  tb <- as.data.frame(gg$gof); print(tb[grepl("Sst", tb$test), ], row.names = FALSE)
}
cur <- pk$gof_stukel(pk$.gof_context(fex))
cat(sprintf("  package gof_stukel (2.7.0): Statistic %.8f p %.8f\n", cur$Statistic, cur$p_value))

cat("\n=== 3. near-degenerate: a sample with 1-3 fitted risks >= 0.5 ===\n")
found <- FALSE
for (sd0 in 1:400) {
  set.seed(sd0); x <- runif(400, -3, 3); y <- rbinom(400, 1, plogis(-2.2 + 0.6 * x))
  f <- suppressWarnings(glm(y ~ x, family = binomial())); k <- sum(fitted(f) >= 0.5)
  if (k >= 1 && k <= 3) { found <- TRUE; break }
}
if (found) {
  cat(sprintf("  seed %d: %d fitted >= 0.5, events %d\n", sd0, k, sum(y)))
  cm <- pk$gof_stukel(pk$.gof_context(f)); pj <- proto_stukel(f)
  cat(sprintf("  current marginal: Statistic %.6f p %.6f | joint: T %.6f df %d p %.6f rho %.4f\n",
              cm$Statistic, cm$p_value, pj$Statistic, pj$df, pj$p_value, pj$rho))
}

cat("\n=== 4. standalone fixtures for the proposed tests ===\n")
set.seed(5); x <- runif(400, -3, 3); y <- rbinom(400, 1, plogis(-4 + 0.5 * x)); flo <- glm(y ~ x, family = binomial())
pj <- proto_stukel(flo); cm <- pk$gof_stukel(pk$.gof_context(flo))
cat(sprintf("  LOW  set.seed(5); x <- runif(400,-3,3); y <- rbinom(400,1,plogis(-4+0.5*x)): joint T %.7f df %d p %.7f [%s] | marginal Statistic %s p %s\n",
            pj$Statistic, pj$df, pj$p_value, pj$Note, format(cm$Statistic), format(cm$p_value)))
set.seed(6); x <- runif(400, -3, 3); y <- rbinom(400, 1, plogis(4 + 0.5 * x)); fhi <- glm(y ~ x, family = binomial())
pj <- proto_stukel(fhi); cm <- pk$gof_stukel(pk$.gof_context(fhi))
cat(sprintf("  HIGH set.seed(6); x <- runif(400,-3,3); y <- rbinom(400,1,plogis(4+0.5*x)): joint T %.7f df %d p %.7f [%s] | marginal Statistic %s p %s | range %.3f-%.3f\n",
            pj$Statistic, pj$df, pj$p_value, pj$Note, format(cm$Statistic), format(cm$p_value), min(fitted(fhi)), max(fitted(fhi))))
rightm <- glm(outcome ~ poly(age, 2) + bmi + sex + treatment, data = pk$gof_demo, family = binomial())
for (bs in c("sym", "poly3", "poly2", "stukel")) for (wt in c("unit", "score")) {
  r <- proto_def(wrong, basis = bs, weights = wt); r2 <- proto_def(rightm, basis = bs, weights = wt)
  cat(sprintf("  gof_demo %-6s %-5s wrong: S %.7f df %.6f p %.7f | right: S %.7f df %.6f p %.7f\n",
              bs, wt, r$Test_Statistic, r$df, r$p_value, r2$Test_Statistic, r2$df, r2$p_value))
}
nai <- suppressWarnings(proto_def(as.numeric(wrong$y), fitted(wrong), basis = "sym", weights = "score"))
cat(sprintf("  gof_demo wrong, (y, ph) with no X, sym score (Omega = I): T %.7f p %.7f\n", nai$Test_Statistic, nai$p_value))

cat("\n=== 5. timing ===\n")
set.seed(9); n <- 20000; x <- runif(n, -3, 3); d <- rbinom(n, 1, .5); y <- rbinom(n, 1, plogis(0.6 * x + 0.5 * d))
fb <- glm(y ~ x + d, family = binomial())
cat(sprintf("  n = 20000: joint Stukel %.3fs; marginal (statmod) %.3fs\n",
            system.time(proto_stukel(fb))[["elapsed"]], system.time(pk$gof_stukel(pk$.gof_context(fb)))[["elapsed"]]))
for (G in c(10, 800)) cat(sprintf("  n = 20000, G = %d: def.gof poly3 unit %.3fs; sym score %.3fs; poly3 score %.3fs\n", G,
    system.time(pk$def.gof(fb, G = G))[["elapsed"]], system.time(proto_def(fb, G = G, basis = "sym", weights = "score"))[["elapsed"]],
    system.time(proto_def(fb, G = G, basis = "poly3", weights = "score"))[["elapsed"]]))
cat("\ndone\n")
