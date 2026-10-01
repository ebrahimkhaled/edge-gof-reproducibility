## pkg280_verify_num_size.R -- null size of the 2.8.0 Stukel forms and of DEF sym / poly3 / stukel (unit, score),
## using the INSTALLED ebrahim.gof. Verification of the package, not a paper result.
##
## Design (the run-L / MAP null-check design): x ~ U(-3,3), d ~ Bern(0.5), eta = c0 + s(0.6x + 0.5d), fit y ~ x + d.
## Cells: base n = 1000 with G = 10 and with G = "auto"; c0 = -2, n = 500, G = "auto"; s = 2, n = 1000, G = "auto".
## Each cell has its own independent draws. B = 2000 per cell, 6 workers.
##
## Stukel is called through run.all.gof(tests = "Stukel", control = list(Stukel = list(form = ...))); every replicate
## also records anova(fit, augmented, test = "Rao") and the deviance drop, so the joint and lr forms are checked
## against base R on all 8000 data sets. DEF is called through def.gof().
##
## Run: Rscript pkg280_verify_num_size.R > ../paper_EDGE/theory/pkg280_verify_num_size.log 2>&1
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

suppressPackageStartupMessages({ library(ebrahim.gof); library(parallel) })
OUT <- edge_path("declarations")
PKGV <- as.character(packageVersion("ebrahim.gof"))
stopifnot(PKGV == "2.8.0")
cat("ebrahim.gof", PKGV, "(installed) |", R.version.string, "|", format(Sys.time()), "\n")

VARS <- c("G", "ev_min", "few_calls",
          "stk.joint", "stk.joint.stat", "stk.joint.df", "rao.stat", "rao.df",
          "stk.lr", "stk.lr.stat", "stk.lr.df", "lrt.stat", "lrt.df",
          "stk.marg", "stk.marg.stat",
          "sym.unit", "sym.score", "sym.score.df",
          "poly3.unit", "poly3.score", "poly3.score.df",
          "stukel.unit", "stukel.score", "stukel.score.df",
          "tight.joint.stat", "tight.rao.stat")

one_rep <- function(i, s, c0, n, G) {
  out <- setNames(rep(NA_real_, length(VARS)), VARS)
  xa <- runif(n, -3, 3); db <- rbinom(n, 1, 0.5)
  resp <- rbinom(n, 1, plogis(c0 + s * (0.6 * xa + 0.5 * db)))
  Gn <- if (identical(G, "auto")) max(10, round(n / 25)) else G
  out["G"] <- Gn
  out["ev_min"] <- min(sum(resp), n - sum(resp))
  fit <- tryCatch(suppressWarnings(glm(resp ~ xa + db, family = binomial())), error = function(e) NULL)
  if (is.null(fit)) return(out)

  ## Stukel, three forms, through the battery
  for (fm in c("joint", "lr", "marginal")) {
    r <- tryCatch(suppressWarnings(run.all.gof(fit, tests = "Stukel", install = "no",
                                               control = list(Stukel = list(form = fm)))),
                  error = function(e) NULL)
    if (is.null(r)) next
    key <- c(joint = "stk.joint", lr = "stk.lr", marginal = "stk.marg")[[fm]]
    out[key] <- r$p_value
    out[paste0(key, ".stat")] <- r$Statistic
    if (fm != "marginal") out[paste0(key, ".df")] <- r$df
  }

  ## base-R reference: augmented fit, anova Rao and the deviance drop
  eta <- predict(fit, type = "link"); ph <- fitted(fit)
  Zs <- cbind(za = 0.5 * eta^2 * (ph >= 0.5), zb = -0.5 * eta^2 * (ph < 0.5))
  Xm <- model.matrix(fit)
  f1 <- tryCatch(suppressWarnings(glm(resp ~ Xm + Zs - 1, family = binomial())), error = function(e) NULL)
  if (!is.null(f1)) {
    a <- anova(fit, f1, test = "Rao")
    out["rao.stat"] <- a$Rao[2]; out["rao.df"] <- a$Df[2]
    out["lrt.stat"] <- fit$deviance - f1$deviance; out["lrt.df"] <- f1$rank - fit$rank
  }
  ## same comparison on a tightly converged refit: anova's Rao regresses the working residuals with the weights
  ## glm.fit returns, which are one IRLS iteration old, so at the default epsilon it drifts from u'I^-1 u
  tc <- glm.control(epsilon = 1e-13, maxit = 100)
  fit_t <- tryCatch(suppressWarnings(glm(resp ~ xa + db, family = binomial(), control = tc)), error = function(e) NULL)
  if (!is.null(fit_t)) {
    rt <- tryCatch(run.all.gof(fit_t, tests = "Stukel", install = "no"), error = function(e) NULL)
    if (!is.null(rt)) out["tight.joint.stat"] <- rt$Statistic
    eta_t <- predict(fit_t, type = "link"); ph_t <- fitted(fit_t)
    Zt <- cbind(za = 0.5 * eta_t^2 * (ph_t >= 0.5), zb = -0.5 * eta_t^2 * (ph_t < 0.5))
    f1t <- tryCatch(suppressWarnings(glm(resp ~ Xm + Zt - 1, family = binomial(), control = tc)),
                    error = function(e) NULL)
    if (!is.null(f1t)) out["tight.rao.stat"] <- anova(fit_t, f1t, test = "Rao")$Rao[2]
  }

  ## DEF sym / poly3 / stukel, unit and score; count the few-events warnings
  few <- 0L
  for (b in c("sym", "poly3", "stukel")) for (wt in c("unit", "score")) {
    r <- tryCatch(withCallingHandlers(def.gof(fit, G = G, basis = b, weights = wt),
                                      def_few_events = function(w) { few <<- few + 1L; invokeRestart("muffleWarning") }),
                  error = function(e) NULL)
    if (is.null(r)) next
    out[paste0(b, ".", wt)] <- r$p_value
    if (wt == "score") out[paste0(b, ".score.df")] <- r$df
  }
  out["few_calls"] <- few
  out
}

CELLS <- data.frame(label = c("base n1000 G10", "base n1000 Gauto", "c0=-2 n500 Gauto", "s=2 n1000 Gauto"),
                    s = c(1, 1, 1, 2), c0 = c(0, 0, -2, 0), n = c(1000, 1000, 500, 1000),
                    G = c("10", "auto", "auto", "auto"), stringsAsFactors = FALSE)
B <- 2000L
cl <- makeCluster(6L)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(ebrahim.gof))))
clusterExport(cl, c("VARS", "one_rep"))

TESTS <- c("stk.joint", "stk.lr", "stk.marg", "sym.unit", "sym.score", "poly3.unit", "poly3.score",
           "stukel.unit", "stukel.score")
t0 <- Sys.time(); ALL <- list(); SUM <- list()
for (k in seq_len(nrow(CELLS))) {
  ce <- CELLS[k, ]
  Garg <- if (ce$G == "auto") "auto" else as.numeric(ce$G)
  clusterSetRNGStream(cl, 20260914L + 100L * k)
  M <- do.call(rbind, parLapply(cl, seq_len(B), one_rep, s = ce$s, c0 = ce$c0, n = ce$n, G = Garg))
  ALL[[k]] <- data.frame(cell = ce$label, rep = seq_len(B), M, pkg = PKGV, check.names = FALSE)
  cat(sprintf("\n== %s  (s = %g, c0 = %g, n = %d, G = %s -> %d; B = %d; %.1f min elapsed) ==\n",
              ce$label, ce$s, ce$c0, ce$n, ce$G, as.integer(M[1, "G"]), B,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  for (v in TESTS) {
    p <- M[, v]; ok <- is.finite(p); m <- sum(ok)
    s05 <- mean(p[ok] <= 0.05); se05 <- sqrt(0.05 * 0.95 / m); z05 <- (s05 - 0.05) / se05
    s01 <- mean(p[ok] <= 0.01); se01 <- sqrt(0.01 * 0.99 / m); z01 <- (s01 - 0.01) / se01
    SUM[[length(SUM) + 1]] <- data.frame(cell = ce$label, test = v, m = m, no_p = B - m,
                                         size05 = s05, se05 = se05, z05 = z05, flag05 = abs(z05) > 3,
                                         size01 = s01, se01 = se01, z01 = z01, flag01 = abs(z01) > 3, pkg = PKGV)
    cat(sprintf("  %-13s m %4d  size05 %.4f (MCSE %.4f, z %+5.2f)%s  size01 %.4f (MCSE %.4f, z %+5.2f)%s\n",
                v, m, s05, se05, z05, if (abs(z05) > 3) " FLAG" else "     ",
                s01, se01, z01, if (abs(z01) > 3) " FLAG" else ""))
  }
  okj <- is.finite(M[, "stk.joint.stat"]) & is.finite(M[, "rao.stat"])
  okl <- is.finite(M[, "stk.lr.stat"]) & is.finite(M[, "lrt.stat"])
  cat(sprintf("  joint vs anova Rao: n %d, max |diff| %.2e, max rel %.2e, df mismatches %d; 1-df fallback %d\n",
              sum(okj), max(abs(M[okj, "stk.joint.stat"] - M[okj, "rao.stat"])),
              max(abs(M[okj, "stk.joint.stat"] - M[okj, "rao.stat"]) / pmax(M[okj, "rao.stat"], 1e-8)),
              sum(M[okj, "stk.joint.df"] != M[okj, "rao.df"]), sum(M[, "stk.joint.df"] %in% 1)))
  okt <- is.finite(M[, "tight.joint.stat"]) & is.finite(M[, "tight.rao.stat"])
  cat(sprintf("  joint vs anova Rao, epsilon 1e-13 refit: n %d, max |diff| %.2e, max rel %.2e\n",
              sum(okt), max(abs(M[okt, "tight.joint.stat"] - M[okt, "tight.rao.stat"])),
              max(abs(M[okt, "tight.joint.stat"] - M[okt, "tight.rao.stat"]) / pmax(M[okt, "tight.rao.stat"], 1e-8))))
  cat(sprintf("  lr vs deviance drop: n %d, max |diff| %.2e, df mismatches %d; 1-df %d\n",
              sum(okl), max(abs(M[okl, "stk.lr.stat"] - M[okl, "lrt.stat"])),
              sum(M[okl, "stk.lr.df"] != M[okl, "lrt.df"]), sum(M[, "stk.lr.df"] %in% 1)))
  cat(sprintf("  score df tables: sym %s | poly3 %s | stukel %s\n",
              paste(names(table(M[, "sym.score.df"])), table(M[, "sym.score.df"]), sep = ":", collapse = " "),
              paste(names(table(M[, "poly3.score.df"])), table(M[, "poly3.score.df"]), sep = ":", collapse = " "),
              paste(names(table(M[, "stukel.score.df"])), table(M[, "stukel.score.df"]), sep = ":", collapse = " ")))
  expect_few <- M[, "ev_min"] < M[, "G"]
  cat(sprintf("  few-events: reps with min(events, non-events) < G: %d; reps warned: %d; disagreements: %d; calls per warned rep: %s\n",
              sum(expect_few), sum(M[, "few_calls"] > 0), sum((M[, "few_calls"] > 0) != expect_few),
              paste(unique(M[M[, "few_calls"] > 0, "few_calls"]), collapse = ",")))
  cat(sprintf("  min(events, non-events): min %d, 1%% %g, median %g\n", as.integer(min(M[, "ev_min"])),
              quantile(M[, "ev_min"], 0.01), median(M[, "ev_min"])))
}
stopCluster(cl)
write.csv(do.call(rbind, ALL), file.path(OUT, "pkg280_verify_num_size_pvalues.csv"), row.names = FALSE)
write.csv(do.call(rbind, SUM), file.path(OUT, "pkg280_verify_num_size_summary.csv"), row.names = FALSE)
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 2)), "\n")
