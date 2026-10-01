## pkg280_verify_fix_size.R -- null size after the review fixes, fixed branch loaded with pkgload::load_all() in each
## worker. Base design x ~ U(-3,3), d ~ Bern(0.5), eta = 0.6x + 0.5d, fit y ~ x + d, n = 1000, G = "auto" (40).
## Tests, all through run.all.gof(G = "auto"): Stukel joint; DEF.sym and DEF.stukel, unit (default) and score (control).
## Each replicate also checks the battery rows against def.gof() called directly. B = 1000, 6 workers.
## Run: Rscript pkg280_verify_fix_size.R [package dir] > ../paper_EDGE/theory/pkg280_verify_fix_size.log 2>&1
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
OUT <- edge_path("declarations")
suppressPackageStartupMessages(library(parallel))
sha <- system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE)
cat("fixed branch HEAD", sha, "via load_all in each worker |", R.version.string, "|", format(Sys.time()), "\n")

one_rep <- function(i, n = 1000) {
  V <- c("ev_min", "stk.joint", "stk.joint.stat", "stk.joint.df", "sym.unit", "sym.score", "sym.score.df",
         "stukel.unit", "stukel.score", "stukel.score.df", "notes_auto_ok", "max_abs_vs_defgof")
  out <- setNames(rep(NA_real_, length(V)), V)
  xa <- runif(n, -3, 3); db <- rbinom(n, 1, 0.5); resp <- rbinom(n, 1, plogis(0.6 * xa + 0.5 * db))
  out["ev_min"] <- min(sum(resp), n - sum(resp))
  fit <- suppressWarnings(glm(resp ~ xa + db, family = binomial()))
  u <- as.data.frame(run.all.gof(fit, G = "auto", tests = c("Stukel", "DEF.sym", "DEF.stukel"), install = "no"))
  s <- as.data.frame(run.all.gof(fit, G = "auto", tests = c("DEF.sym", "DEF.stukel"), install = "no",
                                 control = list(DEF.sym = list(weights = "score"), DEF.stukel = list(weights = "score"))))
  g <- function(b, t, col) b[[col]][b$Test == t]
  out["stk.joint"] <- g(u, "Stukel", "p_value"); out["stk.joint.stat"] <- g(u, "Stukel", "Statistic")
  out["stk.joint.df"] <- g(u, "Stukel", "df")
  out["sym.unit"] <- g(u, "DEF.sym", "p_value"); out["stukel.unit"] <- g(u, "DEF.stukel", "p_value")
  out["sym.score"] <- g(s, "DEF.sym", "p_value"); out["sym.score.df"] <- g(s, "DEF.sym", "df")
  out["stukel.score"] <- g(s, "DEF.stukel", "p_value"); out["stukel.score.df"] <- g(s, "DEF.stukel", "df")
  out["notes_auto_ok"] <- all(grepl("G = 40 \\(auto\\)", c(g(u, "DEF.sym", "Note"), g(u, "DEF.stukel", "Note"),
                                                          g(s, "DEF.sym", "Note"), g(s, "DEF.stukel", "Note"))))
  direct <- c(def.gof(fit, G = 40, basis = "sym")$p_value, def.gof(fit, G = 40, basis = "stukel")$p_value,
              def.gof(fit, G = 40, basis = "sym", weights = "score")$p_value,
              def.gof(fit, G = 40, basis = "stukel", weights = "score")$p_value)
  out["max_abs_vs_defgof"] <- max(abs(direct - out[c("sym.unit", "stukel.unit", "sym.score", "stukel.score")]))
  out
}

B <- 1000L
cl <- makeCluster(6L)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl, "PKG")
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(pkgload::load_all(PKG, export_all = FALSE, quiet = TRUE))))
clusterSetRNGStream(cl, 20260915L)
t0 <- Sys.time()
M <- do.call(rbind, parLapply(cl, seq_len(B), one_rep))
stopCluster(cl)
cat(sprintf("B = %d replicates in %.1f min\n", B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
ALL <- data.frame(rep = seq_len(B), M, head = sha, check.names = FALSE)
write.csv(ALL, file.path(OUT, "pkg280_verify_fix_size_pvalues.csv"), row.names = FALSE)

TESTS <- c("stk.joint", "sym.unit", "sym.score", "stukel.unit", "stukel.score")
SUM <- list()
for (v in TESTS) {
  p <- M[, v]; ok <- is.finite(p); m <- sum(ok)
  s05 <- mean(p[ok] <= 0.05); se05 <- sqrt(0.05 * 0.95 / m); z05 <- (s05 - 0.05) / se05
  s01 <- mean(p[ok] <= 0.01); se01 <- sqrt(0.01 * 0.99 / m); z01 <- (s01 - 0.01) / se01
  SUM[[length(SUM) + 1]] <- data.frame(test = v, m = m, no_p = B - m, size05 = s05, se05 = se05, z05 = z05,
                                       size01 = s01, se01 = se01, z01 = z01, head = sha)
  cat(sprintf("  %-13s m %4d  size05 %.4f (MCSE %.4f, z %+5.2f)%s  size01 %.4f (MCSE %.4f, z %+5.2f)%s\n",
              v, m, s05, se05, z05, if (abs(z05) > 3) " FLAG" else "     ", s01, se01, z01, if (abs(z01) > 3) " FLAG" else ""))
}
write.csv(do.call(rbind, SUM), file.path(OUT, "pkg280_verify_fix_size_summary.csv"), row.names = FALSE)
cat(sprintf("  Stukel joint df table: %s\n", paste(names(table(M[, "stk.joint.df"])), table(M[, "stk.joint.df"]), sep = ":", collapse = " ")))
cat(sprintf("  score df tables: sym %s | stukel %s\n",
            paste(names(table(M[, "sym.score.df"])), table(M[, "sym.score.df"]), sep = ":", collapse = " "),
            paste(names(table(M[, "stukel.score.df"])), table(M[, "stukel.score.df"]), sep = ":", collapse = " ")))
cat(sprintf("  Note 'G = 40 (auto)' on every DEF row: %d/%d | battery vs def.gof(G = 40) max abs p diff: %.2e\n",
            sum(M[, "notes_auto_ok"] == 1), B, max(M[, "max_abs_vs_defgof"], na.rm = TRUE)))
cat(sprintf("  min(events, non-events): min %d, median %g\n", as.integer(min(M[, "ev_min"])), median(M[, "ev_min"])))
cat("done", format(Sys.time()), "\n")
