## pkg280_verify_num_size2.R -- fresh-seed re-check of the two borderline unit-form sizes in pkg280_verify_num_size.R:
## DEF stukel unit at alpha = 0.01 for s = 2, n = 1000, G = auto (0.0170, z +3.15) and at alpha = 0.05 for base
## n = 1000, G = auto (0.0355, z -2.98). The unit form is unchanged from 2.7.0; this asks whether the Satterthwaite
## reference is the cause, by running Satterthwaite and Imhof on the same data. INSTALLED ebrahim.gof 2.8.0,
## 5 workers, B = 4000 per cell. Verification of the package, not a paper result.
## Run: Rscript pkg280_verify_num_size2.R > ../paper_EDGE/theory/pkg280_verify_num_size2.log 2>&1
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
PKGV <- as.character(packageVersion("ebrahim.gof")); stopifnot(PKGV == "2.8.0")
cat("ebrahim.gof", PKGV, "(installed) |", R.version.string, "|", format(Sys.time()), "\n")

VARS <- c("stukel.sat", "stukel.imh", "stukel.score", "sym.sat", "sym.imh", "poly3.sat", "poly3.imh", "stk.joint")
one_rep <- function(i, s, c0, n) {
  out <- setNames(rep(NA_real_, length(VARS)), VARS)
  xa <- runif(n, -3, 3); db <- rbinom(n, 1, 0.5)
  resp <- rbinom(n, 1, plogis(c0 + s * (0.6 * xa + 0.5 * db)))
  fit <- tryCatch(suppressWarnings(glm(resp ~ xa + db, family = binomial())), error = function(e) NULL)
  if (is.null(fit)) return(out)
  pv <- function(...) tryCatch(suppressWarnings(def.gof(fit, G = "auto", ...))$p_value, error = function(e) NA_real_)
  out["stukel.sat"] <- pv(basis = "stukel"); out["stukel.imh"] <- pv(basis = "stukel", method = "imhof")
  out["stukel.score"] <- pv(basis = "stukel", weights = "score")
  out["sym.sat"] <- pv(basis = "sym");       out["sym.imh"] <- pv(basis = "sym", method = "imhof")
  out["poly3.sat"] <- pv(basis = "poly3");   out["poly3.imh"] <- pv(basis = "poly3", method = "imhof")
  out["stk.joint"] <- tryCatch(run.all.gof(fit, tests = "Stukel", install = "no")$p_value, error = function(e) NA_real_)
  out
}
CELLS <- data.frame(label = c("s=2 n1000 Gauto", "base n1000 Gauto"), s = c(2, 1), c0 = c(0, 0), n = c(1000, 1000))
B <- 4000L
cl <- makeCluster(5L)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(ebrahim.gof))))
clusterExport(cl, c("VARS", "one_rep"))
t0 <- Sys.time(); ALL <- list()
for (k in seq_len(nrow(CELLS))) {
  ce <- CELLS[k, ]; clusterSetRNGStream(cl, 99020260L + 7L * k)
  M <- do.call(rbind, parLapply(cl, seq_len(B), one_rep, s = ce$s, c0 = ce$c0, n = ce$n))
  ALL[[k]] <- data.frame(cell = ce$label, rep = seq_len(B), M, pkg = PKGV, check.names = FALSE)
  cat(sprintf("\n== %s (B = %d, fresh seeds; %.1f min) ==\n", ce$label, B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  for (v in VARS) {
    p <- M[, v]; ok <- is.finite(p); m <- sum(ok)
    for (al in c(0.05, 0.01)) {
      sz <- mean(p[ok] <= al); se <- sqrt(al * (1 - al) / m); z <- (sz - al) / se
      cat(sprintf("  %-13s alpha %.2f  m %4d  size %.4f (MCSE %.4f, z %+5.2f)%s\n", v, al, m, sz, se, z,
                  if (abs(z) > 3) " FLAG" else ""))
    }
  }
}
stopCluster(cl)
write.csv(do.call(rbind, ALL), file.path(OUT, "pkg280_verify_num_size2_pvalues.csv"), row.names = FALSE)
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 2)), "\n")
