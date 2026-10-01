## map_pkg_nullcheck.R -- sanity check of the PROTOTYPES in map_pkg_proto.R (NOT a paper result): null size of the
## joint Stukel score with its one-df fallback, of the current marginal sum, and of EDGE poly3 / sym / stukel in unit
## and score form, on the run-L design (x ~ U(-3,3), d ~ Bern(0.5), eta = c0 + s(0.6x + 0.5d)). 5 workers, B = 2000.
## Per-replicate p-values go to the session scratchpad, not to the project.
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
PROTO <- edge_path("code/simulations/map_pkg_proto.R")
source(PROTO)
suppressPackageStartupMessages(library(parallel))
OUTDIR <- file.path(tempdir())
VARS <- c("Stk.marg", "Stk.joint", "fallback", "rho", "poly3.unit", "poly3.score", "sym.unit", "sym.score",
          "stukel.unit", "stukel.score", "stukel2.score", "stukel.score.mineig", "events_lt_G")

one_rep <- function(i, s, c0, n, G) {
  out <- setNames(rep(NA_real_, length(VARS)), VARS)
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); y <- rbinom(n, 1, plogis(c0 + s * (0.6 * x + 0.5 * d)))
  out["events_lt_G"] <- as.numeric(min(sum(y), n - sum(y)) < G)
  fit <- tryCatch(suppressWarnings(glm(y ~ x + d, family = binomial())), error = function(e) NULL)
  if (is.null(fit)) return(out)
  out["Stk.marg"] <- tryCatch(suppressWarnings(pk$gof_stukel(pk$.gof_context(fit))$p_value), error = function(e) NA_real_)
  sj <- tryCatch(proto_stukel(fit, "joint"), error = function(e) NULL)
  if (!is.null(sj)) { out["Stk.joint"] <- sj$p_value; out["fallback"] <- as.numeric(isTRUE(sj$fallback)); out["rho"] <- sj$rho }
  for (b in c("poly3", "sym", "stukel")) for (wt in c("unit", "score")) {
    r <- tryCatch(suppressWarnings(proto_def(fit, G = G, basis = b, weights = wt, diag = TRUE)), error = function(e) NULL)
    if (!is.null(r)) {
      out[paste0(b, ".", wt)] <- r$p_value
      if (b == "stukel" && wt == "score") out["stukel.score.mineig"] <- min(attr(r, "eig_rel"))
    }
  }
  r <- tryCatch(suppressWarnings(proto_def(fit, G = G, basis = "stukel2", weights = "score")), error = function(e) NULL)
  if (!is.null(r)) out["stukel2.score"] <- r$p_value
  out
}

CELLS <- data.frame(label = c("base n500 G10", "base n1000 G40", "c0=-2 n500 G20", "s=2 n1000 G40"),
                    s = c(1, 1, 1, 2), c0 = c(0, 0, -2, 0), n = c(500, 1000, 500, 1000), G = c(10, 40, 20, 40))
B <- 2000L
cl <- makeCluster(5L)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl, c("PROTO", "VARS", "one_rep"))
invisible(clusterEvalQ(cl, { source(PROTO); NULL }))
t0 <- Sys.time(); ALL <- list()
for (k in seq_len(nrow(CELLS))) {
  ce <- CELLS[k, ]; clusterSetRNGStream(cl, 20260914L + k)
  M <- do.call(rbind, parLapply(cl, seq_len(B), one_rep, s = ce$s, c0 = ce$c0, n = ce$n, G = ce$G))
  ALL[[k]] <- data.frame(cell = ce$label, rep = seq_len(B), M, check.names = FALSE)
  cat(sprintf("\n== %s  (B = %d, %.1f min) ==\n", ce$label, B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  for (v in c("Stk.marg", "Stk.joint", "poly3.unit", "poly3.score", "sym.unit", "sym.score",
              "stukel.unit", "stukel.score", "stukel2.score")) {
    p <- M[, v]; ok <- is.finite(p); m <- sum(ok)
    cat(sprintf("  %-14s size05 %.4f (se %.4f)  size01 %.4f (se %.4f)  no-p %.3f\n", v,
                mean(p[ok] <= 0.05), sqrt(.05 * .95 / max(m, 1)), mean(p[ok] <= 0.01), sqrt(.01 * .99 / max(m, 1)), 1 - m / B))
  }
  fb <- M[, "fallback"] %in% 1
  cat(sprintf("  joint fallback rate %.3f; joint size05 within fallback reps %.4f (m = %d); events < G rate %.3f\n",
              mean(fb), if (any(fb)) mean(M[fb, "Stk.joint"] <= .05, na.rm = TRUE) else NA, sum(fb), mean(M[, "events_lt_G"])))
  cat(sprintf("  post-fit corr of the two Stukel columns: median %.4f\n", median(M[, "rho"], na.rm = TRUE)))
  q <- quantile(M[, "stukel.score.mineig"], c(0, .01, .5), na.rm = TRUE)
  cat(sprintf("  stukel (3-col) score form, min relative eigenvalue of I: min %.2e  1%% %.2e  median %.2e\n", q[1], q[2], q[3]))
}
stopCluster(cl)
write.csv(do.call(rbind, ALL), file.path(OUTDIR, "map_pkg_nullcheck_pvalues.csv"), row.names = FALSE)
cat("\nwall", format(round(difftime(Sys.time(), t0, units = "mins"), 2)), "\n")
