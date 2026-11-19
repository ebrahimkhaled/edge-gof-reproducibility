## check_edge_dxd_equality.R -- edge.gof()/def.gof() without the G x G matrix Omega: same answers?
## The new code (branch edge-fast-omega) builds every B' Omega B from U and the basis and never forms
## Omega. This script runs the released 3.1.0 code and the new code on the same data sets, in two
## separate R processes, and compares the statistic, df, p-value and the diagonal of Omega kept for
## the plot. 250 data sets: 50 at each n in {100, 500, 1000, 5000, 20000}; on each, the bases poly3,
## sym and stukel, the weights unit and score, G = "auto" and 10, internal and external mode (the
## external mode takes unit weights only). Every fifth data set is also run without X (Omega = I).
##
##   Rscript check_edge_dxd_equality.R          run both versions, compare, write battery/T/edge_dxd_equality.csv
##   Rscript check_edge_dxd_equality.R --child <pkg> <out.rds>     (internal: one version)
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
SIMDIR <- edge_path("code/simulations")
## The two package sources compared: OLD is ebrahim.gof 3.1.0 from CRAN, NEW is the source in
## code/package/ (the same package with the d x d computation). Unpack both and point these at them.
OLD    <- Sys.getenv("EBGOF_OLD_SRC", edge_path("ebrahim.gof-3.1.0"))
NEW    <- Sys.getenv("EBGOF_NEW_SRC", edge_path("ebrahim.gof-dxd"))
args   <- commandArgs(TRUE)

make_data <- function(i, n) {
  set.seed(20261119 + 1000 * i + n)
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  bend <- c(0, 0.3, 0.6)[i %% 3 + 1]                       # a correct model, and two misspecified ones
  y <- rbinom(n, 1, plogis(-0.5 + 0.8 * x1 + 0.5 * x2 + bend * x1^2))
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  glm(y ~ x1 + x2, family = binomial(), data = dat)
}
NS   <- c(100L, 500L, 1000L, 5000L, 20000L)
PER  <- 50L
CONF <- rbind(expand.grid(basis = c("poly3", "sym", "stukel"), weights = c("unit", "score"),
                          G = c("auto", "10"), external = FALSE, stringsAsFactors = FALSE),
              expand.grid(basis = c("poly3", "sym", "stukel"), weights = "unit",
                          G = c("auto", "10"), external = TRUE, stringsAsFactors = FALSE))

if (length(args) >= 3 && args[1] == "--child") {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  suppressPackageStartupMessages(pkgload::load_all(args[2], quiet = TRUE, export_all = TRUE))
  res <- list()
  for (n in NS) for (i in seq_len(PER)) {
    fit <- make_data(i, n)
    inputs <- list(glm = list(object = fit, X = NULL, pp = NULL))
    if (i %% 2 == 0) inputs$yX <- list(object = fit$y, X = model.matrix(fit), pp = fitted(fit))
    if (i %% 5 == 0) inputs$naive <- list(object = fit$y, X = NULL, pp = fitted(fit))
    for (inp in names(inputs)) for (k in seq_len(nrow(CONF))) {
      cf <- CONF[k, ]
      if (cf$external && inp != "glm") next
      g  <- if (cf$G == "auto") "auto" else 10
      a  <- inputs[[inp]]
      o  <- tryCatch(suppressWarnings(.def_gof(a$object, predicted_probs = a$pp, X = a$X, G = g,
                                             basis = cf$basis, weights = cf$weights,
                                             external = cf$external, keep = TRUE)),
                     error = function(e) NULL)
      res[[length(res) + 1]] <- data.frame(
        n = n, dataset = i, input = inp, basis = cf$basis, weights = cf$weights, G = cf$G,
        external = cf$external,
        stat = if (is.null(o)) NA_real_ else o$Test_Statistic,
        df   = if (is.null(o)) NA_real_ else o$df,
        p    = if (is.null(o)) NA_real_ else o$p_value,
        omega_diag = if (is.null(o)) "" else paste(sprintf("%.17g", attr(o, "details")$omega_diag), collapse = ";"),
        error = is.null(o), stringsAsFactors = FALSE)
    }
  }
  saveRDS(do.call(rbind, res), args[3])
  quit(save = "no")
}

## ---- the parent: run both versions one after the other, then compare ----
rs  <- file.path(R.home("bin"), "Rscript")
me  <- file.path(SIMDIR, "check_edge_dxd_equality.R")
tmp <- tempfile(c("old", "new"), fileext = ".rds")
system2(rs, c(shQuote(me), "--child", shQuote(OLD), shQuote(tmp[1])))
system2(rs, c(shQuote(me), "--child", shQuote(NEW), shQuote(tmp[2])))
A <- readRDS(tmp[1]); B <- readRDS(tmp[2])
stopifnot(nrow(A) == nrow(B), identical(A[, 1:7], B[, 1:7]))
dmax <- function(a, b) {
  a <- as.numeric(strsplit(a, ";")[[1]]); b <- as.numeric(strsplit(b, ";")[[1]])
  if (length(a) != length(b)) return(Inf)
  if (!length(a)) 0 else max(abs(a - b))
}
same_na <- function(a, b) (is.na(a) & is.na(b))
out <- data.frame(A[, 1:7],
                  stat_old = A$stat, stat_new = B$stat, abs_diff_stat = ifelse(same_na(A$stat, B$stat), 0, abs(A$stat - B$stat)),
                  df_old = A$df, df_new = B$df, abs_diff_df = ifelse(same_na(A$df, B$df), 0, abs(A$df - B$df)),
                  p_old = A$p, p_new = B$p, abs_diff_p = ifelse(same_na(A$p, B$p), 0, abs(A$p - B$p)),
                  abs_diff_omega_diag = mapply(dmax, A$omega_diag, B$omega_diag, USE.NAMES = FALSE),
                  error_old = A$error, error_new = B$error)
dir.create(edge_battery("T"), showWarnings = FALSE, recursive = TRUE)
write.csv(out, edge_battery("T", "edge_dxd_equality.csv"), row.names = FALSE)
cat(sprintf("%d comparisons on %d data sets; errors old %d, new %d, mismatched errors %d\n",
            nrow(out), length(unique(paste(out$n, out$dataset))), sum(out$error_old), sum(out$error_new),
            sum(out$error_old != out$error_new)))
cat(sprintf("max |diff|: statistic %.3g, df %.3g, p-value %.3g, diag(Omega) %.3g\n",
            max(out$abs_diff_stat, na.rm = TRUE), max(out$abs_diff_df, na.rm = TRUE),
            max(out$abs_diff_p, na.rm = TRUE), max(out$abs_diff_omega_diag, na.rm = TRUE)))
cat(sprintf("NA mismatches: statistic %d, p-value %d\n",
            sum(xor(is.na(out$stat_old), is.na(out$stat_new))), sum(xor(is.na(out$p_old), is.na(out$p_new)))))
