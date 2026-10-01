## pkg280_verify_corr_size.R -- null size of the installed ebrahim.gof 2.8.0 on the sparse design after the correction.
## Data: gen_sparse_link(n, -4.9) of _dgp_library.R (x = standardised chi-square(4), eta = -4.9 + x, fitted y ~ x),
## n = 200 and 300, B = 2000 each, new seeds 78e6 + 1000 n + r. Tests: the Stukel row (joint form, the default) and
## def.gof with basis stukel and sym in the unit form (Satterthwaite, and Imhof as a side check) and the score form, at
## G = 10 and the rule G = max(10, round(n / 25)). A missing p-value (no fit, guard, error) counts as no rejection.
## Flag: rejection rate at 0.05 above 0.05 + 3 MCSE. Per-replicate p-values are stored.
## Verification of the package, not a paper result.
## Run: Rscript pkg280_verify_corr_size.R <workers <= 8> > ../paper_EDGE/theory/pkg280_verify_corr_size.log 2>&1
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
NW  <- min(8L, as.integer(args[1]))
SIM <- edge_path("code/simulations")
OUT <- edge_path("declarations")
options(width = 200, warn = 1)
library(parallel)
t0 <- Sys.time()
B <- 2000L

size_rep <- function(r, n) {
  set.seed(78000000 + 1000 * n + r)
  g <- gen_sparse_link(n, intercept = -4.9)
  fit <- suppressWarnings(glm(g$f, family = binomial(), data = g$d))
  out <- c(r = r, n = n, events = sum(fit$y))
  st <- tryCatch(suppressWarnings(run.all.gof(fit, tests = "Stukel")), error = function(e) NULL)
  out["Stukel.joint"] <- if (is.null(st)) NA else st$p_value[1]
  out["Stukel.joint.df"] <- if (is.null(st)) NA else st$df[1]
  out["Stukel.note.nofit"] <- as.numeric(!is.null(st) && grepl("no maximum-likelihood", st$Note[1]))
  out["Stukel.note.noinfo"] <- as.numeric(!is.null(st) && grepl("no information", st$Note[1]))
  for (G in unique(c(10, max(10, round(n / 25))))) for (b in c("stukel", "sym")) for (fm in c("unit", "imhof", "score")) {
    cls <- character(0)
    v <- tryCatch(withCallingHandlers(
           switch(fm, unit = def.gof(fit, G = G, basis = b), imhof = def.gof(fit, G = G, basis = b, method = "imhof"),
                      score = def.gof(fit, G = G, basis = b, weights = "score"))$p_value,
           warning = function(w) { cls <<- c(cls, class(w)[1]); invokeRestart("muffleWarning") }),
         error = function(e) -1)
    nm <- sprintf("DEF.%s.%s.G%d", if (b == "stukel") "stk" else b, fm, G)
    out[nm] <- if (identical(v, -1)) NA else v
    out[paste0(nm, ".err")] <- as.numeric(identical(v, -1))
    out[paste0(nm, ".noinfo")] <- as.numeric(any(cls == "def_no_information"))
    out[paste0(nm, ".degenerate")] <- as.numeric(any(cls == "def_degenerate"))
  }
  out
}

cat("pkg280_verify_corr_size | workers", NW, "| B", B, "|", format(t0), "\n")
cl <- makeCluster(NW)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
invisible(clusterCall(cl, function(SIM) { suppressPackageStartupMessages(library(ebrahim.gof)); source(file.path(SIM, "_dgp_library.R")); NULL }, SIM))
clusterExport(cl, "size_rep")
cat("installed ebrahim.gof", clusterCall(cl[1], function() as.character(packageVersion("ebrahim.gof")))[[1]], "\n")
res <- list()
for (n in c(200L, 300L)) {
  m <- do.call(rbind, parLapplyLB(cl, seq_len(B), size_rep, n = n))
  res[[as.character(n)]] <- as.data.frame(m)
  cat(sprintf("  n = %d done at %.1f min\n", n, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
cols <- union(names(res[["200"]]), names(res[["300"]]))
P <- do.call(rbind, lapply(res, function(d) { for (k in setdiff(cols, names(d))) d[[k]] <- NA; d[, cols] }))
rownames(P) <- NULL
write.csv(P, file.path(OUT, "pkg280_verify_corr_size_pvalues.csv"), row.names = FALSE)

mcse <- sqrt(0.05 * 0.95 / B); lim <- 0.05 + 3 * mcse
cat(sprintf("\nMCSE at 0.05 with B = %d: %.5f; flag above %.5f. NA counts as no rejection.\n", B, mcse, lim))
tests <- grep("^(Stukel.joint|DEF\\.[a-z]+\\.[a-z]+\\.G[0-9]+)$", cols, value = TRUE)
S <- do.call(rbind, lapply(split(P, P$n), function(d) do.call(rbind, lapply(tests, function(t) {
  if (all(is.na(d[[t]])) && !grepl("Stukel", t) && !(paste0(t, ".err") %in% names(d) && any(!is.na(d[[paste0(t, ".err")]])))) return(NULL)
  p <- d[[t]]
  rej <- sum(p < 0.05, na.rm = TRUE)
  data.frame(n = d$n[1], test = t, B = nrow(d), rejections = rej, rate = rej / nrow(d), na = sum(is.na(p)),
             rate_given_p = rej / max(1, sum(!is.na(p))),
             noinfo = if (grepl("Stukel", t)) sum(d$Stukel.note.noinfo) else sum(d[[paste0(t, ".noinfo")]], na.rm = TRUE),
             nofit = if (grepl("Stukel", t)) sum(d$Stukel.note.nofit) else sum(d[[paste0(t, ".degenerate")]], na.rm = TRUE),
             errors = if (grepl("Stukel", t)) NA else sum(d[[paste0(t, ".err")]], na.rm = TRUE),
             above_limit = rej / nrow(d) > lim)
}))))
rownames(S) <- NULL
print(S, digits = 4)
write.csv(S, file.path(OUT, "pkg280_verify_corr_size_summary.csv"), row.names = FALSE)
cat(sprintf("\nrows above 0.05 + 3 MCSE: %d of %d\n", sum(S$above_limit), nrow(S)))
cat("\nevents per sample: \n"); print(aggregate(events ~ n, data = P, FUN = function(v) c(mean = mean(v), zero = sum(v == 0), le3 = sum(v <= 3))))
cat("\ndone", format(Sys.time()), sprintf("(%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
