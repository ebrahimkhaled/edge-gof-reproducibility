## map_pkg_eventsG.R -- does the package's def.gof return a p-value when events < groups? (MAP_package.md) Read-only.
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
set.seed(5); x <- runif(400, -3, 3); y <- rbinom(400, 1, plogis(-4 + 0.5 * x)); fit <- glm(y ~ x, family = binomial())
cat(sprintf("events = %d, n = %d\n", sum(y), length(y)))
for (G in c(10, 20, 40)) {
  a <- tryCatch(pk$def.gof(fit, G = G)$p_value, error = function(e) paste("ERROR", conditionMessage(e)))
  b <- tryCatch(proto_def(fit, G = G, basis = "sym")$p_value, error = function(e) paste("ERROR", conditionMessage(e)))
  s <- tryCatch(proto_def(fit, G = G, basis = "sym", weights = "score")$p_value, error = function(e) paste("ERROR", conditionMessage(e)))
  r <- as.data.frame(pk$run.all.gof(fit, tests = c("DEF.poly3", "HL", "EF"), G = G, install = "no"))
  cat(sprintf("G = %2d (events < G: %s): def.gof poly3 p = %s | sym unit p = %s | sym score p = %s | battery HL p = %s, EF p = %s\n",
              G, sum(y) < G, format(a), format(b), format(s), format(r$p_value[r$Test == "HL"]), format(r$p_value[r$Test == "EF"])))
}
