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
suppressMessages({ library(stats); library(ebrahim.gof) })
SIM <- edge_path("code/simulations")
source(file.path(SIM, "_proj_test.R"))
uis <- readRDS(file.path(SIM, "uis_data.rds"))
d <- uis; d$y <- as.integer(d$DFREE == "no")
d$IVn <- as.integer(d$IVHX)                 # numeric IV coding 1/2/3 (Liu's "1=never,2=prev,3=recent")
d$NDRGFP1 <- 10/(d$NDRGTX + 1)

key <- function(g, nm) g$p_value[g$Test==nm]
show <- function(fit, tag){
  X <- model.matrix(fit)
  pj <- proj_pvalue(fit$y, X, B=1000, seed=7)$p_value
  g <- run.all.gof(fit)
  cat(sprintf("\n[%s]  n=%d p=%d\n", tag, length(fit$y), ncol(X)))
  cat(sprintf("  proj=%.4f | EDGE-poly3=%.4f EDGE-stk=%.4f | HL=%.4f Stukel=%.4f | Tsiatis=%.4f Xie=%.4f Xie-GAM=%s\n",
      pj, key(g,"DEF.poly3"), key(g,"DEF.stukel"), key(g,"HL"), key(g,"Stukel"),
      key(g,"Tsiatis"), key(g,"Xie"),
      ifelse(length(key(g,"Xie-GAM")), sprintf("%.4f",key(g,"Xie-GAM")), "NA")))
}

set.seed(20260707)
## (a) simple model, IVHX as NUMERIC (Liu may have coded it this way)
f_num <- glm(y ~ AGE + NDRGTX + IVn + RACE + TREAT + SITE, data=d, family=binomial)
show(f_num, "model17, IVHX numeric")

## (b) interaction-context model with LINEAR NDRGTX (Liu model 19): + AGE:NDRGFP1 + RACE:SITE
f19 <- glm(y ~ AGE + NDRGTX + IVHX + RACE + TREAT + SITE + AGE:NDRGFP1 + RACE:SITE,
           data=d, family=binomial)
show(f19, "model19 (interactions, linear NDRGTX)")

cat("\nDONE2\n")
