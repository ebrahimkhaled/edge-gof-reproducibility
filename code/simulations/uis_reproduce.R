## Persist the EXACT UIS runs behind the §7 vignette (reproducibility fix, Phase D).
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
suppressMessages({ library(stats); library(ebrahim.gof); library(splines) })
SIM <- edge_path("code/simulations")
source(file.path(SIM, "_proj_test.R"))
uis <- readRDS(file.path(SIM, "uis_data.rds"))
d <- uis; d$y <- as.integer(d$DFREE == "no"); d$NDRGFP1 <- 10/(d$NDRGTX + 1)

fits <- list(
  linear_ndrgtx = glm(y ~ AGE + NDRGTX + IVHX + RACE + TREAT + SITE, data=d, family=binomial),
  ## Liu et al. (2024) model 19 = Hosmer's final model with NDRGTX entered linearly
  hosmer_model19 = glm(y ~ AGE + NDRGTX + IVHX + RACE + TREAT + SITE + AGE:NDRGFP1 + RACE:SITE,
                       data=d, family=binomial))

rows <- list()
for (nm in names(fits)) {
  f <- fits[[nm]]
  g <- run.all.gof(f)
  rows[[length(rows)+1]] <- data.frame(model=nm, test=g$Test, statistic=g$Statistic,
                                       df=g$df, pvalue=g$p_value, stringsAsFactors=FALSE)
  pj <- proj_pvalue(f$y, model.matrix(f), B=1000, seed=2026)
  rows[[length(rows)+1]] <- data.frame(model=nm, test="proj(Escanciano-Liu)",
                                       statistic=pj$stat, df=NA, pvalue=pj$p_value, stringsAsFactors=FALSE)
  cat(sprintf("[%s] n=%d  proj p=%.4f (B=1000, min resolution 1/1001)\n", nm, length(f$y), pj$p_value))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(SIM, "uis_gof_reproducible.csv"), row.names=FALSE)

## echo the exact numbers the vignette cites (model 19)
m19 <- res[res$model=="hosmer_model19", ]
key <- function(t) round(m19$pvalue[m19$test==t], 4)
cat(sprintf("\nMODEL 19 vignette numbers: proj=%.4f EDGE-poly3=%.3f EDGE-stk=%.3f HL=%.3f Stukel=%.3f Tsiatis=%.3f Xie=%.3f\n",
    key("proj(Escanciano-Liu)"), key("DEF.poly3"), key("DEF.stukel"), key("HL"), key("Stukel"), key("Tsiatis"), key("Xie")))
cat("WROTE uis_gof_reproducible.csv  rows:", nrow(res), "\nDONE\n")
