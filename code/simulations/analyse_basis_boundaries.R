## analyse_basis_boundaries.R -- at each usage boundary of Section 8, can a different BASIS of the same
## test do what the paper currently says needs a different TEST?
##
## The paper's "when not to" list sends the analyst to a rival in six situations. The directed test has
## three bases, and Table 3 already shows that none leads everywhere, so the question is whether the
## boundary is a boundary of the test or only of its default basis. This reads the census that is
## already deposited; nothing new is computed and no hypothesis was declared about it, so the result is
## a reading of the census in the same sense as the cubic-LR comparison of Section 5.6.
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
suppressPackageStartupMessages(library(data.table))
setwd(SIMDIR)
AN <- edge_battery("analysis")

P <- fread(edge_path("results/analysis/_census_power_paper2.csv"))
MEM <- unique(fread(file.path(AN, "rule_A_membership.csv"))[, .(cell, family, family_name)], by = "cell")
P <- merge(P, MEM, by = "cell", all.x = TRUE)

BASES <- c(`EDGE-poly3` = "EDGE.poly3.u.Grule", `EDGE-stk` = "EDGE.stk.u.Grule",
           `EDGE-sym` = "EDGE.sym.u.Grule")
RIVALS <- c(GiViTI = "GiViTI", `cubic LR` = "Cubic.LR", `Stukel joint` = "Stk.joint",
            `Stukel one-par.` = "Stk.sym1", `HL (rule G)` = "HL.Grule")
keep <- c(BASES, RIVALS)
cat("tests present in the census:\n"); print(intersect(keep, unique(P$test)))

W <- dcast(P[test %in% keep], cell + n + family + family_name ~ test, value.var = "power")
setnames(W, BASES, names(BASES), skip_absent = TRUE)
setnames(W, RIVALS, names(RIVALS), skip_absent = TRUE)

## the six situations of Section 8, as the scenarios that define them
SIT <- list(
  `events are scarce`      = quote(grepl("_e12", cell)),
  `asymmetric links`       = quote(family == 2),
  `crossover`              = quote(grepl("crossover", cell)),
  `rough misfit`           = quote(family == 5),
  `skewed covariate`       = quote(grepl("skew", cell)),
  `symmetric tails`        = quote(family == 1)
)
RIVAL_OF <- c(`events are scarce` = "GiViTI", `asymmetric links` = "cubic LR",
              `crossover` = "GiViTI", `rough misfit` = "GiViTI",
              `skewed covariate` = "Stukel one-par.", `symmetric tails` = "GiViTI")

cat("\n== at each boundary: the best EDGE basis against the leading rival named there ==\n")
out <- list()
for (nm in names(SIT)) {
  X <- W[eval(SIT[[nm]])]
  if (!nrow(X)) { cat(sprintf("\n-- %-20s no scenarios matched\n", nm)); next }
  riv <- RIVAL_OF[[nm]]
  if (!riv %in% names(X)) { cat(sprintf("\n-- %-20s rival %s absent\n", nm, riv)); next }
  m <- sapply(names(BASES), function(b) if (b %in% names(X)) mean(X[[b]], na.rm = TRUE) else NA_real_)
  mr <- mean(X[[riv]], na.rm = TRUE)
  best <- names(which.max(m))
  ## how often each basis leads or ties the rival, at the paper's 0.01 margin
  lt <- sapply(names(BASES), function(b)
    if (b %in% names(X)) sum(X[[b]] - X[[riv]] > -0.01, na.rm = TRUE) else NA_integer_)
  cat(sprintf("\n-- %s (%d scenarios), rival = %s (mean %.3f)\n", nm, nrow(X), riv, mr))
  for (b in names(BASES))
    cat(sprintf("     %-11s mean %.3f   leads or ties in %2d of %d\n", b, m[[b]], lt[[b]], nrow(X)))
  cat(sprintf("     best basis: %s (%.3f), gap to rival %+.3f\n", best, m[[best]], m[[best]] - mr))
  out[[nm]] <- data.table(situation = nm, scenarios = nrow(X), rival = riv, rival_mean = mr,
                          best_basis = best, best_mean = m[[best]],
                          default_mean = m[["EDGE-poly3"]],
                          best_leads_ties = lt[[best]])
}
R <- rbindlist(out)
fwrite(R, file.path(AN, "basis_at_boundaries.csv"))
cat("\nwrote battery/analysis/basis_at_boundaries.csv\n")
