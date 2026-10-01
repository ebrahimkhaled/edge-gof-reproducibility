## analyse_corruption_partitions.R -- false alarms under the two covariate errors for the forms and rivals paper 2 did
## not tabulate: EDGE (cubic and symmetric basis) at ten groups as well as at the rule G, and the grouped rivals the
## battery ran (Hosmer-Lemeshow on equal-width risk groups HL_w, Pigeon-Heyse, Tsiatis, Xie). Read from the deposited per-replicate
## p-values of block 9 (the exaggeration, x -> 4x) and block C1b (the sign error, x -> -4x), logistic truth, so every
## rejection is a false alarm. A replicate without a p-value counts as a non-rejection, as everywhere in the paper.
## The published cells are recomputed first and must agree, or the script stops.
##
##   Rscript analyse_corruption_partitions.R   -> battery/analysis/corruption_partitions.csv
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
suppressMessages(library(data.table))
SIM <- edge_path("code/simulations")
B <- edge_battery()
TESTS <- c("EDGE.poly3.u.Grule", "EDGE.poly3.u.G10", "EDGE.sym.u.Grule", "EDGE.sym.u.G10", "HL.Grule", "HL.G10",
           "HL_w", "PH", "Tsiatis", "Xie", "Stk.joint", "Cubic.LR", "GiViTI")
rate_k <- c(clean = 0, r001 = 0.001, r002 = 0.002, r005 = 0.005, r010 = 0.010)

one <- function(file, error, n, rate) {
  P <- fread(file)
  k <- round(rate_k[[rate]] * n)
  data.table(error = error, n = n, k = k, test = TESTS, reps = nrow(P),
             false_alarm = vapply(TESTS, function(t) mean(!is.na(P[[t]]) & P[[t]] <= 0.05), 0))
}
R <- list()
for (n in c(1000L, 5000L)) for (r in names(rate_k)) {
  f9 <- edge_battery("9", if (r == "clean") sprintf("logit_clean_n%d_pvalues.csv.gz", n) else sprintf("logit_C1_%s_n%d_pvalues.csv.gz", r, n))
  fc <- edge_battery("C1b", if (r == "clean") sprintf("c1b_clean_n%d_pvalues.csv.gz", n) else sprintf("c1b_%s_n%d_pvalues.csv.gz", r, n))
  R[[length(R) + 1]] <- one(f9, "exaggeration", n, r)
  R[[length(R) + 1]] <- one(fc, "sign error", n, r)
}
R <- rbindlist(R)

## the published cells (Table 5 of paper 3, from battery/C1b/_summary.csv) must be reproduced
pub <- fread(edge_battery("C1b", "_summary.csv"))[test %in% c("EDGE", "EDGE-sym", "Stukel joint", "cubic LR")]
pub[, test := c(EDGE = "EDGE.poly3.u.Grule", `EDGE-sym` = "EDGE.sym.u.Grule", `Stukel joint` = "Stk.joint",
                `cubic LR` = "Cubic.LR")[test]]
chk <- merge(melt(pub, id.vars = c("n", "k", "test"), measure.vars = c("sign_error", "exaggerated"),
                  variable.name = "error", value.name = "published")[, error := ifelse(error == "sign_error", "sign error", "exaggeration")],
             R, by = c("error", "n", "k", "test"))
d <- max(abs(chk$published - round(chk$false_alarm, 3)))
cat(sprintf("published cells reproduced: %d, max |difference| = %.3f\n", nrow(chk), d))
if (nrow(chk) < 40 || d > 0.0005) stop("the per-replicate files do not reproduce the published table")

fwrite(R, edge_battery("analysis", "corruption_partitions.csv"))
print(dcast(R, error + n + k ~ test, value.var = "false_alarm")[, c("error", "n", "k", TESTS), with = FALSE], digits = 3)
