## make_tab_external.R -- Online Resource 1, Table S6 (tab_external.tex): external validation.
## The table block of make_tables_paper3.R, with one more column: the cumulative calibration test of
## Sadatsafavi and Petkau (cumulcalib 0.2.0, two-part Brownian-bridge test), block EXT_SP, run on block EXT's
## data sets under PREDECLARATION_blockEXT_SP.md (sha256 7455de2d...).
## Sources, and nothing else:
##   battery/EXT/_summary.csv      EDGE (ten groups, default), Cox, Spiegelhalter, GiViTI, HL, Stukel
##   battery/EXT_SP/_summary.csv   SP_BB (and, as a check, EDGE_G10 / cox repeated from EXT)
## Writes into this folder's tables/ only.
suppressMessages(library(data.table))
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
OUT <- edge_path("manuscript", "tables")
f3 <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
wr <- function(txt, file) { writeLines(txt, file.path(OUT, file)); cat("wrote", file, "\n") }

EX0 <- fread(edge_battery("EXT", "_summary.csv"))
SP  <- fread(edge_battery("EXT_SP", "_summary.csv"))
## the rows EXT_SP repeats from EXT must agree with EXT, or the two files are not the same study
chk <- merge(EX0[test %in% c("EDGE_G10", "EDGE_auto", "cox")], SP[test %in% c("EDGE_G10", "EDGE_auto", "cox")],
             by = c("id", "test"))
stopifnot(nrow(chk) == 63L, isTRUE(all.equal(chk$raw.x, chk$raw.y)), isTRUE(all.equal(chk$adj.x, chk$adj.y)))
EX <- rbind(EX0, SP[test == "SP_BB"], fill = TRUE)

ext_tests <- c("EDGE_G10", "EDGE_auto", "cox", "SP_BB", "spiegelhalter", "giviti", "HL", "stukel")
ext_row <- function(z, label, v = "raw") sprintf("%s & %s \\\\", label,
  paste(vapply(ext_tests, function(t) { x <- z[z$test == t, ][[v]]
    if (v == "raw" && z$role[1] == "corrupted" && x > 0.10) paste0("\\cellcolor{black!12}", f3(x)) else f3(x) }, ""), collapse = " & "))
FAMX <- c(large = "calibration in the large", slope = "calibration slope", ushape = "U-shaped term",
          thresh = "threshold term", inter = "interaction")
ex_ids <- function(r) unique(EX[role == r, .(id, fam, C, n, k)])[order(id)]
ex_lev <- ex_ids("level"); ex_alt <- ex_ids("alternative"); ex_cor <- ex_ids("corrupted")
mean_adj <- EX[role == "alternative", .(m = mean(adj)), by = test]
held <- EX[role == "corrupted", .(h = sum(raw <= 0.10)), by = test]
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{External validation: a published model with frozen coefficients checked on new patients. Five",
     "standard normal covariates; the truth departs from the published linear predictor in five ways, two",
     "strengths each, at $n=500$; in the last block the model is correct and $k$ of $1000$ records carry the",
     "second covariate multiplied by four (exaggeration) or by $-4$ (sign error), with their prediction",
     "recomputed at the recorded value. EDGE is in external mode (cubic basis and constant, four",
     "degrees of freedom) at ten groups and at the default partition. Cox is the recalibration test on",
     "intercept and slope, SP the two-part Brownian-bridge test of Sadatsafavi and Petkau,",
     "Spiegelhalter the $z$ test, GiViTI the belt in external mode, HL the",
     "Hosmer--Lemeshow statistic at ten groups and Stukel his two terms added to the frozen linear predictor.",
     "$1000$ replicates a scenario, the same data sets for every test. Power is size-adjusted at the null",
     "with $n=500$; shaded false-alarm rates exceed $0.10$.}",
     "\\label{tab:external}",
     "\\resizebox{\\textwidth}{!}{\\begin{tabular}{lrrrrrrrr}\\hline",
     "Scenario & EDGE, $G{=}10$ & EDGE, default & Cox & SP & Spiegelhalter & GiViTI & HL & Stukel \\\\ \\hline",
     "\\multicolumn{9}{l}{\\emph{Level: rejection rate of a correct model}} \\\\",
     vapply(seq_len(nrow(ex_lev)), function(i) ext_row(EX[id == ex_lev$id[i]], sprintf("\\quad $n=%d$", ex_lev$n[i])), ""),
     "\\multicolumn{9}{l}{\\emph{Power, size-adjusted, $n=500$}} \\\\",
     vapply(seq_len(nrow(ex_alt)), function(i) ext_row(EX[id == ex_alt$id[i]],
       sprintf("\\quad %s, $C=%s$", FAMX[[ex_alt$fam[i]]], format(ex_alt$C[i])), "adj"), ""),
     sprintf("\\quad mean of the ten & %s \\\\", paste(f3(mean_adj[match(ext_tests, test), m]), collapse = " & ")),
     "\\multicolumn{9}{l}{\\emph{False alarms, correct model, $n=1000$}} \\\\",
     vapply(seq_len(nrow(ex_cor)), function(i) ext_row(EX[id == ex_cor$id[i]],
       sprintf("\\quad %s, $k=%d$", if (ex_cor$fam[i] == "exag") "exaggeration" else "sign error", ex_cor$k[i])), ""),
     sprintf("\\quad settings held, of eight & %s \\\\", paste(held[match(ext_tests, test), h], collapse = " & ")),
     "\\hline\\end{tabular}}\\end{table}"), "tab_external.tex")
