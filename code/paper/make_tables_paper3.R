## make_tables_paper3.R -- every table of EDGE paper 3, generated from the deposits (copied from paper 2; output now paper3/tables).
##
## Sources, and nothing else:
##   battery/analysis/hypotheses.csv   the pre-declared hypotheses and their verdicts (E12)
##   battery/analysis/mcnemar.csv      per-cell paired comparisons against GiViTI and Stukel
##   battery/analysis/rule_A_gate.csv  the size of each EDGE variant over 131 null cells
##   battery/analysis/rule_A_families.csv  mean power by departure family
##   battery/8/_summary.csv, _paired.csv   the projection test and BAGofT (block 8)
##   battery/9/analysis/cell_test.csv      contamination (block 9)
##   battery/9c/analysis/cell_G_test.csv   the partition surface (block 9c)
##   battery/9b/analysis/rates_9b.csv, seconds_9b.csv   the resampling tests under corruption (block 9b)
##
## No number is recomputed here from raw p-values: the analysis scripts are the source of record,
## and a table that disagreed with them would be a second opinion, not a table.
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
BAT <- edge_battery()
OUT <- Sys.getenv("EDGE_TABLES_OUT", edge_path("manuscript/tables"))   # archive: set EDGE_TABLES_OUT to write elsewhere
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
f3 <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
wr <- function(txt, file) { writeLines(txt, file.path(OUT, file)); cat("wrote", file, "\n") }

## ---------------------------------------------------------------------------------------------
## Table 1: the pre-declared hypotheses and what happened to them
## ---------------------------------------------------------------------------------------------
H <- fread(edge_battery("analysis", "hypotheses.csv"))[forms_used == TRUE]
## H1 is declared over family 1 (symmetric tails) and H3 over family 2 (asymmetric links); the labels name
## the families the declaration names, and hypotheses.csv carries the rule each was read against.
lab <- c(H1 = "EDGE-poly3 leads GiViTI on symmetric tails",
         H2 = "EDGE-sym is not behind Stukel's joint score",
         H3 = "EDGE-poly3 is not far behind GiViTI on asymmetric links",
         H4 = "EDGE-poly3 leads or ties the partition tests at $G=10$",
         H4.ruleG = "the same at the default partition",
         H5 = "EDGE-sym is not behind Stukel's one-parameter score",
         H5.skew = "the same under a skewed covariate")
verd <- c(holds = "\\textbf{holds}", against = "\\textbf{against}", reported = "reported")

## The declared verdict is a statement about the declared STATISTIC against its threshold. It says
## nothing about how many individual cells go each way, and for H1 the two readings differ sharply:
## the hypothesis fails on a design condition while the test still leads or ties in 26 of 32 cells.
## Both are reported, in adjacent columns, so neither can stand in for the other.
M <- fread(edge_battery("analysis", "mcnemar.csv"))[form == "unit" & arm == "Grule"]
lead_tie <- function(h) {
  x <- M[hypothesis == h]
  if (!nrow(x)) return(NA_character_)
  d <- x$diff
  if (startsWith(h, "H5")) d <- -d          # H5 declares Stukel - EDGE, so the sign is reversed
  sprintf("%d of %d", sum(d > -0.01), nrow(x))
}
cnt <- vapply(H$hypothesis, function(h) {
  ## H4 counts the detectable, non-saturated scenarios, not all 29: 19 at G = 10 and 21 at the rule G
  ## (hypotheses.csv; the paper-3 referee round found "of 29" printed here in paper 2)
  if (startsWith(h, "H4")) sprintf("%s of %s", format(H$statistic[H$hypothesis == h][1]),
                                   if (h == "H4") "19" else "21")
  else lead_tie(h)
}, character(1))

rows <- sprintf("%s & %s & %d & %s & %s & %s & %s \\\\",
                gsub("\\.", "-", H$hypothesis), lab[H$hypothesis], H$cells,
                sprintf("%.3f", H$statistic),
                ifelse(H$threshold == "" | is.na(H$threshold), "--", H$threshold),
                verd[H$final_verdict], ifelse(is.na(cnt), "--", cnt))
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{The seven hypotheses declared before the study was run, and their outcomes. The",
     "statistic and the threshold are those fixed in the pre-declaration; no threshold was moved",
     "after the results were seen. The last column counts the scenarios in which EDGE",
     "leads or ties, which is a different question from the declared verdict and answers it",
     "differently: H1 fails on a condition about one design while the test leads or ties in $26$",
     "of $32$ scenarios. Section 5.2 of the main text reports the results as findings.}",
     "\\label{tab:hypotheses}",
     "\\begin{tabular}{lp{0.40\\textwidth}rrrll}\\hline",
     "Hypothesis & Statement & Scenarios & Statistic & Threshold & Outcome & Leads or ties \\\\ \\hline",
     rows, "\\hline\\end{tabular}\\end{table}"), "tab_hypotheses.tex")

## ---------------------------------------------------------------------------------------------
## Table 2: size of the EDGE variants over the 131 null cells
## ---------------------------------------------------------------------------------------------
G <- fread(edge_battery("analysis", "rule_A_gate.csv"))
## "failed" is the DECLARED gate, which is checked at both 0.05 and 0.01: counting only the 5%
## column understates it (the Stukel basis fails three cells, all of them at the 1% level).
S <- G[, .(cells = .N, min = min(size_05), median = median(size_05), max = max(size_05),
           fail05 = sum(pass_05 == FALSE), fail01 = sum(pass_01 == FALSE),
           failed = sum(pass == FALSE)), by = .(basis, form)]
S[, form := ifelse(form == "u", "unit", "score")]
setorder(S, basis, form)
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{Size over the $131$ null scenarios of the study, for each basis and weighting at the",
     "default partition. The three columns on the left summarise the realised size at the nominal $5\\%$",
     "level. The level check is made at both $5\\%$ and $1\\%$, and the last three columns",
     "count the scenarios outside a band of three Monte Carlo standard errors at each level and at",
     "either. No form passed in every scenario; under the fallback clause of the rule fixed in advance the",
     "unit form, the specified default, was kept. Power is size-adjusted throughout.}",
     "\\label{tab:size}",
     "\\begin{tabular}{llrrrrrrr}\\hline",
     " & & & \\multicolumn{3}{c}{size at $5\\%$} & \\multicolumn{3}{c}{scenarios failing the gate} \\\\",
     "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}",
     "Basis & Weighting & Scenarios & Min & Median & Max & at $5\\%$ & at $1\\%$ & either \\\\ \\hline",
     sprintf("%s & %s & %d & %s & %s & %s & %d & %d & %d \\\\", S$basis, S$form, S$cells,
             f3(S$min), f3(S$median), f3(S$max), S$fail05, S$fail01, S$failed),
     "\\hline\\end{tabular}\\end{table}"), "tab_size.tex")

## ---------------------------------------------------------------------------------------------
## Table 3: mean power by departure family, the six families of the pre-declaration
## ---------------------------------------------------------------------------------------------
F <- fread(edge_battery("analysis", "rule_A_families.csv"))[form == "u" & family != "macro"]
W <- dcast(F, family + family_name + cells ~ basis, value.var = "mean_power")
setorder(W, family)
MAC <- fread(edge_battery("analysis", "rule_A_families.csv"))[form == "u" & family == "macro"]
mac <- setNames(MAC$mean_power, MAC$basis)
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{Mean size-adjusted power of the three bases, unit weighting at the default partition, over",
     "the six departure families. The last row is the macro-average of",
     "the six family means, on which the cubic basis is the default.}",
     "\\label{tab:families}",
     "\\begin{tabular}{llrrrr}\\hline",
     "Family & Departure & Scenarios & EDGE-poly3 & EDGE-stk & EDGE-sym \\\\ \\hline",
     sprintf("%s & %s & %d & %s & %s & %s \\\\", W$family, W$family_name, W$cells,
             f3(W$poly3), f3(W$stk), f3(W$sym)),
     "\\hline",
     sprintf("\\multicolumn{3}{l}{macro-average of the six} & %s & %s & %s \\\\",
             f3(mac[["poly3"]]), f3(mac[["stk"]]), f3(mac[["sym"]])),
     "\\hline\\end{tabular}\\end{table}"), "tab_families.tex")

## ---------------------------------------------------------------------------------------------
## Table 4: contamination -- the spine of paper 2 (block 9, C1, logistic truth = false alarms)
## ---------------------------------------------------------------------------------------------
B9 <- fread(edge_battery("9", "analysis", "cell_test.csv"))
## The Hosmer-Lemeshow statistic is in this study's own deposit at both partitions and belongs in the table:
## it is the other grouped test, and at the cells where EDGE fails it is the better protected.
## The clean cells are included too, so the section reads its own null off its own table.
keep <- c("EDGE.poly3.u.Grule", "EDGE.sym.u.Grule", "HL.Grule", "HL.G10",
          "Stk.joint", "Stk.sym1", "GiViTI", "Cubic.LR")
C1 <- B9[truth == "logit" & corruption %in% c("C1", "clean") & test %in% keep]
T4 <- dcast(C1, n + k ~ test, value.var = "rejection")
## le Cessie (block 8L Part B, sha256 0c1a1731) on block 9's own data sets at n = 1000; its n x n kernel was not
## run at n = 5000. EDGE recomputed on 8L's data sets must equal block 9's EDGE exactly.
L8B <- fread(edge_battery("8L", "_summary.csv"))[part == "B" & role %in% c("corrupted", "clean"),
                                                   .(n, k, LC = lecessie_rejection, E8L = edge_rejection)]
T4 <- merge(T4, L8B, by = c("n", "k"), all.x = TRUE)
if (!isTRUE(all.equal(T4[!is.na(E8L), E8L], T4[!is.na(E8L), EDGE.poly3.u.Grule])))
  stop("Table 4: EDGE on block 8L's data sets differs from block 9")
setorder(T4, n, k)
## what each column costs (block T, one core, battery/T/timing.csv), so that protection is read beside its price.
## Le Cessie's test is timed in its published O(n^3) form (bench_time_T_lecessie_n3.R), as the author asked
## (2026-10-01); the statistic is the one every table reports.
TT  <- fread(edge_battery("T", "timing.csv"))
TT3 <- fread(edge_battery("T", "timing_lecessie_n3.csv"))
ms  <- function(tname, nn, tab = TT) { v <- tab[tab$test == tname & tab$n == nn, ]$median_sec
  if (!length(v) || is.na(v[1])) "--" else if (v[1] < 1) sprintf("%.0f", 1000 * v[1]) else
    formatC(1000 * v[1], format = "d", big.mark = ",") }
lc_ms <- function(nn) ms("le Cessie, O(n^3) form", nn, TT3)
time_row <- function(nn) sprintf("\\multicolumn{2}{l}{ms, $n=%d$} & %s & -- & %s & %s & %s & %s & -- & %s & %s \\\\", nn,
  ms("EDGE", nn), ms("Hosmer-Lemeshow", nn), ms("Hosmer-Lemeshow", nn), lc_ms(nn), ms("Stukel joint", nn),
  ms("GiViTI", nn), ms("cubic calib. LR", nn))
wr(c("\\begin{table}[tp]\\centering\\small",
     "\\caption{False-alarm rate at the nominal $5\\%$ level when $k$ records carry a corrupted",
     "covariate and the model is correct for every other record; the $k=0$ rows are sizes on the same",
     "design. $1000$ replicates a scenario; the nominal standard error is $0.007$. EDGE is in its",
     "unit form at the default partition, with the cubic (default) and symmetric bases, and the",
     "Hosmer--Lemeshow statistic is shown at both partitions. The pooled tests combine each record's",
     "residual with those of other records --- in groups of predicted risk, or for le Cessie's test",
     "over neighbourhoods in covariate space. Every test built on each record's own",
     "linear predictor rises with $k$. The last two rows give the median time a data set in milliseconds on",
     "one core (the symmetric forms were not timed). Le Cessie's test is protected but holds an $n\\times n$",
     "kernel and grows with the cube of $n$; it was run with corrupted records at $n=1000$ only.}",
     "\\label{tab:contamination}",
     "\\resizebox{\\textwidth}{!}{\\begin{tabular}{rrrrrrrrrrr}\\hline",
     " & & \\multicolumn{5}{c}{pooled} & \\multicolumn{4}{c}{record-level} \\\\",
     "\\cmidrule(lr){3-7}\\cmidrule(lr){8-11}",
     "$n$ & $k$ & EDGE & EDGE-sym & HL ($G$ rule) & HL ($G{=}10$) & le Cessie & Stukel joint & Stukel sym & GiViTI & cubic LR \\\\ \\hline",
     sprintf("%d & %d & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", T4$n, T4$k,
             f3(T4$EDGE.poly3.u.Grule), f3(T4$EDGE.sym.u.Grule), f3(T4$HL.Grule), f3(T4$HL.G10),
             f3(T4$LC), f3(T4$Stk.joint),
             f3(T4$Stk.sym1), f3(T4$GiViTI), f3(T4$Cubic.LR)),
     "\\hline", time_row(1000L), time_row(5000L),
     "\\hline\\end{tabular}}\\end{table}"), "tab_contamination.tex")

## ---------------------------------------------------------------------------------------------
## Table 5: the two resampling rivals (block 8), size-adjusted, paired with EDGE
## ---------------------------------------------------------------------------------------------
P <- fread(edge_battery("8", "_paired.csv"))
P <- P[comparator == "EDGE.poly3.u.Grule" & in_holm == TRUE, ]
R <- dcast(P, cell + n ~ rival, value.var = c("rival_rejection", "comparator_rejection", "holm_reject_05"))
## le Cessie (block 8L Part A, sha256 0c1a1731) on the projection test's replicates 1-500; EDGE on those
## replicates must equal block 8's EDGE exactly
L8A <- fread(edge_battery("8L", "_paired.csv"))[part == "A" & role == "alternative",
                                                  .(cell, lc = lecessie_rejection, lc_edge = edge_rejection, lc_star = holm_reject_05)]
R <- merge(R, L8A, by = "cell", all.x = TRUE)
if (!isTRUE(all.equal(R$comparator_rejection_proj, R$lc_edge))) stop("Table 5: EDGE on block 8L's replicates differs from block 8")
## paper 3: readable rows grouped by family, the best test of each row in bold, and a closing block that states
## what the table shows -- power at least equal to the three tests built without a partition, at a thousandth of
## their cost -- from the table's own rows and block T's timings.
MEM3 <- unique(fread(edge_battery("analysis", "rule_A_membership.csv"))[, .(cell, family_name)], by = "cell")
R <- merge(R, MEM3, by = "cell", all.x = TRUE)
stopifnot(!anyNA(R$family_name))
R[, dep := sub("_n[0-9]+$", "", sub("_(auc|e12|base)_n[0-9]+$", "", cell))]
R[, design := ifelse(grepl("_auc_", cell), "high AUC", ifelse(grepl("_e12_", cell), "$12\\%$ events", ifelse(grepl("_base_", cell), "base", "one covariate")))]
DEP <- c(cauchit = "cauchit link", t4 = "$t_4$ link", loglog = "log--log link", stk_asym = "Stukel, asymmetric",
         stk_long = "Stukel, long tails", stk_short = "Stukel, short tails",
         crossover = "curve crossing the diagonal twice")
stopifnot(all(R$dep %in% names(DEP)))
FAM <- c("symmetric tails" = "Symmetric tails", "asymmetric links" = "Asymmetric links",
         "Stukel symmetric family" = "Stukel symmetric family", "off-index" = "A crossing curve")
R[, fam := factor(family_name, levels = names(FAM))]
R[, dord := match(dep, names(DEP))]
setorder(R, fam, dord, n)
R[, `:=`(e = comparator_rejection_proj, l = lc, pj = rival_rejection_proj, b = rival_rejection_bagoft)]
cellfmt <- function(x, best, star) paste0(ifelse(best, "\\textbf{", ""), f3(x), ifelse(best, "}", ""),
                                          ifelse(isTRUE_v(star), "$^\\star$", ""))
isTRUE_v <- function(x) !is.na(x) & x == TRUE
body <- unlist(lapply(levels(R$fam), function(f) {
  X <- R[fam == f]
  mx <- pmax(X$e, X$l, X$pj, X$b, na.rm = TRUE)
  c(sprintf("\\multicolumn{7}{l}{\\emph{%s}} \\\\", FAM[[f]]),
    sprintf("\\quad %s & %s & %d & %s & %s & %s & %s \\\\", DEP[X$dep], X$design, X$n,
            cellfmt(X$e, X$e == mx, NA), cellfmt(X$l, X$l == mx, X$lc_star),
            cellfmt(X$pj, X$pj == mx, X$holm_reject_05_proj), cellfmt(X$b, X$b == mx, X$holm_reject_05_bagoft)))
}))
## EDGE ahead / behind after the Holm correction, read from the starred differences
ab <- function(r, star) sprintf("%d / %d", sum(isTRUE_v(star) & R$e > r), sum(isTRUE_v(star) & R$e < r))
TT3 <- rbind(fread(edge_battery("T", "timing.csv"))[n == 1000 & test != "le Cessie", .(test, median_sec)],
             fread(edge_battery("T", "timing_lecessie_n3.csv"))[n == 1000, .(test = "le Cessie", median_sec)])
tfmt <- function(tn) { v <- TT3[test == tn, median_sec]
  if (v < 1) sprintf("%.0f ms", 1000 * v) else if (v < 10) sprintf("%.1f s", v) else sprintf("%.0f s", v) }
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{EDGE against the three tests built without a partition, on the same data sets:",
     "le Cessie's smoothed test \\citep{le1991goodness}, the Liu projection test \\citep{liu2024comprehensive} and",
     "BAGofT \\citep{BAGofT2019}, over the $20$ link and tail scenarios at $n\\ge500$. Rejection rates at the",
     "nominal $5\\%$ level on replicates $1$--$500$ ($1$--$200$ for BAGofT); the highest rate in each row is in",
     "bold, and $\\star$ marks a paired difference from EDGE significant after the Holm correction.",
     "EDGE is at its default partition. The closing rows summarise the table: after the Holm correction EDGE",
     "is ahead of le Cessie's test and BAGofT in $14$ of the $20$ scenarios and of the projection test in",
     "$7$, and each rival is ahead in one, the crossing curve; the mean is of raw rates, while the text and",
     "Figure~\\ref{fig:cost} give size-adjusted power. The rivals lead on the families they are built for,",
     "interactions and rough misfit (Section~\\ref{sec:rivals}, Table~\\ref{tab:whennot}). Time is the median",
     "for one data set at $n=1000$ on one core, le Cessie's test in its published form.}",
     "\\label{tab:rivals}",
     "\\begin{tabular}{llrrrrr}\\hline",
     "Departure & Design & $n$ & EDGE & le Cessie & Liu projection & BAGofT \\\\ \\hline",
     body,
     "\\hline",
     sprintf("\\multicolumn{3}{l}{Mean rejection rate (raw)} & %s & %s & %s & %s \\\\",
             f3(mean(R$e)), f3(mean(R$l)), f3(mean(R$pj)), f3(mean(R$b))),
     sprintf("\\multicolumn{3}{l}{EDGE ahead / behind after Holm} & -- & %s & %s & %s \\\\",
             ab(R$l, R$lc_star), ab(R$pj, R$holm_reject_05_proj), ab(R$b, R$holm_reject_05_bagoft)),
     sprintf("\\multicolumn{3}{l}{Time a data set, $n=1000$} & %s & %s & %s & %s \\\\",
             tfmt("EDGE"), tfmt("le Cessie"), tfmt("Liu projection"), tfmt("BAGofT")),
     "\\hline\\end{tabular}\\end{table}"), "tab_rivals.tex")

## ---------------------------------------------------------------------------------------------
## Table 6: does a robust estimator solve the problem instead? (block 9d, sha256 799a85aa)
## ---------------------------------------------------------------------------------------------
R9d <- fread(edge_battery("9d", "analysis_9d.csv"))[truth == "logit"]
## Every form of EDGE is listed, not only the one that survives the robust fit: three of the four
## break their level under it, which is the point of the last paragraph of the section.
shown <- c("EDGE.poly3.u.Grule", "EDGE.sym.u.Grule", "EDGE.poly3.sc.Grule", "EDGE.sym.sc.Grule",
           "HL.Grule", "HLF.Grule", "GiViTI", "Cubic.LR", "Stk.joint", "Stk.sym1")
nice  <- c(EDGE.poly3.u.Grule = "EDGE (default)", EDGE.sym.u.Grule = "EDGE, symmetric basis",
           EDGE.poly3.sc.Grule = "EDGE, score-weighted", EDGE.sym.sc.Grule = "EDGE, symmetric score-weighted",
           HL.Grule = "Hosmer--Lemeshow",
           HLF.Grule = "HL, Farrington-corrected", GiViTI = "GiViTI belt", Cubic.LR = "cubic calibration LR",
           Stk.joint = "Stukel joint score", Stk.sym1 = "Stukel one-parameter")
T6 <- dcast(R9d[test %in% shown], test ~ est + k, value.var = "rate")
T6[, test := factor(test, levels = shown)]; setorder(T6, test)
row6 <- function(i) {
  z <- T6[i]
  sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", nice[as.character(z$test)],
          f3(z$ML_0), f3(z$robust_0), f3(z$ML_10), f3(z$robust_10), f3(z$ML_25), f3(z$robust_25))
}
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{Rejection rate at the nominal $5\\%$ level under a correct logistic model with $k$",
     "corrupted covariates, $n=1000$, when the model is fitted by maximum likelihood (ML) and by the",
     "Huber-type quasi-likelihood estimator of \\citet{cantoni2001robust} --- \\textsf{robustbase}'s",
     "\\texttt{Mqle} at tuning constant $c=1.345$, with no weight on the covariates --- on the same",
     "data. This study was run separately from that of",
     "the main contamination table and on its own data sets, so its maximum-likelihood columns differ",
     "from the corresponding rows of that table by Monte Carlo error alone. The robust fit converged in",
     "all $1000$ replicates of every configuration and no replicate was excluded. With $k=0$ the",
     "rates are sizes. The robust fit leaves every test's false-alarm rate where it was or raises",
     "it, and it breaks the level of five of the ten tests the gate covers, even without contamination",
     "--- Stukel's two score tests and three of the four forms of EDGE --- because their",
     "reference is built on the maximum-likelihood information matrix. Only the default unit form of",
     "EDGE, the two Hosmer--Lemeshow statistics, the belt and the cubic test are read",
     "against a reference that still holds. Eleven tests were run and ten are shown: Stukel's",
     "likelihood-ratio refit is omitted because a robust fit has no deviance, so under the robust arm",
     "that refit compares a maximum-likelihood augmented fit against a robust base and is not a",
     "likelihood-ratio test. $1000$ replicates; the nominal standard error is $0.007$.}",
     "\\label{tab:robust}",
     "\\begin{tabular}{lrrrrrr}\\hline",
     " & \\multicolumn{2}{c}{$k=0$} & \\multicolumn{2}{c}{$k=10$} & \\multicolumn{2}{c}{$k=25$} \\\\",
     "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
     "Test & ML & robust & ML & robust & ML & robust \\\\ \\hline",
     vapply(seq_len(nrow(T6)), row6, character(1)),
     "\\hline\\end{tabular}\\end{table}"), "tab_robust.tex")

## ---------------------------------------------------------------------------------------------
## Table 7: the resampling tests under the same corruption (block 9b, sha256 72d703a9)
## ---------------------------------------------------------------------------------------------
R9b <- fread(edge_battery("9b", "analysis", "rates_9b.csv"))
S9b <- fread(edge_battery("9b", "analysis", "seconds_9b.csv"))
if (any(R9b$reps != 100L)) stop("block 9b: a cell does not have its 100 replicates")
shown9b <- c("EDGE.poly3.u.Grule", "EDGE.sym.u.Grule", "HL.G10", "BAGofT",
             "proj", "Stk.joint", "Stk.sym1", "Stk.LR", "GiViTI", "Cubic.LR")
nice9b <- c(EDGE.poly3.u.Grule = "EDGE (default)", EDGE.sym.u.Grule = "EDGE, symmetric basis",
            HL.G10 = "Hosmer--Lemeshow, $G=10$", BAGofT = "BAGofT", proj = "Liu projection test",
            Stk.joint = "Stukel joint score", Stk.sym1 = "Stukel one-parameter", Stk.LR = "Stukel LR refit",
            GiViTI = "GiViTI belt", Cubic.LR = "cubic calibration LR")
T7 <- dcast(R9b[test %in% shown9b], test ~ cell, value.var = "rate")
## le Cessie on the same 300 data sets (block 9bL, sha256 e3c89be7); EDGE there must equal block 9b's EDGE exactly
L9 <- fread(edge_battery("9bL", "_summary.csv"))
L9 <- L9[match(c("clean", "x4", "x8"), cell)]
if (!isTRUE(all.equal(L9$edge_rejection, unlist(T7[test == "EDGE.poly3.u.Grule", .(clean, x4, x8)]), check.attributes = FALSE)))
  stop("Table 6: EDGE on block 9bL's data sets differs from block 9b")
T7 <- rbind(T7, data.table(test = "lecessie", clean = L9$lecessie_rejection[1], x4 = L9$lecessie_rejection[2],
                           x8 = L9$lecessie_rejection[3]))
nice9b <- c(nice9b, lecessie = "le Cessie--van Houwelingen")
shown9b <- c(shown9b[1:4], "lecessie", shown9b[5:10])
T7 <- T7[match(shown9b, test)]
row7 <- function(i) sprintf("\\quad %s & %s & %s & %s \\\\", nice9b[[T7$test[i]]], f3(T7$clean[i]), f3(T7$x4[i]), f3(T7$x8[i]))
sec7 <- S9b[, .(sec = median(median_sec)), by = test]
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{False-alarm rate at the nominal $5\\%$ level when three records of $n=500$ have their",
     "covariate multiplied by four or by eight and the model is correct for every other record. $100$",
     "replicates a scenario, every test on the same data sets; the nominal standard error is $0.022$. EDGE is",
     "in its unit form at the default partition, which is twenty groups here. The",
     "Liu projection test and BAGofT reject at $p<0.05$, the others at $p\\le0.05$. BAGofT pools the residuals",
     "within groups it learns from the data and le Cessie's test over neighbourhoods in covariate space, so",
     "both are listed with the pooled tests; le Cessie's test was run on the same data sets.",
     sprintf(" The Liu projection test took a median of $%.0f$ seconds a data set and BAGofT $%.0f$.}",
             sec7[test == "projection", sec], sec7[test == "BAGofT", sec]),
     "\\label{tab:rivals-contam}",
     "\\begin{tabular}{lrrr}\\hline",
     "Test & clean & $\\times4$ & $\\times8$ \\\\ \\hline",
     "\\multicolumn{4}{l}{\\emph{pooled}} \\\\",
     vapply(1:5, row7, character(1)),
     "\\multicolumn{4}{l}{\\emph{record-level}} \\\\",
     vapply(6:11, row7, character(1)),
     "\\hline\\end{tabular}\\end{table}"), "tab_rivals_contam.tex")

## ---------------------------------------------------------------------------------------------
## Table 8: the robust setting EDGE-FR (blocks 9R2, sha256 9943e551, and 9R3, sha256 421217c0)
## ---------------------------------------------------------------------------------------------
FR <- fread(edge_battery("analysis", "edgefr.csv"))
if (any(FR$missing > 0)) stop("Table 8: a replicate has no p-value")
## false alarms and level are rejection rates at the nominal level; power is size-adjusted, as everywhere else
FR[, shown := ifelse(grepl("^power", role), adjusted, rejection)]
FW <- dcast(FR, role + cell + n + k ~ variant, value.var = "shown")
FRP <- fread(edge_battery("analysis", "edgefr_paired.csv"))[kind == "power (size-adjusted)"]
fr_row <- function(z, label) sprintf("\\quad %s & %s & %s & %s & %s & %s \\\\", label,
                                     f3(z$V0), f3(z$HYB05), f3(z$HYB10), f3(z$HYB20), f3(z$G10))
alarm <- FW[role == "false alarm"][order(n, k)]
lev   <- FW[role == "level" & cell %in% c("logit_clean_n1000", "logit_clean_n5000")][order(n)]
fresh <- FW[role == "level" & cell %in% c("fr_clean_n1000", "fr_clean_n5000")]
pw    <- FW[cell %in% c("battery mean (6)", "cloglog_clean_n1000", "probit_clean_n5000")]
pw    <- pw[match(c("battery mean (6)", "cloglog_clean_n1000", "probit_clean_n5000"), cell)]
frange <- range(unlist(fresh[, c("V0", "HYB05", "HYB10", "HYB20", "G10"), with = FALSE]))

wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{The tail-pooled setting EDGE-FR, at three tail fractions, against the default partition and",
     "against ten groups. Its default column repeats the EDGE column of",
     "Table S5 for the eight scenarios the two tables share, so that the",
     "comparison can be read in place. EDGE-FR$(c)$ puts the outer $5c$ of the risk scale at each end into one group",
     "and leaves the middle at the default partition, which is used throughout otherwise. $1000$ replicates",
     "a scenario, the same data sets for every setting. The false alarms are rejection rates at the",
     "nominal $5\\%$ level under a correct logistic model with $k$ records carrying a corrupted",
     "covariate; the rows at $k \\ge 20$ for $n=1000$ and $k \\ge 100$ for $n=5000$ were run on fresh",
     "seeds, whose uncontaminated scenarios read between",
     sprintf("$%s$ and $%s$ for every setting. The power rows are size-adjusted, each setting at its own", f3(frange[1]), f3(frange[2])),
     "critical value read from the matched null scenario, and the first is the mean over six link and",
     "tail scenarios at $n=1000$ (cauchit, log--log, three Stukel designs, $t_4$); the paired readings behind them are in the",
     "archive. Widening the tails costs power without adding protection in proportion, and the",
     "protection ends in every setting once the corrupted records pass a few tens.}",
     "\\label{tab:edgefr}",
     "\\begin{tabular}{lrrrrr}\\hline",
     " & & \\multicolumn{3}{c}{EDGE-FR} & \\\\",
     "\\cmidrule(lr){3-5}",
     "Scenario & default & $c=1\\%$ & $c=2\\%$ & $c=4\\%$ & ten groups \\\\ \\hline",
     "\\multicolumn{6}{l}{\\emph{false alarms, $k$ records with a corrupted covariate}} \\\\",
     vapply(seq_len(nrow(alarm)), function(i)
       fr_row(alarm[i], sprintf("$n=%d$, $k=%d$", alarm$n[i], alarm$k[i])), character(1)),
     "\\multicolumn{6}{l}{\\emph{level, no corruption}} \\\\",
     vapply(seq_len(nrow(lev)), function(i)
       fr_row(lev[i], sprintf("$n=%d$", lev$n[i])), character(1)),
     "\\multicolumn{6}{l}{\\emph{power against a real misfit}} \\\\",
     fr_row(pw[1], "six link and tail scenarios, mean"),
     fr_row(pw[2], "complementary log--log, $n=1000$"),
     fr_row(pw[3], "probit, $n=5000$"),
     "\\hline\\end{tabular}\\end{table}"), "tab_edgefr.tex")

## ---------------------------------------------------------------------------------------------
## Table 5 of paper 3: false alarms under the two covariate errors, by partition, basis and grouped rival
## (battery/analysis/corruption_partitions.csv, from analyse_corruption_partitions.R, which first reproduces the
## published cells of blocks 9 and C1b exactly). Paper 2's Table 9 showed the cubic basis at the rule G only; the
## referee round asked for the symmetric basis, ten groups and the grouped rivals the battery also ran.
## ---------------------------------------------------------------------------------------------
CP <- fread(edge_battery("analysis", "corruption_partitions.csv"))
cp_tests <- c("EDGE.poly3.u.Grule", "EDGE.poly3.u.G10", "EDGE.sym.u.Grule", "EDGE.sym.u.G10", "HL.G10", "HL_w", "PH", "Tsiatis", "Stk.joint")
cp_row <- function(err, nn, kk) {
  z <- CP[error == err & n == nn & k == kk]
  ## a rate above 0.10, where Section 6 reads a test as no longer holding, is shaded
  v <- vapply(cp_tests, function(t) { x <- z$false_alarm[z$test == t][1]
    if (!is.na(x) && x > 0.10) paste0("\\cellcolor{black!12}", f3(x)) else f3(x) }, character(1))
  sprintf("%d & %d & %s \\\\", nn, kk, paste(v, collapse = " & "))
}
cp_block <- function(err) unlist(lapply(c(1000L, 5000L), function(nn)
  vapply(sort(unique(CP$k[CP$n == nn & CP$error == err])), function(kk) cp_row(err, nn, kk), character(1))))
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{False-alarm rate at the nominal $5\\%$ level under a correct logistic model when $k$ records",
     "carry a corrupted covariate, for the two errors of Section~\\ref{sec:contamination-design}: the",
     "covariate multiplied by four (exaggeration) and by $-4$ (sign error). EDGE is shown with the cubic",
     "and symmetric bases at the default partition and at ten groups. HL-ew is the Hosmer--Lemeshow statistic with equal-width risk groups; PH is the",
     "Pigeon--Heyse test and Tsiatis the score test on indicators of ten $K$-means clusters in covariate space, all at ten groups.",
     "$1000$ replicates a scenario, every test on the same data sets; the Monte Carlo standard error of",
     "a rate of $0.05$ is $0.007$. A replicate with no $p$-value counts as a non-rejection. Shaded cells exceed",
     "$0.10$, the rate above which a test is read as no longer holding its level.}",
     "\\label{tab:severity}",
     "\\begin{adjustbox}{max width=\\textwidth}",
     "\\begin{tabular}{rrrrrrrrrrr}\\hline",
     " & & \\multicolumn{4}{c}{EDGE} & \\multicolumn{4}{c}{grouped rivals} & \\\\",
     "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
     "$n$ & $k$ & cubic, default & cubic, $G{=}10$ & symmetric, default & symmetric, $G{=}10$ & HL, $G{=}10$ & HL-ew & PH & Tsiatis & Stukel joint \\\\ \\hline",
     "\\multicolumn{11}{l}{\\emph{Exaggeration}, $x\\mapsto4x$} \\\\",
     cp_block("exaggeration"),
     "\\hline",
     "\\multicolumn{11}{l}{\\emph{Sign error}, $x\\mapsto-4x$} \\\\",
     cp_block("sign error"),
     "\\hline\\end{tabular}\\end{adjustbox}\\end{table}"), "tab_severity.tex")

## ---------------------------------------------------------------------------------------------
## Table S8: external validation (block EXT, PREDECLARATION_blockEXT_external_validation.md, sha256 754cf309):
## EDGE in external mode against the usual external-validation tests
## ---------------------------------------------------------------------------------------------
EX <- fread(edge_battery("EXT", "_summary.csv"))
ext_tests <- c("EDGE_G10", "EDGE_auto", "cox", "spiegelhalter", "giviti", "HL", "stukel")
ext_row <- function(z, label, v = "raw") sprintf("%s & %s \\\\", label,
  paste(vapply(ext_tests, function(t) { x <- z[z$test == t, ][[v]]
    if (v == "raw" && z$role[1] == "corrupted" && x > 0.10) paste0("\\cellcolor{black!12}", f3(x)) else f3(x) }, ""), collapse = " & "))
FAMX <- c(large = "calibration in the large", slope = "calibration slope", ushape = "U-shaped term",
          thresh = "threshold term", inter = "interaction")
## EX is long, one row a test and scenario; the table has one row a scenario
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
     "intercept and slope, Spiegelhalter the $z$ test, GiViTI the belt in external mode, HL the",
     "Hosmer--Lemeshow statistic at ten groups and Stukel his two terms added to the frozen linear predictor.",
     "$1000$ replicates a scenario. Power is size-adjusted at the null with $n=500$; shaded false-alarm rates",
     "exceed $0.10$.}",
     "\\label{tab:external}",
     "\\resizebox{\\textwidth}{!}{\\begin{tabular}{lrrrrrrr}\\hline",
     "Scenario & EDGE, $G{=}10$ & EDGE, default & Cox & Spiegelhalter & GiViTI & HL & Stukel \\\\ \\hline",
     "\\multicolumn{8}{l}{\\emph{Level: rejection rate of a correct model}} \\\\",
     vapply(seq_len(nrow(ex_lev)), function(i) ext_row(EX[id == ex_lev$id[i]], sprintf("\\quad $n=%d$", ex_lev$n[i])), ""),
     "\\multicolumn{8}{l}{\\emph{Power, size-adjusted, $n=500$}} \\\\",
     vapply(seq_len(nrow(ex_alt)), function(i) ext_row(EX[id == ex_alt$id[i]],
       sprintf("\\quad %s, $C=%s$", FAMX[[ex_alt$fam[i]]], format(ex_alt$C[i])), "adj"), ""),
     sprintf("\\quad mean of the ten & %s \\\\", paste(f3(mean_adj[match(ext_tests, test), m]), collapse = " & ")),
     "\\multicolumn{8}{l}{\\emph{False alarms, correct model, $n=1000$}} \\\\",
     vapply(seq_len(nrow(ex_cor)), function(i) ext_row(EX[id == ex_cor$id[i]],
       sprintf("\\quad %s, $k=%d$", if (ex_cor$fam[i] == "exag") "exaggeration" else "sign error", ex_cor$k[i])), ""),
     sprintf("\\quad settings held, of eight & %s \\\\", paste(held[match(ext_tests, test), h], collapse = " & ")),
     "\\hline\\end{tabular}}\\end{table}"), "tab_external.tex")

## ---------------------------------------------------------------------------------------------
## Table S10: mean difference in size-adjusted power by departure family, with Monte Carlo intervals
## (analyse_paired_families.R; referee M6.1, 2026-10-01)
## ---------------------------------------------------------------------------------------------
PF <- fread(edge_battery("analysis", "paired_families.csv"))
fam_order <- c("symmetric tails", "asymmetric links", "Stukel symmetric family", "omitted terms", "rough misfit",
               "off-index", "all families")
cmp_order <- c("Hosmer-Lemeshow, ten groups", "Stukel joint score", "Stukel one-parameter", "cubic calibration LR",
               "GiViTI belt", "EDGE at ten groups")
cell_pf <- function(fm, cm) { z <- PF[family_name == fm & comparator == cm]
  if (!nrow(z)) "--" else sprintf("$%+.3f$ {\\footnotesize(%.3f, %.3f)}", z$mean_diff, z$lo, z$hi) }
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{Mean difference in size-adjusted power, EDGE (cubic basis, unit form, default partition) minus",
     "each comparator, by departure family, with a Monte Carlo $95\\%$ interval in brackets. The interval",
     "ignores the positive correlation between two tests read on the same replicates, which makes it",
     "conservative, and the Monte Carlo error of the critical values used for size adjustment. The last",
     "column is the price of ten groups.}",
     "\\label{tab:pairedfam}",
     "\\resizebox{\\textwidth}{!}{\\begin{tabular}{lrllllll}\\hline",
     "Family & Scenarios & HL, ten groups & Stukel joint & Stukel one-par. & cubic LR & GiViTI & EDGE, ten groups \\\\ \\hline",
     vapply(fam_order, function(fm) sprintf("%s & %d & %s \\\\", fm, PF[family_name == fm][1]$cells,
            paste(vapply(cmp_order, function(cm) cell_pf(fm, cm), ""), collapse = " & ")), ""),
     "\\hline\\end{tabular}}\\end{table}"), "tab_pairedfam.tex")

## ---------------------------------------------------------------------------------------------
## Table S11: block TW, what protects a test (PREDECLARATION_blockTW_mechanism.md and addenda 1-2)
## ---------------------------------------------------------------------------------------------
TWS <- fread(edge_battery("TW", "_summary.csv"))
tw_tests <- c("EDGE.default", "EDGE.G10", "HL.G10", "HH.HLnp", "HH.X2np", "TWIN.score", "EDGE.Gn", "SPZ", "Stk.joint",
              "Cubic.LR", "GiViTI")
tw_cells <- c("TWA_logit_clean_k00_n1000", "TWA_logit_C1_k01_n1000", "TWA_logit_C1_k10_n1000", "TWA_logit_C1b_k01_n1000",
              "TWA_logit_C1b_k10_n1000", "TWA_logit_clean_k00_n5000", "TWA_logit_C1_k50_n5000", "TWA_logit_C1b_k25_n5000",
              "TWB_logit_clean_k00_n2000", "TWB_logit_C1_k25_n2000", "TWB_logit_C1b_k10_n2000",
              "TWA_cloglog_clean_k00_n1000", "TWA_tail0.6_clean_k00_n1000", "TWA_logit_C2_k10_n1000",
              "TWA_tail0.8_clean_k00_n5000", "TWB_cloglog_clean_k00_n2000", "TWB_inter_clean_k00_n2000")
tw_lab <- function(cl) { z <- TWS[cell == cl][1]
  what <- switch(z$corruption, clean = if (z$truth == "logit") "correct model" else
                   switch(sub("[0-9.]+$", "", z$truth), cloglog = "cloglog truth", inter = "omitted interaction",
                          tail = sprintf("top 2\\%%, risk %s of prediction", sub("tail", "", z$truth))),
                 C1 = sprintf("$\\times4$, $k=%d$", z$k), C1b = sprintf("$\\times-4$, $k=%d$", z$k),
                 C2 = sprintf("C2, $k=%d$", z$k))
  sprintf("%s & %d & %s", if (z$part == "TWB") "twelve covariates" else "two covariates", z$n, what) }
tw_val <- function(cl, t) { z <- TWS[cell == cl & test == t]
  if (!nrow(z) || (t == "EDGE.Gn" && z$n > 1000)) "--" else {
    v <- if ((TWS[cell == cl][1]$corruption == "clean" && TWS[cell == cl][1]$truth != "logit") || TWS[cell == cl][1]$corruption == "C2") z$adj else z$raw
    if (TWS[cell == cl][1]$corruption %in% c("C1", "C1b") && v > 0.10) paste0("\\cellcolor{black!12}", f3(v)) else f3(v) } }
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{What protects a test. Rejection rate at $5\\%$ on the same data sets, $1000$ replicates a row:",
     "false alarms under a correct model, with $k$ records carrying the covariate multiplied by $4$ or $-4$",
     "(shaded above $0.10$), and size-adjusted power under misfit (last six rows). EDGE at the default",
     "partition and at ten groups; HL, the Hosmer--Lemeshow statistic at ten groups; HH, the weighted grouped",
     "tests of Hosmer and Hjort with their general-purpose weight; Twin, the ungrouped score test for the",
     "cubic directions of EDGE in the predicted risk; EDGE$_n$, EDGE with one record a group ($n\\le1000$);",
     "Spz, Spiegelhalter's $z$, whose raw rate is $0$ throughout and which is shown at its own realised level.}",
     "\\label{tab:tw}",
     "\\resizebox{\\textwidth}{!}{\\begin{tabular}{llrrrrrrrrrrrrr}\\hline",
     "Design & $n$ & Setting & EDGE & EDGE$_{10}$ & HL & HH$_{\\mathrm{HL}}$ & HH$_{X^2}$ & Twin & EDGE$_n$ & Spz & Stukel & cubic & GiViTI \\\\ \\hline",
     vapply(tw_cells, function(cl) sprintf("%s & %s \\\\", tw_lab(cl),
            paste(vapply(tw_tests, function(t) if (t == "SPZ") { v <- TWS[cell == cl & test == t]$adj; if (TWS[cell == cl][1]$corruption %in% c("C1", "C1b") && v > 0.10) paste0("\\cellcolor{black!12}", f3(v)) else f3(v) } else tw_val(cl, t), ""),
                  collapse = " & ")), ""),
     "\\hline\\end{tabular}}\\end{table}"), "tab_tw.tex")

## ---------------------------------------------------------------------------------------------
## Table S12: block TW_ext, the same tests on frozen predictions (ADDENDUM3)
## ---------------------------------------------------------------------------------------------
TX <- fread(edge_battery("TW_ext", "_summary.csv"))
tx_val <- function(v, fam) if (fam != "null" && v > 0.10) paste0("\\cellcolor{black!12}", f3(v)) else f3(v)
wr(c("\\begin{table}[htbp]\\centering\\small",
     "\\caption{The same question on frozen predictions: rejection rate at $5\\%$, $1000$ replicates a row,",
     "on the corrupted data sets of the external-validation study ($n=1000$, five covariates), where nothing",
     "is refitted. Twin, the score test for the cubic directions of EDGE in the frozen prediction, on three",
     "degrees of freedom; Spz, Spiegelhalter's $z$; EDGE in external mode at ten groups and at the default",
     "partition. Shaded above $0.10$.}",
     "\\label{tab:twext}",
     "\\begin{tabular}{lrrrrr}\\hline",
     "Setting & Twin & Spz & EDGE$_{10}$ & EDGE \\\\ \\hline",
     vapply(seq_len(nrow(TX)), function(j) { z <- TX[j]
       lab <- if (z$fam == "null") "correct model" else sprintf("$\\times%s$, $k=%d$", if (z$fam == "exag") "4" else "-4", z$k)
       sprintf("%s & %s & %s & %s & %s \\\\", lab, tx_val(z$TWIN.ext, z$fam), tx_val(z$SPZ.ext, z$fam),
               tx_val(z$EDGE.G10.ext, z$fam), tx_val(z$EDGE.auto.ext, z$fam)) }, ""),
     "\\hline\\end{tabular}\\end{table}"), "tab_twext.tex")

cat("\nall tables written to", OUT, "\n")
