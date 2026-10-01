## analyse_census_extended.R -- the census of paper 3, extended to the comparators the battery ran but paper 2 did not
## report: Hosmer-Lemeshow and its Farrington form at G = 10, the Hosmer-Lemeshow statistic on equal-width risk groups (HL_w), Pigeon-Heyse
## (PH), Tsiatis, Xie and Pulkstenis-Robinson (PR). Read from the block summaries, the same size-adjusted power the
## declared analysis uses (each test at the critical value of its own matched null). Before anything new is computed,
## the tests already in _census_power_paper2.csv are recomputed the same way and must agree exactly, or the script stops.
##
##   Rscript analyse_census_extended.R      -> battery/analysis/census_extended.csv, census_extended_counts.csv
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
AN  <- edge_battery("analysis")
MARGIN <- 0.01                                                    # the paper's lead-or-tie margin

CEN <- fread(edge_path("results/analysis/_census_power_paper2.csv"))
cells <- unique(CEN$cell)
stopifnot(length(cells) == 158L)
MEM <- unique(fread(file.path(AN, "rule_A_membership.csv"))[, .(cell, family_name)], by = "cell")

S <- rbindlist(lapply(c("1a", "1b", "2", "3", "4", "5", "6", "7"), function(b) {
  x <- fread(edge_battery(b, "_summary.csv"), colClasses = list(character = c("subset", "null_cell")))
  x[, block := NULL][, block := b]
}), fill = TRUE)
S <- S[role == "alternative" & abs(alpha - 0.05) < 1e-9 & cell %in% cells & subset == "all"]
S <- unique(S, by = c("cell", "test"))                            # a cell reported by two blocks carries the same data

## 1. reproduce the published census exactly
chk <- merge(CEN, S[, .(cell, test, sap = size_adj_power)], by = c("cell", "test"))
if (nrow(chk) != nrow(CEN)) stop(sprintf("only %d of %d published census values found in the summaries", nrow(chk), nrow(CEN)))
dmax <- max(abs(chk$power - chk$sap), na.rm = TRUE)
nas  <- sum(is.na(chk$power) != is.na(chk$sap))
cat(sprintf("published census reproduced: %d values, max |difference| = %.2e, NA mismatches = %d\n", nrow(chk), dmax, nas))
if (dmax > 1e-12 || nas > 0) stop("the summaries do not reproduce the published census; nothing new is written")

## 2. the extension
NEW <- c("HL.G10", "HLF.G10", "HL_w", "PH", "Tsiatis", "Xie", "PR", "Stk.LR", "EDGE.poly3.u.G10", "EDGE.sym.sc.Grule")
X <- S[test %in% c(unique(CEN$test), NEW), .(cell, test, power = size_adj_power)]
fwrite(X, file.path(AN, "census_extended.csv"))
W <- dcast(X, cell ~ test, value.var = "power")
W <- merge(W, MEM, by = "cell", all.x = TRUE)
W[, best_part_G10 := pmax(HL.G10, HLF.G10, na.rm = TRUE)]
W[, best_part_rule := pmax(HL.Grule, HLF.Grule, na.rm = TRUE)]
W[, best_part_any := pmax(HL.G10, HLF.G10, HL.Grule, HLF.Grule, HL_w, PH, na.rm = TRUE)]
EDGE <- "EDGE.poly3.u.Grule"
cmp <- c("best_part_G10", "best_part_rule", "best_part_any", "HL.G10", "HL_w", "PH", "Tsiatis", "Xie", "PR",
         "GiViTI", "Stk.joint", "Stk.sym1", "Stk.LR", "Cubic.LR")
cnt <- rbindlist(lapply(cmp, function(r) {
  ok <- is.finite(W[[EDGE]]) & is.finite(W[[r]])
  lo <- ok & pmax(W[[EDGE]], W[[r]]) < 0.10                        # both tests below 0.10: a tie that says nothing
  data.table(comparator = r, cells = sum(ok), leads_or_ties = sum(W[[EDGE]][ok] > W[[r]][ok] - MARGIN),
             both_below_0.10 = sum(lo), leads_or_ties_informative = sum((W[[EDGE]] > W[[r]] - MARGIN)[ok & !lo]),
             informative = sum(ok & !lo), mean_edge = mean(W[[EDGE]][ok]), mean_rival = mean(W[[r]][ok]))
}))
fwrite(cnt, file.path(AN, "census_extended_counts.csv"))
print(cnt, digits = 3)

## by family for the grouped rivals, to see where each leads
fam <- rbindlist(lapply(c("HL_w", "PH", "Tsiatis", "Xie", "best_part_G10"), function(r) {
  W[is.finite(get(EDGE)) & is.finite(get(r)), .(comparator = r, cells = .N,
     leads_or_ties = sum(get(EDGE) > get(r) - MARGIN), mean_edge = mean(get(EDGE)), mean_rival = mean(get(r))),
    by = family_name]
}))
print(fam[order(comparator, family_name)], digits = 3)
cat("written: battery/analysis/census_extended.csv, census_extended_counts.csv\n")
