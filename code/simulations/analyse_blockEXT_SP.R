## analyse_blockEXT_SP.R -- the readings fixed in PREDECLARATION_blockEXT_SP.md (sha256 7455de2d...).
##   (a) level: SP-BB rejection at 5% in cells 1-3; holds within 0.029 to 0.071
##   (b) power: size-adjusted at the 5% quantile (type 1) of each test's p-values in cell 2; paired differences
##       EDGE (ten groups) minus SP-BB and Cox minus SP-BB on the same replicates, by cell, family and overall
##   (c) protection: raw false-alarm rate in cells 14-21; holds where at most 0.10
## EDGE and Cox p-values are block EXT's stored values for the same data sets. No p-value = non-rejection.
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
suppressMessages(library(data.table))
IN_EXT <- edge_battery("EXT"); IN <- edge_battery("EXT_SP")
fs <- sort(list.files(IN, pattern = "_pvalues\\.csv\\.gz$", full.names = TRUE))
stopifnot(length(fs) == 21L)
P <- rbindlist(lapply(fs, function(f) {
  sp <- fread(f); ext <- fread(file.path(IN_EXT, basename(f)))
  stopifnot(nrow(sp) == nrow(ext), identical(sp$rep, ext$rep), identical(sp$events, ext$events),
            isTRUE(all.equal(sp$spz_check, ext$spiegelhalter)))
  m <- regmatches(basename(f), regexec("c(\\d+)_([a-z]+)_C([-0-9.]+)_n(\\d+)_k(\\d+)", basename(f)))[[1]]
  data.table(id = as.integer(m[2]), fam = m[3], C = as.numeric(m[4]), n = as.integer(m[5]), k = as.integer(m[6]),
             rep = sp$rep, SP_BB = sp$SP_BB, SP_BM = sp$SP_BM, EDGE_G10 = ext$EDGE_G10, EDGE_auto = ext$EDGE_auto,
             cox = ext$cox, giviti = ext$giviti, spiegelhalter = ext$spiegelhalter, notes = sp$notes)
}))
cat(sprintf("replicates %d; missing SP_BB %d; tie/other notes %d\n", nrow(P), sum(is.na(P$SP_BB)), sum(P$notes)))
TESTS <- c("SP_BB", "SP_BM", "EDGE_G10", "EDGE_auto", "cox", "giviti", "spiegelhalter")
crit <- P[id == 2, lapply(.SD, function(p) as.numeric(quantile(p, 0.05, type = 1, na.rm = TRUE))), .SDcols = TESTS]
rej <- function(p, c) !is.na(p) & p <= c

## ---- rates per cell and test
L <- melt(P, id.vars = c("id", "fam", "C", "n", "k", "rep"), measure.vars = TESTS, variable.name = "test", value.name = "p")
L[, test := as.character(test)]
L[, crit := as.numeric(crit[[test[1]]]), by = test]
S <- L[, .(raw = mean(!is.na(p) & p < 0.05), adj = mean(rej(p, crit)), reps = .N), by = .(id, fam, C, n, k, test)]
S[, role := ifelse(fam == "null", "level", ifelse(k > 0, "corrupted", "alternative"))]
S[, holds := ifelse(role == "level", raw >= 0.029 & raw <= 0.071, ifelse(role == "corrupted", raw <= 0.10, NA))]
setorder(S, test, id)
fwrite(S, file.path(IN, "_summary.csv"))

## ---- paired size-adjusted differences on the same replicates
FAMLAB <- c(large = "shift", slope = "slope", ushape = "U-shape", thresh = "threshold", inter = "interaction")
A <- P[id %in% 4:13]
A[, `:=`(e = rej(EDGE_G10, crit$EDGE_G10), s = rej(SP_BB, crit$SP_BB), x = rej(cox, crit$cox))]
cellD <- A[, .(power_SP = mean(s), power_EDGE_G10 = mean(e), power_cox = mean(x),
               d_EDGE_SP = mean(e - s), se_EDGE_SP = sd(e - s) / sqrt(.N),
               d_cox_SP = mean(x - s), se_cox_SP = sd(x - s) / sqrt(.N), reps = .N), by = .(id, fam, C)]
agg <- function(D, lab) D[, .(level = lab, power_SP = mean(power_SP), power_EDGE_G10 = mean(power_EDGE_G10),
                              power_cox = mean(power_cox),
                              d_EDGE_SP = mean(d_EDGE_SP), se_EDGE_SP = sqrt(sum(se_EDGE_SP^2)) / .N,
                              d_cox_SP = mean(d_cox_SP), se_cox_SP = sqrt(sum(se_cox_SP^2)) / .N, cells = .N)]
fam <- rbindlist(lapply(names(FAMLAB), function(f) agg(cellD[fam == f], paste("family", FAMLAB[[f]]))))
all10 <- agg(cellD, "overall (ten alternatives)")
cl <- cellD[, .(level = sprintf("cell %d %s %s", id, FAMLAB[fam], C), power_SP, power_EDGE_G10, power_cox,
                d_EDGE_SP, se_EDGE_SP, d_cox_SP, se_cox_SP, cells = 1L)]
PD <- rbind(cl, fam, all10)
PD[, `:=`(lo_EDGE_SP = d_EDGE_SP - 1.96 * se_EDGE_SP, hi_EDGE_SP = d_EDGE_SP + 1.96 * se_EDGE_SP,
          lo_cox_SP = d_cox_SP - 1.96 * se_cox_SP, hi_cox_SP = d_cox_SP + 1.96 * se_cox_SP)]
fwrite(PD, file.path(IN, "_paired.csv"))

## ---- printed readings
wide <- function(D, v) dcast(D, id + fam + C + n + k ~ test, value.var = v)[, c("id", "fam", "C", "n", "k", TESTS), with = FALSE]
cat("\n(a) level, raw rejection at 5%\n"); print(wide(S[role == "level"], "raw"), digits = 3)
cat("\n(b) size-adjusted power (critical values from cell 2)\n"); print(wide(S[role == "alternative"], "adj"), digits = 3)
cat("\nmean size-adjusted power over the ten alternatives\n")
print(S[role == "alternative", .(mean_adj = round(mean(adj), 3)), by = test][order(-mean_adj)])
cat("\npaired differences (95% intervals)\n"); print(PD[, .(level, power_SP, power_EDGE_G10, power_cox, d_EDGE_SP,
  lo_EDGE_SP, hi_EDGE_SP, d_cox_SP, lo_cox_SP, hi_cox_SP)], digits = 3)
cat("\n(c) false alarms under corruption, raw at 5%\n"); print(wide(S[role == "corrupted"], "raw"), digits = 3)
H <- S[role == "corrupted", .(holds = sum(holds)), by = test][order(-holds)]
cat("\ncorrupted cells held (<= 0.10), of 8\n"); print(H)

## ---- the declared reading of what counts against EDGE
sp_holds_all <- H[test == "SP_BB", holds] == 8L
o <- all10
verdict <- if (sp_holds_all && o$d_EDGE_SP <= 0 && o$d_EDGE_SP + 1.96 * o$se_EDGE_SP < 0) "STRONGLY AGAINST EDGE (SP-BB dominates)" else
  if (sp_holds_all && o$d_EDGE_SP <= 0) "AGAINST EDGE (SP-BB holds everywhere and is at least as powerful)" else
  if (!sp_holds_all || o$d_EDGE_SP - 1.96 * o$se_EDGE_SP > 0) "FOR EDGE" else
  "NEUTRAL (protection shared; EDGE's power edge not shown)"
cat(sprintf("\nSP-BB holds in all eight corrupted cells: %s; overall EDGE - SP %.3f (95%% %.3f to %.3f)\nDeclared reading: %s\n",
            sp_holds_all, o$d_EDGE_SP, o$d_EDGE_SP - 1.96 * o$se_EDGE_SP, o$d_EDGE_SP + 1.96 * o$se_EDGE_SP, verdict))
