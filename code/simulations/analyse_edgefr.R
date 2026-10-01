## analyse_edgefr.R -- the one analysis file behind the paper's EDGE-FR table.
## Reads blocks 9R2 (sha256 9943e551), 9R3 (sha256 421217c0) and the matched null of 9R2n, and writes
## battery/analysis/edgefr.csv and battery/analysis/edgefr_paired.csv.
##
## The paper's own reading rules (Section 4.3) are applied here to the paper's own new results:
##   - power is SIZE-ADJUSTED, using each variant's critical value read from its matched null scenario
##     (the type-1 5% quantile of that variant's null p-values), not the nominal 0.05;
##   - comparisons between variants are PAIRED inside the scenario, by exact McNemar on the discordant
##     pairs, Holm-corrected within the family of comparisons reported.
## False-alarm rates under contamination are rejection rates at the nominal level, as everywhere else in
## Section 6: they are not power and are not size-adjusted.
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
BAT <- edge_battery()
V <- c("V0", "HYB05", "HYB10", "HYB20", "G10")
BATT <- c("cauchit_n1000", "t4_n1000", "loglog_n1000", "stk_short_n1000", "stk_long_n1000", "stk_asym_n1000")
## every power scenario and the null its size adjustment is read from
NULL_OF <- c(setNames(rep("null_link_n1000", length(BATT)), BATT),
             cloglog_clean_n1000 = "logit_clean_n1000", probit_clean_n5000 = "logit_clean_n5000",
             probit_clean_n1000  = "logit_clean_n1000", cloglog_clean_n5000 = "logit_clean_n5000",
             cloglog_C1_r005_n1000 = "logit_clean_n1000", probit_C1_r005_n1000 = "logit_clean_n1000")
DIRS <- c("9R2", "9R3", "9R2n")

## one read per scenario, kept in memory: every rate below is asked for five variants and several comparisons,
## and re-opening a gzipped file each time costs minutes for no reason
RD <- new.env(parent = emptyenv())
rd <- function(cell) {
  if (!is.null(RD[[cell]])) return(RD[[cell]])
  for (d in DIRS) {
    f <- edge_battery(d, sprintf("%s_variants.csv.gz", cell))
    if (file.exists(f)) { RD[[cell]] <- fread(f); return(RD[[cell]]) }
  }
  stop("no per-replicate file for ", cell)
}
cells_in <- function(d) sub("_variants.csv.gz$", "", basename(Sys.glob(edge_battery(d, "*_variants.csv.gz"))))
null_of <- function(cl) { x <- NULL_OF[cl]; if (is.na(x)) NULL else unname(x) }

## the critical value of one variant in one null scenario: the type-1 5% quantile, a replicate with no
## p-value never counted as a rejection
crit <- function(null_cell, v) {
  p <- rd(null_cell)[[v]]
  p[!is.finite(p)] <- 1
  as.numeric(stats::quantile(p, 0.05, type = 1))
}
CRIT <- new.env()
crit_cached <- function(null_cell, v) {
  k <- paste(null_cell, v)
  if (is.null(CRIT[[k]])) CRIT[[k]] <- crit(null_cell, v)
  CRIT[[k]]
}

## rejections of one variant in one scenario, at the nominal level and size-adjusted where a null is matched
rej <- function(cell, v, adjusted) {
  p <- rd(cell)[[v]]
  p[!is.finite(p)] <- 1
  if (!adjusted) return(p <= 0.05)
  nc <- null_of(cell)
  if (is.null(nc)) return(rep(NA, length(p)))
  p <= crit_cached(nc, v)
}

all_cells <- unique(c(cells_in("9R2"), cells_in("9R3"), cells_in("9R2n")))
R <- rbindlist(lapply(all_cells, function(cl) {
  M <- rd(cl)
  n <- as.integer(sub(".*_n(\\d+)$", "\\1", cl))
  k <- if (grepl("^logit_C1_r", cl)) as.integer(round(as.numeric(sub(".*_r(\\d+)_.*", "\\1", cl)) / 1000 * n))
       else if (grepl("^fr_C1_k", cl)) as.integer(sub(".*_k(\\d+)_.*", "\\1", cl)) else NA_integer_
   role <- if (!is.na(k)) "false alarm"
           else if (grepl("^(logit|fr)_clean", cl) || grepl("^null_", cl)) "level"
           else if (cl %in% BATT) "power (battery)" else "power"
  data.table(cell = cl, role = role, n = n, k = k, variant = V, reps = nrow(M),
             rejection = vapply(V, function(v) mean(rej(cl, v, FALSE), na.rm = TRUE), 0),
             adjusted = vapply(V, function(v) mean(rej(cl, v, TRUE), na.rm = TRUE), 0),
             crit = vapply(V, function(v) if (is.null(null_of(cl))) NA_real_ else crit_cached(null_of(cl), v), 0),
             missing = vapply(V, function(v) sum(!is.finite(M[[v]])), 0L))
}))
mb <- R[cell %in% BATT, .(cell = "battery mean (6)", role = "power (mean)", n = 1000L, k = NA_integer_,
                          reps = 1000L, rejection = mean(rejection), adjusted = mean(adjusted),
                          crit = NA_real_, missing = sum(missing)), by = variant]
R <- rbind(R, mb, fill = TRUE)
setcolorder(R, c("cell", "role", "n", "k", "variant", "rejection", "adjusted", "crit", "reps", "missing"))
setorder(R, role, n, k, cell, variant)

## paired comparisons, size-adjusted, on the six battery scenarios and the two other power scenarios
PAIRS <- list(c("HYB05", "V0"), c("HYB05", "G10"), c("HYB10", "G10"), c("HYB20", "G10"), c("G10", "V0"))
POWER_CELLS <- c(BATT, "cloglog_clean_n1000", "probit_clean_n5000")
P <- rbindlist(lapply(PAIRS, function(ab) rbindlist(lapply(POWER_CELLS, function(cl) {
  ra <- rej(cl, ab[1], TRUE); rb <- rej(cl, ab[2], TRUE)
  n10 <- sum(ra & !rb); n01 <- sum(!ra & rb)
  data.table(a = ab[1], b = ab[2], cell = cl, a_power = mean(ra), b_power = mean(rb),
             diff = mean(ra) - mean(rb), a_only = n10, b_only = n01,
             p = if (n10 + n01 == 0) 1 else stats::binom.test(n10, n10 + n01, 0.5)$p.value)
}))))
P[, p_holm := p.adjust(p, "holm"), by = .(a, b)]

## the same pairing on the contaminated scenarios, where the question is protection, not power
FA <- R[role == "false alarm", unique(cell)]
Q <- rbindlist(lapply(list(c("HYB05", "G10"), c("HYB05", "V0"), c("HYB20", "G10")), function(ab)
  rbindlist(lapply(FA, function(cl) {
    ra <- rej(cl, ab[1], FALSE); rb <- rej(cl, ab[2], FALSE)
    n10 <- sum(ra & !rb); n01 <- sum(!ra & rb)
    data.table(a = ab[1], b = ab[2], cell = cl, a_rate = mean(ra), b_rate = mean(rb), diff = mean(ra) - mean(rb),
               a_only = n10, b_only = n01,
               p = if (n10 + n01 == 0) 1 else stats::binom.test(n10, n10 + n01, 0.5)$p.value)
  }))))
Q[, p_holm := p.adjust(p, "holm"), by = .(a, b)]

dir.create(edge_battery("analysis"), showWarnings = FALSE, recursive = TRUE)
fwrite(R, edge_battery("analysis", "edgefr.csv"))
fwrite(rbind(P[, kind := "power (size-adjusted)"], Q[, kind := "false alarm (nominal)"], fill = TRUE),
       edge_battery("analysis", "edgefr_paired.csv"))

cat("\n== critical values from the matched nulls (type-1 5% quantile) ==\n")
print(unique(R[role %in% c("power (battery)", "power") & is.finite(crit), .(cell, variant, crit)])[order(cell, variant)][1:10], digits = 3)
cat("\n== power: raw against size-adjusted ==\n")
print(dcast(R[role %in% c("power (battery)", "power", "power (mean)")], cell ~ variant,
            value.var = "adjusted")[], digits = 3)
cat("\n== mean over the six battery scenarios ==\n")
mm <- R[cell == "battery mean (6)"]
print(mm[, .(variant, raw = round(rejection, 4), size_adjusted = round(adjusted, 4))])
cat("\n== paired, size-adjusted, over the eight power scenarios ==\n")
print(P[, .(scenarios = .N, mean_diff = round(mean(diff), 4), a_wins = sum(diff > 0), b_wins = sum(diff < 0),
            holm_signif = sum(p_holm < 0.05)), by = .(a, b)])
cat("\n== paired false alarms, HYB05 against G10, by scenario ==\n")
print(Q[a == "HYB05" & b == "G10"][order(cell), .(cell, HYB05 = round(a_rate, 3), G10 = round(b_rate, 3),
                                                  diff = round(diff, 3), holm = signif(p_holm, 2))])
cat(sprintf("\nreplicates with no p-value: %d\n", sum(R$missing, na.rm = TRUE)))
