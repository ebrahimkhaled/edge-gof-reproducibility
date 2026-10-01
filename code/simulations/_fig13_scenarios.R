## _fig13_scenarios.R -- Figure: the default EDGE against every comparator, one point per cell.
##
## What the figure shows is the comparison an analyst would actually run: the DEFAULT directed test
## (cubic basis, unit form, rule G) against each comparator, on every alternative cell that carries a
## declared departure family. Wins and losses alike; no cell is chosen. The declared verdicts stay in
## Table 1 on the declared cells; this is the descriptive census beside them.
##
## Why the default in every panel. EDGE-sym is a one-column specialist for symmetric departures;
## pairing it with Stukel's two-column joint score across the asymmetric-link cells would manufacture
## losses that say nothing about how the test is used. H2 restricted that pairing to symmetric tails
## for exactly that reason, and it is reported there.
##
## Size-adjusted power uses the study's own rule: the critical value is the alpha quantile of the
## test's own matched null cell, and a test with no p-value counts as no rejection. The resampling
## tests ran in block 8 only, so their panels come from block 8's declared summary.
##
## The computed table is saved (_census_power_paper2.csv) so a re-plot never recomputes it.
## Foundation: _ek_theme.R.
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
SIM <- edge_path("code/simulations")
FIG <- edge_out("paper2_figures")
source(file.path(SIM, "_ek_theme.R"))
suppressMessages({ library(ggplot2); library(data.table) })
BAT <- edge_battery(); AN <- edge_battery("analysis")
CACHE <- edge_path("results/analysis/_census_power_paper2.csv")
ALPHA <- 0.05
EDGE  <- "EDGE.poly3.u.Grule"
WANT  <- c(EDGE, "EDGE.sym.u.Grule", "HL.Grule", "HLF.Grule", "GiViTI", "Stk.joint", "Stk.sym1", "Cubic.LR")

MEM <- unique(fread(file.path(AN, "rule_A_membership.csv"))[, .(cell, family_name)], by = "cell")

if (file.exists(CACHE)) {
  P <- fread(CACHE)
  cat("read the saved census:", nrow(P), "rows\n")
} else {
  C <- fread(edge_battery("cells.csv"),
             colClasses = list(character = c("block", "cell", "null_block", "null_cell")))
  C <- C[!is.na(null_cell) & null_cell != "" & cell %in% MEM$cell]
  cat("alternative cells with a declared family:", nrow(C), "\n")

  read_p <- function(block, cell) {
    f <- edge_battery(block, paste0(cell, "_pvalues.csv.gz"))
    if (!file.exists(f)) return(NULL)
    have <- intersect(WANT, names(fread(f, nrows = 0)))
    if (!length(have)) return(NULL)
    fread(f, select = have)
  }
  ## a null is shared by many alternatives: read each once and keep its critical values
  crit_cache <- new.env()
  crit_of <- function(block, cell) {
    key <- paste(block, cell)
    if (!is.null(crit_cache[[key]])) return(crit_cache[[key]])
    N <- read_p(block, cell)
    cr <- if (is.null(N)) NULL else vapply(names(N), function(t) {
      pn <- N[[t]]; pn <- ifelse(is.finite(pn), pn, 1)
      as.numeric(stats::quantile(pn, ALPHA, type = 1))
    }, numeric(1))
    assign(key, cr, envir = crit_cache)
    cr
  }

  rows <- vector("list", nrow(C))
  for (i in seq_len(nrow(C))) {
    ce <- C[i]
    A <- read_p(ce$block, ce$cell); cr <- crit_of(ce$null_block, ce$null_cell)
    if (is.null(A) || is.null(cr)) next
    tt <- intersect(names(A), names(cr))
    rows[[i]] <- data.table(cell = ce$cell, n = ce$n, test = tt,
                            power = vapply(tt, function(t) mean(is.finite(A[[t]]) & A[[t]] <= cr[[t]]),
                                           numeric(1)))
    if (i %% 25 == 0) cat("  ", i, "of", nrow(C), "cells\n")
  }
  P <- rbindlist(rows)
  fwrite(P, CACHE)
  cat("census computed and saved:", nrow(P), "rows\n")
}

W <- merge(dcast(P, cell + n ~ test, value.var = "power"), MEM, by = "cell")
W[, best_partition := pmax(HL.Grule, HLF.Grule, na.rm = TRUE)]

pan <- function(col, lab, grouped, edge_col = EDGE) {
  x <- W[is.finite(get(col)) & is.finite(get(edge_col))]
  data.table(cell = x$cell, family_name = x$family_name, rival = x[[col]], edge = x[[edge_col]],
             panel = lab, grouped = grouped)
}
## Every panel uses the default EDGE except one. Stukel's one-parameter score tests the single
## symmetric direction eta|eta|, and EDGE-sym is its grouped counterpart -- the same column -- so
## that pairing is like-for-like and is shown on EVERY cell, not a chosen subset.
D <- rbindlist(list(
  pan("best_partition", "best partition test",  TRUE),
  pan("GiViTI",         "GiViTI belt",          FALSE),
  pan("Stk.joint",      "Stukel joint score",   FALSE),
  pan("Stk.sym1",       "Stukel one-parameter (EDGE-sym)", FALSE, edge_col = "EDGE.sym.u.Grule"),
  pan("Cubic.LR",       "cubic calibration LR", FALSE)))

## the two resampling tests ran in block 8 only
S8 <- fread(edge_battery("8", "_summary.csv"))[role == "alternative" & alpha == 0.05]
for (r in c("proj", "bagoft")) {
  e <- S8[rival == r & test == EDGE, .(cell, edge = size_adj_power)]
  v <- S8[rival == r & test == ifelse(r == "proj", "proj", "BAGofT"), .(cell, rival = size_adj_power)]
  m <- merge(merge(e, v, by = "cell"), MEM, by = "cell", all.x = TRUE)
  D <- rbind(D, data.table(cell = m$cell, family_name = m$family_name, rival = m$rival,
                           edge = m$edge, panel = ifelse(r == "proj", "projection test", "BAGofT"),
                           grouped = FALSE), use.names = TRUE)
}
D <- D[is.finite(edge) & is.finite(rival)]
D[is.na(family_name), family_name := "other"]

## the same lead-or-tie margin as Table 1, so figure and table cannot disagree
lab <- D[, .(cells = .N, win = sum(edge > rival - 0.01)), by = .(panel, grouped)]
lab[, strip := sprintf("%s\nEDGE leads or ties %d of %d", panel, win, cells)]
D <- merge(D, lab[, .(panel, strip)], by = "panel")
lev <- c("best partition test", "GiViTI belt", "Stukel joint score",
         "Stukel one-parameter (EDGE-sym)", "cubic calibration LR", "projection test", "BAGofT")
D[, strip := factor(strip, levels = lab[match(lev, panel), strip])]

pal <- c("symmetric tails" = "#0072B2", "asymmetric links" = "#D55E00",
         "Stukel symmetric family" = "#CC79A7", "omitted terms" = "#009E73",
         "rough misfit" = "#E69F00", "off-index" = "#7F7F7F", "other" = "#BBBBBB")

## tint the half of each panel in which EDGE is ahead green, the other half red
tri <- data.table(x = c(0, 0, 1, 0, 1, 1), y = c(0, 1, 1, 0, 0, 1),
                  half = rep(c("EDGE ahead", "comparator ahead"), each = 3))

p <- ggplot(D, aes(rival, edge)) +
  geom_polygon(data = tri, inherit.aes = FALSE, aes(x, y, fill = half), alpha = 0.15) +
  scale_fill_manual(values = c("EDGE ahead" = "#009E73", "comparator ahead" = "#D55E00"),
                    guide = "none") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey45", linewidth = 0.4) +
  geom_point(aes(colour = family_name), size = 1.8, alpha = 0.85) +
  facet_wrap(~ strip, nrow = 2) +
  scale_colour_manual(values = pal, name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "size-adjusted power of the comparator", y = "size-adjusted power of EDGE") +
  theme_ek() +
  theme(legend.position = "bottom", panel.spacing = unit(0.9, "lines"),
        strip.text = element_text(size = rel(0.80), lineheight = 1.15))

dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIG, "Fig13_scenarios.pdf"), p, width = EK_W2, height = EK_W2 * 0.68,
       device = cairo_pdf)
ggsave(file.path(SIM, "_preview_Fig13_scenarios.png"), p, width = EK_W2, height = EK_W2 * 0.68,
       dpi = 150)

print(lab[match(lev, panel), .(panel, cells, leads_or_ties = win)])
cat("\ncells in the census:", uniqueN(D$cell), "\n")
cat("written: figures/Fig13_scenarios.pdf\n")
