## _fig13_final.R -- the two figures the authors chose (2026-09-19): v3 and v4 of _fig13_versions.R.
##
##   Fig13_scenarios.pdf   (v3) all 158 cells, the default EDGE (cubic basis) against every comparator
##   Fig14_census.pdf      (v4) the same census as one table, comparator x departure family
##
## Authors' changes to v3: a little padding inside the grey strip, and room above 1.00 on both axes
## so the strip no longer cuts the top tick label or the points that sit at power 1. No in-image
## title: the caption carries it. Numbers from _census_power_paper2.csv, which reproduces the declared
## analysis to 0.00e+00 on the 85 cells where they overlap, and from block 8's declared summary.
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
AN <- edge_battery("analysis")

P   <- fread(edge_path("results/analysis/_census_power_paper2.csv"))
MEM <- unique(fread(file.path(AN, "rule_A_membership.csv"))[, .(cell, family, family_name)], by = "cell")
W   <- merge(dcast(P, cell + n ~ test, value.var = "power"), MEM, by = "cell")
W[, best_partition := pmax(HL.Grule, HLF.Grule, na.rm = TRUE)]
EDGE <- "EDGE.poly3.u.Grule"

S8 <- fread(edge_battery("8", "_summary.csv"))[role == "alternative" & alpha == 0.05]
rival8 <- function(r) {
  tn <- if (r == "proj") "proj" else "BAGofT"
  x  <- dcast(S8[rival == r & test %in% c(tn, EDGE)], cell ~ test, value.var = "size_adj_power")
  setnames(x, tn, "rival")
  merge(x, MEM, by = "cell", all.x = TRUE)
}

mk <- function(col, lab) W[, .(cell, family_name, edge = get(EDGE), rival = get(col), panel = lab)]
D <- rbindlist(list(mk("best_partition", "best partition test"), mk("GiViTI", "GiViTI belt"),
                    mk("Stk.joint", "Stukel joint score"), mk("Stk.sym1", "Stukel one-parameter"),
                    mk("Cubic.LR", "cubic calibration LR")))
for (r in c("proj", "bagoft")) {
  x <- rival8(r)
  D <- rbind(D, data.table(cell = x$cell, family_name = x$family_name, edge = x[[EDGE]], rival = x$rival,
                           panel = ifelse(r == "proj", "Liu projection test", "BAGofT")))
}
## le Cessie-van Houwelingen (block 8L, declaration sha256 0c1a1731): size-adjusted power on block 8's twenty alternative
## cells, replicates 1-500, EDGE on the same replicates; both read from the runner's own summary. Drawn once 8L has run.
L8 <- edge_battery("8L", "_summary.csv")
if (file.exists(L8)) {
  x <- merge(fread(L8)[part == "A" & role == "alternative", .(cell, edge = edge_size_adj, rival = lecessie_size_adj)],
             MEM[, .(cell, family_name)], by = "cell", all.x = TRUE)
  D <- rbind(D, x[, .(cell, family_name, edge, rival, panel = "le Cessie-van Houwelingen")], use.names = TRUE)
}
## block 8R (declaration sha256 586d8f51): the same three rivals on the omitted-term and rough-misfit scenarios the
## rivals study had not covered. Same quantity, size-adjusted power, from that block's own summary.
R8 <- edge_battery("8R", "_summary.csv")
if (file.exists(R8)) {
  y <- fread(R8)[role == "alternative", .(cell, test, edge = edge_size_adj, rival = rival_size_adj)]
  y <- merge(y, MEM[, .(cell, family_name)], by = "cell", all.x = TRUE)
  y[, panel := c(lecessie = "le Cessie-van Houwelingen", proj = "Liu projection test", bagoft = "BAGofT")[test]]
  D <- rbind(D, y[, .(cell, family_name, edge, rival, panel)], use.names = TRUE)
}
## block 8S (declaration sha256 712e944d): the same three rivals at n = 200, the one sample size the
## rivals study had not covered. Same schema and the same quantity as 8R.
S8S <- edge_battery("8S", "_summary.csv")
if (file.exists(S8S)) {
  z <- fread(S8S)[role == "alternative", .(cell, test, edge = edge_size_adj, rival = rival_size_adj)]
  z <- merge(z, MEM[, .(cell, family_name)], by = "cell", all.x = TRUE)
  z[, panel := c(lecessie = "le Cessie-van Houwelingen", proj = "Liu projection test", bagoft = "BAGofT")[test]]
  D <- rbind(D, z[, .(cell, family_name, edge, rival, panel)], use.names = TRUE)
}
D <- D[is.finite(edge) & is.finite(rival)]
D[is.na(family_name), family_name := "other"]
LEV <- c("best partition test", "GiViTI belt", "Stukel joint score", "Stukel one-parameter",
         "cubic calibration LR", "le Cessie-van Houwelingen", "Liu projection test", "BAGofT")
LEV <- LEV[LEV %in% D$panel]

lab <- D[, .(cells = .N, win = sum(edge > rival - 0.01)), by = panel]     # Table 1's lead-or-tie margin
lab[, strip := sprintf("%s\nEDGE leads or ties %d of %d", panel, win, cells)]
D <- merge(D, lab[, .(panel, strip)], by = "panel")
D[, strip := factor(strip, levels = lab[match(LEV, panel), strip])]

pal <- c("symmetric tails" = "#0072B2", "asymmetric links" = "#D55E00",
         "Stukel symmetric family" = "#CC79A7", "omitted terms" = "#009E73",
         "rough misfit" = "#E69F00", "off-index" = "#7F7F7F", "other" = "#BBBBBB")

## the axes run a little past 1 so nothing at power 1 meets the strip; the green and red halves are
## drawn over the whole enlarged square, split on the diagonal
TOP <- 1.05
tri <- data.table(x = c(0, 0, TOP, 0, TOP, TOP), y = c(0, TOP, TOP, 0, 0, TOP),
                  half = rep(c("EDGE ahead", "comparator ahead"), each = 3))

p3 <- ggplot(D, aes(rival, edge)) +
  geom_polygon(data = tri, inherit.aes = FALSE, aes(x, y, fill = half), alpha = 0.15) +
  scale_fill_manual(values = c("EDGE ahead" = "#009E73", "comparator ahead" = "#D55E00"), guide = "none") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey45", linewidth = 0.4) +
  geom_point(aes(colour = family_name), size = 1.7, alpha = 0.8) +
  facet_wrap(~ strip, ncol = 3) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) + scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_equal(xlim = c(0, TOP), ylim = c(0, TOP), expand = FALSE) +
  labs(x = "size-adjusted power of the comparator", y = "size-adjusted power of EDGE") +
  theme_ek() +
  theme(legend.position = "bottom", panel.spacing = unit(1.0, "lines"),
        strip.text = element_text(size = rel(0.80), lineheight = 1.15,
                                  margin = margin(t = 3.5, b = 3.5, l = 2, r = 2)))
h3 <- EK_W2 * (0.16 + 0.29 * 3)
ggsave(file.path(FIG, "Fig13_scenarios.pdf"), p3, width = EK_W2, height = h3, device = cairo_pdf)
ggsave(file.path(SIM, "_preview_Fig13_final.png"), p3, width = EK_W2, height = h3, dpi = 150)
print(lab[match(LEV, panel), .(panel, cells, leads_or_ties = win)])

## ---- v4: the same census as one table ---------------------------------------------------------
D4 <- merge(D[, .(cell, panel, edge, rival)], MEM[, .(cell, family_name)], by = "cell", all.x = TRUE)
T4 <- D4[!is.na(family_name), .(cells = .N, win = sum(edge > rival - 0.01)), by = .(panel, family_name)]
T4[, `:=`(share = win / cells, lab = sprintf("%d/%d", win, cells))]
tot <- D4[, .(cells = .N, win = sum(edge > rival - 0.01)), by = panel]
tot[, `:=`(family_name = "all families", share = win / cells, lab = sprintf("%d/%d", win, cells))]
T4 <- rbind(T4, tot, fill = TRUE)
fam_lev <- c("symmetric tails", "asymmetric links", "Stukel symmetric family", "omitted terms",
             "rough misfit", "off-index", "all families")
T4[, family_name := factor(family_name, levels = fam_lev)]
T4[, panel := factor(panel, levels = rev(LEV))]
p4 <- ggplot(T4, aes(family_name, panel, fill = share)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = lab), size = 2.7, family = EK_FAMILY) +
  scale_fill_gradient2(low = "#F4C7B8", mid = "#F7F7F7", high = "#9FD5B8", midpoint = 0.5,
                       limits = c(0, 1), labels = scales::percent, name = "EDGE leads\nor ties") +
  geom_vline(xintercept = 6.5, colour = "grey30", linewidth = 0.5) +
  labs(x = NULL, y = NULL) +
  theme_ek() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(),
        legend.position = "right")
ggsave(file.path(FIG, "Fig14_census.pdf"), p4, width = EK_W2, height = EK_W2 * 0.50, device = cairo_pdf)
ggsave(file.path(SIM, "_preview_Fig14_final.png"), p4, width = EK_W2, height = EK_W2 * 0.50, dpi = 150)
cat("written: Fig13_scenarios.pdf and Fig14_census.pdf\n")
