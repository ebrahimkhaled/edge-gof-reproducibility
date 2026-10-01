## _fig13_versions.R -- four candidate designs for Figure 13, for the authors to choose between.
##
##   v1  the DECLARED comparisons: each panel is a pairing the pre-declaration named, on its cells
##   v2  full census, 158 cells, each EDGE basis in the domain it is built for
##   v3  full census, 158 cells, the default EDGE (cubic basis) in every panel
##   v4  full census as one compact table: comparator x departure family, share of cells led or tied
##
## "Each basis in its domain" is fixed here, before any version is drawn, from what each basis is
## FOR (Section 2): EDGE-sym is the symmetric-link specialist, so it takes the two symmetric-link
## families; the default EDGE-poly3 takes everything else. It is not chosen cell by cell.
##
## Every number comes from _census_power_paper2.csv, which reproduces the declared analysis to
## 0.00e+00 on all 85 cells where the two overlap, or from block 8's declared summary.
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
OUT <- file.path(SIM, "_fig13_versions"); dir.create(OUT, showWarnings = FALSE)
source(file.path(SIM, "_ek_theme.R"))
suppressMessages({ library(ggplot2); library(data.table) })
AN <- edge_battery("analysis")

P   <- fread(edge_path("results/analysis/_census_power_paper2.csv"))
MEM <- unique(fread(file.path(AN, "rule_A_membership.csv"))[, .(cell, family, family_name)], by = "cell")
W   <- merge(dcast(P, cell + n ~ test, value.var = "power"), MEM, by = "cell")
W[, best_partition := pmax(HL.Grule, HLF.Grule, na.rm = TRUE)]
SYM_FAM <- c(1L, 3L)                                  # the two symmetric-link families
W[, edge_domain := ifelse(family %in% SYM_FAM, EDGE.sym.u.Grule, EDGE.poly3.u.Grule)]

S8 <- fread(edge_battery("8", "_summary.csv"))[role == "alternative" & alpha == 0.05]
s8 <- function(rival, test) S8[rival == get("rival", parent.frame(0)) & test == get("test", parent.frame(0))]
rival8 <- function(r) {
  tn <- if (r == "proj") "proj" else "BAGofT"
  x  <- dcast(S8[rival == r & test %in% c(tn, "EDGE.poly3.u.Grule", "EDGE.sym.u.Grule")],
              cell ~ test, value.var = "size_adj_power")
  setnames(x, tn, "rival")
  merge(x, MEM, by = "cell", all.x = TRUE)
}

pal <- c("symmetric tails" = "#0072B2", "asymmetric links" = "#D55E00",
         "Stukel symmetric family" = "#CC79A7", "omitted terms" = "#009E73",
         "rough misfit" = "#E69F00", "off-index" = "#7F7F7F", "census designs" = "#56B4E9", "other" = "#BBBBBB")
tri <- data.table(x = c(0, 0, 1, 0, 1, 1), y = c(0, 1, 1, 0, 0, 1),
                  half = rep(c("EDGE ahead", "comparator ahead"), each = 3))

scatter <- function(D, lev, file, subtitle, ncol = 3) {
  D <- D[is.finite(edge) & is.finite(rival)]
  D[is.na(family_name) & grepl("^census_", cell), family_name := "census designs"]
  D[is.na(family_name), family_name := "other"]
  lab <- D[, .(cells = .N, win = sum(edge > rival - 0.01)), by = panel]
  lab[, strip := sprintf("%s\nEDGE leads or ties %d of %d", panel, win, cells)]
  D <- merge(D, lab[, .(panel, strip)], by = "panel")
  D[, strip := factor(strip, levels = lab[match(lev, panel), strip])]
  p <- ggplot(D, aes(rival, edge)) +
    geom_polygon(data = tri, inherit.aes = FALSE, aes(x, y, fill = half), alpha = 0.15) +
    scale_fill_manual(values = c("EDGE ahead" = "#009E73", "comparator ahead" = "#D55E00"), guide = "none") +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey45", linewidth = 0.4) +
    geom_point(aes(colour = family_name), size = 1.7, alpha = 0.8) +
    facet_wrap(~ strip, ncol = ncol) +
    scale_colour_manual(values = pal, name = NULL) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(subtitle = subtitle, x = "size-adjusted power of the comparator", y = "size-adjusted power of EDGE") +
    theme_ek() +
    theme(legend.position = "bottom", panel.spacing = unit(0.9, "lines"),
          strip.text = element_text(size = rel(0.80), lineheight = 1.15))
  rows <- ceiling(length(lev) / ncol)
  h <- EK_W2 * (0.16 + 0.29 * rows)
  ggsave(file.path(OUT, paste0(file, ".png")), p, width = EK_W2, height = h, dpi = 150)
  ggsave(file.path(OUT, paste0(file, ".pdf")), p, width = EK_W2, height = h, device = cairo_pdf)
  cat("\n", file, "\n", sep = ""); print(lab[match(lev, panel), .(panel, cells, leads_or_ties = win)])
  invisible(lab)
}

## ---- v1: the declared comparisons ---------------------------------------------------------------
M   <- fread(file.path(AN, "mcnemar.csv"))[form == "unit" & arm == "Grule"]
CEN <- fread(file.path(AN, "h4_census.csv"))[form == "unit" & detectable == TRUE & saturated == FALSE]
v1 <- rbindlist(list(
  data.table(cell = CEN$cell, edge = CEN$edge_poly3, rival = CEN$best_partition, panel = "best partition test"),
  M[hypothesis %in% c("H1", "H3"), .(cell, edge = power_a, rival = power_b, panel = "GiViTI belt")],
  M[hypothesis == "H2", .(cell, edge = power_a, rival = power_b, panel = "Stukel joint score")]))
for (r in c("proj", "bagoft")) {
  x <- rival8(r)
  v1 <- rbind(v1, data.table(cell = x$cell, edge = x$EDGE.poly3.u.Grule, rival = x$rival,
                             panel = ifelse(r == "proj", "projection test", "BAGofT")))
}
v1 <- merge(v1, MEM[, .(cell, family_name)], by = "cell", all.x = TRUE)
scatter(v1, c("best partition test", "GiViTI belt", "Stukel joint score", "projection test", "BAGofT"),
        "v1_declared", "v1 -- the comparisons declared in advance, on their own cells (matches Table 1)")

## ---- v2 and v3: the full census ------------------------------------------------------------------
census <- function(edge_col) {
  mk <- function(col, lab) W[, .(cell, family_name, edge = get(edge_col), rival = get(col), panel = lab)]
  D <- rbindlist(list(mk("best_partition", "best partition test"), mk("GiViTI", "GiViTI belt"),
                      mk("Stk.joint", "Stukel joint score"), mk("Stk.sym1", "Stukel one-parameter"),
                      mk("Cubic.LR", "cubic calibration LR")))
  for (r in c("proj", "bagoft")) {
    x <- rival8(r)
    e <- if (edge_col == "edge_domain")
           ifelse(x$family %in% SYM_FAM, x$EDGE.sym.u.Grule, x$EDGE.poly3.u.Grule) else x$EDGE.poly3.u.Grule
    D <- rbind(D, data.table(cell = x$cell, family_name = x$family_name, edge = e, rival = x$rival,
                             panel = ifelse(r == "proj", "projection test", "BAGofT")))
  }
  D
}
LEV <- c("best partition test", "GiViTI belt", "Stukel joint score", "Stukel one-parameter",
         "cubic calibration LR", "projection test", "BAGofT")
scatter(census("edge_domain"), LEV, "v2_census_domain",
        "v2 -- all 158 cells; EDGE-sym on the symmetric-link families, the default EDGE elsewhere")
scatter(census("EDGE.poly3.u.Grule"), LEV, "v3_census_default",
        "v3 -- all 158 cells; the default EDGE (cubic basis) in every panel")

## ---- v4: the census as one compact table ----------------------------------------------------------
D4 <- census("edge_domain")[is.finite(edge) & is.finite(rival)]
D4 <- merge(D4, MEM[, .(cell, family)], by = "cell", all.x = TRUE)
T4 <- D4[, .(cells = .N, win = sum(edge > rival - 0.01)), by = .(panel, family_name)]
T4 <- T4[!is.na(family_name)]
T4[, share := win / cells][, lab := sprintf("%d/%d", win, cells)]
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
  labs(subtitle = "v4 -- all 158 cells as one table: cells in which EDGE leads or ties, by comparator and family",
       x = NULL, y = NULL) +
  theme_ek() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(),
        legend.position = "right")
ggsave(file.path(OUT, "v4_table.png"), p4, width = EK_W2, height = EK_W2 * 0.52, dpi = 150)
ggsave(file.path(OUT, "v4_table.pdf"), p4, width = EK_W2, height = EK_W2 * 0.52, device = cairo_pdf)
cat("\nv4_table\n"); print(dcast(T4, panel ~ family_name, value.var = "lab"))
cat("\nall four written to", OUT, "\n")
