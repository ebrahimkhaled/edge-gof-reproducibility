## _fig16_8R.R -- block 8R (declaration sha256 586d8f51): the three outside rivals on the two families the rivals study
## had not covered, omitted terms and rough misfit. One panel per rival, one point per scenario, size-adjusted power,
## drawn from battery/8R/_summary.csv (the declared summary). Panels appear only for tests that have finished.
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

S <- fread(edge_battery("8R", "_summary.csv"))[role == "alternative"]
MEM <- unique(fread(edge_battery("analysis", "rule_A_membership.csv"))[, .(cell, family_name)], by = "cell")
S <- merge(S, MEM, by = "cell", all.x = TRUE)
S[, ladder := sub("_[0-9.]+_n[0-9]+$|_n[0-9]+$", "", cell)]
S[, departure := fifelse(family_name == "rough misfit", "rough misfit",
                  fifelse(grepl("^binint", ladder), "omitted interaction (binary)",
                   fifelse(grepl("^contint", ladder), "omitted interaction (continuous)", "omitted quadratic term")))]
nice <- c(lecessie = "le Cessie", proj = "Liu projection", bagoft = "BAGofT")
S[, panel := nice[test]]
S[, lead_or_tie := edge_size_adj > rival_size_adj - 0.01]

lab <- S[, .(scen = .N, win = sum(lead_or_tie)), by = .(test, panel)]
lab[, strip := sprintf("%s\nEDGE leads or ties %d of %d", panel, win, scen)]
S <- merge(S, lab[, .(test, strip)], by = "test")
S[, strip := factor(strip, levels = lab[match(intersect(c("lecessie", "proj", "bagoft"), lab$test), test), strip])]

TOP <- 1.05
tri <- data.table(x = c(0, 0, TOP, 0, TOP, TOP), y = c(0, TOP, TOP, 0, 0, TOP),
                  half = rep(c("EDGE ahead", "rival ahead"), each = 3))
pal <- c("omitted quadratic term" = "#009E73", "omitted interaction (binary)" = "#0072B2",
         "omitted interaction (continuous)" = "#56B4E9", "rough misfit" = "#E69F00")

p <- ggplot(S, aes(rival_size_adj, edge_size_adj)) +
  geom_polygon(data = tri, inherit.aes = FALSE, aes(x, y, fill = half), alpha = 0.15) +
  scale_fill_manual(values = c("EDGE ahead" = "#009E73", "rival ahead" = "#D55E00"), guide = "none") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey45", linewidth = 0.4) +
  geom_point(aes(colour = departure, shape = factor(n)), size = 1.9, alpha = 0.85) +
  facet_wrap(~ strip, nrow = 1) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_shape_manual(values = c("200" = 1, "500" = 16, "1000" = 17), name = "n") +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) + scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_equal(xlim = c(0, TOP), ylim = c(0, TOP), expand = FALSE) +
  labs(x = "size-adjusted power of the rival", y = "size-adjusted power of EDGE") +
  theme_ek() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(1, "pt"),
        panel.spacing = unit(1.0, "lines"),
        strip.text = element_text(size = rel(0.80), lineheight = 1.15, margin = margin(t = 3.5, b = 3.5, l = 2, r = 2)))
np <- length(unique(S$strip))
h <- EK_W2 * (if (np >= 3) 0.48 else 0.66)
ggsave(file.path(FIG, "Fig16_rough_omitted.pdf"), p, width = EK_W2, height = h, device = cairo_pdf)
ggsave(file.path(SIM, "_preview_Fig16_8R.png"), p, width = EK_W2, height = h, dpi = 150)
print(S[, .(scenarios = .N, EDGE_leads_or_ties = sum(lead_or_tie), mean_EDGE = round(mean(edge_size_adj), 3),
            mean_rival = round(mean(rival_size_adj), 3)), by = .(panel, departure)][order(panel, departure)])
cat("written: Fig16_rough_omitted.pdf (", np, "panels )\n")
