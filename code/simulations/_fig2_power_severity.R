## _fig2_power_severity.R — HEADLINE POWER FIGURE (Fig 2)
## Across omitted-curvature and interaction misfit, EDGE gains power fastest and
## dominates the omnibus partition tests (HL, EF).
## CAPTION: EDGE gains power fastest under omitted curvature and interaction.
## Renders: figures/Fig2_power_severity.pdf (vector) + figures/_preview_Fig2_power_severity.png
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
suppressMessages({
  library(ggplot2); library(dplyr); library(readr); library(ggrepel); library(scales)
})

sim_dir <- edge_path("code/simulations")
fig_dir <- edge_out("figures")
source(file.path(sim_dir, "_ek_theme.R"))
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- curated tests ----
curated <- c("DEF.poly2","DEF.poly3","DEF.stukel","EF","HL","HL-equalwidth",
             "Pigeon-Heyse","Stukel","Tsiatis","Xie")

fam_labels <- c(
  quad    = "Omitted curvature (quadratic)",
  binint  = "Omitted interaction (binary)",
  contint = "Omitted interaction (cont.)"
)

raw <- read_csv(file.path(sim_dir, "sim_power_broad.csv"), show_col_types = FALSE)

d <- raw %>%
  filter(alpha == 0.05, G == 10, n == 1000,
         family %in% c("quad","binint","contint"),
         test %in% curated) %>%
  mutate(param = suppressWarnings(as.numeric(param))) %>%
  filter(!is.na(param)) %>%
  mutate(Test   = ek_factor(test),
         family = factor(family, levels = c("quad","binint","contint"),
                         labels = fam_labels[c("quad","binint","contint")]),
         is_edge = grepl("^EDGE", as.character(Test)))

## ---- EDGE >= HL advantage band (per family, per severity) ----
## For each family x param: shaded ribbon between the HL curve and the best EDGE curve.
band <- d %>%
  group_by(family, param) %>%
  summarise(
    hl_y   = reject_rate[Test == "HL"][1],
    edge_y = max(reject_rate[grepl("^EDGE", as.character(Test))], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(hl_y), !is.na(edge_y)) %>%
  mutate(ymin = pmin(hl_y, edge_y), ymax = pmax(hl_y, edge_y))

## (in-plot end-labels removed --- the bottom legend carries test identification)

## midpoint annotation: where does EDGE first clear a big margin over HL?
## pick the family panel for the "EDGE > HL" callout (quad, mid severity)
callout <- band %>%
  filter(grepl("quadratic", family)) %>%
  slice(which.min(abs(param - 0.02))) %>%
  mutate(lab = sprintf("At severity 0.02:\nEDGE %.0f%% power  vs  HL %.0f%%",
                       100 * pmax(hl_y, edge_y), 100 * pmin(hl_y, edge_y)))

p <- ggplot(d, aes(param, reject_rate)) +
  ## advantage band: EDGE-best over HL
  geom_ribbon(data = band,
              aes(x = param, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = ek_pal[["EDGE-poly2"]], alpha = 0.11) +
  ek_nominal(0.05) +
  ## non-EDGE lines (thin, dashed by family via lty)
  geom_line(data = filter(d, !is_edge),
            aes(colour = Test, linetype = Test, group = Test), linewidth = 0.55) +
  geom_point(data = filter(d, !is_edge),
             aes(colour = Test), size = 0.9, alpha = 0.8) +
  ## EDGE lines emphasised: thicker, solid, on top
  geom_line(data = filter(d, is_edge),
            aes(colour = Test, linetype = Test, group = Test), linewidth = 1.15) +
  geom_point(data = filter(d, is_edge),
             aes(colour = Test), size = 1.5) +
  scale_color_ek() + scale_lty_ek() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = label_percent(accuracy = 1),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  facet_wrap(~ family, scales = "free_x", nrow = 1) +
  labs(
    title = "Power versus misspecification severity (n = 1000)",
    subtitle = "Rejection rate vs. misspecification severity  ·  n = 1000, G = 10, α = 0.05  ·  shaded band = EDGE advantage over HL",
    x = "Misspecification severity (family-specific scale)",
    y = "Power  (rejection rate)"
  ) +
  ## Two rows: ten entries on one row overflow the canvas and clip the last
  ## label ("Xie") at the right edge.
  guides(colour = guide_legend(nrow = 2, override.aes = list(
    linewidth = c(1.15,1.15,1.15, rep(0.55,7)), size = 1.2))) +
  theme_ek() +
  theme(
    plot.subtitle = element_text(size = rel(0.82), hjust = 0.5, colour = "grey30", margin = margin(b = 6)),
    legend.key.width = unit(16, "pt"),
    legend.text = element_text(size = rel(0.78)),
    panel.spacing.x = unit(11, "pt")
  )

## boxed key-cell callout on the quad panel (EDGE clears HL early), placed in the
## empty lower-right of the panel with a segment pointing to the severity-0.02 gap.
callout_q <- transform(callout, family = factor(fam_labels["quad"], levels = fam_labels))
p <- p +
  geom_segment(
    data = callout_q,
    aes(x = 0.13, xend = param + 0.004, y = 0.36, yend = (ymin + ymax) / 2),
    inherit.aes = FALSE, colour = ek_pal[["EDGE-poly2"]],
    linewidth = 0.35, alpha = 0.8
  ) +
  geom_label(
    data = callout_q,
    aes(x = 0.13, y = 0.28, label = lab),
    inherit.aes = FALSE, family = EK_FAMILY, size = 2.4, fontface = "bold", lineheight = 0.95,
    colour = ek_pal[["EDGE-poly2"]], fill = "white",
    label.size = 0.3, label.padding = unit(2.2, "pt"), hjust = 0
  )

ggsave(file.path(fig_dir, "Fig2_power_severity.pdf"),
       p, device = cairo_pdf, width = EK_W2, height = 3.75, units = "in")
ggsave(file.path(fig_dir, "_preview_Fig2_power_severity.png"),
       p, width = EK_W2, height = 3.75, units = "in", dpi = 150)

cat("Wrote:\n",
    file.path(fig_dir, "Fig2_power_severity.pdf"), "\n",
    file.path(fig_dir, "_preview_Fig2_power_severity.png"), "\n")
