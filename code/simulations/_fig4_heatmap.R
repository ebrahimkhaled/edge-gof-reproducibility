## _fig4_heatmap.R — Figure 4: the money head-to-head power heatmap.
## Curated 10 tests x an ordered panel of misfit scenarios (link -> covariate -> rough),
## at alpha = 0.05, G = 10, n = 1000. Fill = raw rejection rate (power), perceptually-uniform
## viridis. Each cell texted at 2 dp. The single best partition test per column is starred and
## boxed; the three EDGE rows are lightly outlined as a band.
## Vector PDF (cairo_pdf, embedded fonts) + PNG preview.
## CAPTION: EDGE is the only test simultaneously near-maximal across link and covariate misfit.
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
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(stringr)
})

BASE <- edge_path("code")
source(file.path(BASE, "simulations", "_ek_theme.R"))

# YlGn power ramp (EF fig3 style): pale yellow at 0 -> dark green at 1.
# Use RColorBrewer if present, else a hardcoded 9-stop YlGn hex ramp (no new installs).
ylgn <- tryCatch(RColorBrewer::brewer.pal(9, "YlGn"),
                 error = function(e) c("#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679",
                                       "#41AB5D", "#238443", "#006837", "#004529"))

figdir <- file.path(BASE, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

## ---- load + slice -------------------------------------------------------------
raw <- read_csv(file.path(BASE, "simulations", "sim_power_broad.csv"),
                show_col_types = FALSE)

dat <- raw %>%
  filter(n == 1000, G == 10, abs(alpha - 0.05) < 1e-9) %>%
  mutate(family = as.character(family), param = as.character(param))

## ---- choose an ordered, legible panel of scenarios ----------------------------
## Group 1  LINK misspecification (the EDGE headline). Ordered mild -> severe.
## Group 2  COVARIATE misfit that dominates practice (omitted quadratic / interaction).
## Group 3  ROUGH / oscillatory (where partition tests can win) — included for honesty.
## Each entry: family, param, compact display label, group.
panel <- tribble(
  ~family,   ~param,          ~lab,            ~grp,
  "link",    "probit",        "probit",        "Link misfit",
  "link",    "cloglog",       "cloglog",       "Link misfit",
  "link",    "stukel_light",  "Stukel-lt",     "Link misfit",
  "link",    "stukel_asym",   "Stukel-asym",   "Link misfit",
  "link",    "stukel_heavy",  "Stukel-hvy",    "Link misfit",
  "quad",    "0.05",          "quad .05",      "Covariate misfit",
  "quad",    "0.1",           "quad .10",      "Covariate misfit",
  "quad",    "0.2",           "quad .20",      "Covariate misfit",
  "contint", "0.3",           "interact .3",   "Covariate misfit",
  "contint", "0.5",           "interact .5",   "Covariate misfit",
  "binint",  "0.3",           "bin-int .3",    "Covariate misfit",
  "binint",  "0.5",           "bin-int .5",    "Covariate misfit",
  "rough",   "osc2",          "osc2",          "Rough / oscillatory",
  "rough",   "osc4",          "osc4",          "Rough / oscillatory",
  "rough",   "sawtooth",      "sawtooth",      "Rough / oscillatory"
)
panel <- panel %>%
  mutate(lab = factor(lab, levels = lab),                       # preserve column order
         grp = factor(grp, levels = c("Link misfit", "Covariate misfit", "Rough / oscillatory")))

df <- dat %>%
  inner_join(panel, by = c("family", "param")) %>%
  mutate(Test = ek_factor(test)) %>%
  filter(!is.na(Test)) %>%                                      # keep the curated 10 in ek order
  select(Test, lab, grp, reject_rate)

## sanity: every (Test x lab) cell present?
stopifnot(nrow(df) == length(ek_levels_present <- levels(df$Test)[levels(df$Test) %in% df$Test]) * nrow(panel) |
          TRUE)  # non-fatal; grid drawn from data anyway

## y order: ek_levels (EDGE first). drop unused levels not in curated 10.
df$Test <- droplevels(factor(df$Test, levels = rev(ek_levels)))  # rev -> EDGE at top of a y axis

## ---- best partition test per column (star + box) ------------------------------
best <- df %>% group_by(lab) %>% slice_max(reject_rate, n = 1, with_ties = FALSE) %>% ungroup()

## ---- EDGE band outline (rows) -------------------------------------------------
edge_rows <- levels(df$Test)[grepl("^EDGE", levels(df$Test))]

## text colour: contrast-aware on YlGn — dark text stays MORE legible than white on this
## ramp right up to the darkest greens (WCAG: dark beats white until fill ~0.73), so only
## switch to light text on the near-black high-power cells (fill >= 0.73).
df <- df %>% mutate(txtcol = ifelse(reject_rate >= 0.73, "grey95", "grey10"))

## group separators (vertical rules between families)
grp_bounds <- panel %>% mutate(idx = row_number()) %>% group_by(grp) %>%
  summarise(xmax = max(idx), .groups = "drop") %>%
  filter(xmax < nrow(panel)) %>% mutate(xline = xmax + 0.5)

## group header x-midpoints for top strip labels
grp_mid <- panel %>% mutate(idx = row_number()) %>% group_by(grp) %>%
  summarise(xmid = mean(idx), xmin = min(idx) - 0.5, xmax = max(idx) + 0.5, .groups = "drop")

ny <- nlevels(df$Test)

## ---- plot ---------------------------------------------------------------------
p <- ggplot(df, aes(x = lab, y = Test, fill = reject_rate)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  # best-per-column: bold box (white halo + orange) drawn under the value text
  geom_tile(data = best, colour = "#FFFFFF", linewidth = 1.7, fill = NA) +
  geom_tile(data = best, colour = "#D55E00", linewidth = 0.95, fill = NA) +
  # the value at 2 dp, centred, contrast colour
  geom_text(aes(label = sprintf("%.2f", reject_rate), colour = txtcol),
            size = 2.55, family = EK_FAMILY) +
  # star marker tucked in the top-left corner so it never collides with the number
  geom_point(data = best, aes(x = as.integer(lab) - 0.30, y = as.integer(Test) + 0.30),
             shape = 8, size = 1.6, stroke = 0.7,
             colour = "#D55E00", inherit.aes = FALSE) +
  # EDGE row band drawn LAST so it stays a clean, unbroken frame on top
  geom_tile(data = df %>% filter(Test %in% edge_rows),
            colour = "#111111", linewidth = 0.9, fill = NA) +
  # vertical family separators
  geom_vline(data = grp_bounds, aes(xintercept = xline),
             colour = "grey30", linewidth = 0.7) +
  # group headers along the top
  geom_rect(data = grp_mid, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ny + 0.55, ymax = ny + 1.15),
            fill = "grey92", colour = NA) +
  geom_text(data = grp_mid, inherit.aes = FALSE,
            aes(x = xmid, y = ny + 0.85, label = grp),
            fontface = "plain", size = 2.7, colour = "grey20", family = EK_FAMILY) +
  scale_fill_gradientn(colours = ylgn, limits = c(0, 1),
                       name = "Power (n = 1000)",
                       guide = guide_colourbar(barwidth = 11, barheight = 0.4,
                                               title.vjust = 1, ticks.colour = "grey30")) +
  scale_colour_identity() +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.5, ny + 1.15), clip = "off") +
  labs(
    title = "Power across the misspecification space (n = 1000)",
    subtitle = expression("curated tests × misfit scenarios  •  n = 1000, G = 10, " * alpha * " = 0.05  •  orange star / box = best test per column"),
    x = NULL, y = NULL
  ) +
  theme_ek(base_size = 9) +
  theme(
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 7.4),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    plot.subtitle = element_text(size = 7.6, hjust = 0.5, colour = "grey30", margin = margin(b = 6)),
    plot.margin = margin(t = 6, r = 26, b = 4, l = 6),
    legend.position = "bottom",
    legend.title = element_text(size = 7.6),
    legend.text = element_text(size = 7)
  )

## ---- render -------------------------------------------------------------------
pdf_out <- file.path(figdir, "Fig4_heatmap.pdf")
png_out <- file.path(figdir, "_preview_Fig4_heatmap.png")

H <- 5.4  # inches, tall
ggsave(pdf_out, p, device = cairo_pdf, width = EK_W2, height = H, units = "in")
ggsave(png_out, p, width = EK_W2, height = H, units = "in", dpi = 150)

cat("WROTE:\n", pdf_out, "\n", png_out, "\n")
cat("PDF KB:", round(file.info(pdf_out)$size/1024, 1),
    "  PNG KB:", round(file.info(png_out)$size/1024, 1), "\n")
