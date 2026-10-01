## _fig3_boundary.R — Fig 3: EDGE's two honest blind spots.
## From sim_edge_loses.csv, alpha=0.05, largest n per scenario, five tests:
##   EDGE-poly3 (the directed protagonist), EF, HL (omnibus), Tsiatis, Xie (covariate-space).
## Story: crossover -> covariate-space tests (Tsiatis/Xie) win, EDGE/EF/HL ~null.
##        osc4 / sawtooth -> omnibus EF/HL (and covariate tests) win, directed EDGE loses.
## A horizontal lollipop, faceted by scenario, EDGE-poly3 pulled out as the protagonist,
## shaded null band + winner-family callouts.
## CAPTION: EDGE's two honest blind spots: covariate-space structure and rough high-frequency misfit.
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
  library(ggplot2); library(dplyr); library(tidyr)
})
here <- edge_path("code/simulations")
figdir <- edge_out("figures")
source(file.path(here, "_ek_theme.R"))
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

raw <- read.csv(file.path(here, "sim_edge_loses.csv"), stringsAsFactors = FALSE)

targets_raw <- c("DEF.poly3", "EF", "HL", "Tsiatis", "Xie")

## largest n per scenario, alpha = 0.05, the five target tests
d <- raw %>%
  filter(alpha == 0.05, test %in% targets_raw) %>%
  group_by(scenario) %>% filter(n == max(n)) %>% ungroup() %>%
  mutate(Test = ek_factor(test),
         Test = droplevels(Test),
         reject = reject_rate,
         se = mcse)

## scenario -> readable label + which family wins + the losing/blind-spot cue
scen_meta <- tibble::tibble(
  scenario = c("crossover", "osc4", "sawtooth"),
  scen_lab = c("Crossover misfit", "Oscillation (4 cycles)", "Sawtooth (high-freq.)"),
  blindspot = c("covariate-space structure", "rough, wiggly misfit", "rough, wiggly misfit"),
  winner    = c("Only covariate-space tests win\n(Tsiatis, Xie)",
                "Everything but EDGE wins\n(HL-F, HL, Tsiatis, Xie)",
                "Everything but EDGE wins\n(HL-F, HL, Tsiatis, Xie)")
)
d <- d %>% left_join(scen_meta, by = "scenario")
d$scen_lab <- factor(d$scen_lab, levels = scen_meta$scen_lab)

## family grouping for the annotation logic (protagonist vs winners vs also-ran)
fam_of <- function(t) dplyr::case_when(
  t == "EDGE-poly3"          ~ "directed",
  t %in% c("HL-F", "HL")     ~ "omnibus",
  t %in% c("Tsiatis", "Xie") ~ "covariate")
d$family <- fam_of(as.character(d$Test))

## order rows within each facet by power so the lollipops read as a ranking,
## but keep EDGE-poly3 always visible; use a per-facet y via reorder_within trick
d <- d %>%
  group_by(scen_lab) %>%
  mutate(yorder = rank(reject, ties.method = "first"),
         ylab = paste0(as.character(Test), "___", scen_lab)) %>%
  ungroup()
# stable factor for y
lev <- d %>% arrange(scen_lab, reject) %>% pull(ylab) %>% unique()
d$ylab <- factor(d$ylab, levels = lev)
ylab_names <- setNames(sub("___.*$", "", levels(d$ylab)), levels(d$ylab))

## highlight EDGE-poly3 (the honest protagonist) with a heavier point + label
d$is_edge <- d$Test == "EDGE-poly3"

## per-facet winner callout position (top-left inside panel)
call_df <- d %>% distinct(scen_lab, winner)

nband <- d %>% distinct(scenario, n) %>% arrange(scenario)

p <- ggplot(d, aes(x = reject, y = ylab)) +
  ## null band 0..0.10 to show "near-null" territory
  annotate("rect", xmin = 0, xmax = 0.10, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.16) +
  ## power reference lines
  geom_vline(xintercept = c(0.5), linetype = "13", colour = "grey80", linewidth = 0.35) +
  ek_nominal(0.05) +           # nominal-alpha dotted vertical? ek_nominal is hline; add manual vline below
  ## lollipop stems from 0
  geom_segment(aes(x = 0, xend = reject, yend = ylab, colour = Test, linetype = Test),
               linewidth = 0.7) +
  geom_point(aes(colour = Test, size = is_edge)) +
  ## value labels at the dot
  geom_text(aes(label = ifelse(reject >= 0.9995, "1.00", sprintf("%.2f", reject))),
            hjust = -0.28, size = 2.6, colour = "grey20", family = EK_FAMILY) +
  ## EDGE protagonist callout: left-anchored above the dot (never clips at the left edge;
  ## when the dot sits far right, anchor from a fixed x so the label stays inside the panel)
  geom_text(
    data = transform(subset(d, is_edge), lx = pmin(reject, 0.30)),
    aes(x = lx, label = "EDGE (directed)"),
    hjust = 0, nudge_y = 0.44, size = 2.55, fontface = "bold",
    colour = ek_pal[["EDGE-poly3"]], family = EK_FAMILY) +
  ## winner-family callout, top-left of each facet
  geom_label(data = call_df, aes(x = 0.02, y = Inf, label = winner),
             hjust = 0, vjust = 1.25, size = 2.35, fontface = "bold",
             label.size = 0, fill = "white", alpha = 0.72, lineheight = 0.95,
             colour = "grey15", inherit.aes = FALSE, family = EK_FAMILY) +
  facet_wrap(~ scen_lab, ncol = 3, scales = "free_y") +
  scale_y_discrete(labels = ylab_names) +
  scale_x_continuous(limits = c(0, 1.22), breaks = seq(0, 1, 0.25),
                     labels = c("0", ".25", ".50", ".75", "1"),
                     expand = expansion(mult = c(0, 0))) +
  scale_color_ek() + scale_lty_ek() +
  scale_size_manual(values = c(`FALSE` = 2.1, `TRUE` = 3.6), guide = "none") +
  labs(
    title = "Where EDGE is not the most powerful test",
    subtitle = "Rejection rate at α = 0.05, largest n per scenario (n = 2000). Shaded band = near-null (≤ 0.10).",
    x = "Power  (rejection rate)", y = NULL) +
  theme_ek(base_size = 10) +
  theme(legend.position = "none",
        panel.spacing.x = unit(12, "pt"),
        plot.title = element_text(face = "plain", size = rel(1.05), hjust = 0.5, margin = margin(b = 3)),
        plot.subtitle = element_text(size = rel(0.82), hjust = 0.5, colour = "grey35", margin = margin(b = 6)),
        strip.text = element_text(face = "plain", size = rel(0.82)),
        axis.text.y = element_text(size = rel(0.95)),
        plot.margin = margin(6, 8, 4, 4),
        panel.grid.major.y = element_blank())

## add a proper nominal-alpha vertical marker (ek_nominal draws a horizontal line; we want vertical)
p <- p + geom_vline(xintercept = 0.05, linetype = "13", colour = "grey55", linewidth = 0.4)

ggsave(file.path(figdir, "Fig3_boundary.pdf"), p,
       device = cairo_pdf, width = EK_W2, height = 4.3, units = "in")
ggsave(file.path(figdir, "_preview_Fig3_boundary.png"), p,
       width = EK_W2, height = 4.3, units = "in", dpi = 150)

cat("done: rows used =", nrow(d), "\n")
print(d %>% select(scen_lab, Test, reject) %>% arrange(scen_lab, -reject))
