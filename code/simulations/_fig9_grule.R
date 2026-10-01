## _fig9_grule.R -- Figure 9: "the group count works in opposite directions for a directed and
## an omnibus partition test." Size-adjusted power against observations-per-group m = n/G, for
## two departures (cloglog bow, probit S) at n = 1000 (cloglog) and n = 5000 (probit; at n =
## 1000 nothing has power against the probit), from run F (runF_grule_summary.csv).
## Directed tests RISE as m falls (finer partition); HL FALLS to its size. GiViTI and Stukel do
## not partition and are drawn flat at their own value in the same cell, as the reference the
## reader will ask for. Foundation: _ek_theme.R.
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
FIG <- edge_out("submission_statistics_in_medicine", "figures")
source(file.path(SIM, "_ek_theme.R"))
suppressMessages({ library(dplyr); library(ggplot2) })

S <- read.csv(file.path(SIM, "runF_grule_summary.csv"), stringsAsFactors = FALSE)
S$m_actual <- round(S$n / S$G)          # the n=20000 rows above G=400 are mislabelled by the cap

panels <- bind_rows(
  S %>% filter(scen == "cloglog", n == 1000) %>% mutate(panel = "cloglog bow, n = 1000"),
  S %>% filter(scen == "probit",  n == 5000) %>% mutate(panel = "probit S, n = 5000"))

keep <- c("EDGE.poly3", "EDGE.stk", "HL", "HL_F", "GiViTI", "Stukel")
lab  <- c(EDGE.poly3 = "EDGE-poly3", EDGE.stk = "EDGE-stk", HL = "HL", HL_F = "HL-F",
          GiViTI = "GiViTI", Stukel = "Stukel")
D <- panels %>% filter(test %in% keep) %>%
  mutate(series = factor(lab[test], levels = c("EDGE-poly3", "EDGE-stk", "Stukel", "GiViTI", "HL-F", "HL")),
         partitions = !test %in% c("GiViTI", "Stukel"))

cols <- c(ek_pal[c("EDGE-poly3", "EDGE-stk", "Stukel", "HL-F", "HL")], GiViTI = "#009E73")
ltys <- c("EDGE-poly3" = "solid", "EDGE-stk" = "solid", "Stukel" = "31", "GiViTI" = "31",
          "HL-F" = "22", "HL" = "42")

## direct labels at the finest partition (left edge) for the partition tests, and at the
## coarsest (right edge) for the two non-partition references
lab_left  <- D %>% filter(partitions) %>% group_by(panel, series) %>% slice_min(m_actual, n = 1) %>% ungroup()
lab_right <- D %>% filter(!partitions) %>% group_by(panel, series) %>% slice_max(m_actual, n = 1) %>% ungroup()

p <- ggplot(D, aes(m_actual, power_adj, colour = series, linetype = series, group = series)) +
  ek_nominal() +
  geom_text(data = data.frame(panel = "cloglog bow, n = 1000", m_actual = 200, power_adj = 0.075),
            aes(m_actual, power_adj, label = "nominal 0.05"), inherit.aes = FALSE, hjust = 1,
            size = 2.6, colour = "grey45", family = EK_FAMILY) +
  geom_line(aes(linewidth = partitions)) +
  geom_point(data = filter(D, partitions), size = 1.9) +
  ggrepel::geom_text_repel(data = lab_left, aes(label = series), size = 2.7, family = EK_FAMILY,
                           direction = "y", nudge_x = -0.08, hjust = 1, segment.colour = "grey70",
                           min.segment.length = 0, seed = 2, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = lab_right, aes(label = series), size = 2.7, family = EK_FAMILY,
                           direction = "y", nudge_x = 0.08, hjust = 0, segment.colour = "grey70",
                           min.segment.length = 0, seed = 3, show.legend = FALSE) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = ltys, guide = "none") +
  scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.55), guide = "none") +
  scale_x_log10(breaks = c(15, 25, 50, 100, 200), expand = expansion(mult = c(0.22, 0.20))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.10))) +
  labs(x = "Observations per risk group  m = n / G  (log scale; finer partition to the left)",
       y = "Size-adjusted power") +
  theme_ek() +
  theme(panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92"),
        strip.text = element_text(size = rel(0.95)))

dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIG, "Fig9_grule.pdf"), p, device = cairo_pdf, width = EK_W2, height = 3.6, units = "in")
ggsave(file.path(FIG, "_preview_Fig9_grule.png"), p, width = EK_W2, height = 3.6, units = "in", dpi = 150)
cat("wrote Fig9_grule.pdf\n")
