## _fig8_stukel.R — Figure 8 "the robustness clincher"
## Under sparsity, Stukel's auxiliary refit separates (fails to converge) in 20-28% of
## samples; EDGE never refits, so it is always computable (0% failure).
## CAPTION: Stukel refuses to compute; EDGE never does.
## Data: bench_stukel_failure.csv. Foundation: _ek_theme.R (source()d, not recreated).
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
FIG <- edge_out("figures")
source(file.path(SIM, "_ek_theme.R"))
suppressMessages({ library(dplyr); library(ggplot2) })
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

## ggrepel is nice-to-have; fall back to geom_text with manual nudges if absent.
HAS_REPEL <- requireNamespace("ggrepel", quietly = TRUE)
repel_label <- function(data, mapping, colour, size, fontface, nudge_x = 0, nudge_y = 0,
                        hjust = 0.5, vjust = 0.5, ...) {
  if (HAS_REPEL) {
    ggrepel::geom_text_repel(data = data, mapping = mapping, colour = colour, size = size,
                             fontface = fontface, nudge_x = nudge_x, nudge_y = nudge_y,
                             seed = 1, segment.colour = "grey70",
                             min.segment.length = 0, direction = "y")
  } else {
    geom_text(data = data, mapping = mapping, colour = colour, size = size,
              fontface = fontface, hjust = hjust, vjust = vjust,
              position = position_nudge(x = nudge_x, y = nudge_y))
  }
}

raw <- read.csv(file.path(SIM, "bench_stukel_failure.csv"), stringsAsFactors = FALSE)

## ---- reshape to a long "who fails" frame ----------------------------------------
## Stukel: refit-separation failure rate. In the SPARSE design at tiny n the rate is either
## NA (n=100,150: Stukel could not be attempted) or estimated from a handful of computable
## fits (n=200 -> only 3 of 5000 computable). Both are statistically undefined for a failure-
## rate curve, so we plot the sparse line only where it is estimable (n_computable >= MIN_NC)
## and flag the undefined low-n region with a note. The moderate design is fully defined.
MIN_NC <- 30
stukel <- raw %>%
  mutate(fail = stukel_fail_rate,
         estimable = design == "moderate" | (!is.na(fail) & n_computable >= MIN_NC)) %>%
  transmute(design, n,
            fail = ifelse(estimable, fail, NA_real_),
            event_rate = event_rate_mean,
            series = ifelse(design == "sparse", "Stukel refit (sparse)", "Stukel refit (moderate)"))

## EDGE contrast: never refits -> ef_fail_rate is 0 throughout. One flat reference line.
edge <- raw %>%
  filter(design == "sparse") %>%          # single row-set is enough for a flat 0 line
  transmute(n, fail = ef_fail_rate,
            series = "EDGE (never refits)")

col_sparse   <- ek_pal[["Stukel"]]        # black  — the failing method, sparse
col_moderate <- ek_pal[["Tsiatis"]]       # grey   — the failing method, moderate
col_edge     <- ek_pal[["EDGE-stk"]]      # magenta — EDGE, the hero flat line

series_lvls <- c("Stukel refit (sparse)", "Stukel refit (moderate)", "EDGE (never refits)")
series_cols <- setNames(c(col_sparse, col_moderate, col_edge), series_lvls)
series_lty  <- setNames(c("solid", "solid", "solid"), series_lvls)
series_shp  <- setNames(c(16, 17, NA), series_lvls)

stukel$series <- factor(stukel$series, levels = series_lvls)
edge$series   <- factor(edge$series,   levels = series_lvls)

## ---- the "20-28%" danger zone (sparse Stukel plateau) ---------------------------
zone_lo <- 0.20; zone_hi <- 0.28
xr <- range(raw$n)

## label anchors (direct labels, no legend clutter) — interior points so text stays in-panel
lab_sparse   <- stukel %>% filter(design == "sparse", n == 2000)
lab_moderate <- stukel %>% filter(design == "moderate", n == 500)
lab_edge     <- edge   %>% filter(n == 500)

## sparse "undefined" region: n where Stukel could not be assessed (NA) or was estimated
## from too few computable fits (n_computable < MIN_NC).
na_gap <- raw %>% filter(design == "sparse",
                         is.na(stukel_fail_rate) | n_computable < MIN_NC)
na_note_x <- if (nrow(na_gap)) exp(mean(log(range(na_gap$n)))) else NA
na_note_lab <- if (nrow(na_gap))
  sprintf("Stukel undefined\n(n ≤ %d, sparse)", max(na_gap$n)) else ""

## ---- plot ------------------------------------------------------------------------
p <- ggplot() +
  ## danger band
  annotate("rect", xmin = xr[1], xmax = xr[2], ymin = zone_lo, ymax = zone_hi,
           fill = col_sparse, alpha = 0.07) +
  annotate("segment", x = xr[1], xend = xr[2], y = zone_lo, yend = zone_lo,
           linetype = "13", colour = col_sparse, linewidth = 0.3) +
  annotate("segment", x = xr[1], xend = xr[2], y = zone_hi, yend = zone_hi,
           linetype = "13", colour = col_sparse, linewidth = 0.3) +
  annotate("text", x = xr[1] * 1.05, y = (zone_lo + zone_hi) / 2,
           label = "20–28% of samples\nlost under sparsity", hjust = 0, vjust = 0.5,
           size = 2.9, colour = col_sparse, fontface = "bold", lineheight = 0.95,
           family = EK_FAMILY) +

  ## EDGE hero line at 0 — flat, always computable
  geom_line(data = edge, aes(n, fail, colour = series, linetype = series),
            linewidth = 1.1) +

  ## Stukel failure curves
  geom_line(data = filter(stukel, !is.na(fail)),
            aes(n, fail, colour = series, linetype = series, group = series),
            linewidth = 0.9) +
  geom_point(data = filter(stukel, !is.na(fail)),
             aes(n, fail, colour = series, shape = series), size = 2.2) +

  ## one direct label on the EDGE hero line (the two Stukel series are named in the legend)
  geom_label(data = lab_edge, aes(n, fail, label = "EDGE = 0% failure (never refits)"),
             colour = col_edge, fontface = "bold", size = 3.0, family = EK_FAMILY,
             linewidth = 0, fill = "white", hjust = 0.5, vjust = -0.5, alpha = 0.95) +

  scale_colour_manual(values = series_cols, breaks = series_lvls, name = NULL) +
  scale_linetype_manual(values = series_lty, guide = "none") +
  scale_shape_manual(values = series_shp, guide = "none", na.translate = FALSE) +
  scale_x_log10(breaks = c(100, 200, 500, 1000, 2000, 5000),
                labels = c("100", "200", "500", "1k", "2k", "5k")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.36), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Auxiliary-refit failure under sparsity",
       subtitle = "Sparse design (event rate ≈ 1.5%) vs. moderate design",
       x = "Sample size  n  (log scale)", y = "Test-not-computable rate") +
  theme_ek() +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = rel(0.9), hjust = 0.5, colour = "grey30", margin = margin(b = 6)),
        panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92")) +
  guides(colour = guide_legend(override.aes = list(
           linewidth = c(0.9, 0.9, 1.1), shape = c(16, 17, NA))))

## ---- render both ----------------------------------------------------------------
pdf_path <- file.path(FIG, "Fig8_stukel.pdf")
png_path <- file.path(FIG, "_preview_Fig8_stukel.png")
ggsave(pdf_path, p, device = cairo_pdf, width = EK_W2, height = 4.4, units = "in")
ggsave(png_path, p, width = EK_W2, height = 4.4, units = "in", dpi = 150)

cat("Wrote:\n  ", pdf_path, "\n  ", png_path, "\n")
cat("PDF size (KB): ", round(file.info(pdf_path)$size / 1024, 1), "\n")
cat("PNG size (KB): ", round(file.info(png_path)$size / 1024, 1), "\n")
