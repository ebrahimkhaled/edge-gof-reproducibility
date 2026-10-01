## _fig6_surface.R — Fig6: power surface (reject_rate) over severity x sample size,
## faceted EDGE-poly3 vs HL, with the 0.8-power contour highlighted.
## Story: EDGE reaches 80% power at a smaller effect size / smaller n than Hosmer-Lemeshow.
## CAPTION: EDGE hits 80% power for less: a smaller effect size, or a smaller sample.
## FIX (2026-07-07): the fill + the 0.8 contour are now BOTH built from ONE bilinearly
## interpolated surface over (log10 severity) x (log10 n), so the colour varies in BOTH
## directions and the gold 80% line sits on a single constant shade. (The previous
## geom_raster(interpolate=TRUE) on log10 scales mis-tiled -> vertical stripes.)
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

suppressMessages({ library(ggplot2); library(dplyr) })

sim_dir <- edge_path("code/simulations")
fig_dir <- edge_out("figures")
source(file.path(sim_dir, "_ek_theme.R"))
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- data: quad family, DEF.poly3 vs HL, power = reject_rate ------------------------
raw <- read.csv(file.path(sim_dir, "sim_power_broad.csv"), stringsAsFactors = FALSE)
d <- raw %>%
  filter(family == "quad", alpha == 0.05, G == 10, test %in% c("DEF.poly3", "HL")) %>%
  mutate(severity = as.numeric(param)) %>%
  transmute(test, severity, n, power = reject_rate)

sev  <- sort(unique(d$severity)); nn <- sort(unique(d$n))
lsev <- log10(sev);              ln <- log10(nn)
lsev_fine <- seq(min(lsev), max(lsev), length.out = 160)
ln_fine   <- seq(min(ln),   max(ln),   length.out = 160)

## bilinear interpolation of a grid Z (indexed by xs, ys) onto (xo, yo)
bilin <- function(xs, ys, Z, xo, yo) {
  ix <- findInterval(xo, xs, all.inside = TRUE)
  iy <- findInterval(yo, ys, all.inside = TRUE)
  out <- matrix(NA_real_, length(xo), length(yo))
  for (a in seq_along(xo)) {
    i <- ix[a]; tx <- (xo[a] - xs[i]) / (xs[i + 1] - xs[i])
    for (b in seq_along(yo)) {
      j <- iy[b]; ty <- (yo[b] - ys[j]) / (ys[j + 1] - ys[j])
      out[a, b] <- (1 - tx) * (1 - ty) * Z[i, j]   + tx * (1 - ty) * Z[i + 1, j] +
                   (1 - tx) *      ty  * Z[i, j + 1] + tx *      ty  * Z[i + 1, j + 1]
    }
  }
  out
}

build_surface <- function(tt) {
  sub <- d[d$test == tt, ]
  Z <- matrix(NA_real_, length(sev), length(nn))
  for (i in seq_along(sev)) for (j in seq_along(nn)) {
    v <- sub$power[sub$severity == sev[i] & sub$n == nn[j]]
    if (length(v)) Z[i, j] <- v[1]
  }
  stopifnot(!any(is.na(Z)))                     # grid must be complete
  Ffine <- bilin(lsev, ln, Z, lsev_fine, ln_fine)
  expand.grid(lsev = lsev_fine, ln = ln_fine) |>
    transform(power = as.vector(Ffine), test = tt)
}
surf <- rbind(build_surface("DEF.poly3"), build_surface("HL"))
surf$Test <- factor(surf$test, levels = c("DEF.poly3", "HL"),
                    labels = c("EDGE-poly3  (our test)", "Hosmer–Lemeshow"))

## ---- plot (linear axes in log-space; labels show the original severity/n) -----------
p <- ggplot(surf, aes(lsev, ln, fill = power)) +
  geom_raster() +
  geom_contour(aes(z = power), breaks = c(0.2, 0.4, 0.6),
               colour = "white", linewidth = 0.25, alpha = 0.55) +
  geom_contour(aes(z = power), breaks = 0.8, colour = "black",   linewidth = 1.15) +
  geom_contour(aes(z = power), breaks = 0.8, colour = "#FFD24C", linewidth = 0.5) +
  facet_wrap(~ Test) +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#D0E1F2", "#94C4DF", "#4A98C9", "#1F6FB2", "#08519C", "#083A7A"),
    limits = c(0, 1), oob = scales::squish, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    name = "Power  ",
    guide = guide_colourbar(barwidth = 12, barheight = 0.5, title.vjust = 1, ticks.colour = "grey30")) +
  scale_x_continuous(breaks = lsev, labels = formatC(sev, format = "g"), expand = c(0, 0)) +
  scale_y_continuous(breaks = ln,   labels = scales::comma(nn),          expand = c(0, 0)) +
  labs(title = "Power surface: effect size × sample size",
       subtitle = "Quadratic link departure · α = 0.05 · G = 10 · black/gold = 80%-power frontier (further left & lower = better)",
       x = "Departure severity  (quadratic coefficient, log scale)",
       y = "Sample size  n  (log scale)") +
  coord_cartesian(expand = FALSE) +
  theme_ek(base_size = 10.5) +
  theme(panel.grid = element_blank(), panel.spacing = unit(10, "pt"),
        plot.subtitle = element_text(size = rel(0.8), hjust = 0.5, colour = "grey30", margin = margin(b = 6)),
        strip.text = element_text(face = "plain", size = rel(0.95)),
        legend.position = "bottom", legend.title = element_text(size = rel(0.85), face = "bold"),
        axis.text.x = element_text(size = rel(0.75)))

## annotate the win directly (n = 2000, severity 0.02), in log10 coords -----------------
mark_edge <- data.frame(lsev = log10(0.02), ln = log10(2000),
                        Test = factor("EDGE-poly3  (our test)", levels = levels(surf$Test)))
mark_hl   <- data.frame(lsev = log10(0.02), ln = log10(2000),
                        Test = factor("Hosmer–Lemeshow", levels = levels(surf$Test)))
callout   <- data.frame(lsev = log10(0.02), ln = log10(3600),
                        Test = factor("EDGE-poly3  (our test)", levels = levels(surf$Test)),
                        lab = "at n = 2000, severity 0.02:\nEDGE 82%  vs  HL 69%")
## place each "80% power" tag just BELOW-LEFT of its own frontier (in the lighter,
## sub-80% region) so the black/gold contour stays uncovered.
lab80     <- data.frame(lsev = c(log10(0.135), log10(0.155)),
                        ln   = c(log10(235),   log10(235)),
                        Test = factor(levels(surf$Test), levels = levels(surf$Test)),
                        lab = "80% power")
p <- p +
  geom_point(data = mark_edge, aes(lsev, ln), inherit.aes = FALSE, shape = 21, size = 3.4,
             colour = "black", fill = "#FFD24C", stroke = 0.9) +
  geom_point(data = mark_hl, aes(lsev, ln), inherit.aes = FALSE, shape = 21, size = 3.4,
             colour = "black", fill = "white", stroke = 0.9) +
  geom_label(data = callout, aes(lsev, ln, label = lab), inherit.aes = FALSE, size = 2.55,
             hjust = 0.5, lineheight = 0.95, label.size = 0, fill = "white", alpha = 0.85,
             family = EK_FAMILY) +
  geom_label(data = lab80, aes(lsev, ln, label = lab), inherit.aes = FALSE, size = 2.5,
             colour = "black", label.size = 0, fill = "#FFD24C", alpha = 0.9,
             family = EK_FAMILY)

## ---- render ------------------------------------------------------------------------
pdf_out <- file.path(fig_dir, "Fig6_surface.pdf")
png_out <- file.path(fig_dir, "_preview_Fig6_surface.png")
ggsave(pdf_out, p, device = cairo_pdf, width = EK_W2, height = 4.6, units = "in")
ggsave(png_out, p, width = EK_W2, height = 4.6, units = "in", dpi = 150)
cat("WROTE", pdf_out, "|", png_out, "\n")
cat("PDF bytes:", file.info(pdf_out)$size, " PNG bytes:", file.info(png_out)$size, "\n")
