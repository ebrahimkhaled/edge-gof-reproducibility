## _ek_theme.R — one shared, colourblind-safe (Okabe-Ito), greyscale-legible figure system
## for the EDGE paper. source() at the top of every plotting script. Display keys are the
## paper's names (EDGE-poly2/3, EDGE-stk); plotting scripts recode CSV test names
## (DEF.poly2 -> EDGE-poly2, HL-equalwidth -> HLw, Pigeon-Heyse -> PH, DEF.stukel -> EDGE-stk).
## base_family="sans" => Arial on Windows (Springer Helvetica/Arial requirement, A8).
suppressMessages({ library(ggplot2); library(scales) })

ek_pal <- c("EDGE-poly2" = "#E69F00", "EDGE-poly3" = "#D55E00", "EDGE-stk" = "#CC79A7",
            "HL-F" = "#0072B2", "HL" = "#56B4E9", "HLw" = "#80B1D3", "PH" = "#009E73",
            "Tsiatis" = "#7F7F7F", "Stukel" = "#000000", "Xie" = "#B4A0C7", "PR" = "#BBBBBB")
ek_lty <- c("EDGE-poly2" = "solid", "EDGE-poly3" = "solid", "EDGE-stk" = "solid", "Stukel" = "solid",
            "HL-F" = "22", "HL" = "42", "HLw" = "42", "PH" = "42", "Tsiatis" = "42", "Xie" = "42", "PR" = "42")
ek_levels <- c("EDGE-poly3", "EDGE-poly2", "EDGE-stk", "Stukel", "HL-F", "HL", "HLw", "PH", "Tsiatis", "Xie", "PR")

## map raw CSV Test names -> display keys above
ek_relabel <- function(x) {
  m <- c("DEF.poly2" = "EDGE-poly2", "DEF.poly3" = "EDGE-poly3", "DEF.stukel" = "EDGE-stk",
         "EF" = "HL-F", "HL-equalwidth" = "HLw", "Pigeon-Heyse" = "PH", "Pulkstenis-Robinson" = "PR")
  x <- as.character(x); ifelse(x %in% names(m), m[x], x)
}
ek_factor <- function(x) factor(ek_relabel(x), levels = ek_levels)

scale_color_ek <- function(...) scale_color_manual(values = ek_pal, breaks = ek_levels, name = NULL, ...)
scale_lty_ek   <- function(...) scale_linetype_manual(values = ek_lty, breaks = ek_levels, guide = "none", ...)
scale_fill_alpha_ek <- function(mid = 0.05, lim = c(0, 0.12))
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = mid,
                       limits = lim, oob = scales::squish, name = expression(hat(alpha)))

## Springer requires figure lettering in Helvetica/Arial (sans serif) at 8-12 pt
## final size, and forbids titles/captions inside the artwork. Both are opt-in so
## the Wiley/SiM figure set, which has neither restriction, is untouched:
##   EK_FIG_FAMILY=sans   -> sans-serif lettering (Arial on Windows)
##   EK_FIG_NOTITLE=1     -> strip plot title and subtitle at save time
EK_FAMILY <- { v <- Sys.getenv("EK_FIG_FAMILY", ""); if (nzchar(v)) v else "serif" }

theme_ek <- function(base_size = 10.5, base_family = EK_FAMILY) theme_minimal(base_size, base_family) %+replace%
  theme(plot.title = element_text(face = "plain", size = rel(1.05), hjust = 0.5, margin = margin(b = 4)),
        plot.subtitle = element_text(size = rel(0.82), hjust = 0.5, colour = "grey35"),
        plot.title.position = "plot", panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "plain", size = rel(0.9)), legend.position = "bottom")

ek_nominal <- function(a = 0.05) geom_hline(yintercept = a, linetype = "13", colour = "grey55", linewidth = 0.4)
ek_diag    <- function() geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey55", linewidth = 0.4)

## Springer figure widths (mm): single-column 84 (EK_W1), full/double-column text area 174 (EK_W2) (A8).
EK_W1 <- 84 / 25.4; EK_W2 <- 174 / 25.4; EK_DPI <- 320

## Per-venue width override, in MILLIMETRES, via the environment.
##
## 174 mm is Springer's DOUBLE-column text area. The ADAC manuscript is typeset
## single-column with \textwidth = 372 pt = 131 mm, so a figure authored at 174 mm
## is scaled to 0.75 when included at \textwidth and every label inside it shrinks
## by 25% -- which is what pushed the ADAC figure text to 4.5-6 pt. Re-running with
##   EK_FIG_W2_MM=131
## emits at the true column width so the figure is placed 1:1 and the type prints
## at the size it was authored. Wiley's USG text width is 506 pt = 178 mm, close
## enough to 174 mm that the default needs no override there.
.ek_mm <- function(var, default_in) {
  v <- Sys.getenv(var, "")
  if (nzchar(v)) as.numeric(v) / 25.4 else default_in
}
EK_W1 <- .ek_mm("EK_FIG_W1_MM", EK_W1)
EK_W2 <- .ek_mm("EK_FIG_W2_MM", EK_W2)
message(sprintf("[ek_theme] figure widths: EK_W1=%.2f in, EK_W2=%.2f in", EK_W1, EK_W2))

## Strip titles/subtitles at save time when EK_FIG_NOTITLE is set. Done by masking
## ggsave rather than editing eight plot scripts: the scripts add their own
## theme(plot.subtitle=...) AFTER theme_ek(), so blanking inside the theme would be
## overridden. Masking the writer catches every figure regardless of how it was built.
if (nzchar(Sys.getenv("EK_FIG_NOTITLE"))) {
  .ek_real_ggsave <- ggplot2::ggsave
  ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
    if (inherits(plot, "ggplot")) plot <- plot + ggplot2::labs(title = NULL, subtitle = NULL)
    .ek_real_ggsave(filename, plot, ...)
  }
  message("[ek_theme] EK_FIG_NOTITLE set: plot titles/subtitles will be stripped on save")
}
