## _fig1_calibration.R — Figure 1: the calibration-curve signature.
## Self-generate ONE large draw (n=20000) from a cloglog truth fitted with logit, group into
## deciles, plot group observed rate vs mean predicted prob, smooth line+points, identity diag,
## shaded ribbon for the gap, annotate the largest-gap region. ONE panel, width EK_W1.
## Caption: "Under a misspecified link the calibration curve bends smoothly off the diagonal —
##  a low-dimensional, directed signal that a projection test captures."
## CAPTION: A misspecified link bends calibration smoothly off the diagonal.
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

suppressMessages({ library(ggplot2) })
sim_dir <- edge_path("code/simulations")
fig_dir <- edge_out("figures")
source(file.path(sim_dir, "_ek_theme.R"))
source(file.path(sim_dir, "_dgp_library.R"))
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. one large misspecified draw: cloglog truth, fit logit ------------------------------
set.seed(20260706)
G0  <- dgp_alt("link", "cloglog", 20000)          # y ~ x + d, truth = cloglog(0.6 x + 0.5 d)
fit <- suppressWarnings(glm(G0$f, data = G0$d, family = binomial()))   # fit the WRONG link (logit)
ph  <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6)
y   <- G0$d$y

## ---- 2. decile grouping on predicted prob; observed vs predicted per group -----------------
G   <- 10
grp <- pmin(ceiling(rank(ph, ties.method = "first") / (length(ph) / G)), G)
agg <- data.frame(
  grp   = seq_len(G),
  pred  = as.numeric(tapply(ph, grp, mean)),           # mean predicted prob (x)
  obs   = as.numeric(tapply(y,  grp, mean)),           # observed rate (y)
  n_g   = as.numeric(tapply(y,  grp, length))
)
agg$se   <- sqrt(agg$obs * (1 - agg$obs) / agg$n_g)     # binomial SE of observed rate
agg$gap  <- agg$obs - agg$pred                          # signed calibration gap
## largest-gap decile (by |gap|) — the region we annotate
jmax <- which.max(abs(agg$gap))
lab_pt <- agg[jmax, ]

## smooth monotone-ish calibration curve through the decile points (for the bending line)
sm <- as.data.frame(spline(agg$pred, agg$obs, n = 300, method = "natural"))
names(sm) <- c("pred", "obs")
sm$pred <- pmin(pmax(sm$pred, 0), 1); sm$obs <- pmin(pmax(sm$obs, 0), 1)

## ribbon of the gap: between the identity (y = pred) and the smooth calibration curve
rib <- data.frame(pred = sm$pred, lo = pmin(sm$pred, sm$obs), hi = pmax(sm$pred, sm$obs))

rng <- range(c(agg$pred, agg$obs)); pad <- 0.04 * diff(rng)
lims <- c(max(0, rng[1] - pad), min(1, rng[2] + pad))

## ---- 3. plot -------------------------------------------------------------------------------
edge_blue <- ek_pal[["EF"]]      # anchor accent (deep blue) for the fitted curve
gap_fill  <- "#D55E00"           # EDGE-poly3 vermillion for the gap band

p <- ggplot() +
  ## shaded gap between diagonal and calibration curve
  geom_ribbon(data = rib, aes(x = pred, ymin = lo, ymax = hi),
              fill = gap_fill, alpha = 0.16) +
  ## identity (perfect calibration)
  ek_diag() +
  ## smooth bending calibration curve
  geom_line(data = sm, aes(x = pred, y = obs), colour = edge_blue, linewidth = 1.05) +
  ## decile points with binomial error bars
  geom_errorbar(data = agg, aes(x = pred, ymin = obs - se, ymax = obs + se),
                width = 0, colour = "grey45", linewidth = 0.4) +
  geom_point(data = agg, aes(x = pred, y = obs),
             shape = 21, fill = edge_blue, colour = "white", stroke = 0.6, size = 2.6) +
  ## highlight the largest-gap decile: a vertical connector from curve to diagonal + ring
  geom_segment(data = lab_pt, aes(x = pred, xend = pred, y = pred, yend = obs),
               colour = gap_fill, linewidth = 0.9, linetype = "solid") +
  geom_point(data = lab_pt, aes(x = pred, y = obs),
             shape = 21, fill = NA, colour = gap_fill, stroke = 1.3, size = 5.4) +
  ## leader line + direct label on the largest-gap region (placed lower-right, clear of curve)
  annotate("segment", x = lab_pt$pred, y = lab_pt$obs,
           xend = lab_pt$pred + 0.14, yend = lab_pt$obs - 0.12,
           colour = gap_fill, linewidth = 0.4) +
  annotate("label", x = lab_pt$pred + 0.14, y = lab_pt$obs - 0.12,
           label = sprintf("largest gap: %+.1f pts", 100 * lab_pt$gap),
           hjust = 0, vjust = 0.5, size = 2.6, fill = "white",
           colour = gap_fill, fontface = "bold", family = EK_FAMILY,
           label.padding = unit(1.4, "pt")) +
  ## label the identity line — placed low on the diagonal and offset just BELOW it,
  ## into the empty lower-right triangle, so it never overlaps the calibration curve
  annotate("text", x = lims[1] + 0.34 * diff(lims), y = lims[1] + 0.34 * diff(lims),
           label = "perfect calibration", angle = 45, vjust = 2.5,
           size = 2.5, colour = "grey45", fontface = "italic", family = EK_FAMILY) +
  coord_equal(xlim = lims, ylim = lims, expand = FALSE, clip = "off") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Calibration curve under a misspecified link",
       subtitle = expression("cloglog truth fitted with logit  ("*italic(n)==20000*", deciles)"),
       x = "mean predicted probability", y = "observed event rate") +
  theme_ek(base_size = 9) +
  theme(plot.title = element_text(face = "plain", size = rel(1.05), hjust = 0.5,
                                  lineheight = 0.98, margin = margin(b = 3)),
        plot.subtitle = element_text(size = rel(0.86), hjust = 0.5, colour = "grey35", margin = margin(b = 6)),
        plot.margin = margin(6, 14, 6, 6))

## ---- 4. render vector PDF + PNG preview -----------------------------------------------------
pdf_out <- file.path(fig_dir, "Fig1_calibration.pdf")
png_out <- file.path(fig_dir, "_preview_Fig1_calibration.png")
H <- EK_W1 * 1.12   # near-square panel + headroom for the two-line title
ggsave(pdf_out, p, device = cairo_pdf, width = EK_W1, height = H, units = "in")
ggsave(png_out, p, width = EK_W1, height = H, units = "in", dpi = 150)

cat("gap by decile (obs-pred, pts):\n"); print(round(100 * agg$gap, 2))
cat(sprintf("largest gap: decile %d, pred=%.3f obs=%.3f gap=%+.1f pts\n",
            jmax, lab_pt$pred, lab_pt$obs, 100 * lab_pt$gap))
cat("PDF:", pdf_out, "\nPNG:", png_out, "\n")
cat("PDF KB:", round(file.info(pdf_out)$size / 1024, 1),
    " PNG KB:", round(file.info(png_out)$size / 1024, 1), "\n")
