## _fig11_belt_p3.R -- Figure 6 of EDGE paper 3: _fig11_belt.R on the paper-3 cohort (runI_p3_calcurve.csv), written to paper3/figures.
## The same risk groups Figure 10 draws -- ten of them, hence "deciles", a word this script
## derives from G rather than assumes -- now with the null bands the test itself implies.
## The predictions were frozen on the development half, so nothing was estimated from these
## outcomes and Omega = I is the EXACT reference, not a conservative fallback: the standardized
## residual r_g is a standard normal under a calibrated model, and the only thing left to decide
## is the multiplicity of having looked at ten groups.
##
## Input: runI_p3_calcurve.csv, written by run_I_bigdata.R (columns decile, mean_pred,
## obs_rate, n, r_g, gap_pts). V_g is recovered exactly from r_g = (O_g - E_g)/sqrt(V_g).
## Foundation: _ek_theme.R.
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
FIG <- edge_path("manuscript/figures")
source(file.path(SIM, "_ek_theme.R"))
suppressMessages({ library(ggplot2); library(patchwork) })

cc <- read.csv(edge_path("results/cohort/runI_p3_calcurve.csv"))
G  <- nrow(cc)
LEVEL <- 0.95

## --- the pieces the belt needs, all recovered from what run_I wrote ---------------------
O  <- cc$obs_rate  * cc$n
E  <- cc$mean_pred * cc$n
sV <- (O - E) / cc$r_g                       # sqrt(V_g), exactly
cc$z       <- cc$r_g                         # Omega = I: the null SD is one
cc$rate_sd <- sV / cc$n                      # the same band, on the rate scale

## Pointwise, and simultaneous over the G groups. Under Omega = I the residuals are independent,
## so the simultaneous constant is closed-form and needs no simulation.
crit_p <- qnorm(1 - (1 - LEVEL) / 2)
crit_s <- qnorm(1 - (1 - LEVEL^(1 / G)) / 2)

## The word for one group depends on how many there are: ten of them are deciles, and anything
## else is not. Deriving it from G stops the labels lying if this file is ever regenerated at a
## different partition.
grp_word <- if (G == 10) "decile" else if (G == 4) "quartile" else if (G == 5) "quintile" else
            if (G == 20) "vigintile" else "risk group"

## The direction EDGE tests in external mode: the constant is in the basis (nothing removed it),
## with the three polynomial shapes on the group mean risk -- the four degrees of freedom the
## paper reports for this half.
B   <- cbind(1, stats::poly(cc$mean_pred, 3))
B   <- B / rep(sqrt(colSums(B^2)), each = nrow(B))
cc$direction <- drop(B %*% solve(crossprod(B), crossprod(B, cc$r_g)))
cc$flag <- abs(cc$z) > crit_s

col_e  <- ek_pal[["EDGE-poly3"]]          # the direction curve, as everywhere in this paper
col_lo <- "#0072B2"                        # observed below predicted, as in Figure 10
col_hi <- "#D55E00"

band <- data.frame(x = range(cc$mean_pred) + c(-0.006, 0.006))

## ---- left: the belt -------------------------------------------------------------------
pL <- ggplot(cc, aes(mean_pred, z)) +
  geom_ribbon(data = band, inherit.aes = FALSE,
              aes(x = x, ymin = -crit_s, ymax = crit_s), fill = "grey88") +
  geom_hline(yintercept = c(-crit_p, crit_p), linetype = "22", colour = "grey45", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "grey25", linewidth = 0.4) +
  geom_line(aes(y = direction), colour = col_e, linewidth = 0.9) +
  geom_point(aes(shape = flag, fill = flag), size = 2.6, colour = "white", stroke = 0.5) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), name = NULL,
                     labels = c(`FALSE` = "within the simultaneous band", `TRUE` = "beyond the simultaneous band")) +
  scale_fill_manual(values = c(`FALSE` = "grey25", `TRUE` = col_lo), name = NULL,
                    labels = c(`FALSE` = "within the simultaneous band", `TRUE` = "beyond the simultaneous band")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  annotate("text", x = max(cc$mean_pred), y = -crit_s - 0.55, hjust = 1, size = 2.5,
           family = EK_FAMILY, colour = "grey30",
           label = sprintf("simultaneous 95%% band, %d groups: %.2f", G, crit_s)) +
  annotate("text", x = max(cc$mean_pred), y = crit_p + 0.38, hjust = 1, size = 2.5,
           family = EK_FAMILY, colour = "grey45", label = "pointwise 95% band: 1.96") +
  labs(title = "(a) Residual scale",
       x = sprintf("Mean predicted risk in %s", grp_word),
       y = expression(paste("standardized residual  ", r[g]))) +
  theme_ek() + theme(legend.position = "bottom", legend.text = element_text(size = rel(0.8)))

## ---- right: the same band, read on the readmission-rate scale ---------------------------
cc$lo <- cc$mean_pred - crit_s * cc$rate_sd
cc$hi <- cc$mean_pred + crit_s * cc$rate_sd
pR <- ggplot(cc, aes(mean_pred, obs_rate)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey88") +
  ek_diag() +
  geom_segment(aes(xend = mean_pred, yend = mean_pred,
                   colour = ifelse(obs_rate > mean_pred, "above", "below")),
               linewidth = 1.6, alpha = 0.5) +
  geom_point(aes(shape = flag, fill = flag), size = 2.6, colour = "white", stroke = 0.5) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") +
  scale_fill_manual(values = c(`FALSE` = "grey25", `TRUE` = col_lo), guide = "none") +
  scale_colour_manual(values = c(above = col_hi, below = col_lo), guide = "none") +
  coord_cartesian(xlim = range(cc$mean_pred) + c(-0.008, 0.008)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "(b) Rate scale",
       x = sprintf("Mean predicted risk in %s", grp_word), y = "Observed readmission rate") +
  theme_ek()

p <- pL + pR + plot_layout(widths = c(1, 1))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIG, "Fig11_belt.pdf"), p, width = EK_W2, height = EK_W2 * 0.42, device = cairo_pdf)
ggsave(file.path(SIM, "_preview_Fig11_p3.png"), p, width = EK_W2, height = EK_W2 * 0.42, dpi = 160)

cat(sprintf("G = %d | pointwise +/-%.2f | simultaneous +/-%.2f\n", G, crit_p, crit_s))
print(data.frame(decile = cc$decile, mean_pred = round(cc$mean_pred, 4),
                 obs_rate = round(cc$obs_rate, 4), z = round(cc$z, 2),
                 beyond_pointwise = abs(cc$z) > crit_p, beyond_simultaneous = cc$flag),
      row.names = FALSE)
cat("written: figures/Fig11_belt.pdf\n")
