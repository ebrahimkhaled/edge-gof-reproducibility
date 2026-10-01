## _fig10_diabetes_p3.R -- Figure 5 of EDGE paper 3: the Diabetes-130 cohort, one admission per patient.
## _fig10_diabetes.R with the paper-3 cohort (runI_p3_*.csv from run_I_bigdata_p3.R and run_I_p3_calslope.R).
## Every number printed in the figure is read from those files; nothing is typed in.
## Left: validation-half calibration curve by decile of predicted risk (frozen model), the gap from the diagonal
## drawn as bars. Right: development half, the p-value of the directed default and of the Hosmer-Lemeshow tests as
## the partition is refined from G = 10 to G = 2000.
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
suppressMessages({ library(dplyr); library(ggplot2); library(patchwork) })

cc  <- read.csv(edge_path("results/cohort/runI_p3_calcurve.csv"))
dev <- read.csv(edge_path("results/cohort/runI_p3_dev.csv"))
sm  <- read.csv(edge_path("results/cohort/runI_p3_summary.csv"))
## p-values as a reader writes them: 3 x 10^-8, with a true superscript, not the "3e-08" of R
sup <- function(x) paste(c("0" = "\u2070", "1" = "\u00b9", "2" = "\u00b2", "3" = "\u00b3", "4" = "\u2074", "5" = "\u2075", "6" = "\u2076", "7" = "\u2077", "8" = "\u2078", "9" = "\u2079", "-" = "\u207b")[strsplit(x, "")[[1]]], collapse = "")
fmt_num <- function(p) if (p >= 1e-3) sprintf("%.3f", p) else {
  e <- floor(log10(p)); m <- round(p / 10^e); if (m == 10) { m <- 1; e <- e + 1 }
  sprintf("%d \u00d7 10%s", m, sup(as.character(e))) }
fmt_p <- function(p) paste("p =", fmt_num(p))

col_e <- ek_pal[["EDGE-poly3"]]; col_h <- ek_pal[["HL"]]; col_f <- ek_pal[["HL-F"]]
lim <- range(c(cc$mean_pred, cc$obs_rate)) + c(-0.01, 0.01)

cc$sign <- ifelse(cc$gap_pts > 0, "observed above predicted", "observed below predicted")
pL <- ggplot(cc, aes(mean_pred, obs_rate)) +
  ek_diag() +
  geom_segment(aes(xend = mean_pred, yend = mean_pred, colour = sign), linewidth = 2.2, alpha = 0.55) +
  geom_line(colour = "grey40", linewidth = 0.5) +
  geom_point(aes(colour = sign), size = 2.4) +
  scale_colour_manual(values = c("observed above predicted" = "#D55E00",
                                 "observed below predicted" = "#0072B2"), name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = lim) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = lim) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 2.7, family = EK_FAMILY, colour = "grey30",
           label = sprintf("intercept-and-slope test: %s\nEDGE, 4 df: p = %s to %s\nover all partitions",
                           fmt_p(sm$cox_p), fmt_num(sm$edge_val_p_min), fmt_num(sm$edge_val_p_max))) +
  labs(title = sprintf("(a) Validation half (n = %s), frozen model", format(sm$n_val, big.mark = ",")),
       x = "Mean predicted risk in decile", y = "Observed readmission rate") +
  theme_ek() + theme(legend.position = "inside", legend.position.inside = c(0.98, 0.03),
                     legend.justification = c(1, 0), legend.direction = "vertical",
                     legend.text = element_text(size = rel(0.85)))

long <- bind_rows(
  dev %>% transmute(G, p = EDGE.poly3, test = "EDGE-poly3"),
  dev %>% transmute(G, p = HL, test = "HL"),
  dev %>% transmute(G, p = HL_F, test = "HL-F"))
long$test <- factor(long$test, levels = c("EDGE-poly3", "HL", "HL-F"))
## the legend names the tests in words; the internal codes stay in the data
LAB <- c("EDGE-poly3" = "EDGE", "HL" = "Hosmer\u2013Lemeshow", "HL-F" = "Hosmer\u2013Lemeshow, Farrington")
pR <- ggplot(long, aes(G, p, colour = test, linetype = test)) +
  ek_nominal() +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = c("EDGE-poly3" = col_e, "HL-F" = col_f, "HL" = col_h), labels = LAB, name = NULL) +
  scale_linetype_manual(values = c("EDGE-poly3" = "solid", "HL-F" = "22", "HL" = "42"), labels = LAB, name = NULL) +
  scale_x_log10(breaks = c(10, 50, 200, 1000, 2000), labels = c("10", "50", "200", "1000", "2000")) +
  scale_y_log10(breaks = 10^c(0, -1, -2, -3, -4, -5, -6), labels = c("1", "0.1", "0.01", "0.001", vapply(-4:-6, function(e) paste0("10", sup(as.character(e))), ""))) +
  labs(title = sprintf("(b) Development half (n = %s), same model", format(sm$n_dev, big.mark = ",")),
       x = "Number of risk groups  G  (log scale)", y = "p-value (log scale)") +
  theme_ek() + theme(legend.position = "inside", legend.position.inside = c(0.98, 0.42),
                     legend.justification = c(1, 0.5), legend.direction = "vertical",
                     legend.key.width = unit(1.6, "lines"), legend.text = element_text(size = rel(0.85)))

p <- pL + pR + plot_layout(widths = c(1, 1))
ggsave(file.path(FIG, "Fig10_diabetes.pdf"), p, device = cairo_pdf, width = EK_W2, height = 3.8, units = "in")
ggsave(file.path(SIM, "_preview_Fig10_p3.png"), p, width = EK_W2, height = 3.8, units = "in", dpi = 150)
cat("wrote paper3/figures/Fig10_diabetes.pdf\n")
