## _fig10_diabetes.R -- Figure 10: the Diabetes-130 demonstration in one picture.
## Left: the validation-half calibration curve by decile of predicted risk (frozen model), with
## the gap from the diagonal drawn as bars -- the S that the intercept-and-slope test passes.
## Right: development half, the p-value of the directed default and of the Hosmer-Lemeshow test
## as the partition is refined from G = 10 to G = 2000. Foundation: _ek_theme.R.
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
suppressMessages({ library(dplyr); library(ggplot2); library(patchwork) })

cc  <- read.csv(file.path(SIM, "runI_bigdata_calcurve.csv"))
dev <- read.csv(file.path(SIM, "runI_bigdata_dev.csv"))

col_e <- ek_pal[["EDGE-poly3"]]; col_h <- ek_pal[["HL"]]; col_f <- ek_pal[["HL-F"]]

## ---- left: calibration curve, validation half --------------------------------------
cc$sign <- ifelse(cc$gap_pts > 0, "observed above predicted", "observed below predicted")
pL <- ggplot(cc, aes(mean_pred, obs_rate)) +
  ek_diag() +
  geom_segment(aes(xend = mean_pred, yend = mean_pred, colour = sign), linewidth = 2.2, alpha = 0.55) +
  geom_line(colour = "grey40", linewidth = 0.5) +
  geom_point(aes(colour = sign), size = 2.4) +
  scale_colour_manual(values = c("observed above predicted" = "#D55E00",
                                 "observed below predicted" = "#0072B2"), name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.04, 0.25)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.04, 0.25)) +
  annotate("text", x = 0.245, y = 0.055, hjust = 1, size = 2.7, family = EK_FAMILY, colour = "grey30",
           label = "intercept-and-slope test: p = 0.66\ndirected test, 4 df: p < 1e-12") +
  labs(title = "Validation half (n = 50,882), frozen model",
       x = "Mean predicted risk in decile", y = "Observed readmission rate") +
  theme_ek() + theme(legend.position = "bottom", legend.text = element_text(size = rel(0.85)))

## ---- right: p-value against G, development half ----------------------------------------
long <- bind_rows(
  dev %>% transmute(G, p = EDGE.poly3, test = "EDGE-poly3"),
  dev %>% transmute(G, p = HL, test = "HL"),
  dev %>% transmute(G, p = HL_F, test = "HL-F"))
long$test <- factor(long$test, levels = c("EDGE-poly3", "HL-F", "HL"))
pR <- ggplot(long, aes(G, p, colour = test, linetype = test)) +
  ek_nominal() +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = c("EDGE-poly3" = col_e, "HL-F" = col_f, "HL" = col_h), name = NULL) +
  scale_linetype_manual(values = c("EDGE-poly3" = "solid", "HL-F" = "22", "HL" = "42"), name = NULL) +
  scale_x_log10(breaks = c(10, 50, 200, 1000, 2000), labels = c("10", "50", "200", "1000", "2000")) +
  scale_y_log10(breaks = 10^c(0, -2, -4, -8, -12, -16),
                labels = c("1", "0.01", "1e-4", "1e-8", "1e-12", "1e-16")) +
  annotate("text", x = 12, y = 0.09, hjust = 0, size = 2.6, family = EK_FAMILY, colour = "grey45",
           label = "0.05") +
  labs(title = "Development half (n = 50,881), same model",
       x = "Number of risk groups  G  (log scale)", y = "p-value (log scale)") +
  theme_ek() + theme(legend.position = "bottom")

p <- pL + pR + plot_layout(widths = c(1, 1))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIG, "Fig10_diabetes.pdf"), p, device = cairo_pdf, width = EK_W2, height = 3.8, units = "in")
ggsave(file.path(FIG, "_preview_Fig10_diabetes.png"), p, width = EK_W2, height = 3.8, units = "in", dpi = 150)
cat("wrote Fig10_diabetes.pdf\n")
