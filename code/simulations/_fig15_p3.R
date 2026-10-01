## _fig15_p3.R -- paper-3 Figure 3: power, protection and cost of each test, read together.
##
## Paper 2's version plotted power on clean data against time. It made the cubic calibration test look like the
## best test outright, because the one thing that separates the tests in this paper, what a few wrong records do
## to them, was not on the page. This version puts it there:
##   (a) mean size-adjusted power over the twenty scenarios of the rivals study (clean data) against the false-
##       alarm rate with ten corrupted records in 1000 (corruption C1, correct model); the label gives the time a
##       data set at n = 1000. The two resampling tests were not run with corrupted records and stand apart.
##   (b) median seconds a data set against n, log-log, on one core; le Cessie's test twice, as ebrahim.gof computes
##       it (O(n^2 p)) and in the O(n^3) form of the public code it was adapted from (same statistic).
## Power: the extended census (classical tests), block 8 (the resampling tests), block 8L (le Cessie).
## Protection: block 9, cell_test.csv, C1 at n = 1000, k = 10; le Cessie from block 8L Part B, on block 9's data.
## Time: battery/T/timing.csv and battery/T/timing_lecessie_n3.csv.
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
suppressMessages({ library(ggplot2); library(data.table); library(ggrepel) })

TT  <- fread(edge_battery("T", "timing.csv"))
TT3 <- fread(edge_battery("T", "timing_lecessie_n3.csv"))
S8  <- fread(edge_battery("8", "_summary.csv"))[role == "alternative" & alpha == 0.05]
L8  <- fread(edge_battery("8L", "_summary.csv"))
L8A <- L8[part == "A" & role == "alternative"]
CEN <- fread(edge_battery("analysis", "census_extended.csv"))
B9  <- fread(edge_battery("9", "analysis", "cell_test.csv"))
cells20 <- sort(unique(L8A$cell)); stopifnot(length(cells20) == 20L)

pw_census <- function(tname) mean(CEN[cell %in% cells20 & test == tname, power], na.rm = TRUE)
fa9 <- function(tname) B9[truth == "logit" & corruption == "C1" & n == 1000 & k == 10 & test == tname, rejection]
D <- data.table(
  test  = c("EDGE", "Hosmer-Lemeshow", "Stukel joint", "cubic calib. LR", "GiViTI", "le Cessie", "Liu projection", "BAGofT"),
  power = c(pw_census("EDGE.poly3.u.Grule"), pw_census("HL.G10"), pw_census("Stk.joint"), pw_census("Cubic.LR"),
            pw_census("GiViTI"), mean(L8A$lecessie_size_adj, na.rm = TRUE),
            mean(S8[rival == "proj" & test == "proj", size_adj_power], na.rm = TRUE),
            mean(S8[rival == "bagoft" & test == "BAGofT", size_adj_power], na.rm = TRUE)),
  fa    = c(fa9("EDGE.poly3.u.Grule"), fa9("HL.G10"), fa9("Stk.joint"), fa9("Cubic.LR"), fa9("GiViTI"),
            L8[part == "B" & role == "corrupted" & n == 1000 & k == 10, lecessie_rejection], NA, NA))
stopifnot(all(is.finite(D$power)), sum(is.finite(D$fa)) == 6L)
## le Cessie's test in its published O(n^3) form, not the O(n^2 p) rewrite in ebrahim.gof (the author's choice,
## 2026-10-01): same statistic, the time a reader would meet
TT <- rbind(TT[test != "le Cessie"], TT3[, .(test = "le Cessie", n, median_sec)], fill = TRUE)
D <- merge(D, TT[n == 1000, .(test, sec = median_sec)], by = "test")

## EDGE at ten groups, the setting recommended when records may be wrong (referee round, 2026-10-01)
CPT <- fread(edge_battery("analysis", "corruption_partitions.csv"))
g10 <- data.table(test = "EDGE, 10 groups", power = pw_census("EDGE.poly3.u.G10"),
                  fa = CPT[error == "exaggeration" & n == 1000 & k == 10 & test == "EDGE.poly3.u.G10", false_alarm],
                  sec = D[test == "EDGE", sec])
stopifnot(nrow(g10) == 1L, is.finite(g10$fa), is.finite(g10$power))
D <- rbind(D, g10, fill = TRUE)
fmt_t <- function(s) ifelse(s < 1, sprintf("%.0f ms", 1000 * s), sprintf("%.1f s", s))
D[, label := sprintf("%s\n%s", test, fmt_t(sec))]
D[test == "Hosmer-Lemeshow", label := sprintf("Hosmer\u2013Lemeshow, G = 10\n%s", fmt_t(sec))]
D[test == "EDGE", label := sprintf("EDGE, default, %s", fmt_t(sec))]
D[test == "EDGE, 10 groups", label := sprintf("EDGE, 10 groups, %s", fmt_t(sec))]

XNA <- 0.80                                                        # where the two untested resampling tests stand
D[, x := ifelse(is.finite(fa), fa, XNA)]
D[, measured := is.finite(fa)]
## option B, the author's choice (2026-10-01): the colour gives the cost class, the label the exact time
D[, cost := cut(sec, c(0, 0.05, 20, Inf), labels = c("milliseconds", "seconds",
                                                     "seconds to minutes"))]

## labels are placed by hand, not repelled: the PDF device measures text differently from the PNG preview, and
## repelled labels that clear each other in one collided in the other (author's review, 2026-10-01)
LAB <- data.table(test = c("EDGE", "EDGE, 10 groups", "Hosmer-Lemeshow", "le Cessie", "GiViTI", "cubic calib. LR",
                           "Stukel joint", "Liu projection", "BAGofT"),
                  dx = c(0.025, 0.025, 0.025, 0.025, 0.025, -0.025, 0, 0, 0),
                  dy = c(0.006, -0.008, 0, 0, 0, 0, -0.038, -0.038, -0.038),
                  hj = c(0, 0, 0, 0, 0, 1, 0.5, 0.5, 0.5))
D <- merge(D, LAB, by = "test")
D[test == "Hosmer-Lemeshow", label := sprintf("Hosmer\u2013Lemeshow, G = 10\n%s", fmt_t(sec))]
D[test == "cubic calib. LR", label := sprintf("cubic calibration LR\n%s", fmt_t(sec))]
D[test == "le Cessie", label := sprintf("le Cessie\n%s", fmt_t(sec))]

pa <- ggplot(D, aes(x, power)) +
  annotate("rect", xmin = -Inf, xmax = 0.10, ymin = -Inf, ymax = Inf, fill = "#009E73", alpha = 0.08) +
  annotate("text", x = 0.012, y = 0.015, hjust = 0, vjust = 1, angle = 90, size = 2.4, family = EK_FAMILY,
           colour = "#007A5A", label = "false alarms below 0.10") +
  geom_vline(xintercept = 0.72, colour = "grey75", linewidth = 0.4) +
  annotate("text", x = XNA, y = 0.015, vjust = 0, size = 2.4, family = EK_FAMILY, colour = "grey40",
           label = "not run\nwith corrupted\nrecords") +
  annotate("segment", x = D[test == "EDGE", x], y = D[test == "EDGE", power],
           xend = D[test == "EDGE, 10 groups", x], yend = D[test == "EDGE, 10 groups", power],
           colour = "#0072B2", linewidth = 0.4) +
  geom_point(aes(shape = measured, colour = cost), size = 2.7, stroke = 1) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = c("#0072B2", "#8E44AD", "#B22222"), name = "time at n = 1000:", drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.2))) +
  geom_text(aes(x = x + dx, y = power + dy, label = label, hjust = hj), vjust = 0.5, size = 2.4,
            family = EK_FAMILY, lineheight = 0.88) +
  scale_x_continuous(limits = c(0, 0.90), breaks = c(0.1, 0.2, 0.4, 0.6),
                     labels = c("0.10", "0.20", "0.40", "0.60"), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, 0.52), breaks = seq(0, 0.5, 0.1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "false-alarm rate, 10 corrupted records in 1000",
       y = "mean size-adjusted power, clean data",
       title = "(a) Power, protection and cost") +
  ## the legend sits under the panel in one row, clear of the data (author's review, 2026-10-01)
  theme_ek() + theme(legend.position = "bottom", legend.direction = "horizontal",
                     legend.title = element_text(size = rel(0.72)), legend.text = element_text(size = rel(0.7)),
                     legend.key.height = unit(0.6, "lines"), legend.key.width = unit(0.6, "lines"),
                     legend.spacing.x = unit(2, "pt"), legend.margin = margin(t = -6),
                     plot.title = element_text(size = rel(0.95), hjust = 0), plot.title.position = "plot")

## (b): each line is named at its right end instead of in a legend, so the panel keeps the height of (a)
TB <- TT[is.finite(median_sec), .(test, n, median_sec)]
cols <- c("EDGE" = "#D55E00", "Hosmer-Lemeshow" = "#009E73", "Stukel joint" = "#CC79A7",
          "cubic calib. LR" = "#E69F00", "GiViTI" = "#7F7F7F", "le Cessie" = "#882255",
          "Liu projection" = "#56B4E9", "BAGofT" = "#000000")
NAMES <- c("EDGE" = "EDGE", "Hosmer-Lemeshow" = "Hosmer\u2013Lemeshow", "Stukel joint" = "Stukel joint",
           "cubic calib. LR" = "cubic LR", "GiViTI" = "GiViTI", "le Cessie" = "le Cessie",
           "Liu projection" = "Liu projection", "BAGofT" = "BAGofT")
ENDS <- TB[, .SD[which.max(n)], by = test][, name := NAMES[test]]
pb <- ggplot(TB, aes(n, median_sec, colour = test)) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.3) +
  ggrepel::geom_text_repel(data = ENDS, aes(label = name), hjust = 0, nudge_x = 0.12, direction = "y",
                           size = 2.4, family = EK_FAMILY, segment.size = 0.2, segment.colour = "grey60",
                           min.segment.length = 0.3, box.padding = 0.12, seed = 1) +
  scale_x_log10(limits = c(450, 900000), breaks = c(500, 2000, 10000, 100000), labels = c("500", "2k", "10k", "100k")) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("1 ms", "10 ms", "0.1 s", "1 s", "10 s", "100 s", "1000 s")) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = "sample size (log scale)", y = "time a data set, one core (log scale)", title = "(b) Cost against sample size") +
  theme_ek() + theme(plot.margin = margin(t = 5.5, r = 4, b = 5.5, l = 8),
                     plot.title = element_text(size = rel(0.95), hjust = 0), plot.title.position = "plot")

## patchwork aligns the two plotting regions although only (a) carries a legend
g <- patchwork::wrap_plots(pa, pb, ncol = 2, widths = c(1.1, 1))
ggsave(file.path(FIG, "Fig15_cost.pdf"), g, width = EK_W2, height = EK_W2 * 0.46, device = cairo_pdf)
ggsave(file.path(SIM, "_preview_Fig15_p3.png"), g, width = EK_W2, height = EK_W2 * 0.46, dpi = 150)
print(D[order(-power), .(test, power = round(power, 3), false_alarm = fa, sec)])
cat("written: Fig15_cost.pdf\n")
