## _fig15_options.R -- mock-ups for the author's choice of how Figure 3 shows cost (A, B, C, D). Not a paper figure.
## Reuses the data preparation of _fig15_p3.R (its first 93 lines build D, pa and pb).
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
src <- readLines(file.path(SIM, "_fig15_p3.R"), encoding = "UTF-8")
eval(parse(text = src[1:93], encoding = "UTF-8"))
OUT <- file.path(SIM, "_fig15_options"); dir.create(OUT, showWarnings = FALSE)
H <- EK_W2 * 0.56

## short names for the second panel of A, where the label need not repeat the time
D[, short := c("EDGE" = "EDGE", "Hosmer-Lemeshow" = "HL, G = 10", "Stukel joint" = "Stukel joint",
               "cubic calib. LR" = "cubic LR", "GiViTI" = "GiViTI", "le Cessie" = "le Cessie",
               "Liu projection" = "Liu projection", "BAGofT" = "BAGofT")[test]]
D[, cost := cut(sec, c(0, 0.05, 30, Inf), labels = c("milliseconds (closed form)", "seconds (n \u00d7 n kernel)",
                                                  "minutes (resampling)"))]
shade <- list(annotate("rect", xmin = -Inf, xmax = 0.10, ymin = -Inf, ymax = Inf, fill = "#009E73", alpha = 0.08),
              geom_vline(xintercept = 0.72, colour = "grey75", linewidth = 0.4),
              annotate("text", x = XNA, y = 0.02, vjust = 0, size = 2.5, family = EK_FAMILY, colour = "grey40",
                       label = "not run\nwith corrupted\nrecords"))
xs <- scale_x_continuous(limits = c(0, 0.90), breaks = c(0.1, 0.2, 0.4, 0.6), labels = c("0.10", "0.20", "0.40", "0.60"))
ys <- scale_y_continuous(limits = c(0, 0.52), breaks = seq(0, 0.5, 0.1), expand = expansion(mult = c(0, 0.02)))
xl <- "false-alarm rate with 10 corrupted records in 1000"; yl <- "mean size-adjusted power, clean data"
tag <- function(p, t) p + labs(tag = t) + theme(plot.tag = element_text(family = EK_FAMILY, size = 10))

## ---- A: power vs protection | power vs time, one power axis
a1 <- ggplot(D, aes(x, power)) + shade +
  geom_point(aes(shape = measured), size = 2.6, colour = "#0072B2", stroke = 0.9) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  ggrepel::geom_text_repel(aes(label = short), size = 2.5, family = EK_FAMILY, seed = 3, box.padding = 0.4,
                           min.segment.length = 0.1, max.overlaps = Inf) +
  xs + ys + labs(x = xl, y = yl, title = "(a) Power and protection") + theme_ek()
a2 <- ggplot(D, aes(sec, power)) +
  annotate("rect", xmin = 0, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "#0072B2", alpha = 0.06) +
  geom_point(aes(shape = measured), size = 2.6, colour = "#0072B2", stroke = 0.9) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  ggrepel::geom_text_repel(aes(label = short), size = 2.5, family = EK_FAMILY, seed = 3, box.padding = 0.4,
                           min.segment.length = 0.1, max.overlaps = Inf) +
  scale_x_log10(limits = c(5e-4, 2000), breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("1 ms", "10 ms", "0.1 s", "1 s", "10 s", "100 s", "1000 s")) + ys +
  labs(x = "median time a data set at n = 1000 (log scale)", y = NULL, title = "(b) Power and cost") +
  theme_ek() + theme(axis.text.y = element_blank())
gA <- gridExtra::arrangeGrob(a1, a2, ncol = 2, widths = c(1.1, 1))
ggsave(file.path(OUT, "A_two_panels_shared_power.png"), gA, width = EK_W2, height = H * 0.9, dpi = 150)

## ---- B: one panel, colour = cost class, exact time in the label
b <- ggplot(D, aes(x, power)) + shade +
  geom_point(aes(colour = cost, shape = measured), size = 2.9, stroke = 1) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#D55E00"), name = "cost a data set, n = 1000") +
  ggrepel::geom_text_repel(aes(label = label), size = 2.45, family = EK_FAMILY, lineheight = 0.9, seed = 3,
                           box.padding = 0.45, min.segment.length = 0.1, force = 4, max.overlaps = Inf) +
  xs + ys + labs(x = xl, y = yl) + theme_ek() + theme(legend.position = "right")
gB <- gridExtra::arrangeGrob(b + theme(legend.position = "bottom", legend.direction = "vertical"), pb, ncol = 2, widths = c(1.15, 1))
ggsave(file.path(OUT, "B_colour_cost_class.png"), gB, width = EK_W2, height = H, dpi = 150)

## ---- C: one panel, point size = time (log)
cc <- ggplot(D, aes(x, power)) + shade +
  geom_point(aes(size = log10(sec * 1000), shape = measured), colour = "#0072B2", alpha = 0.75, stroke = 1) +
  scale_size_continuous(range = c(2, 11), breaks = c(0, 2, 4, 5), labels = c("1 ms", "0.1 s", "10 s", "100 s"),
                        name = "time a data set") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  ggrepel::geom_text_repel(aes(label = short), size = 2.5, family = EK_FAMILY, seed = 3, box.padding = 0.7,
                           min.segment.length = 0.1, max.overlaps = Inf) +
  xs + ys + labs(x = xl, y = yl) + theme_ek() + theme(legend.position = "bottom")
gC <- gridExtra::arrangeGrob(cc, pb, ncol = 2, widths = c(1.15, 1))
ggsave(file.path(OUT, "C_size_is_time.png"), gC, width = EK_W2, height = H, dpi = 150)

## ---- D: the current figure
file.copy(file.path(SIM, "_preview_Fig15_p3.png"), file.path(OUT, "D_current.png"), overwrite = TRUE)
cat("written:", paste(list.files(OUT), collapse = ", "), "\n")
