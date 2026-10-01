## _fig15_cost.R -- what each test costs, and what it buys. Two panels, drawn after block T (the timing benchmark,
## design note sha256 85441ae0) has run:
##   (a) mean size-adjusted power over the twenty scenarios of the rivals study against median seconds a data set
##       at n = 1000, one point per test: the compute-power frontier of the EDGE papers, with paper 2's numbers;
##   (b) median seconds a data set against n, log-log, which is where the two resampling tests and le Cessie's
##       kernel separate from the closed-form tests.
## Power comes from the declared analyses only (the census for the classical tests, block 8 for the two resampling
## tests, block 8L for le Cessie); time comes from battery/T/timing.csv.
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
FIG <- edge_out("paper2_figures")
source(file.path(SIM, "_ek_theme.R"))
suppressMessages({ library(ggplot2); library(data.table); library(ggrepel) })

TT <- fread(edge_battery("T", "timing.csv"))
S8 <- fread(edge_battery("8", "_summary.csv"))[role == "alternative" & alpha == 0.05]
L8 <- fread(edge_battery("8L", "_summary.csv"))[part == "A" & role == "alternative"]
CEN <- fread(edge_path("results/analysis/_census_power_paper2.csv"))
cells20 <- sort(unique(L8$cell))                                   # the twenty scenarios of the rivals study
stopifnot(length(cells20) == 20L)

## mean size-adjusted power over those twenty scenarios
pw_census <- function(tname) mean(CEN[cell %in% cells20 & test == tname, power], na.rm = TRUE)
pw <- data.table(
  test  = c("EDGE", "Hosmer-Lemeshow", "Stukel joint", "cubic calib. LR", "GiViTI", "le Cessie", "Liu projection", "BAGofT"),
  power = c(pw_census("EDGE.poly3.u.Grule"), pw_census("HL.Grule"), pw_census("Stk.joint"), pw_census("Cubic.LR"),
            pw_census("GiViTI"), mean(L8$lecessie_size_adj, na.rm = TRUE),
            mean(S8[rival == "proj" & test == "proj", size_adj_power], na.rm = TRUE),
            mean(S8[rival == "bagoft" & test == "BAGofT", size_adj_power], na.rm = TRUE)))

## time at n = 1000, seconds a data set
t1000 <- TT[n == 1000, .(test, sec = median_sec)]
D <- merge(pw, t1000, by = "test")
if (nrow(D) < nrow(pw)) message("no timing at n = 1000 for: ", paste(setdiff(pw$test, t1000$test), collapse = ", "))

pa <- ggplot(D, aes(sec, power)) +
  geom_point(size = 2.6, colour = "#0072B2") +
  ggrepel::geom_text_repel(aes(label = test), size = 2.7, family = EK_FAMILY, min.segment.length = 0.2, seed = 1) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("1 ms", "10 ms", "0.1 s", "1 s", "10 s", "100 s", "1000 s")) +
  ## the vertical axis ranks the tests on ONE family, the twenty link and tail scenarios every test
  ## here was run on, and a reader who takes it for a general ranking reads the figure backwards. The
  ## subtitle says so inside the panel rather than leaving it to the caption alone.
  labs(x = "median time a data set at n = 1000 (log scale)",
       y = "mean size-adjusted power on the twenty\nscenarios all eight tests share",
       subtitle = "one departure family, not a general ranking") +
  theme_ek() +
  theme(plot.subtitle = element_text(size = rel(0.78), colour = "grey35",
                                     margin = margin(b = 5)))

pb <- ggplot(TT[is.finite(median_sec)], aes(n, median_sec, colour = test)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.5) +
  scale_x_log10(breaks = c(500, 1000, 5000, 20000, 100000),
                labels = c("500", "1k", "5k", "20k", "100k")) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c("1 ms", "10 ms", "0.1 s", "1 s", "10 s", "100 s", "1000 s")) +
  scale_colour_manual(values = c("EDGE" = "#0072B2", "Hosmer-Lemeshow" = "#009E73", "Stukel joint" = "#CC79A7",
                                 "cubic calib. LR" = "#E69F00", "GiViTI" = "#7F7F7F", "le Cessie" = "#D55E00",
                                 "Liu projection" = "#56B4E9", "BAGofT" = "#000000"), name = NULL) +
  labs(x = "sample size (log scale)", y = "median time a data set (log scale)") +
  theme_ek() + theme(legend.position = "right")

g <- gridExtra::arrangeGrob(pa, pb, ncol = 2, widths = c(1.05, 1.15))
ggsave(file.path(FIG, "Fig15_cost.pdf"), g, width = EK_W2, height = EK_W2 * 0.42, device = cairo_pdf)
ggsave(file.path(SIM, "_preview_Fig15_cost.png"), g, width = EK_W2, height = EK_W2 * 0.42, dpi = 150)
print(D[order(sec)])
cat("written: Fig15_cost.pdf\n")
