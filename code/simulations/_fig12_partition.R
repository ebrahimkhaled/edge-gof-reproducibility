## _fig12_partition.R -- Figure: the partition is the lever, not the contamination rate.
## Block 9c (PREDECLARATION_block9c_groupsize.md, sha256 6320ff35), EDGE-poly3 unit form.
## Left panel n = 1000, right panel n = 5000; one line per number of corrupted records k; x is the
## number of groups on a log scale, y the false-alarm rate under a CORRECT model.
## The k = 0 line is the control: if the partition itself broke the test, that line would rise too.
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
FIG <- edge_out("paper2_figures")
source(file.path(SIM, "_ek_theme.R"))
suppressMessages({ library(ggplot2); library(data.table) })

L <- fread(edge_battery("9c", "analysis", "cell_G_test.csv"))
E <- L[test == "EDGE.poly3.u"]
E[, k_lab := factor(k, levels = c(0, 5, 10, 25, 50),
                    labels = c("0 (control)", "5", "10", "25", "50"))]
E[, panel := factor(n, levels = c(1000, 5000),
                    labels = c("n = 1000", "n = 5000"))]

## Okabe-Ito, ordered so the control is grey and the damage deepens with k
pal <- c("0 (control)" = "grey55", "5" = "#56B4E9", "10" = "#009E73",
         "25" = "#E69F00", "50" = "#D55E00")

## the control is dashed so that it stays visible where it coincides with k = 5 and k = 10,
## which is the point: at those k the contaminated cells behave like the uncontaminated one
p <- ggplot(E, aes(G, rejection, colour = k_lab, group = k_lab, linetype = k_lab)) +
  geom_hline(yintercept = 0.05, linetype = "13", colour = "grey45", linewidth = 0.4) +
  geom_line(linewidth = 0.85) +
  scale_linetype_manual(values = c("0 (control)" = "22", "5" = "solid", "10" = "solid",
                                   "25" = "solid", "50" = "solid"), name = "corrupted records  k") +
  geom_point(size = 2.1) +
  facet_wrap(~ panel, scales = "free_x") +
  scale_colour_manual(values = pal, name = "corrupted records  k") +
  scale_x_log10(breaks = c(10, 20, 25, 40, 50, 100, 200)) +
  scale_y_continuous(limits = c(0, 0.92), breaks = seq(0, 0.9, 0.15)) +
  annotate("text", x = 11, y = 0.085, hjust = 0, size = 2.6, family = EK_FAMILY,
           colour = "grey35", label = "nominal 0.05") +
  labs(x = "Number of risk groups  G  (log scale);  group size is n/G",
       y = "False-alarm rate") +
  theme_ek() +
  theme(legend.position = "bottom", panel.spacing = unit(1.1, "lines"))

dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIG, "Fig12_partition.pdf"), p, width = EK_W2, height = EK_W2 * 0.44,
       device = cairo_pdf)
ggsave(file.path(SIM, "_preview_Fig12_partition.png"), p, width = EK_W2, height = EK_W2 * 0.44,
       dpi = 160)

cat("the two readings the figure has to make obvious:\n")
cat(sprintf("  same rate 0.5%%: n=1000 k=5 at rule G=40 -> %.3f | n=5000 k=25 at rule G=200 -> %.3f\n",
            E[n == 1000 & k == 5 & G == 40, rejection], E[n == 5000 & k == 25 & G == 200, rejection]))
cat(sprintf("  same data n=5000 k=25: G=10 -> %.3f | G=200 -> %.3f\n",
            E[n == 5000 & k == 25 & G == 10, rejection], E[n == 5000 & k == 25 & G == 200, rejection]))
cat(sprintf("  control k=0 across every G: %.3f to %.3f\n",
            min(E[k == 0, rejection]), max(E[k == 0, rejection])))
cat("written: figures/Fig12_partition.pdf\n")
