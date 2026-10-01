## _fig7_size.R — Type-I / size heatmap (EDGE paper, Fig 7)
## Story: EDGE bases + EF sit at ~nominal 0.05; Pigeon-Heyse is conservative (~0.03);
## covariate-space tests (Stukel/Xie) drift. Fill encodes deviation from 0.05.
## CAPTION: Empirical type-I error at the nominal 5% level, by test and sample size.
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
suppressMessages({ library(ggplot2); library(dplyr); library(readr); library(stringr) })

HERE <- edge_path("code/simulations")
FIGDIR <- edge_out("figures")
source(file.path(HERE, "_ek_theme.R"))
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

curated <- c("DEF.poly2","DEF.poly3","DEF.stukel","EF","HL",
             "HL-equalwidth","Pigeon-Heyse","Stukel","Tsiatis","Xie")

d <- read_csv(file.path(HERE, "sim_null.csv"), show_col_types = FALSE) %>%
  filter(abs(alpha - 0.05) < 1e-9, G == 10, test %in% curated) %>%
  mutate(
    Test = ek_factor(test),
    Test = factor(Test, levels = rev(ek_levels)),   # top-to-bottom reading order
    nlab = ifelse(n >= 1000, paste0(n/1000, "k"), as.character(n)),
    nfac = factor(nlab, levels = { u <- sort(unique(n)); ifelse(u >= 1000, paste0(u/1000, "k"), as.character(u)) }),
    lab  = sprintf("%.3f", null_size),
    # dark cells (strong deviation) get white text for contrast
    txtcol = ifelse(abs(null_size - 0.05) > 0.028, "white", "grey15")
  )

# nice family facet labels
fam_lab <- c(quad = "Quadratic", binint = "Binary interaction",
             contint = "Continuous interaction", link = "Link misspecification",
             rough = "Rough / non-smooth")
d$family_lab <- ifelse(d$family %in% names(fam_lab), fam_lab[d$family], d$family)

p <- ggplot(d, aes(x = nfac, y = Test, fill = null_size)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lab, colour = txtcol), size = 2.15, family = EK_FAMILY) +
  scale_colour_identity() +
  scale_fill_alpha_ek(mid = 0.05) +
  facet_wrap(~ family_lab, nrow = 1, labeller = label_wrap_gen(width = 16)) +
  labs(
    title = "Empirical size at the nominal 5% level",
    subtitle = expression("Empirical type–I error "*hat(alpha)*" at nominal "*alpha == 0.05~"(G = 10)  •  white ≈ nominal, blue = conservative, red = liberal"),
    x = "Sample size  n", y = NULL
  ) +
  theme_ek(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = rel(0.85)),
    axis.text.y = element_text(size = rel(0.9)),
    plot.subtitle = element_text(size = rel(0.92), hjust = 0.5, colour = "grey35", margin = margin(b = 6)),
    legend.position = "right",
    legend.key.height = unit(28, "pt"),
    legend.key.width  = unit(9, "pt"),
    panel.spacing = unit(5, "pt")
  )

ggsave(file.path(FIGDIR, "Fig7_size.pdf"), p, device = cairo_pdf,
       width = EK_W2, height = 4.5, units = "in")
ggsave(file.path(FIGDIR, "_preview_Fig7_size.png"), p,
       width = EK_W2, height = 4.5, units = "in", dpi = 150)

message("done: ", nrow(d), " cells, families = ", paste(unique(d$family), collapse=", "),
        ", n = ", paste(levels(d$nfac), collapse=","))
