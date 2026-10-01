## _fig5_frontier.R — Compute-vs-Power Pareto frontier ("the standout figure")
## x = median time per call (log10, seconds); y = mean power over LINK+COVARIATE scenarios
## (link + binint + contint families) at n=1000, alpha=0.05, per test.
## Fast curated tests: bench_compute_time_summary.csv (n=1000,p=4).
## Slow rivals (le-Cessie, Stute-Zhu, BAGofT, McCullagh, proj): bench_slow_timing_summary.csv (n=1000).
## proj is the Escanciano/Liu bootstrap-based projection test (model-based bootstrap, refit per
## replicate); like BAGofT and McCullagh it has no power point in our grid and appears only as a
## time-only reference tick on the slow end (~42 s per call at n=1000, B=1000).
## Pareto frontier drawn; EDGE annotated as "directed power at HL-level cost".
## CAPTION: EDGE sits on the compute-power frontier: directed power at a few milliseconds per call, orders of magnitude below the kernel, bootstrap-projection, and adaptive tests.
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
suppressMessages({ library(ggplot2); library(dplyr); library(ggrepel); library(scales) })

SIM <- edge_path("code/simulations")
FIG <- edge_out("figures")
source(file.path(SIM, "_ek_theme.R"))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

## ---- POWER: mean reject_rate over link+binint+contint at n=1000, alpha=0.05 ----
pw <- read.csv(file.path(SIM, "sim_power_broad.csv"), stringsAsFactors = FALSE)
pw <- pw[pw$n == 1000 & pw$alpha == 0.05 &
         pw$family %in% c("link", "binint", "contint"), ]
## Stukel has 2 NA reject_rate cells (degenerate link fits) -> na.rm. [QUIRK]
power <- pw %>% group_by(test) %>%
  summarise(power = mean(reject_rate, na.rm = TRUE), .groups = "drop")

## ---- TIME: seconds per call at n=1000 ----
ct <- read.csv(file.path(SIM, "bench_compute_time_summary.csv"), stringsAsFactors = FALSE)
ct <- ct[ct$n == 1000 & ct$p == 4, c("test", "time_median")]
slow <- read.csv(file.path(SIM, "bench_slow_timing_summary.csv"), stringsAsFactors = FALSE)
slow <- slow[slow$n == 1000, c("test", "time_median")]

## Curated fast tests keep their compute-bench time.
## Slow rivals take their slow-bench time (authoritative for the heavy kernels).
## DEF.poly2 has no timing row -> approximate by DEF.poly3 (poly2 is a strictly
## smaller design, i.e. this is a conservative upper bound on its cost). [QUIRK]
## HL-equalwidth (HLw) shares HL's O(n) binning cost -> use HL time. [QUIRK]
## Curated tests -> compute bench (the authoritative apples-to-apples grid, p=4).
## Slow rivals -> slow bench. DEF.poly3 appears in BOTH; the compute-bench value
## is the one used everywhere else in the paper, so prefer it for curated tests.
time_of <- function(t) {
  if (t %in% ct$test)   return(ct$time_median[ct$test == t])
  if (t %in% slow$test) return(slow$time_median[slow$test == t])
  if (t == "DEF.poly2") return(ct$time_median[ct$test == "DEF.poly3"])
  if (t == "HL-equalwidth") return(ct$time_median[ct$test == "HL"])
  NA_real_
}

## Assemble: all curated tests (from power) + slow rivals (power NA, time only).
tests <- unique(power$test)
df <- power %>% mutate(time = vapply(test, time_of, numeric(1)))

## Slow rivals: no power in sim_power_broad; plot them on the time axis only,
## at the mean power of the *closest generic omnibus* (le-Cessie) as a proxy where
## available, else drop the y. We keep only those with a power estimate for the
## frontier; rivals with NO power are shown as time-only reference ticks.
rivals <- c("le-Cessie", "Stute-Zhu", "BAGofT", "McCullagh", "proj")
time_slow <- function(t) slow$time_median[slow$test == t]  # rivals: slow bench
riv <- data.frame(test = rivals,
                  time = vapply(rivals, time_slow, numeric(1)),
                  stringsAsFactors = FALSE)
## le-Cessie IS in the curated power set? check:
riv$power <- vapply(riv$test, function(t)
  if (t %in% power$test) power$power[power$test == t] else NA_real_, numeric(1))

## Build the plotting frame. Points with power get plotted in the scatter.
main <- df %>% filter(!is.na(time)) %>% mutate(kind = "curated")
riv_pts <- riv %>% filter(!is.na(power) & !is.na(time)) %>%
  mutate(kind = "rival") %>% select(test, power, time, kind)
riv_only <- riv %>% filter(is.na(power) & !is.na(time)) %>%
  mutate(kind = "rival_time")

plotdf <- bind_rows(main, riv_pts) %>%
  mutate(Test = ek_factor(test)) %>% filter(!is.na(Test)) %>%
  ## display label: shorten "Tsiatis" -> "Ts" for the scatter callouts (Xie unchanged)
  mutate(disp_label = ifelse(as.character(Test) == "Tsiatis", "Tsi", as.character(Test)))

## ---- Pareto frontier: upper-left (fast & powerful). A point is on the frontier
## if no other point is both faster (smaller time) AND more powerful. ----
po <- plotdf %>% arrange(time)
front_idx <- logical(nrow(po)); best <- -Inf
for (i in seq_len(nrow(po))) { if (po$power[i] > best) { front_idx[i] <- TRUE; best <- po$power[i] } }
frontier <- po[front_idx, ]

## ---- gap annotation: EDGE (poly3) vs the heaviest rival (BAGofT) ----
edge_t  <- df$time[df$test == "DEF.poly3"]
bag_t   <- riv$time[riv$test == "BAGofT"]
gap_x   <- round(bag_t / edge_t)

## ---- x breaks: 1ms .. 500s ----
xbrk <- c(0.001, 0.01, 0.1, 1, 10, 100, 500)
xlab <- c("1 ms", "10 ms", "100 ms", "1 s", "10 s", "100 s", "500 s")

## time-only rivals -> vertical reference ticks with horizontal (de-rotated),
## enlarged labels near the top. le-Cessie ~1.2s and McCullagh ~1.3s nearly
## coincide on the log axis -> stagger le-Cessie one step lower so the two
## horizontal labels don't overprint. BAGofT sits near the right panel edge ->
## right-align its label so it stays inside the panel. [QUIRK]
riv_only <- riv_only[order(riv_only$time), ]
riv_only$y  <- 0.995
riv_only$hj <- 0.5
riv_only$y[riv_only$test == "le-Cessie"] <- 0.955
riv_only$hj[riv_only$test == "BAGofT"]  <- 1

## ---- per-point label nudges: pull the buried "Ts" (Tsiatis) point out of the
## tight cluster into the open gap below-right so it gets a clear leader line.
## nudge_* is an ADDITIVE offset to the point's data coords; on the log10 x-axis
## we compute it as (target - current) so the label lands at an absolute spot. ----
plotdf$nx <- 0; plotdf$ny <- 0
is_ts <- as.character(plotdf$Test) == "Tsiatis"
ts_x <- plotdf$time[is_ts]; ts_y <- plotdf$power[is_ts]
plotdf$nx[is_ts] <- 0.018 - ts_x    # target ~18 ms (open gap, right of the cluster)
plotdf$ny[is_ts] <- 0.305 - ts_y    # target ~30.5% (below the cluster, clear of HL)

p <- ggplot(plotdf, aes(time, power)) +
  ## shaded "efficient region" hugging the fast & powerful EDGE cluster
  annotate("rect", xmin = 0.0008, xmax = 0.011, ymin = 0.32, ymax = 0.58,
           fill = ek_pal[["EDGE-poly3"]], alpha = 0.05) +
  annotate("text", x = 0.00088, y = 0.655, hjust = 0, vjust = 0, size = 2.4,
           fontface = "italic", colour = ek_pal[["EDGE-poly3"]], family = EK_FAMILY,
           label = "efficient region\n(fast & powerful)", lineheight = 0.9) +
  ## Pareto frontier line
  geom_step(data = frontier, aes(time, power), direction = "vh",
            colour = "grey45", linewidth = 0.5, linetype = "22") +
  ## time-only rival reference ticks (le-Cessie w/o power, Stute-Zhu, BAGofT, McCullagh)
  geom_vline(data = riv_only, aes(xintercept = time),
             colour = "grey75", linewidth = 0.35, linetype = "13") +
  geom_text(data = riv_only, aes(x = time, y = y, label = test, hjust = hj),
            vjust = 1, size = 2.9, colour = "grey40",
            fontface = "italic", family = EK_FAMILY) +
  ## the scatter
  geom_point(aes(colour = Test, shape = kind), size = 2.9, stroke = 0.9) +
  ggrepel::geom_text_repel(aes(label = disp_label, colour = Test),
                           nudge_x = plotdf$nx, nudge_y = plotdf$ny,
                           size = 2.7, fontface = "bold", seed = 3, family = EK_FAMILY,
                           min.segment.length = 0, box.padding = 0.7, point.padding = 0.35,
                           force = 2.5, segment.size = 0.3, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_ek() +
  scale_shape_manual(values = c(curated = 16, rival = 17), guide = "none") +
  scale_x_log10(breaks = xbrk, labels = xlab,
                limits = c(0.0008, 600),
                minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = percent_format(accuracy = 1)) +
  annotation_logticks(sides = "b", colour = "grey70",
                      short = unit(0.04, "cm"), mid = unit(0.07, "cm"),
                      long = unit(0.11, "cm")) +
  labs(
    title = "Computational cost versus power",
    subtitle = "Mean power over link + covariate misspecification (n = 1000, α = 0.05) vs. median time per call",
    x = "Median time per test call (log scale)",
    y = "Mean power (link + covariate scenarios)") +
  theme_ek(base_size = 10) +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = rel(0.85), hjust = 0.5, colour = "grey30", margin = margin(b = 6)),
        panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92"))

## ---- key annotations ----
edge_pt <- plotdf[plotdf$test == "DEF.poly3", ]
p <- p +
  ## the big cost gap arrow, EDGE -> BAGofT
  annotate("segment", x = edge_t * 1.7, xend = bag_t * 0.55,
           y = 0.14, yend = 0.14,
           arrow = arrow(length = unit(0.18, "cm"), ends = "both", type = "closed"),
           colour = "grey35", linewidth = 0.5) +
  annotate("label", x = sqrt(edge_t * bag_t), y = 0.14,
           label = sprintf("~%s× slower", format(gap_x, big.mark = ",")),
           size = 2.7, fontface = "bold", colour = "grey20", family = EK_FAMILY,
           fill = "white", label.padding = unit(0.12, "cm")) +
  annotate("text", x = 0.0011, y = 0.035, hjust = 0, size = 2.4, colour = "grey50",
           family = EK_FAMILY,
           label = "← faster                                                           slower →")

ggsave(file.path(FIG, "Fig5_frontier.pdf"), p, device = cairo_pdf,
       width = EK_W2, height = 5.0, units = "in")
ggsave(file.path(FIG, "_preview_Fig5_frontier.png"), p,
       width = EK_W2, height = 5.0, units = "in", dpi = 150)

cat("power/time table:\n")
print(plotdf %>% arrange(desc(power)) %>%
        mutate(time_ms = round(time * 1000, 2), power = round(power, 3)) %>%
        select(test, power, time_ms, kind))
cat("\nframtier tests:", paste(as.character(frontier$Test), collapse = ", "), "\n")
cat("EDGE(poly3) time:", round(edge_t, 4), "s  BAGofT:", round(bag_t, 1), "s  gap:", gap_x, "x\n")
cat("time-only rivals:", paste(riv_only$test, collapse = ", "), "\n")
