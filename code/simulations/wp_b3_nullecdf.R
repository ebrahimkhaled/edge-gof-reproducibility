# WP-B3 — combiner NULL ECDF/QQ (Fig 2 panel b). Closed-form recompute from stored components.
# Shows Fisher/Stouffer null p-values are stochastically too small (liberal) while CCT/HMP/minP sit on
# the diagonal. No new simulation. Combiner math = verbatim from _ensemble_combine.py:15-25.
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
setwd(edge_path("code/simulations"))
suppressMessages({ library(ggplot2) })
ok <- tryCatch({ source("_ek_theme.R"); TRUE }, error=function(e) FALSE)   # theme_ek, EK_W2, EK_DPI, ek_diag
if (!ok || !exists("theme_ek")) theme_ek <- function(...) theme_minimal(base_size=10)
if (!exists("EK_W2"))  EK_W2  <- 174/25.4   # Springer double-column width (in)
if (!exists("EK_DPI")) EK_DPI <- 600

df <- read.csv("_archive/Ensemble_pvalues_full.csv", stringsAsFactors=FALSE)
# NULL scenario is stored with a BLANK label (SPECKIT §0.2b) -> match blank/NA/"null".
sc <- as.character(df$scenario)
nullmask <- is.na(sc) | trimws(sc)=="" | tolower(trimws(sc))=="null"
N <- df[nullmask, c("poly2","poly3","stk")]          # 3DEF component matrix (the default set)
stopifnot(nrow(N) >= 4000)                            # ~5000 null reps expected
P <- as.matrix(N)

# ---- combiners on a p-matrix P (columns = the set) : verbatim R translation of _ensemble_combine.py ----
kf   <- rowSums(is.finite(P))
minp <- 1 - (1 - apply(P,1,min,na.rm=TRUE))^kf
fish <- pchisq(-2*rowSums(log(pmin(pmax(P,1e-300),1)),na.rm=TRUE), 2*kf, lower.tail=FALSE)
stou <- pnorm(rowSums(qnorm(pmin(pmax(P,1e-12),1-1e-12),lower.tail=FALSE),na.rm=TRUE)/sqrt(kf), lower.tail=FALSE)
hmp  <- kf / rowSums(1/pmin(pmax(P,1e-300),1),na.rm=TRUE)
cct  <- 0.5 - atan(rowMeans(tan((0.5-pmin(pmax(P,1e-15),1-1e-15))*pi),na.rm=TRUE))/pi
COMB <- list(minP=minp, Fisher=fish, Stouffer=stou, HMP=hmp, CCT=cct)

# ---- calibration-check table (gate #2: reproduce the released 3DEF size column) ----
sizes <- sapply(COMB, function(p) mean(p <= 0.05))
write.csv(data.frame(combiner=names(sizes), size=round(as.numeric(sizes),4)),
          "wp_b3_null_size.csv", row.names=FALSE)
cat("=== WP-B3 null size (target: minP .033 / Fisher .108 / Stouffer .139 / HMP .050 / CCT .050) ===\n")
print(round(sizes,4))

# ---- ECDF on a grid ----
u <- seq(0, 1, by=0.005)
ecdf_long <- do.call(rbind, lapply(names(COMB), function(nm){
  e <- ecdf(COMB[[nm]]); data.frame(combiner=nm, u=u, ecdf=e(u))
}))
write.csv(ecdf_long, "wp_b3_null_ecdf.csv", row.names=FALSE)

# ---- plot ----
pal <- c(minP="#999999", Fisher="#D55E00", Stouffer="#B2182B", HMP="#66C2A4", CCT="#238B45")
ecdf_long$combiner <- factor(ecdf_long$combiner, levels=c("Fisher","Stouffer","minP","HMP","CCT"))
p <- ggplot(ecdf_long, aes(u, ecdf, colour=combiner)) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey55", linewidth=0.4) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=pal, name=NULL) +
  coord_equal(xlim=c(0,0.2), ylim=c(0,0.35)) +   # zoom on the small-p tail where liberality shows
  labs(x="nominal level u", y="empirical CDF of the null p-value under H0") +
  theme_ek() + theme(legend.position=c(0.82,0.28))
ggsave("Fig2b_null_ecdf.pdf", p, width=EK_W2, height=EK_W2*0.62, device=cairo_pdf)
ggsave("Fig2b_null_ecdf_preview.png", p, width=EK_W2, height=EK_W2*0.62, dpi=150)

# ---- gate checks ----
at05 <- sapply(COMB, function(pv) mean(pv <= 0.05))
cat(sprintf("\nGATE #3: ECDF at u=0.05  CCT=%.3f HMP=%.3f (near .05) | Fisher=%.3f Stouffer=%.3f (above .05)\n",
            at05["CCT"], at05["HMP"], at05["Fisher"], at05["Stouffer"]))
cat("saved: Fig2b_null_ecdf.pdf, wp_b3_null_ecdf.csv, wp_b3_null_size.csv\n")
