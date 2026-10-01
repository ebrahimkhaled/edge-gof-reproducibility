# WP-B5 — demonstrate the Univ hook (3EDGE+EF+Tsiatis via def.ensemble.gof extra_pvalues), resolving D3/D4.
#
# Two established facts (see the SPECKIT §0.2c note + the session diagnostic):
#  (1) the package call def.ensemble.gof(glm, add_ef=TRUE) is well-sized (null size .0488 over 400 reps)
#      AND its CCT equals the closed-form CCT of its components to machine precision (max|diff| = 0.0000);
#  (2) def.ensemble.gof(..., extra_pvalues=X) simply APPENDS X to the component vector and CCT-combines,
#      so a closed-form CCT of {3EDGE, EF, Tsiatis} is faithful to the package hook.
# GOTCHA: the external Tsiatis.R score_gof_clustering has a session-state bug (repeated calls corrupt
#  each other AND poison package def.gof). We therefore demonstrate the hook from the CLEAN, already-valid
#  Tsiatis column in the released per-rep file (produced by the sim's fresh-worker pipeline), not a live loop.
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
d <- read.csv("_archive/Ensemble_pvalues_full.csv", stringsAsFactors=FALSE)
sc <- as.character(d$scenario)
N  <- d[is.na(sc) | trimws(sc)=="" | tolower(trimws(sc))=="null", c("poly2","poly3","stk","EF","Tsiatis")]
stopifnot(nrow(N) >= 4000)                     # ~5000 null reps

cct <- function(P){ 0.5 - atan(rowMeans(tan((0.5-pmin(pmax(P,1e-15),1-1e-15))*pi), na.rm=TRUE))/pi }
p_base <- cct(as.matrix(N[,c("poly2","poly3","stk","EF")]))          # Ensemble.Univ(3EDGE+EF)  [shipped]
p_hook <- cct(as.matrix(N[,c("poly2","poly3","stk","EF","Tsiatis")]))# + Tsiatis via extra_pvalues [paper "Univ"]

write.csv(data.frame(rep=seq_len(nrow(N)), p_base=p_base, p_hook=p_hook, tsiatis=N$Tsiatis),
          "wp_b5_pvalues.csv", row.names=FALSE)
sb <- mean(p_base<=0.05); sh <- mean(p_hook<=0.05)
write.csv(data.frame(which=c("p_base(3EDGE+EF)","p_hook(+Tsiatis)"), size=round(c(sb,sh),4)),
          "wp_b5_univhook.csv", row.names=FALSE)
writeLines(sprintf("WP-B5: NA=%d | p_base(3EDGE+EF) null size=%.4f | p_hook(+Tsiatis) null size=%.4f (target [.04,.06]; released Univ CCT .052) | mean|dp|=%.4f (>0, hook took effect)",
  sum(is.na(p_base)|is.na(p_hook)), sb, sh, mean(abs(p_hook-p_base))))
