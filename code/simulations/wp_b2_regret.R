# WP-B2 — headline regret (deterministic, no RNG). Single provenance = the WP-B1 (1000,10) run.
# <<REGRET_MAX>>     = max over detectable scenarios of  bestDEF - (3DEF:CCT)   (EDGES trails oracle)
# <<FIXED_LOSS_MAX>> = max over (basis, detectable scenario) of  bestDEF - power(basis)  (cost of a fixed choice)
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
scns_order <- c("null","cloglog","loglog","cauchit","stukel_heavy","stukel_asym","quad","cubic","sawtooth","crossover")
REPS <- 2000L

POW <- read.csv("wp_b1_power.csv", check.names=FALSE)      # keep "3DEF:CCT" literal
POW <- POW[POW$n==1000 & POW$G==10, ]
CCTCOL <- "3DEF:CCT"
if (!(CCTCOL %in% names(POW))) { cat("COLS:", paste(names(POW), collapse=" | "), "\n"); stop("3DEF:CCT column not found in POW") }

PV <- read.csv("wp_b1_pvalues.csv")
PV <- PV[PV$n==1000 & PV$G==10, ]                          # per-basis raw single-test powers

scns_alt <- setdiff(scns_order, "null")
H <- do.call(rbind, lapply(scns_alt, function(s){
  bestDEF <- POW[POW$scenario==s, "bestDEF"]
  cct     <- POW[POW$scenario==s, CCTCOL]
  m <- PV$scenario==s & !is.na(PV$scenario)
  pw <- c(poly2=mean(PV$poly2[m] < 0.05),
          poly3=mean(PV$poly3[m] < 0.05),
          stk  =mean(PV$stk[m]   < 0.05))
  data.frame(scenario=s, detectable=(bestDEF >= 0.15), bestDEF=bestDEF, cct=cct,
             edges_regret=bestDEF - cct,
             poly2=pw["poly2"], poly3=pw["poly3"], stk=pw["stk"],
             loss_poly2=bestDEF - pw["poly2"], loss_poly3=bestDEF - pw["poly3"], loss_stk=bestDEF - pw["stk"],
             row.names=NULL, stringsAsFactors=FALSE)
}))

det <- H$detectable
if (any(is.na(H$bestDEF[det]) | is.na(H$cct[det]))) stop("NA in bestDEF/CCT for a detectable scenario")

REGRET_MAX  <- max(H$edges_regret[det])
regret_scn  <- H$scenario[det][which.max(H$edges_regret[det])]
loss_long   <- data.frame(basis=rep(c("poly2","poly3","stk"), each=sum(det)),
                          scenario=rep(H$scenario[det], 3),
                          loss=c(H$loss_poly2[det], H$loss_poly3[det], H$loss_stk[det]))
FIXED_LOSS_MAX <- max(loss_long$loss); fl_i <- which.max(loss_long$loss)

# __MAX__ summary row so both headline numbers are literally readable from the CSV
Hout <- rbind(H, data.frame(scenario="__MAX__", detectable=NA, bestDEF=NA, cct=NA,
                            edges_regret=REGRET_MAX, poly2=NA, poly3=NA, stk=NA,
                            loss_poly2=max(H$loss_poly2[det]), loss_poly3=max(H$loss_poly3[det]),
                            loss_stk=max(H$loss_stk[det]), row.names=NULL))
write.csv(Hout, "wp_b2_headline.csv", row.names=FALSE)

se <- function(p) sqrt(p*(1-p)/REPS)
cat(sprintf("\n<<REGRET_MAX>>     = %.4f  (worst detectable scenario: %s;  SE~%.4f)\n", REGRET_MAX, regret_scn, se(REGRET_MAX)))
cat(sprintf("<<FIXED_LOSS_MAX>> = %.4f  (basis %s on %s;  SE~%.4f)\n", FIXED_LOSS_MAX, loss_long$basis[fl_i], loss_long$scenario[fl_i], se(FIXED_LOSS_MAX)))
cat("\n=== per-scenario (detectable only) ===\n")
print(H[det, c("scenario","bestDEF","cct","edges_regret","loss_poly2","loss_poly3","loss_stk")], row.names=FALSE, digits=3)
cat("\n=== non-detectable (excluded) ===\n"); print(H$scenario[!det])
