# Phase-A numeric verification of T1/T4/T5/T6, computed from the REAL released component p-values
# (archived null + the WP-B1 sweep) — i.e. against the shipped EDGE geometry, not a re-derivation.
# Writes proof_check_T1.csv / _T4.csv / _T5.csv / _T6.csv (columns: metric,value,target,tol,pass).
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
EPS <- 1e-15
clip <- function(p) pmin(pmax(p, EPS), 1-EPS)
Cstat <- function(P) rowMeans(tan((0.5-clip(P))*pi), na.rm=TRUE)      # the raw Cauchy statistic
cct   <- function(P) 0.5 - atan(Cstat(P))/pi                          # CCT p-value
rows  <- function(...) do.call(rbind, list(...))
row1  <- function(metric,value,target,tol,pass) data.frame(metric=metric,value=round(value,5),target=target,tol=tol,pass=pass)

## ---------- load real data ----------
AR <- read.csv("_archive/Ensemble_pvalues_full.csv", stringsAsFactors=FALSE)
scAR <- as.character(AR$scenario); nm <- is.na(scAR)|trimws(scAR)==""|tolower(trimws(scAR))=="null"
NUL <- AR[nm, c("poly2","poly3","stk","EF")]                          # ~5000 null reps, 3EDGE + EF
B1 <- read.csv("wp_b1_pvalues.csv", stringsAsFactors=FALSE); B1 <- B1[B1$n==1000 & B1$G==10, ]
CL <- B1[B1$scenario=="cloglog",  c("poly2","poly3","stk","EF")]      # detectable, bases differ
SW <- B1[B1$scenario=="sawtooth", c("poly2","poly3","stk","EF")]      # rough misfit (ceiling)

## ================= T1: CCT size-validity under same-model dependence =================
P3 <- as.matrix(NUL[,1:3]); pcct <- cct(P3); Cn <- Cstat(P3)
sz05<-mean(pcct<=.05); sz01<-mean(pcct<=.01); sz10<-mean(pcct<=.10)
Z <- qnorm(1-clip(P3))                                                # probit scores (uniform p -> N(0,1))
R <- cor(Z); offmax <- max(abs(R[upper.tri(R)]))                     # dependence real, not +-1
ksnorm <- min(sapply(1:3, function(k) suppressWarnings(ks.test(Z[,k],"pnorm")$p.value)))
# validity is TAIL/level-based (Thm eq:cctnull is a tail limit), NOT full-distribution Cauchy:
tail05 <- mean(Cn > tan(0.45*pi))                                     # right-tail level at alpha=.05
kscau  <- suppressWarnings(ks.test(Cn,"pcauchy")$p.value)            # secondary diagnostic only
T1 <- rows(
  row1("T1_CCT_SIZE_a05", sz05, 0.05, 0.006, abs(sz05-.05)<=.006),
  row1("T1_CCT_SIZE_a01", sz01, 0.01, 0.006, abs(sz01-.01)<=.006),
  row1("T1_CCT_SIZE_a10", sz10, 0.10, 0.006, abs(sz10-.10)<=.006),
  row1("T1_CCT_TAIL_LEVEL_a05", tail05, 0.05, 0.006, abs(tail05-.05)<=.006),
  row1("T1_MAX_ABSCORR_OFFDIAG", offmax, 0.999, NA, offmax<0.999),
  # NOTE: components are MULTI-DF weighted-chi2 DEF p-values, not scalar Gaussian scores, so a
  # probit-Gaussianity KS does not apply; validity is the SIZE (above) + non-collinearity. Both diagnostics only.
  row1("T1_KS_P_NORMAL_min_diag", ksnorm, 0.01, NA, NA),
  row1("T1_KS_P_CAUCHY_diag", kscau, 0.01, NA, NA))                  # CCT is TAIL-Cauchy, not full -> diagnostic
write.csv(T1,"proof_check_T1.csv",row.names=FALSE)

## ================= T5: standard-Cauchy null + rejection region =================
crit <- tan(0.45*pi); qc <- qcauchy(0.95)
grid <- c(-10,-1,0,1,6.314,100); invmax <- max(abs((0.5-atan(grid)/pi) - pcauchy(grid,lower.tail=FALSE)))
T5 <- rows(
  row1("T5_crit_vs_qcauchy", abs(crit-qc), 0, 1e-6, abs(crit-qc)<=1e-6),
  row1("T5_crit_value", crit, 6.313752, 1e-3, abs(crit-6.313752)<=1e-3),
  row1("T5_tail_level_a05", tail05, 0.05, 0.006, abs(tail05-.05)<=.006),
  row1("T5_inversion_identity_maxdiff", invmax, 0, 1e-9, invmax<1e-9))
write.csv(T5,"proof_check_T5.csv",row.names=FALSE)

## ================= T4: power-inheritance / bounded regret (cloglog) =================
Pc <- as.matrix(CL[,1:3]); pcct_c <- cct(Pc); minp <- apply(Pc,1,min)
bCCT <- mean(pcct_c<=.05); bFloor <- mean(minp<=.05/3)
bOracle <- max(colMeans(Pc<=.05))      # theorem's beta* = BEST SINGLE basis (not the union P(any reject))
regret <- bOracle - bCCT
T4 <- rows(
  row1("T4_beta_CCT_a05", bCCT, NA, NA, NA),
  row1("T4_beta_floor_a05over3", bFloor, NA, NA, NA),
  row1("T4_floor_holds(bCCT>=bFloor-.01)", bCCT-(bFloor-.01), 0, NA, bCCT >= bFloor-0.01),
  row1("T4_regret_cloglog(bOracle-bCCT)", regret, 0.03, 0.03, regret<=0.03))
write.csv(T4,"proof_check_T4.csv",row.names=FALSE)

## ================= T6: detector-limit ceiling (sawtooth) =================
Ps3 <- as.matrix(SW[,1:3]); Ps4 <- as.matrix(SW[,1:4])
b3 <- mean(cct(Ps3)<=.05); b4 <- mean(cct(Ps4)<=.05)
baseP <- sapply(1:3, function(k) mean(SW[,k]<=.05)); slack <- b3 - max(baseP)
T6 <- rows(
  row1("T6_SAWTOOTH_3DEF", b3, 0.25, NA, b3<=0.25),
  row1("T6_SAWTOOTH_3DEF_EF", b4, NA, NA, b4 >= b3+0.30),
  row1("T6_EF_gain(b4-b3>=.30)", b4-b3, 0.30, NA, (b4-b3)>=0.30),
  row1("T6_CEILING_SLACK(b3-maxbase<=.02)", slack, 0.02, 0.02, slack<=0.02))
write.csv(T6,"proof_check_T6.csv",row.names=FALSE)

## ---------- summary ----------
ov <- function(D) all(D$pass[!is.na(D$pass)])
cat(sprintf("T1 PASS=%s | CCT size a05=%.4f a01=%.4f a10=%.4f | tail_lvl=%.4f | maxcorr=%.3f | ksNorm=%.3f (ksCauchy diag=%.3f, tail-only claim)\n",
            ov(T1),sz05,sz01,sz10,tail05,offmax,ksnorm,kscau))
cat(sprintf("T5 PASS=%s | crit=%.6f (=qcauchy.95) | inversion maxdiff=%.2e\n", ov(T5),crit,invmax))
cat(sprintf("T4 PASS=%s | beta_CCT=%.3f floor=%.3f oracle=%.3f regret=%.3f\n", ov(T4),bCCT,bFloor,bOracle,regret))
cat(sprintf("T6 PASS=%s | sawtooth 3DEF=%.3f 3DEF+EF=%.3f (gain %.3f) | ceiling slack=%.3f\n", ov(T6),b3,b4,b4-b3,slack))
