## WP-0.1 + WP-0.2 + WP-0.4 ground truth for the Fable revision (read-only analyses)
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
suppressMessages({library(stats)})
SIM <- edge_path("code/simulations")
setwd(SIM)

pb <- read.csv("sim_power_broad.csv", stringsAsFactors=FALSE)
pb <- pb[pb$n==1000 & pb$alpha==0.05,]
cat("== sim_power_broad columns:", paste(names(pb),collapse=","), "\n")
cat("== families/params at n=1000:\n")
print(unique(pb[,c("family","param")]))

PART <- c("EF","HL","HL-equalwidth","Pigeon-Heyse")
ALL6 <- c(PART,"Tsiatis","Xie")
TIE <- 0.014
get <- function(fam,par,test,col="power_size_adj"){
  v <- pb[pb$family==fam & pb$param==par & pb$test==test, col]; if(length(v)) v[1] else NA }

cells <- unique(pb[pb$family!="null" ,c("family","param")])
res <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i){
  fam<-cells$family[i]; par<-cells$param[i]
  p3 <- get(fam,par,"DEF.poly3")
  rivP <- sapply(PART, get, fam=fam, par=par); riv6 <- sapply(ALL6, get, fam=fam, par=par)
  allv <- sapply(unique(pb$test[pb$family==fam & pb$param==par]), get, fam=fam, par=par)
  detectable <- max(allv,na.rm=TRUE) > 0.07              # all <= size+0.02 -> undetectable
  saturated  <- min(allv,na.rm=TRUE) >= 0.98
  data.frame(family=fam,param=par,poly3=round(p3,3),
             best_part=round(max(rivP,na.rm=TRUE),3), best_part_name=names(which.max(rivP)),
             best_all6=round(max(riv6,na.rm=TRUE),3), best_all6_name=names(which.max(riv6)),
             detectable=detectable, saturated=saturated,
             beats_part = p3 > max(rivP,na.rm=TRUE)-TIE,
             beats_all6 = p3 > max(riv6,na.rm=TRUE)-TIE)
}))
cat("\n== FULL CELL TABLE (size-adjusted, n=1000) ==\n"); print(res, row.names=FALSE)

el <- read.csv("sim_edge_loses.csv", stringsAsFactors=FALSE)
cat("\n== sim_edge_loses cols:", paste(names(el),collapse=","),"\n")
el1 <- el[el$n==1000 & el$alpha==0.05,]
cat("rough/crossover cells (raw reject rates):\n")
print(unique(el1[,c("scenario","test", grep("reject|power",names(el1),value=TRUE)[1])])[1:40,])

## band endpoints: poly3 gain over max(HL,EF), detectable non-saturated cells where poly3 leads
res2 <- res[res$detectable & !res$saturated,]
gain <- sapply(seq_len(nrow(res2)), function(i){
  fam<-res2$family[i]; par<-res2$param[i]
  base <- max(get(fam,par,"HL"), get(fam,par,"EF"))
  100*(res2$poly3[i]/base - 1)})
res2$gain_pct <- round(gain,1)
cat("\n== detectable non-saturated cells with poly3 gain% over max(HL,EF):\n")
print(res2[order(res2$gain_pct),c("family","param","poly3","best_part","best_part_name","beats_part","beats_all6","gain_pct")], row.names=FALSE)

## closest-basis max gain (poly2/stk)
for (b in c("DEF.poly2","DEF.stukel")){
  g <- sapply(seq_len(nrow(res2)), function(i){
    fam<-res2$family[i]; par<-res2$param[i]
    100*(get(fam,par,b)/max(get(fam,par,"HL"),get(fam,par,"EF")) - 1)})
  cat(sprintf("max gain%% %s over max(HL,EF): %.1f (cell %s/%s)\n", b, max(g,na.rm=TRUE),
      res2$family[which.max(g)], res2$param[which.max(g)]))
}

## == WP-0.2 speedups, single baseline ==
ct <- read.csv("bench_compute_time_summary.csv", stringsAsFactors=FALSE)
base <- ct$time_median[ct$test=="DEF.poly3" & ct$n==1000 & ct$p==4][1]
cat(sprintf("\n== BASELINE EDGE-poly3 n=1000 p=4: %.5f s ==\n", base))
sl <- read.csv("bench_slow_timing_summary.csv", stringsAsFactors=FALSE)
sl1 <- sl[sl$n==1000,]
sl1$ratio_new <- round(sl1$time_median/base)
print(sl1[,c("test","n","time_median","ratio_new")], row.names=FALSE)
## also fast rivals per tab:whentouse
for(t in c("HL","EF","Stukel","Tsiatis","Xie","DEF.poly3")) {
  v<-ct$time_median[ct$test==t & ct$n==1000 & ct$p==4]; if(length(v)) cat(sprintf("%s: %.4f s\n",t,v[1])) }

## imhof tail max
im <- read.csv("imhof_satterthwaite.csv", stringsAsFactors=FALSE)
cat("\nimhof cols:", paste(names(im),collapse=","),"\n")
dcol <- grep("delta", names(im), value=TRUE)[1]; tcol <- grep("tail",names(im),value=TRUE)[1]
cat(sprintf("max|dp| overall=%.5f  tail=%.5f\n", max(abs(im[[dcol]]),na.rm=TRUE),
    if(!is.na(tcol)) max(abs(im[[dcol]][im[[tcol]]==TRUE]),na.rm=TRUE) else NA))

## GLOW G-sweep search
for(f in c("Enhancement_A_grouping.csv","RealPower_power.csv","sim_g_sensitivity.csv")){
  if(file.exists(f)){ d<-read.csv(f,nrows=5); cat("\n--",f,":",paste(names(d),collapse=","),"\n") }}

## edge_degenerate semantics quick look
bs <- read.csv("bench_stukel_failure.csv"); cat("\nbench_stukel_failure cols:",paste(names(bs),collapse=","),"\n")
print(bs[bs$design=="sparse",c("n","stukel_fail_rate","ef_fail_rate", grep("edge",names(bs),value=TRUE))], row.names=FALSE)
