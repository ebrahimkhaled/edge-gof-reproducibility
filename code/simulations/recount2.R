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
suppressMessages(library(stats))
SIM <- edge_path("code/simulations"); setwd(SIM)
pb <- read.csv("sim_power_broad.csv"); pb <- pb[pb$n==1000 & pb$alpha==0.05,]
g <- function(fam,par,t) { v<-pb$power_size_adj[pb$family==fam&pb$param==par&pb$test==t]; if(length(v)) round(v[1],3) else NA }
PART<-c("EF","HL","HL-equalwidth","Pigeon-Heyse")
## saturation check cells
for(cc in list(c("quad","0.05"), c("binint","0.7"))) {
  cat(cc[1],cc[2],": poly3",g(cc[1],cc[2],"DEF.poly3"),"|",paste(PART,sapply(PART,g,fam=cc[1],par=cc[2]),collapse=" | "),"\n")}
## gains poly3 vs max(HL,EF) on the winning detectable cells
cells <- list(c("link","cloglog"),c("link","stukel_heavy"),c("link","stukel_light"),c("link","stukel_asym"),
  c("quad","0.02"),c("quad","0.03"),c("quad","0.05"),c("binint","0.2"),c("binint","0.3"),c("binint","0.5"),
  c("contint","0.3"),c("contint","0.5"),c("contint","0.7"),c("rough","bump"))
for(cc in cells){ b<-max(g(cc[1],cc[2],"HL"),g(cc[1],cc[2],"EF"))
  cat(sprintf("%s/%s: poly3=%.3f base=%.3f gain=%+.0f%%  poly2=%.3f stk=%.3f\n", cc[1],cc[2],
      g(cc[1],cc[2],"DEF.poly3"), b, 100*(g(cc[1],cc[2],"DEF.poly3")/b-1), g(cc[1],cc[2],"DEF.poly2"), g(cc[1],cc[2],"DEF.stukel")))}
## speedups single baseline
ct <- read.csv("bench_compute_time_summary.csv"); base <- ct$time_median[ct$test=="DEF.poly3"&ct$n==1000&ct$p==4][1]
cat(sprintf("\nBASELINE poly3 = %.5f s\n", base))
sl <- read.csv("bench_slow_timing_summary.csv"); sl1<-sl[sl$n==1000,]
for(i in seq_len(nrow(sl1))) cat(sprintf("%-12s %8.2f s  ratio=%s\n", sl1$test[i], sl1$time_median[i], format(round(sl1$time_median[i]/base), big.mark=",")))
for(t in c("HL","EF","Stukel","Tsiatis","Xie")) cat(sprintf("%s=%.4f s  ", t, ct$time_median[ct$test==t&ct$n==1000&ct$p==4][1])); cat("\n")
## imhof
im <- read.csv("imhof_satterthwaite.csv")
cat("imhof cols:", paste(names(im),collapse=","),"\n")
dc <- grep("delta",names(im),value=TRUE)[1]; tc <- grep("tail",names(im),value=TRUE)[1]
cat(sprintf("max|dp| overall=%.5f", max(abs(im[[dc]]),na.rm=TRUE)))
if(!is.na(tc)) cat(sprintf("  tailTRUE=%.5f", max(abs(im[[dc]][im[[tc]]=="TRUE"|im[[tc]]==TRUE]),na.rm=TRUE))); cat("\n")
## glow sweep search
for(f in c("Enhancement_A_grouping.csv","RealPower_power.csv")) if(file.exists(f)){d<-read.csv(f); cat("\n",f," cols:",paste(names(d),collapse=","),"; G vals:",paste(unique(d$G),collapse=","),"\n")}
## stukel failure cols
bs<-read.csv("bench_stukel_failure.csv"); cat("\nstukelfail cols:",paste(names(bs),collapse=","),"\n")
