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
suppressMessages({library(parallel); library(ebrahim.gof); library(ResourceSelection)})
PROJ<-edge_path("code", "legacy_not_deposited"); sw<-function(e) suppressWarnings(suppressMessages(e))
load_tests<-function(){suppressMessages({library(ebrahim.gof);library(ResourceSelection);library(MASS);library(dplyr)})
  invisible(sw(source(file.path(PROJ,"pigeonheyse.R")))); invisible(sw(source(file.path(PROJ,"Hosmer (H) (equal width interval).R"))))
  invisible(sw(source(file.path(PROJ,"Tsiatis.R")))); invisible(sw(source(file.path(PROJ,"Xie.R")))); invisible(sw(source(file.path(PROJ,"PR_test_only.R")))); invisible(NULL)}
load_tests(); G<-10; alpha<-0.05
sc<-function(expr){v<-tryCatch(suppressWarnings(expr),error=function(e) NA_real_); if(is.null(v)||length(v)!=1||!is.finite(v)) NA_real_ else as.numeric(v)}
inv_stukel<-function(ev,a1,a2){z<-numeric(length(ev));pos<-ev>=0
  if(any(pos))  z[pos]<- if(abs(a1)<1e-12) ev[pos]  else (-1+sqrt(pmax(0,1+2*a1*ev[pos])))/a1
  if(any(!pos)) z[!pos]<-if(abs(a2)<1e-12) ev[!pos] else (-1+sqrt(pmax(0,1+2*a2*ev[!pos])))/a2; plogis(z)}
gen<-function(scn,n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d
  p<-switch(scn, cauchit=pcauchy(eta), robit_t4=pt(eta,df=4), scobit2=plogis(eta)^2, scobit_half=plogis(eta)^0.5,
            stukel_heavy=inv_stukel(eta,-1,-1), stukel_light=inv_stukel(eta,1,1), stukel_asym=inv_stukel(eta,-1,1), stop("bad"))
  data.frame(x=x,d=d,y=rbinom(n,1,p))}
one_rep<-function(scn,n){dat<-gen(scn,n); saved<-get(".Random.seed",envir=.GlobalEnv)
  fit<-suppressWarnings(glm(y~x+d,data=dat,family=binomial())); ph<-as.numeric(fitted(fit))
  z<-ef.gof(y=dat$y,predicted_probs=ph,G=G)$Test_Statistic; f2<-fit; f2$predicted_probs<-ph
  res<-c(EF=sc(1-pchisq(z*sqrt(2*(G-2))+(G-2),G-2)),HL=sc(hoslem.test(dat$y,ph,g=G)$p.value),
    HLeqw=sc(hosmer_lemeshow_equal_intervals(ph,dat$y,num_groups=G)$p.value),PH=sc(pigeon_heyse_test(data.frame(y=dat$y),fit,g=G)$p_value),
    Tsiatis=sc(score_gof_clustering(fit,num_groups=G,y=dat$y)$p_value),Xie=sc(as.numeric(XieGoodnessOfFitTest(dat,f2))),PR=sc(pr_test(dat,"y","d",ph)$p_value))
  assign(".Random.seed",saved,envir=.GlobalEnv); res}
scns<-c("cauchit","robit_t4","scobit2","scobit_half","stukel_heavy","stukel_light","stukel_asym"); ns<-c(500,1000,2000); REPS<-2000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L)); clusterSetRNGStream(cl,424242)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl,c("PROJ","sw","load_tests","G","alpha","sc","inv_stukel","gen","one_rep")); invisible(clusterEvalQ(cl,load_tests())); on.exit(stopCluster(cl))
cat(sprintf("%-12s %5s | %6s %6s %6s %6s %7s %6s %6s\n","scenario","n","EF","HL","HLeqw","PH","Tsiatis","Xie","PR")); out<-list()
for(scn in scns) for(n in ns){clusterExport(cl,c("scn","n"),envir=environment()); m<-parSapply(cl,1:REPS,function(i) one_rep(scn,n))
  pw<-rowMeans(m<alpha,na.rm=TRUE); out[[length(out)+1]]<-data.frame(scenario=scn,n=n,t(round(pw,4)))
  cat(sprintf("%-12s %5d | %6.3f %6.3f %6.3f %6.3f %7.3f %6.3f %6.3f\n",scn,n,pw["EF"],pw["HL"],pw["HLeqw"],pw["PH"],pw["Tsiatis"],pw["Xie"],pw["PR"]))}
write.csv(do.call(rbind,out),"Benchmark_partition_extra.csv",row.names=FALSE); cat("\nSaved Benchmark_partition_extra.csv\n")
