suppressMessages({library(parallel); library(BAGofT)})
gen <- function(scn,n){x<-runif(n,-3,3); d<-rbinom(n,1,.5); eta<-0.6*x+0.5*d
  p<-switch(scn, null=plogis(eta), cloglog=1-exp(-exp(eta)), loglog=exp(-exp(-eta)), probit=pnorm(eta))
  data.frame(x=x,d=d,y=rbinom(n,1,p))}
one_bag <- function(scn,n){dat<-gen(scn,n)
  cap<-capture.output(b<-suppressMessages(BAGofT(testModel=testGlmBi(formula=y~x+d,link="logit"),data=dat,nsplits=1)))
  tryCatch(b$p.value, error=function(e) NA_real_)}
ns<-c(500,1000); scns<-c("null","cloglog","loglog","probit"); REPS<-as.integer(Sys.getenv("REPS","300"))
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L)); clusterSetRNGStream(cl,20250911)
invisible(clusterEvalQ(cl, suppressMessages(library(BAGofT)))); clusterExport(cl,c("gen","one_bag"))
cat(sprintf("=== BAGofT power (nsplits=1, REPS=%d) ===\n",REPS))
out<-list()
for(scn in scns) for(n in ns){clusterExport(cl,c("scn","n"),envir=environment())
  p<-parSapply(cl,1:REPS,function(i) one_bag(scn,n)); pw<-mean(p<0.05,na.rm=TRUE)
  out[[length(out)+1]]<-data.frame(scenario=scn,n=n,BAGofT=round(pw,4))
  cat(sprintf("%-8s n=%-5d BAGofT=%.3f\n",scn,n,pw))}
stopCluster(cl); write.csv(do.call(rbind,out),"Benchmark_BAGofT.csv",row.names=FALSE)
cat("Saved Benchmark_BAGofT.csv\n")
