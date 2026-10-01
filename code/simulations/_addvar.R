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
suppressMessages(library(parallel)); alpha<-0.05; lg<-function(p) log(p/(1-p))
PROJ<-edge_path("code", "legacy_not_deposited"); sw<-function(x)suppressWarnings(suppressMessages(x))
load_tests<-function(){invisible(sw(source(file.path(PROJ,"Tsiatis.R"))));invisible(NULL)}
inv_stukel<-function(ev,a1,a2){z<-numeric(length(ev));pos<-ev>=0
  if(any(pos)) z[pos]<- if(abs(a1)<1e-12) ev[pos] else (-1+sqrt(pmax(0,1+2*a1*ev[pos])))/a1
  if(any(!pos))z[!pos]<-if(abs(a2)<1e-12) ev[!pos] else (-1+sqrt(pmax(0,1+2*a2*ev[!pos])))/a2; plogis(z)}
sqb<-function(J) solve(rbind(c(1,-1.5,2.25),c(1,3,9),c(1,-3,9)),c(lg(.05),lg(.95),lg(J)))
linkp<-function(scn,eta) switch(scn, cloglog=1-exp(-exp(eta)), loglog=exp(-exp(-eta)), cauchit=pcauchy(eta), stukel_heavy=inv_stukel(eta,-1,-1))
linkscns<-c("cloglog","loglog","cauchit","stukel_heavy")
gen<-function(scn,n){
  if(scn=="null"){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d))),f=y~x+d)}
  else if(scn%in%linkscns){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,linkp(scn,0.6*x+0.5*d))),f=y~x+d)}
  else if(scn=="quad"){b<-sqb(0.02);x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*x^2))),f=y~x)}
  else if(scn=="cubic"){x<-runif(n,-2.5,2.5);list(d=data.frame(x=x,y=rbinom(n,1,plogis(0.3+0.7*x+0.12*(x^3-3*x)))),f=y~x)}
  else if(scn=="sawtooth"){x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(0.8*x+1.2*(2*(x/1.5-floor(x/1.5+0.5)))))),f=y~x)}
  else if(scn=="osc4"){x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(0.8*x+1.5*sin(4*x)))),f=y~x)}
  else {x<-runif(n,-2.5,2.5);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.9*x-0.5*x*d))),f=y~x+d)}}
ctr<-function(v) v-mean(v)
defcal<-function(y,ph,fit,Zbuild){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-pmin(pmax(sapply(idx,function(I)mean(ph[I])),1e-4),1-1e-4);r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])))/sqrt(Vg);Om<-diag(G)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  Z<-Zbuild(pbar,qlogis(pbar)); Z<-as.matrix(Z); Z<-Z[,apply(Z,2,function(c)sd(c)>1e-8),drop=FALSE]
  Zi<-tryCatch(solve(crossprod(Z)),error=function(e)return(NULL)); if(is.null(Zi))return(NA_real_)
  Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA_real_)
  cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}
Bpoly<-function(k)function(pb,e)poly(pb,k)
Bstk<-function(pb,e)cbind(e,e^2*(e>=0),-e^2*(e<0))
Bev<-function(pb,e)ctr(-log(1-pb)/pb)                  # Prentice/Aranda extreme-value carrier
Bpreg<-function(pb,e)cbind(ctr(log(pb)^2-log(1-pb)^2),ctr(log(pb)^2+log(1-pb)^2))  # Pregibon
ef_omni<-function(y,ph){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);g<-pmin(ceiling(rank(ph,ties.method="first")/(n/10)),10)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),8)}
stuk<-function(fit){e<-predict(fit);d<-fit$data;d$za<-0.5*e^2*(e>=0);d$zb<- -0.5*e^2*(e<0);fa<-suppressWarnings(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial()));pchisq(deviance(fit)-deviance(fa),2,lower.tail=FALSE)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_);if(is.null(v)||length(v)!=1||!is.finite(v))NA_real_ else as.numeric(v)}
nm<-c("poly2","poly3","stk","poly4","evcar","pregibon","EF","Stukel","Tsiatis")
one_rep<-function(scn,n){G0<-gen(scn,n);dat<-G0$d;fit<-suppressWarnings(glm(G0$f,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  saved<-get(".Random.seed",envir=.GlobalEnv)
  res<-c(poly2=sc(defcal(y,ph,fit,Bpoly(2))),poly3=sc(defcal(y,ph,fit,Bpoly(3))),stk=sc(defcal(y,ph,fit,Bstk)),
    poly4=sc(defcal(y,ph,fit,Bpoly(4))),evcar=sc(defcal(y,ph,fit,Bev)),pregibon=sc(defcal(y,ph,fit,Bpreg)),
    EF=sc(ef_omni(y,ph)),Stukel=sc(stuk(fit)),Tsiatis=sc(score_gof_clustering(fit,num_groups=10,y=y)$p_value))
  assign(".Random.seed",saved,envir=.GlobalEnv);res}
scns<-c("null","cloglog","loglog","cauchit","stukel_heavy","quad","cubic","sawtooth","osc4","crossover");n<-1000;REPS<-4000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L));clusterSetRNGStream(cl,55512)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl,c("PROJ","sw","load_tests","alpha","lg","inv_stukel","sqb","linkp","linkscns","gen","ctr","defcal","Bpoly","Bstk","Bev","Bpreg","ef_omni","stuk","sc","one_rep","n","nm"));invisible(clusterEvalQ(cl,load_tests()));on.exit(stopCluster(cl))
pall<-list()
for(scn in scns){clusterExport(cl,"scn",envir=environment());m<-parSapply(cl,1:REPS,function(i)one_rep(scn,n))
  df<-as.data.frame(t(m[nm,]));names(df)<-nm;df<-cbind(scenario=scn,rep=1:REPS,df);pall[[length(pall)+1]]<-df;cat(scn,"done\n")}
write.csv(do.call(rbind,pall),"AddVar_pvalues.csv",row.names=FALSE);cat("saved AddVar_pvalues.csv\n")
