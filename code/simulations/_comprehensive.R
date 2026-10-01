# Comprehensive: 3 DEF variants + EF-omni + all partition tests (no BAGofT) + Stukel, all scenarios, n=1000.
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
suppressMessages({library(parallel); library(ResourceSelection)})
PROJ<-edge_path("code", "legacy_not_deposited"); sw<-function(e) suppressWarnings(suppressMessages(e))
load_tests<-function(){suppressMessages({library(ResourceSelection);library(MASS);library(dplyr)})
  invisible(sw(source(file.path(PROJ,"pigeonheyse.R")))); invisible(sw(source(file.path(PROJ,"Hosmer (H) (equal width interval).R"))))
  invisible(sw(source(file.path(PROJ,"Tsiatis.R")))); invisible(sw(source(file.path(PROJ,"Xie.R")))); invisible(sw(source(file.path(PROJ,"PR_test_only.R")))); invisible(NULL)}
load_tests(); G<-10; alpha<-0.05; lg<-function(p) log(p/(1-p))
inv_stukel<-function(ev,a1,a2){z<-numeric(length(ev));pos<-ev>=0
  if(any(pos)) z[pos]<- if(abs(a1)<1e-12) ev[pos] else (-1+sqrt(pmax(0,1+2*a1*ev[pos])))/a1
  if(any(!pos))z[!pos]<-if(abs(a2)<1e-12) ev[!pos] else (-1+sqrt(pmax(0,1+2*a2*ev[!pos])))/a2; plogis(z)}
sqb<-function(J) solve(rbind(c(1,-1.5,2.25),c(1,3,9),c(1,-3,9)),c(lg(.05),lg(.95),lg(J)))
sib<-function(I){e0<-lg(.1);e1<-lg(.2);e2<-lg(.2+I);b0<-(e0+e1)/2;b1<-(e1-e0)/6;b3<-(e2-(b0+3*b1))/6;c(b0,b1,3*b3,b3)}
scb<-function(K) solve(rbind(c(1,-3,-2,6),c(1,-3,0,0),c(1,3,0,0),c(1,3,2,6)),c(lg(.1),lg(.1),lg(.2),lg(.2+K)))
linkp<-function(scn,eta) switch(scn, cloglog=1-exp(-exp(eta)), loglog=exp(-exp(-eta)), probit=pnorm(eta), cauchit=pcauchy(eta),
  robit_t4=pt(eta,df=4), scobit2=plogis(eta)^2, scobit_half=plogis(eta)^0.5,
  stukel_heavy=inv_stukel(eta,-1,-1), stukel_light=inv_stukel(eta,1,1), stukel_asym=inv_stukel(eta,-1,1))
linkscns<-c("cloglog","loglog","probit","cauchit","robit_t4","scobit2","scobit_half","stukel_heavy","stukel_light","stukel_asym")
gen<-function(scn,n){
  if(scn=="null"){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d))),f=y~x+d,cat="d")}
  else if(scn%in%linkscns){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,linkp(scn,0.6*x+0.5*d))),f=y~x+d,cat="d")}
  else if(scn=="quad"){b<-sqb(0.02);x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*x^2))),f=y~x,cat=NA)}
  else if(scn=="int_bin"){b<-sib(0.5);x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*d+b[4]*x*d))),f=y~x+d,cat="d")}
  else {b<-scb(0.5);x<-runif(n,-3,3);z<-rnorm(n);list(d=data.frame(x=x,z=z,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*z+b[4]*x*z))),f=y~x+z,cat=NA)}}
ef_dir_cal<-function(y,ph,fit,basis,k){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y)
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])));U<-U/sqrt(Vg)
  Om<-diag(G)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  if(basis=="poly")Z<-as.matrix(poly(pbar,k)) else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  if(ncol(Z)<1)return(NA_real_); Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA_real_)
  cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}
ef_omni<-function(y,ph){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);g<-pmin(ceiling(rank(ph,ties.method="first")/(n/10)),10)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),8)}
hleqw<-function(y,ph){g<-cut(ph,seq(0,1,length.out=11),include.lowest=TRUE,labels=FALSE);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);nn<-tapply(y,g,length);k<-!is.na(o);o<-o[k];e<-e[k];nn<-nn[k];st<-sum((o-e)^2/(e+1e-10)+((nn-o)-(nn-e))^2/((nn-e)+1e-10));if(length(o)>2)1-pchisq(st,length(o)-2) else NA}
stuk<-function(fit){e<-predict(fit);d<-fit$data;d$za<-0.5*e^2*(e>=0);d$zb<- -0.5*e^2*(e<0);fa<-suppressWarnings(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial()));pchisq(deviance(fit)-deviance(fa),2,lower.tail=FALSE)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_);if(is.null(v)||length(v)!=1||!is.finite(v))NA_real_ else as.numeric(v)}
one_rep<-function(scn,n){G0<-gen(scn,n);dat<-G0$d;fit<-suppressWarnings(glm(G0$f,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  saved<-get(".Random.seed",envir=.GlobalEnv); f2<-fit; f2$predicted_probs<-ph
  res<-c(DEF.poly2=sc(ef_dir_cal(y,ph,fit,"poly",2)),DEF.poly3=sc(ef_dir_cal(y,ph,fit,"poly",3)),DEF.stk3=sc(ef_dir_cal(y,ph,fit,"stukel",3)),
    EF.omni=sc(ef_omni(y,ph)),HL=sc(hoslem.test(y,ph,g=10)$p.value),HLeqw=sc(hleqw(y,ph)),
    PH=sc(pigeon_heyse_test(data.frame(y=y),fit,g=10)$p_value),Tsiatis=sc(score_gof_clustering(fit,num_groups=10,y=y)$p_value),
    Xie=sc(as.numeric(XieGoodnessOfFitTest(dat,f2))),PR=if(is.na(G0$cat))NA_real_ else sc(pr_test(dat,"y",G0$cat,ph)$p_value),Stukel=sc(stuk(fit)))
  assign(".Random.seed",saved,envir=.GlobalEnv); res}
scns<-c("null",linkscns,"quad","int_bin","int_cont"); REPS<-2000; n<-1000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L));clusterSetRNGStream(cl,20250911)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl,c("PROJ","sw","load_tests","G","alpha","lg","inv_stukel","sqb","sib","scb","linkp","linkscns","gen","ef_dir_cal","ef_omni","hleqw","stuk","sc","one_rep","n"))
invisible(clusterEvalQ(cl,load_tests())); on.exit(stopCluster(cl))
nm<-c("DEF.poly2","DEF.poly3","DEF.stk3","EF.omni","HL","HLeqw","PH","Tsiatis","Xie","PR","Stukel"); out<-list()
for(scn in scns){clusterExport(cl,"scn",envir=environment());m<-parSapply(cl,1:REPS,function(i)one_rep(scn,n))
  pw<-rowMeans(m<alpha,na.rm=TRUE); out[[length(out)+1]]<-data.frame(scenario=scn,t(round(pw[nm],3)))
  cat(sprintf("%-12s done\n",scn))}
res<-do.call(rbind,out); write.csv(res,"Comprehensive_n1000.csv",row.names=FALSE); print(res, row.names=FALSE)
