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
suppressMessages({library(parallel)}); alpha<-0.05
PROJ<-edge_path("code", "legacy_not_deposited"); sw<-function(x)suppressWarnings(suppressMessages(x))
load_tests<-function(){suppressMessages({library(ResourceSelection)})
  invisible(sw(source(file.path(PROJ,"pigeonheyse.R"))));invisible(sw(source(file.path(PROJ,"Hosmer (H) (equal width interval).R"))));invisible(sw(source(file.path(PROJ,"Tsiatis.R"))));invisible(sw(source(file.path(PROJ,"Xie.R"))));invisible(sw(source(file.path(PROJ,"PR_test_only.R"))));invisible(NULL)}
ef_dir_cal<-function(y,ph,fit,basis,k){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])));U<-U/sqrt(Vg);Om<-diag(G)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  if(basis=="poly")Z<-as.matrix(poly(pbar,k)) else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  if(ncol(Z)<1)return(NA_real_);Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA_real_)
  cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}
ef_omni<-function(y,ph){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);g<-pmin(ceiling(rank(ph,ties.method="first")/(n/10)),10)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),8)}
hleqw<-function(y,ph){g<-cut(ph,seq(0,1,length.out=11),include.lowest=TRUE,labels=FALSE);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);nn<-tapply(y,g,length);k<-!is.na(o);o<-o[k];e<-e[k];nn<-nn[k];st<-sum((o-e)^2/(e+1e-10)+((nn-o)-(nn-e))^2/((nn-e)+1e-10));if(length(o)>2)1-pchisq(st,length(o)-2) else NA_real_}
stuk<-function(fit){e<-predict(fit);d<-fit$data;d$za<-0.5*e^2*(e>=0);d$zb<- -0.5*e^2*(e<0);fa<-suppressWarnings(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial()));pchisq(deviance(fit)-deviance(fa),2,lower.tail=FALSE)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_);if(is.null(v)||length(v)!=1||!is.finite(v))NA_real_ else as.numeric(v)}
genX<-function(scn,n){
 if(scn=="omit_z2"){x<-runif(n,-2.5,2.5);z<-rnorm(n);eta<-0.3+0.6*x+0.4*z+0.45*z^2;list(d=data.frame(x=x,z=z,y=rbinom(n,1,plogis(eta))),f=y~x+z,cat=NA)}
 else if(scn=="omit_x2_strong"){x<-runif(n,-2.5,2.5);eta<-0.3+0.7*x+0.30*x^2;list(d=data.frame(x=x,y=rbinom(n,1,plogis(eta))),f=y~x,cat=NA)}
 else if(scn=="omit_x2z2"){x<-runif(n,-2.5,2.5);z<-rnorm(n);eta<-0.2+0.5*x+0.28*x^2+0.4*z+0.30*z^2;list(d=data.frame(x=x,z=z,y=rbinom(n,1,plogis(eta))),f=y~x+z,cat=NA)}
 else if(scn=="omit_intxz"){x<-runif(n,-2.5,2.5);z<-rnorm(n);eta<-0.4*x+0.4*z+0.8*x*z;list(d=data.frame(x=x,z=z,y=rbinom(n,1,plogis(eta))),f=y~x+z,cat=NA)}
 else if(scn=="omit_2int"){x<-runif(n,-2.5,2.5);d<-rbinom(n,1,.5);z<-rnorm(n);eta<-0.4*x+0.4*d+0.4*z+0.6*x*d+0.6*x*z;list(d=data.frame(x=x,d=d,z=z,y=rbinom(n,1,plogis(eta))),f=y~x+d+z,cat="d")}
 else {x<-runif(n,-3,3);eta<- -0.3+0.5*x+0.7*abs(x);list(d=data.frame(x=x,y=rbinom(n,1,plogis(eta))),f=y~x,cat=NA)}}  # vshape
nm<-c("DEF.poly2","DEF.poly3","DEF.stk3","EF","HL","HLeqw","PH","Tsiatis","Xie","PR","Stukel")
one_rep<-function(scn,n){G0<-genX(scn,n);dat<-G0$d;fit<-suppressWarnings(glm(G0$f,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  saved<-get(".Random.seed",envir=.GlobalEnv);f2<-fit;f2$predicted_probs<-ph
  res<-c(DEF.poly2=sc(ef_dir_cal(y,ph,fit,"poly",2)),DEF.poly3=sc(ef_dir_cal(y,ph,fit,"poly",3)),DEF.stk3=sc(ef_dir_cal(y,ph,fit,"stukel",3)),
    EF=sc(ef_omni(y,ph)),HL=sc(hoslem.test(y,ph,g=10)$p.value),HLeqw=sc(hleqw(y,ph)),
    PH=sc(pigeon_heyse_test(data.frame(y=y),fit,g=10)$p_value),Tsiatis=sc(score_gof_clustering(fit,num_groups=10,y=y)$p_value),
    Xie=sc(as.numeric(XieGoodnessOfFitTest(dat,f2))),PR=if(is.na(G0$cat))NA_real_ else sc(pr_test(dat,"y",G0$cat,ph)$p_value),Stukel=sc(stuk(fit)))
  assign(".Random.seed",saved,envir=.GlobalEnv);res}
scns<-c("omit_z2","omit_x2_strong","omit_x2z2","omit_intxz","omit_2int","vshape");n<-1000;REPS<-5000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L));clusterSetRNGStream(cl,5150607)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl,c("PROJ","sw","load_tests","alpha","ef_dir_cal","ef_omni","hleqw","stuk","sc","genX","one_rep","n","nm"));invisible(clusterEvalQ(cl,load_tests()));on.exit(stopCluster(cl))
pall<-list();summ<-list()
for(scn in scns){clusterExport(cl,"scn",envir=environment());m<-parSapply(cl,1:REPS,function(i)one_rep(scn,n))
  df<-as.data.frame(t(m[nm,]));names(df)<-nm;df<-cbind(scenario=scn,rep=1:REPS,df);pall[[length(pall)+1]]<-df
  pw<-rowMeans(m[nm,,drop=FALSE]<alpha,na.rm=TRUE);summ[[length(summ)+1]]<-data.frame(scenario=scn,t(round(pw,4)));cat(scn,"done\n")}
write.csv(do.call(rbind,pall),"MoreDEF_pvalues.csv",row.names=FALSE)
S<-do.call(rbind,summ);write.csv(S,"MoreDEF_power.csv",row.names=FALSE);print(S,row.names=FALSE)
