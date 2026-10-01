# Null-distribution validation: under H0 the DEF p-values must be Uniform(0,1).
# Stores per-rep null p-values at several n for a PP-plot and a size-by-level table.
suppressMessages(library(parallel))
ef_dir_cal<-function(y,ph,fit,basis,k){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])))/sqrt(Vg);Om<-diag(length(idx))-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  if(basis=="poly")Z<-as.matrix(poly(pbar,k)) else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  if(ncol(Z)<1)return(NA_real_);Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA_real_)
  cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_);if(is.null(v)||length(v)!=1||!is.finite(v))NA_real_ else as.numeric(v)}
one<-function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);dat<-data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d)))
  fit<-suppressWarnings(glm(y~x+d,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  c(poly2=sc(ef_dir_cal(y,ph,fit,"poly",2)),poly3=sc(ef_dir_cal(y,ph,fit,"poly",3)),stk=sc(ef_dir_cal(y,ph,fit,"stukel",3)))}
ns<-c(200,500,1000);REPS<-10000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L));clusterSetRNGStream(cl,20260620)
clusterExport(cl,c("ef_dir_cal","sc","one"));on.exit(stopCluster(cl))
out<-list()
for(nn in ns){clusterExport(cl,"nn",envir=environment());m<-parSapply(cl,1:REPS,function(i)one(nn))
  df<-as.data.frame(t(m));names(df)<-c("poly2","poly3","stk");df<-cbind(n=nn,rep=1:REPS,df);out[[length(out)+1]]<-df;cat("n=",nn,"done\n")}
write.csv(do.call(rbind,out),"Null_pvalues.csv",row.names=FALSE);cat("saved Null_pvalues.csv\n")
