# Imhof (exact) vs Satterthwaite calibration: confirm size at n in {200,500,1000} for DEF-poly3 & DEF-stk3.
suppressMessages({library(parallel); library(CompQuadForm)}); alpha<-0.05
gen_null<-function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d)))}
defp<-function(y,ph,fit,basis,k,method){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])));U<-U/sqrt(Vg);Om<-diag(G)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  if(basis=="poly")Z<-as.matrix(poly(pbar,k)) else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  if(ncol(Z)<1)return(NA);Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA)
  if(method=="imhof"){p<-tryCatch(CompQuadForm::imhof(S,lam)$Qq,error=function(e)NA);min(max(p,0),1)}
  else{cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}}
sc<-function(e){v<-tryCatch(e,error=function(x)NA);if(is.null(v)||length(v)!=1||!is.finite(v))NA else as.numeric(v)}
one<-function(n){dat<-gen_null(n);fit<-suppressWarnings(glm(y~x+d,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  c(p3.imhof=sc(defp(y,ph,fit,"poly",3,"imhof")),p3.satt=sc(defp(y,ph,fit,"poly",3,"satt")),
    stk.imhof=sc(defp(y,ph,fit,"stukel",3,"imhof")),stk.satt=sc(defp(y,ph,fit,"stukel",3,"satt")))}
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L));clusterSetRNGStream(cl,20250911)
clusterEvalQ(cl,suppressMessages(library(CompQuadForm)));clusterExport(cl,c("gen_null","defp","sc","one"));on.exit(stopCluster(cl))
cat("Empirical SIZE (target 0.05):\n"); cat(sprintf("%5s | %9s %8s | %9s %8s\n","n","poly3.imhof","poly3.satt","stk3.imhof","stk3.satt"))
for(n in c(200,500,1000)){clusterExport(cl,"n",envir=environment());m<-parSapply(cl,1:4000,function(i)one(n))
  s<-rowMeans(m<0.05,na.rm=TRUE);cat(sprintf("%5d | %9.4f %8.4f | %9.4f %8.4f\n",n,s["p3.imhof"],s["p3.satt"],s["stk.imhof"],s["stk.satt"]))}
