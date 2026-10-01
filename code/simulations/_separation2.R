suppressMessages(library(CompQuadForm))
ef_dir_cal<-function(y,ph,fit,basis,k){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])));U<-U/sqrt(Vg);Om<-diag(G)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  if(basis=="poly")Z<-as.matrix(poly(pbar,k)) else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  if(ncol(Z)<1)return(NA);Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA)
  p<-tryCatch(CompQuadForm::imhof(S,lam)$Qq,error=function(e)NA);min(max(p,0),1)}
stuk_diag<-function(fit){e<-predict(fit);d<-fit$data;d$za<-0.5*e^2*(e>=0);d$zb<- -0.5*e^2*(e<0);sep<-FALSE
  fa<-tryCatch(withCallingHandlers(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial()),
    warning=function(w){if(grepl("did not converge|numerically 0 or 1",conditionMessage(w)))sep<<-TRUE;invokeRestart("muffleWarning")}),
    error=function(er)NULL)
  if(is.null(fa))return(c(p=NA,fail=1)); if(!isTRUE(fa$converged))sep<-TRUE; pf<-fitted(fa); if(any(pf<1e-7|pf>1-1e-7))sep<-TRUE
  c(p=pchisq(deviance(fit)-deviance(fa),2,lower.tail=FALSE),fail=as.numeric(sep))}
sim<-function(n,b0=-2,b1=1.8,REPS=4000,seed=7){set.seed(seed);R<-matrix(NA,REPS,5,dimnames=list(NULL,c("p2","p3","ps","St","Sf")));er<-numeric(REPS)
  for(i in 1:REPS){x<-runif(n,-4,4);y<-rbinom(n,1,plogis(b0+b1*x));er[i]<-mean(y); if(length(unique(y))<2)next
    fit<-tryCatch(suppressWarnings(glm(y~x,family=binomial())),error=function(e)NULL); if(is.null(fit))next
    ph<-as.numeric(fitted(fit)); R[i,"p2"]<-tryCatch(ef_dir_cal(y,ph,fit,"poly",2),error=function(e)NA)
    R[i,"p3"]<-tryCatch(ef_dir_cal(y,ph,fit,"poly",3),error=function(e)NA); R[i,"ps"]<-tryCatch(ef_dir_cal(y,ph,fit,"stukel",3),error=function(e)NA)
    sd<-stuk_diag(fit);R[i,"St"]<-sd["p"];R[i,"Sf"]<-sd["fail"]}
  a<-.05;data.frame(n=n,evrate=round(mean(er),3),DEFp2.size=mean(R[,"p2"]<a,na.rm=T),DEFp3.size=mean(R[,"p3"]<a,na.rm=T),DEFstk.size=mean(R[,"ps"]<a,na.rm=T),
    DEF.NArate=round(mean(is.na(R[,"p3"])),4),Stukel.size=mean(R[,"St"]<a,na.rm=T),Stukel.failrate=round(mean(R[,"Sf"],na.rm=T),3))}
out<-do.call(rbind,lapply(c(60,100,150,250,400),sim))
write.csv(out,"Sep_results.csv",row.names=FALSE); print(out,row.names=FALSE,digits=3)
