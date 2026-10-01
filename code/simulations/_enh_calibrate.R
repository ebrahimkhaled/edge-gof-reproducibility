# Closed-form calibrated directed-EF across ALL link scenarios. Null row = SIZE (gate).
suppressMessages(library(parallel)); alpha<-0.05
inv_stukel<-function(ev,a1,a2){z<-numeric(length(ev));pos<-ev>=0
  if(any(pos))  z[pos]<- if(abs(a1)<1e-12) ev[pos]  else (-1+sqrt(pmax(0,1+2*a1*ev[pos])))/a1
  if(any(!pos)) z[!pos]<-if(abs(a2)<1e-12) ev[!pos] else (-1+sqrt(pmax(0,1+2*a2*ev[!pos])))/a2; plogis(z)}
gen<-function(scn,n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d
  p<-switch(scn, null=plogis(eta), cloglog=1-exp(-exp(eta)), loglog=exp(-exp(-eta)), probit=pnorm(eta),
            cauchit=pcauchy(eta), robit_t4=pt(eta,df=4), scobit2=plogis(eta)^2, scobit_half=plogis(eta)^0.5,
            stukel_heavy=inv_stukel(eta,-1,-1), stukel_light=inv_stukel(eta,1,1), stukel_asym=inv_stukel(eta,-1,1), stop("bad"))
  data.frame(x=x,d=d,y=rbinom(n,1,p))}
ef_dir_cal<-function(y,ph,fit,basis,k){
  ph<-pmin(pmax(ph,1e-6),1-1e-6); n<-length(y); G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G); idx<-split(seq_len(n),grp); Wii<-ph*(1-ph)
  og<-sapply(idx,function(I) sum(y[I])); eg<-sapply(idx,function(I) sum(ph[I]))
  Vg<-sapply(idx,function(I) sum(Wii[I])); pbar<-sapply(idx,function(I) mean(ph[I])); r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit); U<-t(sapply(idx,function(I) colSums(Wii[I]*X[I,,drop=FALSE]))); U<-U/sqrt(Vg)
  Omega<-diag(G) - U %*% solve(crossprod(X, Wii*X)) %*% t(U)
  if(basis=="poly") Z<-poly(pbar,k) else { eta<-qlogis(pbar); Z<-cbind(eta,eta^2*(eta>=0),-eta^2*(eta<0)) }
  Z<-as.matrix(Z); ZtZi<-solve(crossprod(Z)); Ztr<-crossprod(Z,r); S<-as.numeric(t(Ztr)%*%ZtZi%*%Ztr)
  lam<-Re(eigen(ZtZi %*% (t(Z)%*%Omega%*%Z), only.values=TRUE)$values); lam<-lam[lam>1e-9]
  if(length(lam)==0) return(NA_real_); cc<-sum(lam^2)/sum(lam); nu<-sum(lam)^2/sum(lam^2); 1-pchisq(S/cc, nu)}
ef_omni<-function(y,ph){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/10)),10)
  og<-tapply(y,grp,sum);eg<-tapply(ph,grp,sum);ng<-tapply(y,grp,length);pb<-as.numeric(tapply(ph,grp,mean));Vg<-ng*pb*(1-pb);oe<-as.numeric(og-eg)
  1-pchisq(sum(oe^2/Vg)-sum((1-2*pb)*oe/Vg),8)}
hleqw_p<-function(y,ph){gw<-cut(ph,breaks=seq(0,1,length.out=11),include.lowest=TRUE,labels=FALSE)
  o<-tapply(y,gw,sum);e<-tapply(ph,gw,sum);nn<-tapply(y,gw,length);k<-!is.na(o);o<-o[k];e<-e[k];nn<-nn[k]
  st<-sum((o-e)^2/(e+1e-10)+((nn-o)-(nn-e))^2/((nn-e)+1e-10)); if(length(o)>2) 1-pchisq(st,length(o)-2) else NA}
stukel_p<-function(fit){eta<-predict(fit);d<-fit$data;d$za<-0.5*eta^2*(eta>=0);d$zb<- -0.5*eta^2*(eta<0)
  fa<-suppressWarnings(glm(y~x+d+za+zb,data=d,family=binomial()));pchisq(deviance(fit)-deviance(fa),df=2,lower.tail=FALSE)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_); if(is.null(v)||length(v)!=1||!is.finite(v)) NA_real_ else as.numeric(v)}
one_rep<-function(scn,n){dat<-gen(scn,n); fit<-suppressWarnings(glm(y~x+d,data=dat,family=binomial())); ph<-as.numeric(fitted(fit)); y<-dat$y
  c("EFd.poly2"=sc(ef_dir_cal(y,ph,fit,"poly",2)),"EFd.poly3"=sc(ef_dir_cal(y,ph,fit,"poly",3)),"EFd.stk3"=sc(ef_dir_cal(y,ph,fit,"stukel",3)),
    "EF.omni"=sc(ef_omni(y,ph)),"HLeqw"=sc(hleqw_p(y,ph)),"Stukel"=sc(stukel_p(fit)))}
scns<-c("null","cloglog","loglog","probit","cauchit","robit_t4","scobit2","scobit_half","stukel_heavy","stukel_light","stukel_asym")
ns<-c(500,1000); REPS<-3000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L)); clusterSetRNGStream(cl,20250911)
clusterExport(cl,c("alpha","inv_stukel","gen","ef_dir_cal","ef_omni","hleqw_p","stukel_p","sc","one_rep")); on.exit(stopCluster(cl))
cat(sprintf("%-13s %5s | %9s %9s %8s %8s %6s %7s\n","scenario","n","EFd.poly2","EFd.poly3","EFd.stk3","EF.omni","HLeqw","Stukel"))
out<-list()
for(scn in scns) for(n in ns){clusterExport(cl,c("scn","n"),envir=environment()); m<-parSapply(cl,1:REPS,function(i) one_rep(scn,n))
  pw<-rowMeans(m<alpha,na.rm=TRUE); out[[length(out)+1]]<-data.frame(scenario=scn,n=n,t(round(pw,4)))
  cat(sprintf("%-13s %5d | %9.3f %9.3f %8.3f %8.3f %6.3f %7.3f\n",scn,n,pw["EFd.poly2"],pw["EFd.poly3"],pw["EFd.stk3"],pw["EF.omni"],pw["HLeqw"],pw["Stukel"]))}
write.csv(do.call(rbind,out),"Enhancement_calibrated_all.csv",row.names=FALSE); cat("\n(null row = SIZE, must be ~0.05; others = power. Calibrated directed-EF = closed-form weighted-chi2.)\n")
