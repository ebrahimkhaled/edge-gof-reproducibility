# Enhancement B2: SIZE-ADJUSTED power. Calibrate every test to exactly 0.05 via its null
# distribution (same covariate design under null & alt), then compare intrinsic power.
suppressMessages({library(parallel)})
alpha<-0.05
inv_stukel<-function(ev,a1,a2){z<-numeric(length(ev));pos<-ev>=0
  if(any(pos))  z[pos]<- if(abs(a1)<1e-12) ev[pos]  else (-1+sqrt(pmax(0,1+2*a1*ev[pos])))/a1
  if(any(!pos)) z[!pos]<-if(abs(a2)<1e-12) ev[!pos] else (-1+sqrt(pmax(0,1+2*a2*ev[!pos])))/a2; plogis(z)}
gen<-function(scn,n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d
  p<-switch(scn, null=plogis(eta), cloglog=1-exp(-exp(eta)), loglog=exp(-exp(-eta)),
            stukel_heavy=inv_stukel(eta,-1,-1), stukel_asym=inv_stukel(eta,-1,1), cauchit=pcauchy(eta), stop("bad"))
  data.frame(x=x,d=d,y=rbinom(n,1,p))}
# all statistics (reject for LARGE values)
stats_vec<-function(dat){
  fit<-suppressWarnings(glm(y~x+d,data=dat,family=binomial())); ph<-pmin(pmax(as.numeric(fitted(fit)),1e-6),1-1e-6); y<-dat$y; n<-length(y)
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/10)),10)
  og<-tapply(y,grp,sum); eg<-tapply(ph,grp,sum); ng<-tapply(y,grp,length); pbar<-as.numeric(tapply(ph,grp,mean))
  Vg<-ng*pbar*(1-pbar); oe<-as.numeric(og-eg); r<-oe/sqrt(Vg); eta_g<-qlogis(pbar)
  proj<-function(Z){Z<-as.matrix(Z); Ztr<-crossprod(Z,r); as.numeric(t(Ztr)%*%solve(crossprod(Z))%*%Ztr)}
  T_EF<-sum(oe^2/Vg)-sum((1-2*pbar)*oe/Vg); HLdec<-sum(oe^2/Vg)
  S2<-tryCatch(proj(poly(pbar,2)),error=function(e)NA_real_); S3<-tryCatch(proj(poly(pbar,3)),error=function(e)NA_real_)
  S4<-tryCatch(proj(poly(pbar,4)),error=function(e)NA_real_)
  Sstk3<-tryCatch(proj(cbind(eta_g,eta_g^2*(eta_g>=0),-eta_g^2*(eta_g<0))),error=function(e)NA_real_)
  gw<-cut(ph,breaks=seq(0,1,length.out=11),include.lowest=TRUE,labels=FALSE)
  ogw<-tapply(y,gw,sum); egw<-tapply(ph,gw,sum); ngw<-tapply(y,gw,length); k<-!is.na(ogw); ogw<-ogw[k];egw<-egw[k];ngw<-ngw[k]
  HLeqw<-sum((ogw-egw)^2/(egw+1e-10)+((ngw-ogw)-(ngw-egw))^2/((ngw-egw)+1e-10))
  eta_i<-predict(fit); d2<-dat; d2$za<-0.5*eta_i^2*(eta_i>=0); d2$zb<- -0.5*eta_i^2*(eta_i<0)
  fa<-suppressWarnings(glm(y~x+d+za+zb,data=d2,family=binomial())); StukelLR<-deviance(fit)-deviance(fa)
  c(EF.omni=T_EF, EFd.poly2=S2, EFd.poly3=S3, EFd.poly4=S4, EFd.stk3=Sstk3, HL=HLdec, HLeqw=HLeqw, Stukel=StukelLR)}
alts<-c("cloglog","loglog","stukel_heavy","stukel_asym","cauchit"); ns<-c(500,1000); RN<-6000; RA<-3000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L)); clusterSetRNGStream(cl,20250911)
clusterExport(cl,c("alpha","inv_stukel","gen","stats_vec")); on.exit(stopCluster(cl))
q95<-list()
for(n in ns){clusterExport(cl,"n",envir=environment())
  M<-parSapply(cl,1:RN,function(i) stats_vec(gen("null",n))); q95[[as.character(n)]]<-apply(M,1,function(v) quantile(v,0.95,na.rm=TRUE))}
cat(sprintf("%-13s %5s | %7s %8s %8s %8s %8s %5s %6s %6s\n","scenario","n","EF.omni","EFd.pol2","EFd.pol3","EFd.pol4","EFd.stk3","HL","HLeqw","Stukel"))
nm<-c("EF.omni","EFd.poly2","EFd.poly3","EFd.poly4","EFd.stk3","HL","HLeqw","Stukel"); out<-list()
for(scn in alts) for(n in ns){clusterExport(cl,c("scn","n"),envir=environment())
  A<-parSapply(cl,1:RA,function(i) stats_vec(gen(scn,n))); q<-q95[[as.character(n)]]
  pw<-sapply(nm,function(s) mean(A[s,] > q[s], na.rm=TRUE)); out[[length(out)+1]]<-data.frame(scenario=scn,n=n,t(round(pw,4)))
  cat(sprintf("%-13s %5d | %7.3f %8.3f %8.3f %8.3f %8.3f %5.3f %6.3f %6.3f\n",scn,n,pw["EF.omni"],pw["EFd.poly2"],pw["EFd.poly3"],pw["EFd.poly4"],pw["EFd.stk3"],pw["HL"],pw["HLeqw"],pw["Stukel"]))}
write.csv(do.call(rbind,out),"Enhancement_B2_sizeadjusted.csv",row.names=FALSE); cat("\nSaved (size-adjusted power; all tests calibrated to 0.05)\n")
