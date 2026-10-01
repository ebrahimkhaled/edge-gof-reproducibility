# DEF vs field on COVARIATE-SPACE misspecification (omitted quadratic / interaction), calibrated.
suppressMessages(library(parallel)); alpha<-0.05; lg<-function(p) log(p/(1-p))
sqb<-function(J){solve(rbind(c(1,-1.5,2.25),c(1,3,9),c(1,-3,9)),c(lg(.05),lg(.95),lg(J)))}
sib<-function(I){e0<-lg(.1);e1<-lg(.2);e2<-lg(.2+I);b0<-(e0+e1)/2;b1<-(e1-e0)/6;b3<-(e2-(b0+3*b1))/6;c(b0,b1,3*b3,b3)}
scb<-function(K){solve(rbind(c(1,-3,-2,6),c(1,-3,0,0),c(1,3,0,0),c(1,3,2,6)),c(lg(.1),lg(.1),lg(.2),lg(.2+K)))}
gen<-function(scn,n){
  if(scn=="null_x"){x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(0.8*x))),f=y~x)}
  else if(grepl("^quad",scn)){J<-as.numeric(sub("quad","",scn))/1000;b<-sqb(J);x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*x^2))),f=y~x)}
  else if(scn=="null_xd"){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(-1+0.5*x+0.8*d))),f=y~x+d)}
  else if(grepl("^int",scn)&!grepl("^intc",scn)){I<-as.numeric(sub("int","",scn))/100;b<-sib(I);x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*d+b[4]*x*d))),f=y~x+d)}
  else if(scn=="null_xz"){x<-runif(n,-3,3);z<-rnorm(n);list(d=data.frame(x=x,z=z,y=rbinom(n,1,plogis(-1+0.5*x+0.5*z))),f=y~x+z)}
  else {K<-as.numeric(sub("intc","",scn))/100;b<-scb(K);x<-runif(n,-3,3);z<-rnorm(n);list(d=data.frame(x=x,z=z,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*z+b[4]*x*z))),f=y~x+z)}}
ef_dir_cal<-function(y,ph,fit,basis,k){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])));U<-U/sqrt(Vg)
  Om<-diag(G)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  if(basis=="poly")Z<-poly(pbar,k) else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0))};Z<-as.matrix(Z)
  Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA_real_)
  cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}
ef_omni<-function(y,ph){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);g<-pmin(ceiling(rank(ph,ties.method="first")/(n/10)),10)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),8)}
hleqw<-function(y,ph){g<-cut(ph,seq(0,1,length.out=11),include.lowest=TRUE,labels=FALSE);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);nn<-tapply(y,g,length);k<-!is.na(o);o<-o[k];e<-e[k];nn<-nn[k];st<-sum((o-e)^2/(e+1e-10)+((nn-o)-(nn-e))^2/((nn-e)+1e-10));if(length(o)>2)1-pchisq(st,length(o)-2) else NA}
stuk<-function(fit){e<-predict(fit);d<-fit$data;d$za<-0.5*e^2*(e>=0);d$zb<- -0.5*e^2*(e<0);fa<-suppressWarnings(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial()));pchisq(deviance(fit)-deviance(fa),2,lower.tail=FALSE)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_);if(is.null(v)||length(v)!=1||!is.finite(v))NA_real_ else as.numeric(v)}
one_rep<-function(scn,n){G<-gen(scn,n);dat<-G$d;fit<-suppressWarnings(glm(G$f,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  c("DEF.stk"=sc(ef_dir_cal(y,ph,fit,"stukel",3)),"DEF.poly2"=sc(ef_dir_cal(y,ph,fit,"poly",2)),"EF.omni"=sc(ef_omni(y,ph)),"HL"=sc(hoslem<-ResourceSelection::hoslem.test(y,ph,g=10)$p.value),"HLeqw"=sc(hleqw(y,ph)),"Stukel"=sc(stuk(fit)))}
suppressMessages(library(ResourceSelection))
scns<-c("null_x","quad20","quad50","null_xd","int30","int50","null_xz","intc30","intc50");ns<-c(500,1000);REPS<-3000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L));clusterSetRNGStream(cl,20250911)
clusterEvalQ(cl,suppressMessages(library(ResourceSelection)));clusterExport(cl,c("alpha","lg","sqb","sib","scb","gen","ef_dir_cal","ef_omni","hleqw","stuk","sc","one_rep"));on.exit(stopCluster(cl))
cat(sprintf("%-9s %5s | %7s %8s %7s %5s %6s %7s\n","scenario","n","DEF.stk","DEF.poly2","EF.omni","HL","HLeqw","Stukel"));out<-list()
for(scn in scns)for(n in ns){clusterExport(cl,c("scn","n"),envir=environment());m<-parSapply(cl,1:REPS,function(i)one_rep(scn,n))
  pw<-rowMeans(m<alpha,na.rm=TRUE);out[[length(out)+1]]<-data.frame(scenario=scn,n=n,t(round(pw,4)))
  cat(sprintf("%-9s %5d | %7.3f %8.3f %7.3f %5.3f %6.3f %7.3f\n",scn,n,pw["DEF.stk"],pw["DEF.poly2"],pw["EF.omni"],pw["HL"],pw["HLeqw"],pw["Stukel"]))}
write.csv(do.call(rbind,out),"Enhancement_covariate.csv",row.names=FALSE);cat("\n(null_* rows = SIZE; quad/int/intc = power. quad J in /1000, int I & intc K in /100)\n")
