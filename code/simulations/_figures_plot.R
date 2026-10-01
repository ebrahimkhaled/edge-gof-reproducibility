suppressMessages({library(ggplot2); library(tidyr); library(dplyr); library(cowplot)})
cols<-c("DEF.poly3"="#D55E00","DEF.stk3"="#CC79A7","EF.omni"="#0072B2","HL"="#56B4E9","Stukel"="#009E73")
lab<-c(DEF.poly3="DEF-poly3",DEF.stk3="DEF-stk",EF.omni="EF (omnibus)",HL="HL",Stukel="Stukel")

## Fig 1: calibration-curve schematic (cloglog truth, logit fit)
set.seed(1); n<-30000; x<-runif(n,-3,3); y<-rbinom(n,1,1-exp(-exp(0.8*x)))
ph<-as.numeric(fitted(glm(y~x,family=binomial())))
g<-cut(ph,quantile(ph,0:20/20),include.lowest=TRUE,labels=FALSE)
cal<-data.frame(phat=tapply(ph,g,mean),obs=tapply(y,g,mean))
f1<-ggplot(cal,aes(phat,obs))+geom_abline(slope=1,intercept=0,linetype=2,color="grey50")+
  geom_line(color="#0072B2",linewidth=1)+geom_point(color="#0072B2",size=1.6)+
  labs(x="Mean predicted probability (per group)",y="Observed event rate",
       title="Calibration curve: true cloglog link, fitted logit")+
  annotate("text",x=0.62,y=0.30,label="smooth, systematic\ndeviation from the diagonal\n= directed signal",size=3,hjust=0,color="grey25")+
  theme_minimal(base_size=11)+theme(plot.title=element_text(face="bold",size=11))
ggsave("Fig_calibration.png",f1,width=6,height=4.2,dpi=200)

## Fig 2: power vs severity (quad-J and Stukel-alpha)
sev<-read.csv("Fig_severity.csv") %>% pivot_longer(c(DEF.poly3,DEF.stk3,EF.omni,HL,Stukel),names_to="Test",values_to="Power")
sev$Test<-factor(sev$Test,names(cols),lab[names(cols)])
f2<-ggplot(sev,aes(sev,Power,color=Test))+geom_line(linewidth=.9)+geom_point(size=1.6)+
  facet_wrap(~family,scales="free_x")+scale_color_manual(values=setNames(cols,lab[names(cols)]))+
  scale_y_continuous(limits=c(0,1))+labs(x="Misspecification severity",y="Power",title="Power increases with severity; DEF tracks Stukel and dominates the omnibus tests")+
  theme_minimal(base_size=11)+theme(legend.position="bottom",plot.title=element_text(face="bold",size=10.5),strip.text=element_text(face="bold"))
ggsave("Fig_severity.png",f2,width=8,height=4,dpi=200)

## Fig 3: omega crossover (boundary)
om<-read.csv("Fig_omega.csv") %>% pivot_longer(c(DEF.poly3,DEF.stk3,EF.omni,HL,Stukel),names_to="Test",values_to="Power")
om$Test<-factor(om$Test,names(cols),lab[names(cols)])
f3<-ggplot(om,aes(sev,Power,color=Test))+geom_line(linewidth=.9)+geom_point(size=1.6)+
  scale_color_manual(values=setNames(cols,lab[names(cols)]))+scale_y_continuous(limits=c(0,1))+
  annotate("rect",xmin=2.5,xmax=8.5,ymin=0,ymax=1,alpha=.06,fill="red")+
  annotate("text",x=5.5,y=0.5,label="high-frequency regime:\nomnibus EF/HL win,\nall directed tests fail",size=3,color="grey25")+
  labs(x=expression(paste("Oscillation frequency  ",omega,"   (logit = 0.8x + 1.5 sin(",omega,"x))")),y="Power",
       title="Boundary: directed tests (DEF, Stukel) lose to the omnibus EF/HL under rough, high-frequency misfit")+
  theme_minimal(base_size=11)+theme(legend.position="bottom",plot.title=element_text(face="bold",size=10))
ggsave("Fig_omega.png",f3,width=7.5,height=4.2,dpi=200)
cat("Saved Fig_calibration.png, Fig_severity.png, Fig_omega.png\n")
