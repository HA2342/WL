options(stringsAsFactors=FALSE)
library('lme4'); library('plyr'); library('ggplot2'); library('boot')

# create dataset x plus some quick charts

x=read.csv('WL.csv',skip=1)
x=rename(x,c('mouse'='id')) #id = id for each mouse
x$ts=paste(x$id,x$assay,sep='.') #ts = id for each time series
x$ia=paste(x$infection,x$assay,sep='.')
x$day=pmax(0,x$dayPI-x$dayTreat) #day = day post treatment; since antibody levels pre treatment are fairly stable, day -1 is treated as day 0
x$day2=x$dayPI-x$stopTreatDay #day2 = day post stop of treatment
x$titer=log(x$titer)

## quick charts; charts in paper were not made with R
ggplot(data=subset(x,treatment=='treated'&assay=='ELISA'))+ggtitle('ELISA, treated mice only')+geom_point()+geom_line()+aes(x=day,y=titer,col=infection,line=ts)
ggplot(data=subset(x,treatment=='untreated'&assay=='ELISA'))+ggtitle('ELISA, untreated mice only')+geom_point()+geom_line()+aes(x=dayPI,y=titer,col=infection,line=ts)
ggplot(data=subset(x,treatment=='stop'&day2>5))+ggtitle('ELISA, VSV, stop treatment group')+geom_point()+geom_line()+aes(x=day2,y=titer,line=ts)
ggplot(data=subset(x,assay=='HAI'))+ggtitle('HAI')+geom_point()+geom_line()+aes(x=day,y=titer,col=treatment,line=ts)
ggplot(data=subset(x,assay=='VSVNT'))+ggtitle('VSVNT')+geom_point()+geom_line()+aes(x=day,y=titer,col=treatment,line=ts)

# calculate half lives

'ΔAIC>4 was considered statistically significant. For scientific & statistical reasons, only overall half lives were shown in the paper.'

fboot=function(mod){ #function for non parametric bootstrap
 dat=mod@frame
 dat$id=sub('[.][^.]*$','',dat$ts)
 ids=unique(dat$id)
 boot(ids,function(ids,i){
  ids=ids[i]
  ndat=ldply(1:length(ids),function(j){
   o=subset(dat,id==ids[j]); o$id=j; o #each resample needs a unique id
  })
  ndat$ts=paste(ndat$id,ndat$ts,sep='.')
  fixef(lmer(formula(mod),ndat,REML=0))
 },strata=factor(substr(ids,1,3)),R=500) #stratification is by experiment & treatment group
}
SD=function(w){1.06857519*sd(sort(w)[6:495])} #calculates standard error using trimmed distribution of bootstrap resamples; only valid when number of resamples is 500; not used for CIs in paper

set.seed(228,'Mersenne-Twister','Inversion')

x$tv=(x$dayPI-59)/365 #tv = temporary variable
tmp<-lmer(titer~tv*ia+(tv||ts),subset(x,treatment=='untreated'),REML=0)
anova(tmp) #but trend for influenza ELISA to decline in untreated
tmp<-lmer(titer~tv+ia+(tv||ts),subset(x,treatment=='untreated'),REML=0)
bout=fboot(tmp)
365*log(0.5)/(bout$t0[[2]]+c(-1.96,0,1.96)*sd(bout$t[,2])) #HL (half life) for untreated groups

x$tv=log(x$day+13)-log(13) #offset of 13 days was chosen by looking at the data & should be considered an extra parameter
pfa<-lmer(titer~tv*ia+(tv|ts),subset(x,treatment=='treated'),REML=0) #power function model
bout=fboot(pfa)
est=(bout$t0[[7]]+bout$t0[[10]])/2-(bout$t0[[8]]+bout$t0[[9]])/3
se=sd((bout$t[,7]+bout$t[,10])/2-(bout$t[,8]+bout$t[,9])/3)
est/se #t value comparing decline of HAI & VSVNT with ELISA

ea<-lmer(titer~I(day/365)*ia+(I(day/365)|ts),subset(x,treatment=='treated'),REML=0) #exponential model
AIC(ea)-(AIC(pfa)+2) #power function model has much lower AIC

#load('seed1.RData') #loads RNG state used for calculations in paper

## models for early (first 17 days post rituximab) decline
ba<-lmer(titer~day*ia+(1|ts),subset(x,treatment=='treated'&day<18),REML=0)
bout=fboot(ba)
log(0.5)/c(bout$t0[[2]]+c(-1.96,0,1.96)*sd(bout$t[,2])) #HL for influenza ELISA
log(0.5)/c(bout$t0[[2]]+bout$t0[[7]]+c(-1.96,0,1.96)*sd(bout$t[,2]+bout$t[,7])) #HL for influenza HAI
log(0.5)/c(bout$t0[[2]]+bout$t0[[8]]+c(-1.96,0,1.96)*sd(bout$t[,2]+bout$t[,8])) #HL for LCMV ELISA
log(0.5)/c(bout$t0[[2]]+bout$t0[[9]]+c(-1.96,0,1.96)*sd(bout$t[,2]+bout$t[,9])) #HL for VSV ELISA
log(0.5)/c(bout$t0[[2]]+bout$t0[[10]]+c(-1.96,0,1.96)*sd(bout$t[,2]+bout$t[,10])) #HL for VSV NT
((bout$t0[[7]]+bout$t0[[10]])/2-(bout$t0[[8]]+bout$t0[[9]])/3)/sd((bout$t[,7]+bout$t[,10])/2-(bout$t[,8]+bout$t[,9])/3) #t value for HAI & VSV NT vs ELISA
bb<-lmer(titer~day+ia+(1|ts),subset(x,treatment=='treated'&day<18),REML=0)
bout=fboot(bb)
log(0.5)/(bout$t0[[2]]+c(-1.96,0,1.96)*sd(bout$t[,2])) #HL for overall early decline

## models for late (starting day 85 post rituximab) decline
x$tv=(x$day-85)/365
ca<-lmer(titer~tv*ia+(tv||ts),subset(x,treatment=='treated'&day>84),REML=0)
bout=fboot(ca)
365*log(0.5)/c(bout$t0[[2]]+c(-1.96,0,1.96)*sd(bout$t[,2])) #HL for influenza ELISA
365*log(0.5)/c(bout$t0[[2]]+bout$t0[[7]]+c(-1.96,0,1.96)*sd(bout$t[,2]+bout$t[,7])) #HL for influenza HAI
365*log(0.5)/c(bout$t0[[2]]+bout$t0[[8]]+c(-1.96,0,1.96)*sd(bout$t[,2]+bout$t[,8])) #HL for LCMV ELISA
365*log(0.5)/c(bout$t0[[2]]+bout$t0[[9]]+c(-1.96,0,1.96)*sd(bout$t[,2]+bout$t[,9])) #HL for VSV ELISA
365*log(0.5)/c(bout$t0[[2]]+bout$t0[[10]]+c(-1.96,0,1.96)*sd(bout$t[,2]+bout$t[,10])) #HL for VSV NT
((bout$t0[[7]]+bout$t0[[10]])/2-(bout$t0[[8]]+bout$t0[[9]])/3)/sd((bout$t[,7]+bout$t[,10])/2-(bout$t[,8]+bout$t[,9])/3) #t value for HAI & VSV NT vs ELISA
cb<-lmer(titer~tv+ia+(tv||ts),subset(x,treatment=='treated'&day>84),REML=0)
bout=fboot(cb)
365*log(0.5)/(bout$t0[[2]]+c(-1.96,0,1.96)*sd(bout$t[,2])) #HL for overall late decline

x$tv=x$day2-29
tmp<-lmer(titer~tv+(1|ts),subset(x,treatment=='stop'&day2>5),REML=0)
log(0.5)/(fixef(tmp)[[2]]+c(-1.96,0,1.96)*sqrt(vcov(tmp)[2,2])) #HL for stop treatment group, starting 29 days post stopping rituximab
