#-----------------------------------------------------------------------------------------------------------------------

# Title: Summarize camera microhabitat assumption test
# This is non-AIC version - simple bootstrapped estimation of effect sizes
# Author: David J. Huggard

#-----------------------------------------------------------------------------------------------------------------------

# Summarize abundances by stand and opening type, with bootstrap for CI's
d<-read.csv("Microsite and abundance dataset for MS Mar 2021.csv")  # Can start here
# Table of stand type and open type
write.table(table(d$StandType,d$OpenType),file="Sample size by stand type and site type for MS.csv",sep=",",col.names=NA)
SpTable<-c("BlackBear","Coyote","Fisher","GrayWolf","Moose","SnowshoeHare","WhitetailedDeer")
SpTable1<-paste(rep(SpTable,each=3),c("","Summer","Winter"),sep="")
sp.to.omit<-c("BlackBearWinter","BlackBearSummer","FisherSummer","FisherWinter","GrayWolfSummer","GrayWolfWinter","CoyoteSummer","CoyoteWinter")
SpTable1<-SpTable1[SpTable1 %in% sp.to.omit == FALSE]  # No winter black bear so no point doing seasonally.  Others too rare in at least one season.
SpName<-c("Black bear","Coyote","Fisher","Wolf","Moose","Moose - Summer","Moose - Winter",
          "Snowshoe Hare","Snowshoe Hare - Summer","Snowshoe hare - Winter","White-tailed deer","White-tailed deer - Summer","White-tailed deer - Winter")
d$StandOpen<-paste(d$StandType,d$OpenType,sep="_")
niter<-10000  # Bootstrap iterations
bs<-array(NA,c(length(SpTable1),6,niter))  # Save each iteration's means for 6 standXopen types
for (iter in 1:niter) {
  s<-unlist(by(1:nrow(d),d$StandOpen,function(x) sample(x,replace=TRUE)))  # Bootstrap sample, by stand*open type
  if (iter==1) s<-1:nrow(d)  # Use data directly for first iteration
  d1<-d[s,]
  for (sp in 1:length(SpTable1)) {
    # Factor out lure first
    lure<-mean(d1[d1$Lured=="n",SpTable1[sp]],na.rm=TRUE)/mean(d1[d1$Lured=="y",SpTable1[sp]],na.rm=TRUE)  # NA's in some seasonal estimates when that season not sampled
    y<-ifelse(d1$Lured=="y",d1[,SpTable1[sp]]*lure,d1[,SpTable1[sp]])
    bs[sp,c(2,1,3,5,4,6),iter]<-as.numeric(by(y,d1$StandOpen,function(x) mean(x,na.rm=TRUE)))  # Sort as none/low/prod for decid/upcon
  }  # Next sp
}  # Next BS iter
# Summarize quantiles for each type and for abundance relative to treed for open sites
abund.bs.sum<-array(NA,c(length(SpTable1),6,3))  # Abundance for each species-season in each stand typeXopen type, {direct data, 5% quantile, 95% quantile}
relabund.bs.sum<-array(NA,c(length(SpTable1),4,3))  # Abundance for each species-season relative to abundance in treed sites for the 4 stand typesXopening types, {direct data, 5% quantile, 95% quantile}
dimnames(abund.bs.sum)<-list(SpTable1,c("Decid.Treed","Decid.Low","Decid.Prod","Conif.Treed","Conif.Low","Conif.Prod"),c("Direct","q5","q95"))
dimnames(relabund.bs.sum)<-list(SpTable1,c("Decid.Low","Decid.Prod","Conif.Low","Conif.Prod"),c("Direct","q5","q95"))
for (sp in 1:length(SpTable1)) {
  for (i in 1:6) {
    abund.bs.sum[sp,i,]<-c(bs[sp,i,1],quantile(bs[sp,i,],c(0.05,0.95)))
  }
  for (i in 1:4) {
    relabund.bs.sum[sp,i,]<-c(bs[sp,c(2,3,5,6)[i],1]/bs[sp,c(1,1,4,4)[i],1],quantile(bs[sp,c(2,3,5,6)[i],]/bs[sp,c(1,1,4,4)[i],],c(0.05,0.95),na.rm=TRUE))
  }
}

# Summary tables and figures
# Figures
for (sp in 1:length(SpTable1)) {
  # 1. Abundances for each stand X site type
  fname<-paste("Stand and site effects for MS ",SpTable1[sp],".png",sep="")
  png(file=fname,width=450,height=450)
  opi<-par(mai=c(1.0,0.9,0.8,0.4))
  x<-c(1,2,3,4.5,5.5,6.5)
  dplot(x,abund.bs.sum[sp,,1],xaxt="n",xlab="",yaxs="i",ylab="Abundance",ylim=c(0,max(abund.bs.sum[sp,,])),pch=18,cex=2,xlim=c(0.7,6.8))
  for (i in 1:6) lines(rep(x[i],2),c(abund.bs.sum[sp,i,2],abund.bs.sum[sp,i,3]),col="grey60",lwd=2)
  points(x,abund.bs.sum[sp,,1],pch=18,cex=2.5)
  axis(side=1,at=x,lab=c("Treed","Low","Prod","Treed","Low","Prod"),tck=0.015,cex.axis=1.2)
  mtext(side=1,line=2.0,at=c(2.5,6),c("Opening","Opening"),cex=1.2)
  mtext(side=1,line=3.2,at=c(2,5.5),c("Deciduous","Conifer"),cex=1.4)
  mtext(side=3,SpName[sp],at=1,adj=0,cex=1.5)
  graphics.off()
  # 2. Abundance relative to treed for each stand X opening type
  fname<-paste("Figures/Abundance relative to treed for MS ",SpTable1[sp],".png",sep="")
  png(file=fname,width=450,height=450)
  opi<-par(mai=c(1.0,0.9,0.8,0.4))
  x<-c(1,2,3.5,4.5)
  ymax<-log(min(c(30,max(relabund.bs.sum[sp,,]))))
  ymin<-log(max(c(0.1,min(relabund.bs.sum[sp,,]))))
  dplot(x,log(relabund.bs.sum[sp,,1]),xaxt="n",xlab="",yaxs="i",yaxt="n",ylab="Relative Abundance (Treed=1)",ylim=c(ymin,ymax),pch=18,cex=2,xlim=c(0.8,4.7))
  axis(side=2,at=log(c(0.1,0.3,0.5,1,1.5,2,3,5,10)),lab=rep("",9),tck=1,col="grey80")
  mtext(side=2,line=1,at=log(c(0.1,0.3,0.5,1,1.5,2,3,5,10)),c(0.1,0.3,0.5,1,1.5,2,3,5,10),srt=2,cex=1.3)
  axis(side=2,at=log(c(0.1,0.3,0.5,1,1.5,2,3,5,10)),lab=rep("",9),tcl=0.15,cex.axis=1.3,col="grey60")
  abline(0,0,col="grey30")
  uci<-ifelse(relabund.bs.sum[sp,,3]>100,100,relabund.bs.sum[sp,,3])  # Deal with infinity when treed abundance=0 in a BS iteration
  for (i in 1:4) lines(rep(x[i],2),log(c(relabund.bs.sum[sp,i,2],uci[i])),col="grey60",lwd=2)
  points(x,log(relabund.bs.sum[sp,,1]),pch=18,cex=2.5)
  axis(side=1,at=x,lab=c("Low","Productive","Low","Productive"),tck=0.015,cex.axis=1.3)
  box(bty="l")
  mtext(side=1,line=2.2,at=c(1.5,4),c("Deciduous","Conifer"),cex=1.4)
  mtext(side=3,SpName[sp],at=1,adj=0,cex=1.5)
  graphics.off()
}  # Next sp

# Summary table - abundance
q<-data.frame(Sp=SpName)
for (i in 1:6) {
  x<-paste(round(abund.bs.sum[,i,1],3)," (",round(abund.bs.sum[,i,2],3),"-",round(abund.bs.sum[,i,3],3),")",sep="")
  q<-cbind(q,x)
}
names(q)<-c("Species","Decid.Treed","Decid.Low","Decid.Prod","Conif.Treed","Conif.Low","Conif.Prod")
write.table(q,file="Summary table Abundances stand and site types all species.csv",sep=",",row.names=FALSE)
# Summary table - abundance relative to treed sites
q<-data.frame(Sp=SpName)
for (i in 1:4) {
  x<-paste(round(relabund.bs.sum[,i,1],3)," (",round(relabund.bs.sum[,i,2],3),"-",round(relabund.bs.sum[,i,3],3),")",sep="")
  q<-cbind(q,x)
}
names(q)<-c("Species","Decid.Low","Decid.Prod","Conif.Low","Conif.Prod")
write.table(q,file="Summary table Abundances relative to treed stand and opening types all species.csv",sep=",",row.names=FALSE)

# Example calculation
sp<-which(SpTable1=="Moose")
p.rep.Low<-0.05  # Proportion of each stand type assumed to be in low-productivity openings
p.rep.Prod<-0.05  # Proportion of each stand type assumed to be in productive openings
p.Decid.Low<-0.135  # Proportion of deciduous actually sampled in low-productivity openings
p.Decid.Prod<-0.509  # Proportion of deciduous actually sampled in low-productivity openings
p.Conif.Low<-0.419  # Proportion of coniferous actually sampled in low-productivity openings
p.Conif.Prod<-0.325  # Proportion of coniferous actually sampled in low-productivity openings
d.rep.decid<- (1-p.rep.Low-p.rep.Prod)*abund.bs.sum[sp,1,1] + p.rep.Low*abund.bs.sum[sp,2,1] + p.rep.Prod*abund.bs.sum[sp,3,1]
d.rep.conif<- (1-p.rep.Low-p.rep.Prod)*abund.bs.sum[sp,4,1] + p.rep.Low*abund.bs.sum[sp,5,1] + p.rep.Prod*abund.bs.sum[sp,6,1]
d.samp.decid<- (1-p.Decid.Low-p.Decid.Prod)*abund.bs.sum[sp,1,1] + p.Decid.Low*abund.bs.sum[sp,2,1] + p.Decid.Prod*abund.bs.sum[sp,3,1]
d.samp.conif<- (1-p.Conif.Low-p.Conif.Prod)*abund.bs.sum[sp,4,1] + p.Conif.Low*abund.bs.sum[sp,5,1] + p.Conif.Prod*abund.bs.sum[sp,6,1]

