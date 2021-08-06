#-----------------------------------------------------------------------------------------------------------------------

# Title: Summarize pole/camera investigation assumption
# Author: David J. Huggard

#-----------------------------------------------------------------------------------------------------------------------

# Make simplified data frame for MS
d1<-read.csv("all-species-behaviour_wide_2020-10-12.csv")
d<-data.frame(Deployment=d1$deployment,Sp=d1$common_name,VegHF=d1$VegForDetectionDistance,VegHF1=d1$VegHF1,tTotal=d1$total_time,
              tInvest=d1$total_time_IP+d1$total_time_IC,
              tInvestAll=d1$total_time_IP+d1$total_time_IC+d1$total_time_DP+d1$total_time_DC+d1$total_time_L,
              tOther=d1$total_time_O)
d$pInvest<-d$tInvest/d$tTotal
d$pInvestAll<-d$tInvestAll/(d$tTotal+0.0000001)  # One rounding error
d$VegHF<-as.character(d$VegHF)
d$VegHF<-ifelse(d$VegHF=="Water" | d$VegHF=="WetGrass","WetOpen",d$VegHF)

SpTable<-c("Black Bear","Canada Lynx","Coyote","Elk (wapiti)","Gray Wolf","Marten","Moose","Mule deer","Pronghorn","Red fox","Snowshoe Hare","White-tailed Deer","White-tailed Jack Rabbit","Woodland Caribou")
SpName<-c("Black bear","Lynx","Coyote","Elk","Wolf","Marten","Moose","Mule deer","Pronghorn","Red fox","Snowshoe hare","White-tailed deer","White-tailed jack rabbit","Woodland caribou")
VegHF.all<-sort(unique(d$VegHF))
invest.sum<-investall.sum<-array(NA,c(length(SpTable),length(VegHF.all)+1,3))  # Predicted proportion of time investigating or investigating+associated for each species, overall and for each veg type, {mean, lci, uci}
dimnames(invest.sum)<-dimnames(investall.sum)<-list(SpTable,c("All",VegHF.all),c("Mean","LCI","UCI"))
for (sp in 1:length(SpTable)) {
  d1<-d[d$Sp==SpTable[sp],]
  q<-table(d1$VegHF)  # Select only VegHF types with >10 series
  i<-names(q)[q>=10]
  d1<-d1[d1$VegHF %in% i == TRUE,]
  wt<-d1$tTotal/sum(d1$tTotal)*nrow(d1)  # Weights proportional to total time of series, but to sum to number of series - to give correct (time-weighted) mean proportion and (binomial) CI's
  m.invest<-m.investall<-list(NULL)
  m.invest[[1]]<-glm(pInvest~1,weights=wt,data=d1,family="binomial")
  m.invest[[2]]<-glm(pInvest~VegHF,weights=wt,data=d1,family="binomial")
  m.investall[[1]]<-glm(pInvestAll~1,weights=wt,data=d1,family="binomial")
  m.investall[[2]]<-glm(pInvestAll~VegHF,weights=wt,data=d1,family="binomial")
  VegHF.list<-sort(unique(d1$VegHF))
  # Figures - Invest
  p1<-predict(m.invest[[1]],newdata=data.frame(x=1),se.fit=TRUE)
  p2<-predict(m.invest[[2]],newdata=data.frame(VegHF=VegHF.list),se.fit=TRUE)
  x<-c(1,seq(2.5,length(VegHF.list)+1.5))
  y<-plogis(c(p1$fit,p2$fit))
  y.lci<-plogis(c(p1$fit,p2$fit)-1.65*c(p1$se.fit,p2$se.fit))
  y.uci<-plogis(c(p1$fit,p2$fit)+1.65*c(p1$se.fit,p2$se.fit))
  i<-c(1,match(VegHF.list,VegHF.all)+1)
  invest.sum[sp,i,1]<-y
  invest.sum[sp,i,2]<-y.lci
  invest.sum[sp,i,3]<-y.uci
  fname<-paste("Time investigating by VegHF ",SpTable[sp],".png",sep="")
  png(file=fname,height=500,width=500)
  opi<-par(mai=c(1.5,0.9,0.8,0.9))
  dplot(x,y,xlab="",ylab="Time Investigating (proportion)",yaxs="i",xaxt="n",ylim=c(0,1))
  for (i in 1:length(x)) lines(rep(x[i],2),c(y.lci[i],y.uci[i]))
  points(x,y,pch=18,cex=2.5)
  axis(side=1,at=x,lab=c("All",VegHF.list),tck=0.015,cex.axis=1.35,las=2)
  axis(side=4,at=1-1/c(1,1.2,1.5,2,3,5,10),lab=c(1,1.2,1.5,2,3,5,10),tck=0.015,cex.axis=1.3)
  box(bty="U",lwd=2)
  mtext(side=4,at=0.5,adj=0.5,line=2.5,"Density factor",cex=1.5)
  mtext(side=3,at=1,line=2,SpName[sp],adj=0,cex=1.5)
  graphics.off()
  # Figures - InvestAll
  p1<-predict(m.investall[[1]],newdata=data.frame(x=1),se.fit=TRUE)
  p2<-predict(m.investall[[2]],newdata=data.frame(VegHF=VegHF.list),se.fit=TRUE)
  x<-c(1,seq(2.5,length(VegHF.list)+1.5))
  y<-plogis(c(p1$fit,p2$fit))
  y.lci<-plogis(c(p1$fit,p2$fit)-1.65*c(p1$se.fit,p2$se.fit))
  y.uci<-plogis(c(p1$fit,p2$fit)+1.65*c(p1$se.fit,p2$se.fit))
  i<-c(1,match(VegHF.list,VegHF.all)+1)
  investall.sum[sp,i,1]<-y
  investall.sum[sp,i,2]<-y.lci
  investall.sum[sp,i,3]<-y.uci
  fname<-paste("Time investigating and associated by VegHF ",SpTable[sp],".png",sep="")
  png(file=fname,height=500,width=500)
  opi<-par(mai=c(1.5,0.9,0.8,0.9))
  dplot(x,y,xlab="",ylab="Time Invest+Assoc (proportion)",yaxs="i",xaxt="n",ylim=c(0,1))
  for (i in 1:length(x)) lines(rep(x[i],2),c(y.lci[i],y.uci[i]))
  points(x,y,pch=18,cex=2.5)
  axis(side=1,at=x,lab=c("All",VegHF.list),tck=0.015,cex.axis=1.35,las=2)
  axis(side=4,at=1-1/c(1,1.2,1.5,2,3,5,10),lab=c(1,1.2,1.5,2,3,5,10),tck=0.015,cex.axis=1.3)
  box(bty="U",lwd=2)
  mtext(side=4,at=0.5,adj=0.5,line=2.5,"Density factor",cex=1.5)
  mtext(side=3,at=1,line=2,SpName[sp],adj=0,cex=1.5)
  graphics.off()
}  # Next sp

# Tables for appendix 4
q<-data.frame(Sp=SpName)
for (i in 1:8) {
  x<-paste(round(invest.sum[,i,1],2)," (",round(invest.sum[,i,2],2),"-",round(invest.sum[,i,3],2),")",sep="")
  q<-cbind(q,x)
}
names(q)<-c("Species",dimnames(invest.sum)[[2]])
write.table(q,file="Time investigating by VegHF and species.csv",sep=",",row.names=FALSE)
df<-1/(1-invest.sum)  # And Density factor
q<-data.frame(Sp=SpName)
for (i in 1:8) {
  x<-paste(round(df[,i,1],2)," (",round(df[,i,2],2),"-",round(df[,i,3],2),")",sep="")
  q<-cbind(q,x)
}
names(q)<-c("Species",dimnames(invest.sum)[[2]])
write.table(q,file="Density factor for time investigating by VegHF and species.csv",sep=",",row.names=FALSE)

# And for InvestAll
q<-data.frame(Sp=SpName)
for (i in 1:8) {
  x<-paste(round(investall.sum[,i,1],2)," (",round(investall.sum[,i,2],2),"-",round(investall.sum[,i,3],2),")",sep="")
  q<-cbind(q,x)
}
names(q)<-c("Species",dimnames(invest.sum)[[2]])
write.table(q,file="Time investigating and associated behaviours by VegHF and species.csv",sep=",",row.names=FALSE)
df<-1/(1-investall.sum)  # And Density factor
q<-data.frame(Sp=SpName)
for (i in 1:8) {
  x<-paste(round(df[,i,1],2)," (",round(df[,i,2],2),"-",round(df[,i,3],2),")",sep="")
  q<-cbind(q,x)
}
names(q)<-c("Species",dimnames(investall.sum)[[2]])
write.table(q,file="Density factor for time investigating and associated behaviours by VegHF and species.csv",sep=",",row.names=FALSE)










