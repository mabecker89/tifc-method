#-----------------------------------------------------------------------------------------------------------------------

# Title: 02_behaviour
# Description: Summarise pole/camera investigation behaviour assumption
# Author: David J. Huggard

#-----------------------------------------------------------------------------------------------------------------------

# Make simplified data frame for MS
d1 <- readr::read_csv("./data/lookup/all-species-behaviour_assumption2.csv")

#-----------------------------------------------------------------------------------------------------------------------

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
  # Invest
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
  # InvestAll
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
}  # Next sp

#-----------------------------------------------------------------------------------------------------------------------

# Write results for downstream use
# Use `data/processed` - user can change if they wish. This folder is in the .gitignore.

processed <- "./data/processed/"

q<-data.frame(Sp=SpName)
for (i in 1:8) {
  x<-paste(round(invest.sum[,i,1],2)," (",round(invest.sum[,i,2],2),"-",round(invest.sum[,i,3],2),")",sep="")
  q<-cbind(q,x)
}
names(q)<-c("Species",dimnames(invest.sum)[[2]])
write.table(q, file = paste0(processed, "time-investigating-by-habitat-species.csv"), sep=",", row.names=FALSE)
df<-1/(1-invest.sum)  # And Density factor
q<-data.frame(Sp=SpName)
for (i in 1:8) {
  x<-paste(round(df[,i,1],2)," (",round(df[,i,2],2),"-",round(df[,i,3],2),")",sep="")
  q<-cbind(q,x)
}
names(q)<-c("Species",dimnames(invest.sum)[[2]])
write.table(q, file = paste0(processed, "density-factor-for-time-investigating-habitat-species.csv"), sep=",", row.names=FALSE)

# And for InvestAll
q<-data.frame(Sp=SpName)
for (i in 1:8) {
  x<-paste(round(investall.sum[,i,1],2)," (",round(investall.sum[,i,2],2),"-",round(investall.sum[,i,3],2),")",sep="")
  q<-cbind(q,x)
}

names(q)<-c("Species",dimnames(invest.sum)[[2]])

write.table(q, file = paste0(processed, "time-investigating-and-associated-by-habitat-species.csv"), sep=",", row.names=FALSE)

df<-1/(1-investall.sum)  # And Density factor
q<-data.frame(Sp=SpName)

for (i in 1:8) {
  x<-paste(round(df[,i,1],2)," (",round(df[,i,2],2),"-",round(df[,i,3],2),")",sep="")
  q<-cbind(q,x)
}

names(q)<-c("Species",dimnames(investall.sum)[[2]])

write.table(q, file = paste0(processed, "density-factor-for-time-investigating-and-associated-habitat-species.csv"), sep=",", row.names=FALSE)

#-----------------------------------------------------------------------------------------------------------------------










