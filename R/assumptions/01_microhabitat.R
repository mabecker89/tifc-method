#-----------------------------------------------------------------------------------------------------------------------

# Title: 01-microhabitat
# Description: Summarize camera microhabitat assumption test
# This is non-AIC version - simple bootstrapped estimation of effect sizes
# Author: David J. Huggard

#-----------------------------------------------------------------------------------------------------------------------

# Load data

d <- readr::read_csv("./data/lookup/microhabitat-density_assumption1.csv")

#-----------------------------------------------------------------------------------------------------------------------

# Sample sizes
d %>%
  group_by(StandType, OpenType) %>%
  tally()

# Summarize abundances by stand and opening type, with bootstrap for CI's
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

#-----------------------------------------------------------------------------------------------------------------------

# Write results for downstream use
# Use `data/processed` - user can change if they wish. This folder is in the .gitignore.

processed <- "./data/processed/"

q<-data.frame(Sp=SpName)
for (i in 1:6) {
  x<-paste(round(abund.bs.sum[,i,1],3)," (",round(abund.bs.sum[,i,2],3),"-",round(abund.bs.sum[,i,3],3),")",sep="")
  q<-cbind(q,x)
}

names(q)<-c("Species","Decid.Treed","Decid.Low","Decid.Prod","Conif.Treed","Conif.Low","Conif.Prod")

write.table(q, file= paste0(processed, "summary-abundances-stand-site-types-all-species.csv") , sep=",", row.names=FALSE)

# Summary table - abundance relative to treed sites
q<-data.frame(Sp=SpName)

for (i in 1:4) {
  x<-paste(round(relabund.bs.sum[,i,1],3)," (",round(relabund.bs.sum[,i,2],3),"-",round(relabund.bs.sum[,i,3],3),")",sep="")
  q<-cbind(q,x)
}

names(q)<-c("Species","Decid.Low","Decid.Prod","Conif.Low","Conif.Prod")

write.table(q, file = paste0(processed, "summary-abundances-stand-opening-all-species.csv"), sep = ",", row.names=FALSE)

#-----------------------------------------------------------------------------------------------------------------------
