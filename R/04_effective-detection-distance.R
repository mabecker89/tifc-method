#-----------------------------------------------------------------------------------------------------------------------

# Title: Effective detection distance (EDD) modeling
# Author: David J. Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-data

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(mgcv)

# Set path to data folders; note that this will be updated once data has been stored.
root <- ""

# Lure lookup
df_lure <- read_csv(paste0(root, "data/lookup/lure/abmi-cmu_all-years_lure_2021-06-21.csv"))

# Veg/HF lookup
df_veghf <- read_csv(paste0(root,"data/lookup/veghf/abmi-cmu_all-years_veghf-soilhf-detdistveg_2021-06-21.csv"))

# EDD modeling species groups
df_edd_groups <- read_csv(paste0(root, "data/lookup/species-distance-groups.csv"))

# Veg/HF groupings to try in modeling
df_veghf_groups <- read_csv(paste0(root,"data/lookup/veg-pole-distance.csv"))

# Season date cutoffs (julian day)
summer.start.j <- 106
summer.end.j <- 288

#-----------------------------------------------------------------------------------------------------------------------

# Step 1: Retrieve pole position information from previous data (years 2013 through 2018)

# Note: these tags were added with the previous ABMI tagging system (not WildTrax)

# Extra pole tagging done in September 2020
df_pole_extra <- read_csv(paste0(root, "data/processed/detection-distance/clean-data/extra-pole-position_2021-06-21.csv"))

# Pole position data
df_pole <- read_csv(
  paste0(root,"data/base/raw/previous/ALL_native-mammals_2019-09-01.csv"),
  col_types = cols(distance = col_character(),
                   number_during_gap = col_number(),
                   number_individuals = col_character()),
  na = "") %>%
  # Only observations with pole information
  filter(!is.na(distance)) %>%
  mutate(project = paste0("ABMI Ecosystem Health ", Year),
         # Remove 'ABMI-' prefix
         location = str_remove(deployment, "^ABMI-")) %>%
  # Standardize variable names
  select(location, project, date_detected = date_time_taken, common_name, distance) %>%
  mutate(date_detected = ymd_hms(date_detected)) %>%
  # Join extra tagging
  bind_rows(df_pole_extra) %>%
  # Make columns of number of individuals at each pole position
  mutate(at_pole = str_count(distance, "A"),
         behind_pole = str_count(distance, "B"),
         front_pole = str_count(distance, "F"),
         ic_pole = str_count(distance, "IC"),
         ip_pole = str_count(distance, "IP"),
         na_pole = str_count(distance, "NA")) %>%
  # Join lure information
  left_join(df_lure, by = c("location", "project")) %>%
  # Filter out lured deployments; only unlured are used in the EDD modeling
  filter(lure == "No") %>%
  # Create variable for julian date
  mutate(julian = as.numeric(format(ymd_hms(date_detected), "%j")))

# Retrieve Veg and HF information for each deployment
df_pole_veghf <- df_veghf %>%
  select(location, project, VegForDetectionDistance) %>%
  # Get rid of WetShrub - not used for modeling.
  mutate(VegHF = ifelse(VegForDetectionDistance == "WetShrub",
                        "Shrub", VegForDetectionDistance)) %>%
  select(-VegForDetectionDistance) %>%
  # Join back to pole information
  right_join(df_pole, by = c("location", "project")) %>%
  # Join in EDD species grouping lookup
  left_join(df_edd_groups, by = "common_name") %>%
  # Join alternative VegHF categories to try in the modeling
  left_join(df_veghf_groups, by = "VegHF") %>%
  filter(!is.na(VegHF)) %>%
  # Sum total individuals in each of at_pole through ip_pole
  mutate(n = rowSums(select(., at_pole:ip_pole))) %>%
  # Only use records where number_individuals = behind_pole or front_pole
  filter(n == behind_pole | n == front_pole) %>%
  # Create new season variable based on julian day, and calculate the proportion of individuals behind pole
  mutate(season = as.factor(ifelse(julian  >= summer.start.j & julian <= summer.end.j, "summer", "winter")),
         prop_behind = behind_pole / n) %>%
  select(location, project, common_name, dist_group, n, prop_behind, VegHF, VegHF1:VegHF5, season)

# Save data
write_csv(df_pole_veghf, paste0(root, "data/processed/detection-distance/clean-data/pole_position_veghf_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Step 2. EDD modeling.

# Species groups
SpTable <- unique(df_pole_veghf$dist_group)
SpTable <- "Bear"

# Tibbles wreak havoc; back to data.frame.
df_veghf_groups <- as.data.frame(df_veghf_groups)
df_pole_veghf <- as.data.frame(df_pole_veghf)

# Loop through each species group, try each model variation, make predictions, output figures
for (sp in 1:length(SpTable)) {

  print(paste(sp, length(SpTable), SpTable[sp], date()))

  d.sp <- df_pole_veghf[df_pole_veghf$dist_group == SpTable[sp], ]

  m <- list(NULL)
  m[[1]]<-try(gam(prop_behind~1,weights=d.sp$n,data=d.sp,family="binomial"))
  if (SpTable[sp]!="Pronghorn" & SpTable[sp]!="Bighorn sheep") {
    m[[2]]<-try(gam(prop_behind~as.factor(VegHF),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[3]]<-try(gam(prop_behind~as.factor(VegHF1),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[4]]<-try(gam(prop_behind~as.factor(VegHF2),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[5]]<-try(gam(prop_behind~as.factor(VegHF3),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[6]]<-try(gam(prop_behind~as.factor(VegHF4),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[7]]<-try(gam(prop_behind~as.factor(VegHF5),weights=d.sp$n,data=d.sp,family="binomial"))

    m[[8]]<-try(gam(prop_behind~as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))

    m[[9]]<-try(gam(prop_behind~as.factor(VegHF)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[10]]<-try(gam(prop_behind~as.factor(VegHF1)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[11]]<-try(gam(prop_behind~as.factor(VegHF2)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[12]]<-try(gam(prop_behind~as.factor(VegHF3)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[13]]<-try(gam(prop_behind~as.factor(VegHF4)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[14]]<-try(gam(prop_behind~as.factor(VegHF5)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
  }

  nModels<-length(m)

  # Additional modification to not count models where the minimum number of pole positions in a veg-HF type is <20
  min.n<-c(min(by(d.sp$n,d.sp$VegHF,sum)),
           min(by(d.sp$n,d.sp$VegHF1,sum)),
           min(by(d.sp$n,d.sp$VegHF2,sum)),
           min(by(d.sp$n,d.sp$VegHF3,sum)),
           min(by(d.sp$n,d.sp$VegHF4,sum)),
           min(by(d.sp$n,d.sp$VegHF5,sum)))
  min.n<-ifelse(is.na(min.n),0,min.n)

  # BIC calculation
  bic.ta<-rep(999999999,(nModels))
  for (i in 1:(nModels)) {
    if (!is.null(m[[i]]) & class(m[[i]])[1]!="try-error") {  # last part is to not used non-converged models, unless none converged
      bic.ta[i]<-BIC(m[[i]])
    }
  }

  if (SpTable[sp]!="Pronghorn" & SpTable[sp]!="Bighorn sheep") {
    bic.ta[c(2:7)]<-bic.ta[c(2:7)]+(min.n<20)*999999999  # Inflate BIC if <20 records in any one vegHF type
    bic.ta[c(9:14)]<-bic.ta[c(9:14)]+(min.n<20)*999999999  # Inflate BIC if <20 records in any one vegHF type
  }

  if (SpTable[sp]=="Bear" | SpTable[sp]=="Bighorn sheep" | SpTable[sp]=="Elk") bic.ta[8:14]<-999999999  # Too few winter data to fit properly
  bic.delta<-bic.ta-min(bic.ta)
  bic.exp<-exp(-1/2*bic.delta)
  bic.wt<-bic.exp/sum(bic.exp)
  best.model<-which.max(bic.wt)

  # Figures predicting for each veg type and season
  hab.list <- levels(as.factor(d.sp$VegHF))
  v.l1<-df_veghf_groups[df_veghf_groups$VegHF %in% hab.list,]
  if (sum(d.sp$VegHF=="Water")==0) v.l1[7,1]<-"WetGrass"  # For species that have no data in water
  p.veg<-array(0,c(length(hab.list),nModels))  # veg+HF types, for each model
  p.season<-array(0,c(2,nModels))  # Two seasons
  p.veg[,1]<-p.season[,1]<-mean(predict(m[[1]]))  # Null model

  if (SpTable[sp]!="Pronghorn" & SpTable[sp]!="Bighorn sheep") {
    for (i in 2:nModels) {
      if(class(m[[i]])[1]!="try-error" & bic.wt[i]>0.001 ) {
        p.veg[,i]<-predict(m[[i]],newdata=data.frame(v.l1, season="summer"))
        p.season[,i]<-predict(m[[i]],newdata=data.frame(v.l1[1,], season=c("summer","winter")))
      }
    }
    p1.veg<-colSums(bic.wt*t(p.veg))
    p1.season<-colSums(bic.wt*t(p.season))
  } else {
    p1.veg<-t(p.veg)
    p1.season<-p.season
  }
  dist.veg<-5/sqrt(1-plogis(p1.veg))  # Time-between-images effect would be added in here, but no meaningful differences found
  dist.season<-5/sqrt(1-plogis(p1.season))  # Time-between-images effect would be added in here, but no meaningful differences found
  dist.veg<-ifelse(dist.veg>20,20,dist.veg)  # Can be inf when all locations are behind pole
  dist.season<-ifelse(dist.season>20,20,dist.season)  # Can be inf when all locations are behind pole
  ymax<-ifelse(max(c(dist.veg,dist.season))>10,20,10)
  hab.n<-by(d.sp$n,d.sp$VegHF,sum)
  fname<-paste0(root, "data/processed/detection-distance/figures/Detection distance ",SpTable[sp],".jpg",sep="")
  jpeg(width=500,height=900,file=fname)
  par(mfrow=c(2,1),mai=c(1.8,0.8,0.3,0.3))
  # Figures - veg types
  x1<-barplot(dist.veg,ylim=c(0,ymax),xlab="",ylab="Effective distance (m)",cex.axis=1.3,cex.lab=1.4)
  mtext(side=1,at=x1,hab.list,las=2,line=1,cex=1.3)
  box(bty="l")
  text(x1[1],ymax*0.96,SpTable[sp],cex=1.5,adj=0)
  text(x1,rep(ymax*0.89,length(x1)),hab.n[match(v.l1$VegHF,names(hab.n))],cex=1.1)
  # Figures - season
  par(mai=c(1,1.8,0.1,1.3))
  x1<-barplot(dist.season,ylim=c(0,ymax),xlab="",ylab="Effective distance (m)",cex.axis=1.3,cex.lab=1.4)
  mtext(side=1,at=x1,c("Summer","Winter"),las=1,line=1,cex=1.3)
  box(bty="l")
  graphics.off()

  # Save models and bic wt
  fname<-paste0(root, "data/processed/detection-distance/group-models/Detection distance models ",SpTable[sp],".rdata",sep="")
  save(file=fname,m,bic.wt,hab.list)

}  # Next SpGroup

#-----------------------------------------------------------------------------------------------------------------------

# Step 3: Make predictions for each deployment, season, and species combination

df_veghf_detdist <- df_veghf %>%
  select(location, project, VegHF = VegForDetectionDistance) %>%
  mutate(location_project = paste0(location, "_", project)) %>%
  mutate(VegHF = ifelse(VegHF == "WetShrub", "Shrub", VegHF)) %>%
  left_join(df_veghf_groups, by = "VegHF")

pred.season <- c("summer", "winter")

site.list <- unique(df_veghf_detdist$location_project)

fname.models <- paste0(root, "data/processed/detection-distance/group-models/Detection distance models ")

dd <- dd.lci <- dd.uci <- array(NA, c(length(site.list), length(SpTable), length(pred.season)))

for (sp in c(1:(length(SpTable) - 2))) {  # Bighorn sheep and Pronghorn treated separately.

  fname<-paste(fname.models, SpTable[sp], ".rdata", sep="")
  load(fname)

  # Model m, bic.wt
  hab.list<-as.character(m[[2]]$xlevels[[1]])  # The shared subset of habitat types in both VegHF models
  s1 <- df_veghf_detdist
  # Missing data for some veg types for some species, so need to collapse vegHF types in s
  if (("Water" %in% hab.list) == FALSE) s1$VegHF <- ifelse(s1$VegHF == "Water", "WetGrass", s1$VegHF)
  if (("WetTreed" %in% hab.list) == FALSE) s1$VegHF <- ifelse(s1$VegHF == "WetTreed", "WetGrass", s1$VegHF)
  if (("WetGrass" %in% hab.list) == FALSE) s1$VegHF <- ifelse(s1$VegHF == "WetGrass", "Grass", s1$VegHF)
  if (("Grass" %in% hab.list) == FALSE) s1$VegHF <- ifelse(s1$VegHF == "Grass", "WetGrass", s1$VegHF)
  if (("Shrub" %in% hab.list)==FALSE) {
    s1$VegHF <- ifelse(s1$VegHF == "Shrub", "Grass", s1$VegHF)
    s1$VegHF2 <- ifelse(s1$VegHF2 == "Shrub", "GrassWater", s1$VegHF2)
    s1$VegHF4 <- ifelse(s1$VegHF4 == "Shrub", "GrassWaterHF", s1$VegHF4)
  }
  if ("wet" %in% m[[10]]$xlevels == FALSE) {
    s1$VegHF1 <- ifelse(s1$VegHF1 == "Wet", "GrassShrub", s1$VegHF1)
    s1$VegHF3 <- ifelse(s1$VegHF3 == "Wet", "GrassShrubHF", s1$VegHF3)
  }
  for (j in 1:2) {  # Predictions for the two seasons
    p <- p.se<-array(0,c(length(m), length(site.list)))
    p1 <- predict(m[[1]], se.fit=TRUE)
    p[1, ] <- rep(p1$fit[1], length(site.list))
    p.se[1, ] <- rep(p1$se.fit[1], length(site.list))
    for (k in 2:length(m)) {

      p1<-predict(m[[k]], newdata = data.frame(
        s1[,c("VegHF","VegHF1","VegHF2","VegHF3","VegHF4","VegHF5")],
        season=pred.season[j]),
        se.fit=TRUE)

      # The SE produces CI's on the detection distance, but not currently used (would be used if we did full CI's on estimates)
      p[k, ] <- p[k, ] + p1$fit
      p.se[k, ] <- p.se[k, ] + p1$se.fit
    }

    p.all <- colSums(p * bic.wt)  # BIC weighted average
    p.se.all <- colSums(bic.wt * sqrt(p.se ^ 2 + (p - rep(p.all, each = length(m))) ^ 2))
    # Convert probability of being behind 5m pole into effective detection distance
    dd[, sp, j] <- 5 / sqrt(1 - plogis(p.all))
    dd.lci[, sp, j] <- 5 / sqrt(1 - plogis(p.all - 2 * p.se.all))
    dd.uci[, sp, j] <- 5 / sqrt(1 - plogis(p.all + 2 * p.se.all))
  }  # Next time period j
}

# Separate predictions for pronghorn and bighorn, because only have null model
for (sp in (length(SpTable) - 1):length(SpTable)) {

  fname <- paste(fname.models, SpTable[sp], ".rdata", sep = "")
  load(fname)  # Model m, bic.wt, hab.list
  s1 <- df_veghf_detdist
  for (j in 1:2) {
    p1 <- predict(m[[1]], se.fit=TRUE)  # Only null model
    dd[, sp, j] <- rep(5 / sqrt(1 - plogis(p1$fit[1])), length(site.list))  # Only null model
    dd.lci[, sp, j] < -rep(5 / sqrt(1 - plogis(p1$fit[1] - 2 * p1$se.fit[1])), length(site.list))  # Only null model
    dd.uci[, sp, j] <- rep(5 / sqrt(1 - plogis(p1$fit[1] + 2 * p1$se.fit[1])), length(site.list))  # Only null model
  }  # Next time period j
}  # Next of these two species

dimnames(dd)[[1]]<-dimnames(dd.lci)[[1]]<-dimnames(dd.uci)[[1]]<-site.list
dimnames(dd)[[2]]<-dimnames(dd.lci)[[2]]<-dimnames(dd.uci)[[2]]<-SpTable
dimnames(dd)[[3]]<-dimnames(dd.lci)[[3]]<-dimnames(dd.uci)[[3]]<-pred.season

#-----------------------------------------------------------------------------------------------------------------------

# Write results

save(dd, dd.lci, dd.uci, SpTable, site.list,
     file = paste0(root,
                   "data/processed/detection-distance/predictions/Detection distances by site species and season_",
                   Sys.Date(),
                   ".rdata"))

#-----------------------------------------------------------------------------------------------------------------------
