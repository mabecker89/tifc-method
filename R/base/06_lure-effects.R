#-----------------------------------------------------------------------------------------------------------------------

# Title: 06_lure-effects
# Description: Calculate summary of lure effects for each species
# Author: David J. Huggard, Marcus Becker

# Previous scripts: 05_calculate-density

#-----------------------------------------------------------------------------------------------------------------------

# Load data - density for lure at matched (site-level) deployments for 13 common species
d <- readr::read_csv("./data/lookup/abmi_density-for-lure.csv")

#-----------------------------------------------------------------------------------------------------------------------

# Summarise lure effects

sitelist<-sort(unique(d$site))  # Sites as unit for resampling
niter<-10000
pa.bs<-agp.bs<-ta.bs<-array(NA,c(length(SpTable),2,niter))
dimnames(pa.bs)[[1]]<-dimnames(agp.bs)[[1]]<-dimnames(ta.bs)[[1]]<-SpTable
dimnames(pa.bs)[[2]]<-dimnames(agp.bs)[[2]]<-dimnames(ta.bs)[[2]]<-c("Unlured","Lured")

# Bootstrapping
for (iter in 1:niter) {
  if (iter/100==round(iter/100)) print(paste(iter,niter,date()))
  s<-sample(1:length(sitelist),length(sitelist),replace=TRUE)
  if (iter==1) s<-1:length(sitelist)  # Use actual data in first iteration
  i<-unlist(lapply(sitelist[s],function(x) which(d$site %in% x == TRUE)))  # Records to use - each deployment in the resampled sites
  d1<-d[i,]
  for (sp in 1:length(SpTable_new)) {
    pa.bs[sp,,iter]<-as.numeric(by(sign(d1[,SpTable_new[sp]]),d1$lure,mean))
    ta.bs[sp,,iter]<-as.numeric(by(d1[,SpTable_new[sp]],d1$lure,mean))
    agp.bs[sp,,iter]<-ta.bs[sp,,iter]/pa.bs[sp,,iter]
  }  # Next sp
}  # Next iter
pa.bs.sum<-agp.bs.sum<-ta.bs.sum<-array(NA,c(length(SpTable),3))  # Lure effect ratios for each species, {direct data, 5%'ile, 95%'ile}
for (sp in 1:length(SpTable)) {
  pa.bs.sum[sp,]<-c(pa.bs[sp,2,1]/pa.bs[sp,1,1],quantile(pa.bs[sp,2,]/pa.bs[sp,1,],c(0.05,0.95)))
  agp.bs.sum[sp,]<-c(agp.bs[sp,2,1]/agp.bs[sp,1,1],quantile(agp.bs[sp,2,]/agp.bs[sp,1,],c(0.05,0.95)))
  ta.bs.sum[sp,]<-c(ta.bs[sp,2,1]/ta.bs[sp,1,1],quantile(ta.bs[sp,2,]/ta.bs[sp,1,],c(0.05,0.95)))
}

# Naming:

# Presence/Absence (PA)
dimnames(pa.bs.sum)[[1]] <- SpTable_new
dimnames(pa.bs.sum)[[2]] <- c("PA", "PA.lci", "PA.uci")
pa <- as.data.frame(pa.bs.sum)
# Density Given Presence (DGP)
dimnames(agp.bs.sum)[[1]] <- SpTable_new
dimnames(agp.bs.sum)[[2]] <- c("DGP", "DGP.lci", "DGP.uci")
agp <- as.data.frame(agp.bs.sum)
# Total Density (TD)
dimnames(ta.bs.sum)[[1]] <- SpTable_new
dimnames(ta.bs.sum)[[2]] <- c("TD", "TD.lci", "TD.uci")
ta <- as.data.frame(ta.bs.sum)

df_lure_effects <- bind_cols(pa, agp, ta) %>% tibble::rownames_to_column(var = "common_name")

#-----------------------------------------------------------------------------------------------------------------------

# Write results for downstream use
# Use `data/processed` - user can change if they wish. This folder is in the .gitignore.

processed <- "./data/processed/"

readr::write_csv(df_lure_effects, paste0(processed, "lure-effect-summary.csv"))

#-----------------------------------------------------------------------------------------------------------------------

