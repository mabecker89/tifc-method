#-----------------------------------------------------------------------------------------------------------------------

# Title: Lure adjustment
# Author: David J. Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-abmi-cmu, 02_probabilistic-gaps, 03_deployment-timeframes,
#                   04_effective-detection-distance, 05_time-in-field-of-view, 06_calculate-density

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# Set path to data folders; note that this will be updated once data has been stored.
root <- ""

#-----------------------------------------------------------------------------------------------------------------------

# Time by day
df_tbd <- read_csv(paste0(root, "data/processed/time-by-day/abmi-cmu_all-years_tbd-summary_2021-07-14.csv")) %>%
  select(1:3) %>%
  separate(location, into = c("site", "station", "extra"), sep = "-", remove = FALSE) %>%
  # Just ABMI core
  filter(str_detect(project, "ABMI"),
         !str_detect(location, "OG"),
         !str_detect(location, "CUDDE|HF2|HP2|^W")) %>%
  separate(location, into = c("site", "station", "extra"), sep = "-", remove = FALSE) %>%
  filter(!str_detect(site, "B"),
         is.na(extra)) %>%
  select(-extra)

# Lure lookup
df_lure <- read_csv(paste0(root, "data/lookup/lure/abmi-cmu_all-years_lure_2021-06-21.csv"))

df_lure_1 <- df_lure %>%
  filter(!is.na(lure)) %>%
  # Just ABMI core
  filter(str_detect(project, "ABMI"),
         !str_detect(location, "OG"),
         !str_detect(location, "CUDDE|HF2|HP2|^W")) %>%
  separate(location, into = c("site", "station", "extra"), sep = "-", remove = FALSE) %>%
  filter(!str_detect(site, "B"),
         is.na(extra)) %>%
  select(-extra) %>%
  group_by(project, site) %>%
  summarise(lured = sum(lure == "Yes"),
            unlured = sum(lure == "No")) %>%
  mutate(check = lured - unlured) %>%
  # Exclude sites missing one of lured or unlured
  filter(lured > 0, unlured > 0)

df_lure_include <- df_lure_1 %>%
  filter(check == "0") %>%
  select(project, site)

df_lure_include_2 <- df_lure_1 %>%
  filter(!check == "0")

df_tbd_subset <- df_tbd %>%
  inner_join(df_lure_include_2, by = c("site", "project")) %>%
  left_join(df_lure, by = c("project", "location"))

# First do those that have more lured
d_1 <- df_tbd_subset %>%
  filter(check > 0) %>%
  group_by(site) %>%
  mutate(days_check = abs(total_days - total_days[lure == "No"]))

d_2 <- d_1 %>%
  filter(lure == "Yes") %>%
  slice(which.min(days_check))

d <- d_1 %>%
  filter(lure == "No") %>%
  bind_rows(d_2) %>%
  arrange(site) %>%
  select(project, location, site, station, total_days, lure)

# Next do those that have more unlured
t_1 <- df_tbd_subset %>%
  filter(check < 0) %>%
  group_by(site) %>%
  mutate(days_check = abs(total_days - total_days[lure == "Yes"]))

t_2 <- t_1 %>%
  filter(lure == "No") %>%
  slice(which.min(days_check))

t <- t_1 %>%
  filter(lure == "Yes") %>%
  bind_rows(t_2) %>%
  arrange(site) %>%
  select(project, location, site, station, total_days, lure)

together <- bind_rows(t, d)

# Site information
north_sites <- read_csv(paste0(root, "data/lookup/climate/site-climate-summary_v2020.csv")) %>%
  filter(!NATURAL_REGIONS == "Grassland") %>%
  select(site = SITE_ID) %>%
  pull()

# Density data
df_density_1 <- read_csv(paste0(root, "results/density/abmi-all-years_density_long_2021-07-22.csv")) %>%
  separate(location, into = c("site", "station"), sep = "-", remove = FALSE) %>%
  semi_join(df_lure_include, by = c("project", "site"))

df_density_2 <- read_csv(paste0(root, "results/density/abmi-all-years_density_long_2021-07-22.csv")) %>%
  semi_join(together, by = c("project", "location"))

df_density <- bind_rows(df_density_1, df_density_2) %>%
  filter(site %in% north_sites) %>%
  group_by(location, project, common_name) %>%
  summarise(density_km2 = mean(density_km2)) %>%
  ungroup()

# Number of deployments each species occurs at:
SpTable <- df_density %>%
  filter(density_km2 > 0) %>%
  group_by(common_name) %>%
  tally() %>%
  filter(!common_name == "Deer",
         n > 100) %>%
  select(common_name) %>%
  pull()

df_density_for_lure <- df_density %>%
  filter(common_name %in% SpTable,
         !is.na(density_km2),
         !density_km2 == "Inf") %>%
  mutate(common_name = str_remove(common_name, " ")) %>%
  pivot_wider(id_cols = c(project, location), names_from = common_name, values_from = density_km2) %>%
  separate(location, into = c("site", "station"), sep = "-", remove = FALSE) %>%
  left_join(df_lure, by = c("project", "location")) %>%
  mutate(lure = factor(lure),
         site = factor(site)) %>%
  as.data.frame()

SpTable_new <- str_remove(SpTable, " ") %>%
  as.data.frame() %>%
  select(species = 1) %>%
  mutate(species = ifelse(species == "Elk", "Elk(wapiti)", species)) %>%
  pull()

# R Summarize lure effects for MS

d <- df_density_for_lure

sitelist<-sort(unique(d$site))  # Sites as unit for resampling
niter<-10000
pa.bs<-agp.bs<-ta.bs<-array(NA,c(length(SpTable),2,niter))
dimnames(pa.bs)[[1]]<-dimnames(agp.bs)[[1]]<-dimnames(ta.bs)[[1]]<-SpTable
dimnames(pa.bs)[[2]]<-dimnames(agp.bs)[[2]]<-dimnames(ta.bs)[[2]]<-c("Unlured","Lured")

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

dimnames(pa.bs.sum)[[1]] <- SpTable_new
dimnames(pa.bs.sum)[[2]] <- c("PA", "PA.lci", "PA.uci")
pa <- as.data.frame(pa.bs.sum)

dimnames(agp.bs.sum)[[1]] <- SpTable_new
dimnames(agp.bs.sum)[[2]] <- c("AGP", "AGP.lci", "AGP.uci")
agp <- as.data.frame(agp.bs.sum)

dimnames(ta.bs.sum)[[1]] <- SpTable_new
dimnames(ta.bs.sum)[[2]] <- c("TA", "TA.lci", "TA.uci")
ta <- as.data.frame(ta.bs.sum)

check <- bind_cols(pa, agp, ta) %>% tibble::rownames_to_column(var = "common_name")

# Write results
readr::write_csv(check, paste0(root, "data/processed/lure/lure-effect-summary_", Sys.Date(), ".csv"))

