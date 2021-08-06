#-----------------------------------------------------------------------------------------------------------------------

# Title: Calculate density for each camera deployment
# Author: Marcus Becker, David J. Huggard

# Previous scripts: 01_clean-raw-abmi-cmu, 02_probabilistic-gaps, 03_deployment-timeframes,
#                   04_effective-detection-distance, 05_time-in-field-of-view

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# Set path to data folders; note that this will be updated once data has been stored.
root <- ""

#-----------------------------------------------------------------------------------------------------------------------

# Load time in front of camera data:
df_tt_full <- read_csv(paste0(root, "data/processed/time-in-cam-fov/abmi-cmu_all-years_fov-time_long_2021-07-22.csv"))

# Load effective detection distance (EDD) modeling:
load(paste0(root,"data/processed/detection-distance/predictions/Detection distances by site species and season_2021-06-21.rdata"))

# Lookup:
df_dist_groups <- read_csv(paste0(root, "data/lookup/species-distance-groups.csv"))

# Set parameters:

# Camera field of view angle
cam_fov_angle <- 40

#-----------------------------------------------------------------------------------------------------------------------

# Step 1. Append EDD information

df_detdist <- dd %>%
  as.data.frame() %>%
  rownames_to_column(var = "location_project") %>%
  gather(key = "SpGroupSeason", value = "detdist", WTDeer.summer:`Pronghorn.winter`) %>%
  mutate(SpGroupSeason = str_replace_all(SpGroupSeason, "[//(//)// ]", "")) %>%
  # Create two new columns: Detection Distance Group and Season, sep by "."
  separate(SpGroupSeason, into = c("dist_group", "season")) %>%
  mutate(dist_group = str_replace(dist_group, "wapiti", ""),
         season = tolower(season))

df_dens_ing <- df_tt_full %>%
  left_join(df_dist_groups, by = "common_name") %>%
  left_join(df_detdist, by = c("location_project", "dist_group", "season")) %>%
  # Remove random species (mostly birds) <- something to check on though.
  filter(!is.na(detdist))

#-----------------------------------------------------------------------------------------------------------------------

# Step 2. Calculate density

df_density <- df_dens_ing %>%
  mutate(effort = total_season_days * (detdist ^ 2 * pi * (cam_fov_angle / 360)) / 100,
         # Catch per unit effort
         cpue = total_duration / effort,
         # Catch per unit effort in km2
         cpue_km2 = cpue / 60 / 60 / 24 * 10000) %>%
  select(location_project:total_season_days, density_km2 = cpue_km2) %>%
  separate(location_project, into = c("location", "project"), sep = "_", remove = TRUE)

#-----------------------------------------------------------------------------------------------------------------------

# Write results

# Whole file long:
write_csv(df_density, paste0(root, "results/density/abmi-cmu_all-years_density_long_", Sys.Date(), ".csv"))

# Whole file wide:
d <- df_density %>%
  select(location, project, season, total_season_days) %>%
  distinct() %>%
  pivot_wider(id_cols = c(location, project), names_from = season, values_from = total_season_days)

df_density %>%
  pivot_wider(id_cols = c(location, project), names_from = c(common_name, season), values_from = density_km2) %>%
  left_join(d, by = c("location", "project")) %>%
  select(location, project, summer, winter, everything()) %>%
  write_csv(paste0(root, "results/density/abmi-cmu_all-years_density_wide_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------
