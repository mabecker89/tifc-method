#-----------------------------------------------------------------------------------------------------------------------

# Title: 05_calculate-density
# Description: Calculate density for each camera deployment by species and season
# Author: Marcus Becker, David J. Huggard

# Previous scripts required: 04_time-in-field-of-view

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

#-----------------------------------------------------------------------------------------------------------------------

# Load previously processed data:
# Use `data/processed` - user can change if they wish. This folder is in the .gitignore.

processed <- "./data/processed/"

# Load time in front of camera data:
df_tt_full <- read_csv(paste0(processed, "abmi-cmu_all-years_fov-time.csv"))

# Load effective detection distance (EDD) modeling:
load(paste0(processed, "Detection distances by site species and season.rdata"))

# Distance groups lookup:
df_dist_groups <- read_csv("./data/lookup/species-distance-groups.csv")

# Set parameters:

# Camera field of view angle (per Reconyx documentation)
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

# Join all ingredients ('ing') required to calculate density
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

# Write results for downstream use
# Use `data/processed` - user can change if they wish. This folder is in the .gitignore.

processed <- "./data/processed/"

# Whole file long:
write_csv(df_density, paste0(processed, "abmi-cmu_all-years_density.csv"))

#-----------------------------------------------------------------------------------------------------------------------
