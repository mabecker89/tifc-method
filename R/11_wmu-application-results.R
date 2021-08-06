#-----------------------------------------------------------------------------------------------------------------------

# Title: Extract data and calculate results for methods paper
# Author: Marcus Becker, David J. Huggard

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(sf)

# Set path to data folders; note that this will be updated once data has been stored.
root <- ""

#-----------------------------------------------------------------------------------------------------------------------

# Comparison of density estimates in Wildlife Management Units (WMUs) between aerial surveys and cameras

# 1. Available aerial survey results from Alberta Environment and Parks (AEP)
df_aerial <- read_csv(paste0(root, "data/supplemental/aerial-comparison/aerial-moose-estimates_2014-2020.csv")) %>%
  # Update WMU Code for later joining
  mutate(WMUNIT_COD = paste0("00", WMUNIT_COD)) %>%
  # Distance sampling only
  filter(method == "distance",
         # Remove duplicate years that were further removed chronologically from ABMI sampling
         !(WMUNIT_COD == "00356" & year == "2014"),
         !(WMUNIT_COD == "00515" & year == "2019"),
         # Remove WMus located further south in the foothills w/ minimal camera coverage
         !WMUNIT_COD == "00214",
         !WMUNIT_COD == "00314",
         !WMUNIT_COD == "00340") %>%
  select(WMUNIT_COD, aerial_avg = density_avg,
         aerial_lci_0.9 = density_lci_0.9, aerial_uci_0.9 = density_uci_0.9)

# 2. WMU locations
sf_subset_wmus <- st_read(paste0(root, "data/spatial/ab_wmus_all.shp"), quiet = TRUE) %>%
  right_join(df_aerial, by = "WMUNIT_COD") %>%
  select(WMUNIT_COD) %>%
  st_transform(4326)

# 3. Match camera locations to WMU

# ABMI core camera WMU locations (provided by Eric via spatial join)
df_abmi_cam_locations <- read_csv(paste0(root, "data/lookup/locations/OffGridCamARU_SPJoin_WMUsel_rev01.csv")) %>%
  # Remove off-grid locations
  filter(ProjectNam == "Core Terrestrial",
         deployment == "CAM" | deployment == "BOTH") %>%
  # Create project filed
  mutate(project = paste0("ABMI Ecosystem Health ", Year)) %>%
  # Relevant fields
  select(location = Site, project, WMUNIT_COD) %>%
  # Update WMU Code for later joining
  mutate(WMUNIT_COD = paste0("00", WMUNIT_COD))

# CMU locations (obtained directly from WildTrax)
paths_cmu <- fs::dir_ls(path = paste0(root, "data/base/raw/from_WildTrax/CMU"), glob = "*.csv")

df_cmu_cam_locations <- map_df(.x = paths_cmu, .f = vroom::vroom) %>% # Vroom-vroom!!
  # Remove extraneous rows
  select(project, location, latitude, longitude) %>%
  distinct() %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  # Keep CMU locations that are within one of the WMUs
  st_join(sf_subset_wmus, left = FALSE) %>%
  # Remove trail (T) locations
  filter(!str_detect(location, "T$")) %>%
  # Remove geometry
  st_set_geometry(NULL)

# Bind locations together
df_cam_locations <- bind_rows(df_cmu_cam_locations, df_abmi_cam_locations)

# 4. Join density estimates to this subset of cameras

df_density <- read_csv(paste0(root, "results/density/abmi-cmu_all-years_density_long_2021-07-22.csv")) %>%
  # Moose only
  filter(common_name == "Moose",
         # Season-days cut off set at 20 days
         total_season_days >= 20) %>%
  group_by(location, project, common_name) %>%
  # Summarise density by deployment, across the two seasons
  summarise(total_days = sum(total_season_days),
            density_km2 = mean(density_km2)) %>%
  # Total days cut off set at 50 days
  filter(total_days >= 50)

df_results <- df_cam_locations %>%
  left_join(df_density, by = c("project", "location")) %>%
  # Remove deployments without a density (b/c it was filtered out above)
  filter(!is.na(density_km2))

# 5. Lure adjustment

df_lure <- read_csv(paste0(root, "data/lookup/lure/abmi-cmu_all-years_lure_2021-06-21.csv")) %>%
  distinct()

# Moose lure effect
moose_lure_adj <- read_csv(paste0(root, "data/processed/lure/lure-effect-summary_2021-07-22.csv")) %>%
  filter(common_name == "Moose") %>%
  select(TA) %>%
  pull()

df_updated_results <- df_results %>%
  left_join(df_lure, by = c("project", "location")) %>%
  mutate(lure = if_else(is.na(lure), "No", lure)) %>%
  # Update density values that are lured
  mutate(density_km2_adj = if_else(lure == "Yes", density_km2 / moose_lure_adj, density_km2)) %>%
  select(project, location, WMUNIT_COD, common_name, density_km2_adj)

# 6. Investigating time adjustment

# VegHF information for each deployment
df_veghf <- read_csv(paste0(root, "data/lookup/veghf/abmi-cmu_all-years_veghf-soilhf-detdistveg_2021-06-21.csv")) %>%
  select(2, 1, veghf = 5)
# Moose calibration for investigating time
df_invest <- read_csv(paste0(root, "data/supplemental/assumptions-tests/investigation-behaviour-by-veghf/results/tables/moose-calibration-factors.csv")) %>%
  rename(common_name = Species, correction = mean) %>%
  filter(analysis == "invest")

df_updated_results_corr <- df_updated_results %>%
  left_join(df_veghf, by = c("project", "location")) %>%
  mutate(veghf = ifelse(veghf == "Water", "WetOpen", veghf),
         veghf = ifelse(veghf == "WetGrass", "WetOpen", veghf)) %>%
  left_join(df_invest, by = c("common_name", "veghf")) %>%
  mutate(density_km2_adj_corr = density_km2_adj * (1 - correction)) %>%
  select(1:4, 11)

# 7. Summarise density by WMU, re-join aerial estimates

# Load function(s)
source(paste0(root, "src/R/summarise-density_2021-06-23.R"))

df_sum_dens <- df_updated_results %>%
  # Use `summarise_density` function:
  summarise_density(group_id = WMUNIT_COD,
                    agg_samp_per = TRUE,
                    species_col = common_name,
                    dens_col = density_km2_adj,
                    conflevel = 0.9) %>%
  # Join aerial densities back in
  left_join(df_aerial, by = "WMUNIT_COD") %>%
  # Remove WMUs that were deemed unsuitable
  filter(!is.na(aerial_avg),
         n_deployments > 10) %>%
  # Add weight for regression
  mutate(weight = sqrt(n_deployments))

df_sum_dens_corr <- df_updated_results_corr %>%
  # Use `summarise_density` function:
  summarise_density(group_id = WMUNIT_COD,
                    agg_samp_per = TRUE,
                    species_col = common_name,
                    dens_col = density_km2_adj_corr,
                    conflevel = 0.9) %>%
  # Join aerial densities back in
  left_join(df_aerial, by = "WMUNIT_COD") %>%
  # Remove WMUs that were deemed unsuitable
  filter(!is.na(aerial_avg),
         n_deployments > 10) %>%
  # Add weight for regression
  mutate(weight = sqrt(n_deployments))

# 8. Construct linear model

m.lm <- lm(density_avg ~ aerial_avg - 1, data = df_sum_dens, weights = weight)
m.lm.corr <- lm(density_avg ~ aerial_avg - 1, data = df_sum_dens_corr, weights = weight)

#-----------------------------------------------------------------------------------------------------------------------

# Save Results

save(m.lm, m.lm.corr, df_sum_dens, df_sum_dens_corr, df_updated_results, df_updated_results_corr,
     file = paste0(root, "results/methods-paper/WMU-analysis/model_raw_data_", Sys.Date(), ".rdata"))

#-----------------------------------------------------------------------------------------------------------------------
