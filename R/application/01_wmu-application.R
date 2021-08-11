#-----------------------------------------------------------------------------------------------------------------------

# Title: 01_wmu-application
# Description: Comparison of aerial survey and TIFC camera moose density estimates in WMUs of Alberta
# Author: Marcus Becker, David J. Huggard

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(sf)

# Load functions - Monte Carlo simulation function to obtain confidence intervals
source("R/functions/summarise-density.R")

# Load camera density data, with Wildlife Management Unit (WMU) as an attribute
# Note: Includes both lure-adjusted density estimate ('_adj') and lure-adjusted & habitat-corrected estimate ('_corr')
df_cam_results <- read_csv("./data/lookup/abmi-cmu_camera-moose-density-estimates-wmu.csv")

# Aerial moose estimates from AEP
df_aerial <- read_csv("./data/lookup/aep-aerial-moose-density-estimates-wmu.csv")

#-----------------------------------------------------------------------------------------------------------------------

# Summarise density by WMU, re-join aerial estimates

# Summarise density (original)
df_sum_dens <- df_cam_results %>%
  # Use `summarise_density` function:
  summarise_density(group_id = WMUNIT_COD,
                    agg_samp_per = TRUE,
                    species_col = common_name,
                    dens_col = density_km2_adj,
                    # 90% confidence interval
                    conflevel = 0.9) %>%
  # Join aerial densities back in
  left_join(df_aerial, by = "WMUNIT_COD") %>%
  # Remove WMUs that were deemed unsuitable
  filter(!is.na(aerial_avg),
         n_deployments > 10) %>%
  # Add weight for regression
  mutate(weight = sqrt(n_deployments))

# Summarise density (habitat corrected)
df_sum_dens_corr <- df_cam_results %>%
  # Use `summarise_density` function:
  summarise_density(group_id = WMUNIT_COD,
                    agg_samp_per = TRUE,
                    species_col = common_name,
                    dens_col = density_km2_adj_corr,
                    # 90% confidence interval
                    conflevel = 0.9) %>%
  # Join aerial densities back in
  left_join(df_aerial, by = "WMUNIT_COD") %>%
  # Remove WMUs that were deemed unsuitable
  filter(!is.na(aerial_avg),
         n_deployments > 10) %>%
  # Add weight for regression
  mutate(weight = sqrt(n_deployments))

# Construct linear model

m.lm <- lm(density_avg ~ aerial_avg - 1, data = df_sum_dens, weights = weight)
m.lm.corr <- lm(density_avg ~ aerial_avg - 1, data = df_sum_dens_corr, weights = weight)

#-----------------------------------------------------------------------------------------------------------------------

# Save Results

save(m.lm, m.lm.corr, df_sum_dens, df_sum_dens_corr, file = "./results/models/wmu_application_model_output.rdata")

#-----------------------------------------------------------------------------------------------------------------------
