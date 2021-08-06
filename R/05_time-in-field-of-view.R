#-----------------------------------------------------------------------------------------------------------------------

# Title: Identify series, add probabilistic gap assignment, and calculate time in front of camera
# Author: Marcus Becker, David J. Huggard

# Previous scripts: 01_clean-raw-data, 02_probabilistic-gaps, 03_deployment-timeframes,
#                   04_effective-detection-distance

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# Set path to data folders; note that this will be updated once data has been stored.
root <- ""

#-----------------------------------------------------------------------------------------------------------------------

# Previously processed data:

# 1. Probabilistic gaps
df_leave_prob_pred <- read_csv(paste0(root, "data/processed/probabilistic-gaps/gap-leave-prob_predictions_2021-07-14.csv"))
# 2. Time between photos
df_tbp <- read_csv(paste0(root, "data/processed/time-btwn-images/Table of time between photos within series by species incl 2019 May 2020.csv")) %>%
  rename(common_name = Species)
# 3. Time by day summary
df_tbd <- read_csv(paste0(root, "data/processed/time-by-day/abmi-cmu_all-years_tbd-summary_2021-07-14.csv")) %>%
  # Unite location and project into one column
  unite(col = "location_project", location, project, sep = "_", remove = TRUE)
# 4. Gap classes
df_gap <- read_csv(paste0(root, "data/processed/probabilistic-gaps/abmi-cmu_all-years_gap-class-raw_2021-07-14.csv"))

# Lookup:

# Native species
load(paste0(root, "data/lookup/wt_native_sp.RData"))
# Gap groups
df_gap_groups <- read_csv(paste0(root, "data/lookup/species-gap-groups.csv"))

# Parameters:

# Seasonal start/end dates (julian day):
summer.start.j <- 106 # April 16
summer.end.j <- 288 # October 15

# Load tag data
df_all <- read_csv(paste0(root, "data/base/clean/abmi-cmu_all-years_all-data_clean_2021-07-14.csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Quick clean of a few location names (hopefully will be fixed in WildTrax soon)
df_all <- df_all %>%
  mutate(location = case_when(
    location == "1300-NW" & str_detect(project, "2015") ~ "13-NW",
    location == "1300-SE" & str_detect(project, "2015") ~ "13-SE",
    location == "930-NW" & str_detect(project, "2016") ~ "93-NW",
    location == "930-SE" & str_detect(project, "2016") ~ "93-SE",
    location == "930-SW" & str_detect(project, "2016") ~ "93-SW",
    TRUE ~ location
  ))

#-----------------------------------------------------------------------------------------------------------------------

# Step 1. Identify Series

df_series <- df_all %>%
  # Only consider images within the field of view and of native species
  filter(field_of_view == "WITHIN",
         common_name %in% native_sp) %>%
  # Join gap class
  left_join(df_gap, by = c("location", "project", "date_detected", "common_name")) %>%
  # Order observations
  arrange(project, location, date_detected, common_name) %>%
  # Identify series and gaps requiring probabilistic time assignment
  mutate(series_num = 0,
         # Lagged time stamp
         date_detected_previous = lag(date_detected),
         # Lead time stamp
         date_detected_next = lead(date_detected),
         # Calculate difference in time between ordered images
         diff_time_previous = as.numeric(date_detected - date_detected_previous),
         diff_time_next = as.numeric(abs(date_detected - date_detected_next)),
         # Lagged species
         common_name_previous = lag(common_name),
         # Was is a different species?
         diff_sp = ifelse(common_name != common_name_previous, TRUE, FALSE),
         # Lagged deployment
         location_previous = lag(location),
         # Was is a different deployment?
         diff_location = ifelse(location != location_previous, TRUE, FALSE),
         # Flag gaps that will need checking
         gap_check = ifelse(diff_location == FALSE & diff_sp == FALSE & (diff_time_previous <= 120 & diff_time_previous >= 20), 1, 0),
         # Lagged gap class
         gap_class_previous = replace_na(lag(gap_class), ""),
         # Identify new series, based on being different deployment, species, greater than 120 seconds, and approp gaps
         diff_series = ifelse(diff_location == TRUE | diff_sp == TRUE | diff_time_previous > 120 | (gap_class_previous == "L" | gap_class_previous == "N"), 1, 0),
         # Number series
         series_num = c(0, cumsum(diff_series[-1])),
         # Flag gaps that require probabilistic time assignment
         gap_prob = replace_na(ifelse(gap_check == 1 & (gap_class_previous == "" | gap_class_previous == "U"), 1, 0), 0)) %>%
  group_by(series_num) %>%
  mutate(diff_time_previous = ifelse(row_number() == 1, 0, diff_time_previous),
         diff_time_next = ifelse(row_number() == n(), 0, diff_time_next)) %>%
  ungroup() %>%
  # Join gap group lookup table
  left_join(df_gap_groups, by = "common_name") %>%
  # Join gap leaving predictions
  left_join(df_leave_prob_pred, by = c("gap_group", "diff_time_previous" = "diff_time")) %>%
  # Adjust time difference between ordered images that require probabilistic time assignment
  mutate(pred = replace_na(pred, 1),
         diff_time_previous_adj = ifelse(gap_prob == 1, diff_time_previous * (1 - pred), diff_time_previous),
         diff_time_next_adj = ifelse(lead(gap_prob == 1), diff_time_next * (1 - lead(pred)), diff_time_next))

#-----------------------------------------------------------------------------------------------------------------------

# Step 2. Calculate time between photos (tbp), by species. This is ABMI and CMU only.

df_tbp <- df_series %>%
  mutate(series_num_previous = lag(series_num)) %>%
  # Remove first image from each series
  filter(series_num == series_num_previous) %>%
  group_by(common_name) %>%
  # Calculate average tbp and number of images from each species
  summarise(tbp = mean(diff_time_previous),
            sample_size = n())

# Write results
# write_csv(df_tbp, paste0(root, "data/processed/time-btwn-images/abmi-cmu_-all-years_tbp_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Step 3. Calculate total time in front of the camera, by series (tts = Total Time by Series)
df_tts <- df_series %>%
  left_join(df_tbp, by = "common_name") %>%
  group_by(series_num) %>%
  mutate(# Check whether the image was first or last in a series
    bookend = ifelse(row_number() == 1 | row_number() == n(), 1, 0),
    # Calculate time for each individual image
    image_time = ifelse(bookend == 1,
                        ((diff_time_previous_adj + diff_time_next_adj) / 2) + (tbp / 2),
                        (diff_time_previous_adj + diff_time_next_adj) / 2),
    # Multiply image time by the number of animals present
    image_time_ni = image_time * number_individuals) %>%
  # Group by common name as well to add it as a variable to output
  group_by(common_name, .add = TRUE) %>%
  # Calculate total time and number of images for each series
  summarise(n_images = n(),
            series_total_time = sum(image_time_ni)) %>%
  ungroup()

#-----------------------------------------------------------------------------------------------------------------------

# Step 4. Calculate total time in front of camera, by deployment, project, and species (tt = total time)

df_tt <- df_series %>%
  group_by(series_num) %>%
  arrange(date_detected, .by_group = TRUE) %>%
  filter(row_number() == 1) %>%
  left_join(df_tts, by = c("series_num", "common_name")) %>%
  select(project, location, date_detected, common_name, series_num, series_total_time) %>%
  ungroup() %>%
  mutate(julian = as.numeric(format(date_detected, "%j")),
         season = ifelse(julian >= summer.start.j & julian <= summer.end.j, "summer", "winter")) %>%
  unite(location, project, col = "location_project", sep = "_", remove = TRUE) %>%
  mutate_at(c("location_project", "common_name", "season"), factor) %>%
  group_by(location_project, common_name, season, .drop = FALSE) %>%
  summarise(total_duration = sum(series_total_time)) %>%
  ungroup() %>%
  mutate_if(is.factor, as.character) %>%
  left_join(df_tbd, by = "location_project")

# For deployments with no images of native animals (nn = no natives):

# Vector of all deployments in the ABMI projects:
dep <- df_all %>%
  select(location, project) %>%
  distinct() %>%
  unite(col = "location_project", location, project, sep = "_", remove = TRUE) %>%
  pull()

# Unique species seen
sp <- as.character(sort(unique(df_tt$common_name)))

df_tt_nn <- df_tbd %>%
  # Filter out non-ABMI deployments in df_tbd
  filter(location_project %in% dep) %>%
  # Retrieve only those that had no images of native species
  anti_join(df_tt, by = "location_project") %>%
  expand(location_project, season = c("summer", "winter"), common_name = sp) %>%
  # Re-join time-by-day information
  left_join(df_tbd, by = "location_project") %>%
  # Add total_duration column, which is zero in these cases
  mutate(total_duration = 0)

df_tt_full <- df_tt %>%
  bind_rows(df_tt_nn) %>%
  arrange(location_project, common_name, season) %>%
  mutate(total_season_days = ifelse(season == "summer", total_summer_days, total_winter_days)) %>%
  select(location_project:season, total_season_days, total_duration)

#-----------------------------------------------------------------------------------------------------------------------

# Write results

# Full results long:
write_csv(df_tt_full, paste0(root, "data/processed/time-in-cam-fov/abmi-cmu_all-years_fov-time_long_", Sys.Date(), ".csv"))

# Full results wide:
df_tt_full %>%
  pivot_wider(id_cols = location_project, names_from = c(common_name, season), values_from = total_duration) %>%
  write_csv(paste0(root, "data/processed/time-in-cam-fov/abmi-cmu_all-years_fov-time_wide_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------
