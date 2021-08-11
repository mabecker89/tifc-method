#-----------------------------------------------------------------------------------------------------------------------

# Title: 01_probabilistic-gaps
# Description: Modeling the probability of an animal leaving the camera field of view based on gap length between images
# Author: David J. Huggard, Marcus Becker

# Previous scripts required: None

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(tidyr)
library(dplyr)
library(purrr)

#-----------------------------------------------------------------------------------------------------------------------

# Load data:

# Assigned gap classes for sample of images:
#   - P (Present)
#   - L (Left)
#   - U (Uncertain)
#   - N (None)
df_gap <- read_csv("./data/lookup/abmi-cmu_all-years_gap-class-raw.csv")

# Lookup table for species gap groupings
df_gap_groups <- read_csv("./data/lookup/species-gap-groups.csv")

# Load tag data
df_all <- read_csv("./data/base/abmi-cmu_all-years_native-sp.csv")

#-----------------------------------------------------------------------------------------------------------------------

# Probabilistic Gaps

# Set up data for modeling leaving probability using observations from gap analysis and time difference between images:
df_leave_prob <- df_all %>%
  # Join gap class
  left_join(df_gap, by = c("location", "project", "date_detected", "common_name")) %>%
  # Order observations
  arrange(location, date_detected, common_name, project) %>%
  # Identify series and gaps requiring probabilistic time assignment
  # Lagged time stamp
  mutate(date_detected_previous = lag(date_detected),
         # Calculate difference in time between ordered images
         diff_time = as.numeric(date_detected - date_detected_previous),
         # Lagged species
         common_name_previous = lag(common_name),
         # Was it a different species?
         diff_sp = ifelse(common_name != common_name_previous, TRUE, FALSE),
         # Lagged deployment
         location_previous = lag(location),
         # Was is a different deployment?
         diff_location = ifelse(location != location_previous, TRUE, FALSE),
         # Flag gaps that will need checking
         gap_check = ifelse(diff_location == FALSE & diff_sp == FALSE & (diff_time <= 120 & diff_time >= 20), 1, 0),
         # Lagged gap class
         gap_class_previous = lag(gap_class)) %>%
  # Join gap group lookup table
  left_join(df_gap_groups, by = "common_name") %>%
  # Isolate observations requiring a gap check
  filter(!is.na(gap_class_previous) & gap_check == "1") %>%
  select(common_name, gap_group, gap_class_previous, diff_time) %>%
  # Determine whether the animal left the field of view or not (i.e. tagged as 'P' (present) by gap analysis)
  mutate(left = ifelse(gap_class_previous == "P", 0, 1)) %>%
  select(-gap_class_previous)

# Modeling by gap group, followed by predictions using those models:
df_leave_prob_pred <- df_leave_prob %>%
  filter(!is.na(gap_group)) %>%
  group_by(gap_group) %>%
  nest() %>%
  # Define model for each gap group, and then make predictions
  mutate(model = map(.x = data, ~ smooth.spline(x = .$diff_time, y = .$left, df = 3)),
         pred = map(.x = model, ~ predict(., x = 20:120))) %>%
  select(gap_group, pred) %>%
  unnest_wider(pred) %>%
  unnest(cols = c(x, y)) %>%
  rename(diff_time = x, pred = y) %>%
  ungroup()

#-----------------------------------------------------------------------------------------------------------------------

# Write results for downstream use
# Use `data/processed` - user can change if they wish. This folder is in the .gitignore.

processed <- "./data/processed/"

# Predictions
write_csv(df_leave_prob_pred, paste0(processed, "gap-leave-prob_predictions.csv"))
# Data used in predictions
write_csv(df_leave_prob, paste0(processed, "gap-leave-prob_raw-data.csv"))

#-----------------------------------------------------------------------------------------------------------------------
