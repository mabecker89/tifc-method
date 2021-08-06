#-----------------------------------------------------------------------------------------------------------------------

# Title: Detection distance testing using images from the Edmonton Valley Zoo
# Subtitle: Reindeer sub-selection
# Author: Marcus Becker

#-----------------------------------------------------------------------------------------------------------------------

# Set up:

# Load packages
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(purrr)

# Set path to data folders; note that this will be updated once data has been stored.
root <- ""
# Project name
proj <- "Edmonton_Valley_Zoo_Detection_Distance_Testing_WILDTRAX_REPORT"

# Previously processed data:

# 1. Probabilistic gaps
df_leave_prob_pred <- read_csv(paste0(root, "processed/probabilistic-gaps/gap-leave-prob_predictions_2020-06-25.csv"))
# 2. Time between photos
wc_tbp <- read_csv(paste0(root, "processed/time-btwn-images/Table of time between photos within series by species incl 2019 May 2020.csv")) %>%
  rename(common_name = Species) %>%
  filter(common_name == "Woodland Caribou") %>%
  select(TBP) %>%
  pull()

# Adjust times to filter out periods in motion detection camera (Z2A) when the timelapse camera (Z2B_PRESENT) wasn't operating

r_one_start <- as.POSIXct("2020-10-09 07:00:00", tz = "MST")
r_one_end <- as.POSIXct("2020-10-09 10:53:43", tz = "MST")
r_one_interval <- r_one_start %--% r_one_end

r_two_start <- as.POSIXct("2020-10-14 07:00:00", tz = "MST")
r_two_end <- as.POSIXct("2020-10-17 07:19:21", tz = "MST")
r_two_interval <- r_two_start %--% r_two_end

r_three_start <- as.POSIXct("2020-10-18 07:00:00", tz = "MST")
r_three_end <- as.POSIXct("2020-10-18 10:30:00", tz = "MST")
r_three_interval <- r_three_start %--% r_three_end

r_four_start <- as.POSIXct("2020-10-23 07:00:00", tz = "MST")
r_four_end <- as.POSIXct("2020-10-26 07:16:34", tz = "MST")
r_four_interval <- r_four_start %--% r_four_end

all_intervals <- list(r_one_interval, r_two_interval, r_three_interval, r_four_interval)

# Import and quick clean of data
df_reindeer <- read_csv(paste0(root, "supplemental/assumptions-tests/detection-distance/data/", proj, ".csv")) %>%
  # Remove out of range images
  filter(field_of_view == "WITHIN") %>%
  # Remove non-reindeer images
  filter(common_name == "Woodland Caribou") %>%
  # Keep is simple, stupid!
  select(name = location, date_detected, common_name, pole_position = image_comments) %>%
  # Make sure date_detected is a date, rename species ;)
  mutate(date_detected = ymd_hms(date_detected, tz = "MST"),
         common_name = "Reindeer") %>%
  # Subtract 35 seconds from Z2A (motion) camera to get the timing as close as possible
  mutate(date_detected = if_else(name == "Z2A", date_detected %m-% seconds(35), date_detected))

#-----------------------------------------------------------------------------------------------------------------------

# Step 1. Identify series for each of motion and timelapse cameras

# Motion Detection camera

# Series summary
df_motion_series <- df_reindeer %>%
  # Z2A first
  filter(name == "Z2A") %>%
  # Order observations chronologically
  arrange(date_detected) %>%
  # Remove observations which occurred when the timelapse camera was not working
  filter(date_detected %within% all_intervals) %>%
  # Identify series
  mutate(# Lagged time stamp
    date_detected_previous = lag(date_detected),
    # Calculate difference in time between ordered images
    diff_time = as.numeric(date_detected - date_detected_previous),
    # New series?
    diff_series = ifelse(diff_time <= 120, 0, 1),
    # Number different series
    series_num = (c(0, cumsum(diff_series[-1]))) + 1,
    series_num_previous = lag(series_num),
    gap_group = "Most ungulates") %>%
  group_by(series_num) %>%
  mutate(diff_time = ifelse(row_number() == 1, 0, diff_time)) %>%
  ungroup() %>%
  left_join(df_leave_prob_pred, by = c("gap_group", "diff_time")) %>%
  mutate(diff_time_adj = round(ifelse(!is.na(pred), (1 - pred) * diff_time, diff_time), digits = 2),
         # Let's do something special for series 11 ...
         diff_time_adj_excep = round(ifelse(!is.na(pred) & !series_num == 11, (1 - pred) * diff_time, diff_time), digits = 2))

# Series time ranges (bookends)
df_motion_bookend <- df_motion_series %>%
  group_by(name, series_num) %>%
  summarise(date_detected_first = min(date_detected),
            date_detected_last = max(date_detected))

# Pole position and series time summary
df_motion_pole <- df_motion_series %>%
  # New numeric column for whether the image was in front of the pole or not.
  mutate(pole_position = ifelse(pole_position == "F", 1, 0)) %>%
  group_by(series_num) %>%
  # Summarise proportion of images front/behind
  summarise(n_images = n(),
            n_images_front = sum(pole_position),
            prop_front = n_images_front / n_images,
            prop_behind = 1 - prop_front,
            total_time = sum(diff_time_adj),
            total_time_excep = sum(diff_time_adj_excep)) %>%
  ungroup() %>%
  left_join(df_motion_bookend, by = "series_num") %>%
  select(name, series_num, date_detected_first, date_detected_last, n_images,
         prop_front, prop_behind, total_time, total_time_excep) %>%
  # Add time between photos (used woodland caribou)
  mutate(total_time = round(total_time + wc_tbp, digits = 2),
         total_time_excep = round(total_time_excep + wc_tbp, digits = 2))

# Timelapse camera

# Series summary
df_timelapse_series <- df_reindeer %>%
  # Z2A first
  filter(name == "Z2B_PRESENT") %>%
  # Order observations chronologically
  arrange(date_detected) %>%
  # Identify series
  mutate(# Lagged time stamp
    date_detected_previous = lag(date_detected),
    # Calculate difference in time between ordered images
    diff_time = as.numeric(date_detected - date_detected_previous),
    # New series? Only 5 second difference here.
    diff_series = ifelse(diff_time <= 30, 0, 1),
    # Number different series
    series_num = (c(0, cumsum(diff_series[-1]))) + 1,
    series_num_previous = lag(series_num)) %>%
  group_by(series_num) %>%
  mutate(diff_time = ifelse(row_number() == 1, 0, diff_time)) %>%
  ungroup()

# Series time ranges (bookends)
df_timelapse_bookend <- df_timelapse_series %>%
  group_by(name, series_num) %>%
  summarise(date_detected_first = min(date_detected),
            date_detected_last = max(date_detected))

# Pole position and series time summary
df_timelapse_pole <- df_timelapse_series %>%
  mutate(pole_position = ifelse(pole_position == "F", 1, 0)) %>%
  group_by(series_num) %>%
  summarise(n_images = n(),
            n_images_front = sum(pole_position),
            prop_front = n_images_front / n_images,
            prop_behind = 1 - prop_front,
            total_time = sum(diff_time)) %>%
  ungroup() %>%
  left_join(df_timelapse_bookend, by = "series_num") %>%
  select(name, series_num, date_detected_first, date_detected_last, n_images,
         prop_front, prop_behind, total_time) %>%
  # Add 3 seconds to each series, per Dave's old way.
  mutate(total_time = total_time + 3)

# Alternative Series time summary - let's redo the series calculations -> more restrictive time differences.
df_timelapse_series_alt <- df_reindeer %>%
  # Z2A first
  filter(name == "Z2B_PRESENT") %>%
  # Order observations chronologically
  arrange(date_detected) %>%
  # Identify series
  mutate(# Lagged time stamp
    date_detected_previous = lag(date_detected),
    # Calculate difference in time between ordered images
    diff_time = as.numeric(date_detected - date_detected_previous),
    # New series? Only 5 second difference here.
    diff_series = ifelse(diff_time <= 10, 0, 1),
    # Number different series
    series_num = (c(0, cumsum(diff_series[-1]))) + 1,
    series_num_previous = lag(series_num)) %>%
  group_by(series_num) %>%
  mutate(diff_time = ifelse(row_number() == 1, 0, diff_time)) %>%
  ungroup()

# Alternative time ranges
df_timelapse_bookend_alt <- df_timelapse_series_alt %>%
  group_by(name, series_num) %>%
  summarise(date_detected_first = min(date_detected),
            date_detected_last = max(date_detected))

# Alternative pole position and series time summary
df_timelapse_pole_alt <- df_timelapse_series_alt %>%
  mutate(pole_position = ifelse(pole_position == "F", 1, 0)) %>%
  group_by(series_num) %>%
  summarise(n_images = n(),
            n_images_front = sum(pole_position),
            prop_front = n_images_front / n_images,
            prop_behind = 1 - prop_front,
            total_time = sum(diff_time)) %>%
  ungroup() %>%
  left_join(df_timelapse_bookend_alt, by = "series_num") %>%
  # Add 3 seconds to each series, per Dave's old way.
  mutate(total_time = total_time + 3) %>%
  select(name, series_num, date_detected_first, date_detected_last, n_images,
         prop_front, prop_behind, total_time)

#-----------------------------------------------------------------------------------------------------------------------

# Compare motion and timelapse results

# Motion camera intervals
motion_intervals <- df_motion_pole %>%
  # Calculate intervals - add 10 seconds buffer time onto each end.
  mutate(interval = interval(start = date_detected_first %m-% seconds(5),
                             end = date_detected_last %m+% seconds(5))) %>%
  pull(interval)

# Compare with timelapse intervals - essentially asking the question "was there a motion interval that overlapped?"
compare_to_timelapse <- df_timelapse_pole %>%
  mutate(tl_interval = interval(start = date_detected_first %m-% seconds(5), end = date_detected_last %m+% seconds(5))) %>%
  mutate(check_overlap = map(.x = tl_interval, .f = ~ int_overlaps(.x, motion_intervals)),
         overlap_motion = map_dbl(.x = check_overlap, .f = ~ if_else(any(.x), 1, 0))) %>%
  select(name:total_time, overlap_motion)

# For series that included at least one image in front:
compare_to_timelapse %>%
  filter(prop_front > 0) %>%
  summarise(check = sum(overlap_motion) / n()) %>%
  pull() # 94%

# Motion camera intervals

timelapse_intervals <- df_timelapse_pole %>%
  # Calculate intervals - add 5 seconds buffer time onto each end.
  mutate(interval = interval(start = date_detected_first %m-% seconds(5),
                             end = date_detected_last %m+% seconds(5))) %>%
  pull(interval)

# Compare with motion intervals - essentially asking the question "was there a timelapse interval that overlapped?"
compare_to_motion <- df_motion_pole %>%
  mutate(tl_interval = interval(start = date_detected_first %m-% seconds(5), end = date_detected_last %m+% seconds(5))) %>%
  mutate(check_overlap = map(.x = tl_interval, .f = ~ int_overlaps(.x, timelapse_intervals)),
         overlap_timelapse = map_dbl(.x = check_overlap, .f = ~ if_else(any(.x), 1, 0))) %>%
  select(name:total_time, overlap_timelapse)

# For series that included at least one image in front:
check <- compare_to_motion %>%
  mutate(hour_of_day = hour(date_detected_first)) %>%
  filter(prop_front > 0,
         (hour_of_day > 6 & hour_of_day < 21)) %>%
  summarise(check = sum(overlap_timelapse) / n()) %>%
  pull()

#-----------------------------------------------------------------------------------------------------------------------

# Total time - for series that include a front image
motion <- df_motion_pole %>%
  mutate(hour_of_day = hour(date_detected_first)) %>%
  # Remove series with majority images behind pole
  filter(prop_front > 0,
         hour_of_day > 6,
         hour_of_day < 21) %>%
  select(total_time) %>% # Adjusting that weird series 11.
  pull() %>%
  sum() # 9730 seconds

timelapse <- df_timelapse_pole %>%
  # Remove series with majority images behind pole
  filter(prop_front > 0) %>%
  select(total_time) %>%
  pull() %>%
  sum() # 11,640 seconds

motion / timelapse # 84%

#-----------------------------------------------------------------------------------------------------------------------
