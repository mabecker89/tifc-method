#-----------------------------------------------------------------------------------------------------------------------

# Title: Calculate camera deployment operating time ranges by day
# Authors: David J. Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-data

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(tibble)

# Set path to data folders; note that this will be updated once data has been stored.
root <- ""

#-----------------------------------------------------------------------------------------------------------------------

# Load tag data
df_all <- read_csv(paste0(root, "data/base/clean/abmi-cmu_all-years_all-data_clean_2021-07-14.csv"))

# Camera start and end times (from Metadata Corrina used to provide; everything from 2018 and earlier)
df_old_cam_ranges <- read_csv(paste0(root, "data/lookup/start-end/abmi-cmu-other_2013-2018_startend_2021-06-22.csv"))

# CMU grid abbreviations + ABMI and OG
cmu <- read_csv(paste0(root, "data/lookup/cmu-grid-abb.csv")) %>%
  pull()

#-----------------------------------------------------------------------------------------------------------------------

# Obtain date ranges from newer data (ABMI 2019 and 2020, CMU 2019 and 2020)

df_new_cam_ranges <- df_all %>%
  filter(str_detect(project, "2019|2020")) %>%
  group_by(project, location) %>%
  summarise(start_date_time = min(date_detected),
            end_date_time = max(date_detected)) %>%
  ungroup() %>%
  mutate(start_date_time = as.character(start_date_time)) %>%
  # Fix a couple obvious mistakes
  mutate(start_date_time = case_when(
    location == "CHR-102" & str_detect(project, "2020") ~ "2020-02-08 00:00:00",
    location == "CHR-103" & str_detect(project, "2020") ~ "2020-02-08 00:00:00",
    location == "CHR-104" & str_detect(project, "2020") ~ "2020-02-09 00:00:00",
    location == "CHR-113" & str_detect(project, "2020") ~ "2020-02-09 00:00:00",
    location == "CHR-114" & str_detect(project, "2020") ~ "2020-02-09 00:00:00",
    location == "LLB-144" & str_detect(project, "2020") ~ "2020-01-10 00:00:00",
    location == "LRN-68" & str_detect(project, "2020") ~ "2020-02-01 00:00:00",
    location == "MRG-7" & str_detect(project, "2020") ~ "2020-02-20 00:00:00",
    TRUE ~ start_date_time
  )) %>%
  mutate(start_date_time = ymd_hms(start_date_time))

# Bind together
df_cam_ranges <- bind_rows(df_new_cam_ranges, df_old_cam_ranges)

#-----------------------------------------------------------------------------------------------------------------------

# Truncate time ranges if there is a field of view issue (i.e. END without a subsequent START):
df_cam_range_trunc <- df_all %>%
  select(project, location, date_detected, common_name, field_of_view) %>%
  filter(field_of_view == "END - Last Good Image in FOV" | field_of_view == "START - First Good Image in FOV") %>%
  arrange(project, location, date_detected) %>%
  group_by(project, location) %>%
  distinct() %>%
  mutate(starts_again = ifelse(lead(field_of_view) == "START - First Good Image in FOV", 1, 0)) %>%
  filter(field_of_view == "END - Last Good Image in FOV", is.na(starts_again)) %>%
  select(project, location, updated_latest = date_detected)

# Update deployment operating time ranges with new end_date_time, if applicable:
df_cam_ranges_upd <- df_cam_ranges %>%
  left_join(df_cam_range_trunc, by = c("project", "location")) %>%
  mutate(end_date_time = if_else(is.na(updated_latest), end_date_time, updated_latest)) %>%
  select(project, location, start_date_time, end_date_time) %>%
  mutate(project_location = paste(project, location, sep = "_"))

#-----------------------------------------------------------------------------------------------------------------------

# Create intermediate End-Start pairs dataframe, formatted as END / START in subsequent rows:

# Locations which started operation again after an intermediate pause
inter <- df_all %>%
  filter(field_of_view == "START - First Good Image in FOV" | field_of_view == "END - Last Good Image in FOV") %>%
  group_by(project, location) %>%
  tally() %>%
  filter(n > 1) %>%
  select(project, location)

df_inter_pairs <- df_all %>%
  filter(field_of_view == "START - First Good Image in FOV" | field_of_view == "END - Last Good Image in FOV") %>%
  # Reurn all rows with a match in inter
  semi_join(inter, by = c("project", "location")) %>%
  arrange(project, location, date_detected) %>%
  select(project, location, date_detected, field_of_view) %>%
  group_by(project, location) %>%
  # Issue is that cameras <2019 are formatted slightly differently wrt START / END tagging
  mutate(starts_again = ifelse(lead(field_of_view) == "START - First Good Image in FOV" & field_of_view == "END - Last Good Image in FOV", 1, NA),
         restart = ifelse(lag(starts_again) == "1" & lag(field_of_view) == "END - Last Good Image in FOV", 1, NA)) %>%
  filter(starts_again == "1" | restart == "1") %>%
  select(-c(starts_again, restart)) %>%
  ungroup() %>%
  group_split(field_of_view) %>%
  bind_cols() %>%
  mutate(time_diff = difftime(`date_detected...7`, `date_detected...3`, units = "hours")) %>%
  filter(time_diff > 12) %>%
  select(1:4, 7, 8)

ends <- df_inter_pairs %>%
  select(1, 2, date_detected = `date_detected...3`, field_of_view = `field_of_view...4`)

# Final intermediate pairs dataframe
df_inter_pairs <- df_inter_pairs %>%
  select(1, 2, date_detected = `date_detected...7`, field_of_view = `field_of_view...8`) %>%
  bind_rows(ends) %>%
  arrange(project, location, date_detected) %>%
  mutate(project_location = paste(project, location, sep = "_"))

# Temporary - add several END/START pairs from CMU 2018, where refreshes happened after the camera had konked out.
df_inter_pairs_extra <- read_csv(paste0(root, "data/lookup/start-end/cmu-2018_intermediate-end-start-pairs.csv")) %>%
  mutate(field_of_view = ifelse(field_of_view == "END",
                                "END - Last Good Image in FOV",
                                "START - First Good Image in FOV"),
         project_location = paste0(project, "_", location))

df_inter_pairs <- bind_rows(df_inter_pairs, df_inter_pairs_extra)

#-----------------------------------------------------------------------------------------------------------------------

# Summarise time-by-day for each camera deployment

start <- as.Date("2009-01-01")
end <- df_cam_ranges_upd %>% filter(!is.na(end_date_time)) %>% pull(end_date_time) %>% max()
interval <- start %--% end
days <- ceiling(as.duration(interval) / ddays(1))

# Create vector of deployments
dep <- df_cam_ranges_upd %>% filter(!is.na(end_date_time)) %>% pull(project_location)

# Date ranges, no NAs
ranges <- df_cam_ranges_upd %>% filter(!is.na(end_date_time))

# Build arrays
time.by.day <- array(0, c(length(dep), 366))
time.since.2009 <- array(0, c(length(dep), days))

# Populate arrays
for (i in 1:length(dep)) {
  df <- ranges[ranges$project_location == dep[i],]
  for (k in 1:nrow(df)) {
    # For days since 2009. floor() is used so that first & last day of operation is included.
    j1 <- floor(julian(df$start_date_time[k], origin = strptime("2009-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")))
    j2 <- floor(julian(df$end_date_time[k], origin = strptime("2009-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")))
    if (!is.na(j1) & !is.na(j2)) {
      time.since.2009[i, (j1:j2)] <- 1
    }
  }
  # To take off time(s) when camera wasn't working.
  if (dep[i] %in% df_inter_pairs$project_location) {
    df1<- df_inter_pairs[as.character(df_inter_pairs$project_location) == as.character(dep[i]),]
    for (j in seq(1, (nrow(df1) - 1), 2)) { # Assumes all extra times are formatted as END/START pairs
      # Use ceiling() so that day of failure is excluded from operating days
      j1 <- ceiling(julian(df1$date_detected[j], origin = strptime("2009-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")))
      j2 <- floor(julian(df1$date_detected[j + 1], origin = strptime("2009-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")))
      if (j2 > j1) time.since.2009[i, j1:(j2-1)] <- 0
    }
  }
}

# Add to each year (note whether it was a leap year). Currently up to 2021 (from 2009).
days.per.year <- c(365, 365, 365, 366, 365, 365, 365, 366, 365, 365, 365, 366, 365)
Jan1 <- cumsum(c(1, days.per.year))

yrnum <- julday <- NULL

# Row sums
for (i in 1:ncol(time.since.2009)) {
  yrnum[i] <- sum(i >= Jan1)
  julday[i] <- i - Jan1[yrnum[i]] + 1
}
for (i in 1:ncol(time.by.day)) {
  time.by.day[,i] <- rowSums(time.since.2009[,which(julday == i)])
}

rownames(time.by.day) <- rownames(time.since.2009) <- dep
columns <- as.character(1:366)

# Summarise time-by-day for each camera deployment
df_tbd_summary <- time.by.day %>%
  as_tibble(rownames = "project_location", .name_repair = ~ columns) %>%
  mutate(total_days = rowSums(select(., -project_location)),
         total_summer_days = rowSums(select(., -c(project_location:106, 289:366))),
         total_winter_days = rowSums(select(., -c(project_location, 107:288)))) %>%
  select(project_location, total_days, total_summer_days, total_winter_days) %>%
  separate(project_location, into = c("project", "location"), sep = "_")

#-----------------------------------------------------------------------------------------------------------------------

# Write results
write_csv(df_tbd_summary, paste0(root, "data/processed/time-by-day/abmi-cmu_all-years_tbd-summary_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------
