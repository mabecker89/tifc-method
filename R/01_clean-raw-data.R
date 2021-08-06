#-----------------------------------------------------------------------------------------------------------------------

# Title: Clean Raw Data
# Description: Clean raw image data from WildTrax in preparation for downstream density estimation.
#              Also extract gap class information.
# Author: Marcus Becker, David J. Huggard

# Previous scripts: None

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(fs)
library(lubridate)

# Set path to data folders; note that this will be updated once data has been stored.
root <- ""

# Native species strings
load(paste0(root, "data/lookup/wt_native_sp.RData"))

# CMU grid abbreviations
cmu <- read_csv(paste0(root, "data/lookup/cmu-grid-abb.csv")) %>%
  pull()

#-----------------------------------------------------------------------------------------------------------------------

# Import data

# Absolute paths to each data file from ABMI (Alberta Biodiversity Monitoring Institute) and CMU (Caribou Monitoring Unit)
paths_abmi <- fs::dir_ls(path = paste0(root, "data/base/raw/from_WildTrax/ABMI"), glob = "*.csv")
paths_cmu <- fs::dir_ls(path = paste0(root, "data/base/raw/from_WildTrax/CMU"), glob = "*.csv")

paths <- c(paths_abmi, paths_cmu)

# Read data directly as dataframe
df_all <- purrr::map_df(.x = paths, .f = vroom::vroom) %>% # Vroom-vroom!!
  # Remove extraneous rows
  select(project, organization, location, date_detected, field_of_view,
         common_name, age_class, sex, number_individuals)

#-----------------------------------------------------------------------------------------------------------------------

# Clean data

df_native_clean <- df_all %>%
  # Only keep native species tags for this step, takes way too long otherwise
  filter(common_name %in% native_sp) %>%
  # Amalgamate tags of same species in the same image
  group_by(project, location, date_detected, common_name) %>%
  mutate(number_individuals = sum(number_individuals),
         age_class = paste0(age_class, collapse = ", "),
         sex = paste0(sex, collapse = ", ")) %>%
  distinct(project, location, date_detected, common_name, number_individuals, .keep_all = TRUE) %>%
  ungroup()

# Join the remainder non-native mammal images back in
df_all_clean <- df_all %>%
  filter(!common_name %in% native_sp) %>%
  # Change VNA to 1 so that it can be numeric for bind_rows to follow
  mutate(number_individuals = as.numeric(ifelse(number_individuals == "VNA", 1, number_individuals))) %>%
  bind_rows(df_native_clean)

#-----------------------------------------------------------------------------------------------------------------------

# Extract gap class information

# Retrieve gap class designations from previous data (years 2013 through 2018), and include:
#   - P (Present)
#   - L (Left)
#   - U (Uncertain)
#   - N (None)
# Note: These tags were added with the previous ABMI tagging system (not WildTrax)

df_gap_201318 <-
  # Read in previous data (years 2013 through 2018)
  read_csv(paste0(root, "data/base/raw/previous/ALL_native-mammals_2019-09-01.csv"),
           col_types = cols(distance = col_character(),
                            number_during_gap = col_number(),
                            number_individuals = col_character()),
           na = "") %>%
  filter(!is.na(gap_class)) %>%
  mutate(date_detected = ymd_hms(date_time_taken)) %>%
  # Camera sites are now called 'location' in WildTrax data
  select(location = deployment, Year, date_detected, common_name, gap_class) %>%
  # Adjust location names to be same as current WildTrax versions
  # Remove 'ABMI' from core site names
  mutate(location = str_remove(location, "^ABMI-")) %>%
  # Remove 'N' gap classes - we can re-estimate below with all the data together
  filter(!gap_class == "N") %>%
  # Add project variable
  mutate(project = paste0("ABMI Ecosystem Health ", Year)) %>%
  select(-Year)

# Add 'N' gap classes for images following a 'NONE' image
# Note: a 'NONE' image is used to demarcate when a series should be truncated because an animal left the field of view.
df_gap_other <- df_all_clean %>%
  select(project, location, date_detected, common_name) %>%
  arrange(project, location, date_detected) %>%
  # Create gap class column
  mutate(common_name_next = lead(common_name),
         gap_class = ifelse(common_name != "NONE" & common_name_next == "NONE", "N", NA)) %>%
  # Include only N gap class for native mammals
  filter(gap_class == "N",
         common_name %in% native_sp) %>%
  select(-c(common_name_next))

# Combine with previous years data, write to csv
bind_rows(df_gap_201318, df_gap_other) %>%
  write_csv(paste0(root, "data/processed/probabilistic-gaps/abmi-cmu_all-years_gap-class-raw_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Save cleaned data of just native species
df_all_clean %>%
  # Remove all non-native mammal images
  filter(common_name %in% native_sp) %>%
  write_csv(paste0(root, "data/base/clean/abmi-cmu_all-years_native-sp_clean_", Sys.Date(), ".csv"))

# Save (all) cleaned data
df_all_clean %>% write_csv(paste0(root, "data/base/clean/abmi-cmu_all-years_all-data_clean_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------
