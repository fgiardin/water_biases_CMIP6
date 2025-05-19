# script to calculate critical soil moisture thresholds for CMIP6 models
# for every pixel using monthly data
# generalized to calculate thresholds for different models and scenarios

# Load required libraries
library(terra)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(segmented)
library(data.table)
library(lubridate)
library(R.matlab) # to read .mat files
library(ncdf4)
library(ncdf4.helpers)
library(future)
library(furrr)

# load function to extract and process data
source("R/process_moisture_thresholds.R")

# Set start and end dates based on availability of observations
start_date <- as.Date("2007-01-01")
end_date <- as.Date("2014-12-31")

# Variables needed
variables <- c("mrso", "hfls", "rsds", "rsus", "rlds", "rlus")

# Define model names and scenario
model_names <- c(
  "CESM2",
  "EC-Earth3-Veg",
  "IPSL-CM6A-LR",
  "UKESM1-0-LL",
  "CNRM-ESM2-1",
  "E3SM-1-1",
  "CMCC-ESM2",
  "CNRM-CM6-1",
  "MIROC6"
)

scenario_types <- c("land-hist", "historical")

# Use Future for parallel processing
num_cores <- detectCores() - 1
plan(multisession, workers = num_cores)

# Create a list of model-scenario combinations
model_scenario_list <- expand.grid(
  model = model_names,
  scenario = scenario_types,
  stringsAsFactors = FALSE
)

# Extract and process data in parallel using future_map2
results_list <- future_map2(
  model_scenario_list$model,
  model_scenario_list$scenario,
  function(model_name, scenario) {
    result <- process_moisture_thresholds(model_name, scenario)
    return(result)
  }
)

# Filter out NULL results (in case of missing files)
results_list <- results_list[!sapply(results_list, is.null)]

# Combine results from all models and scenarios
dt_allresults <- rbindlist(results_list, use.names = TRUE, fill = TRUE)


# Fit bilinear regression per location ------------------------------------
# Normalize mrso by model, scenario, and location
dt_final <- dt_allresults %>%
  group_by(model, scenario, lon, lat) %>%
  mutate(mrso_norm = (mrso - min(mrso, na.rm = TRUE)) / (max(mrso, na.rm = TRUE) - min(mrso, na.rm = TRUE))) %>%
  ungroup() %>%
  setDT()

# Load necessary functions for bilinear fitting
source("R/fit_bilinear_cmip6.R")
source("R/fit_bilinear_from_combination_cmip6.R") # version that works for this application

# Create combinations of model, scenario, lon, and lat
combinations_dt <- unique(dt_final[, .(model, scenario, lon, lat)])

# Split the data into a list by model and scenario
split_data <- split(dt_final, list(dt_final$model, dt_final$scenario), drop = TRUE)

# Apply the bilinear fitting function in parallel using future_map
results_list_mrso <- future_map(1:nrow(combinations_dt), function(i) {
  # select combination
  combo <- combinations_dt[i]

  # Apply the bilinear fitting function
  result <- fit_bilinear_from_combination(combo, split_data, "EF", "mrso_norm")

  return(result)
})

# Filter out list elements that are NULL or not data frames (due to errors)
filtered_list_mrso <- results_list_mrso[
  sapply(results_list_mrso, function(x) is.data.frame(x) || is.data.table(x))
]

# Combine the results
final_results_mrso <- rbindlist(filtered_list_mrso, fill = TRUE)

# Save the final results
saveRDS(final_results_mrso, "theta_crit_mrso_allmodels.rds", compress = "xz")

# Calculate the percentage of time under soil moisture limitation
df_count <- dt_final %>%
  left_join(final_results_mrso, by = c("lon", "lat", "model", "scenario")) %>%
  group_by(lon, lat, model, scenario) %>%
  mutate(
    count = ifelse(
      all(is.na(theta_crit)), NA,
      sum(mrso_norm < theta_crit, na.rm = TRUE) / n()
    )
  ) %>%
  ungroup()

# Save the count data
saveRDS(df_count, "df_count_mrso.rds", compress = "xz")
