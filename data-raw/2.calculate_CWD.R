# Script to calculate CWD in every grid cell of CMIP6 data and extract maximum CWD per grid cell
# generalized to work for both scenarios land-hist and historical

# Load libraries
library(terra)
library(tidyverse)
library(lubridate)
library(data.table)
library(furrr)
library(future)
library(ncdf4)
library(ncdf4.helpers)
library(R.matlab)

# Load functions
devtools::load_all(".")
source("R/process_model_cwd.R")
source("R/mct.R")

# Set model names
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

# Set scenario types
scenario_types <- c("land-hist", "historical")

# Create a list of model-scenario combinations
model_scenario_list <- expand.grid(model_name = model_names, scenario_type = scenario_types, stringsAsFactors = FALSE)

# # Focus on selection of combinations only
# model_scenario_list <- model_scenario_list %>%
#   dplyr::slice(-c(2, 3, 5, 7))

# Process each model-scenario combination sequentially!! (way faster to parallelize at the grid-cell level)
results_list <- vector("list", nrow(model_scenario_list))

for (i in seq_len(nrow(model_scenario_list))) {
  holga <- model_scenario_list$model_name[i]
  conga <- model_scenario_list$scenario_type[i]

  result <- process_model_cwd(holga, conga)
  results_list[[i]] <- result
}

# Filter out NULL results
results_list <- results_list[!sapply(results_list, is.null)]

# Combine all results into a single data.table
final_results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

# Save the combined results
saveRDS(final_results, "max_cwd_all_models_scenarios.rds", compress = "xz")
