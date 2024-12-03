# Generalized script to calculate CWD in every grid cell of CMIP6 data
# and extract maximum CWD per grid cell

# Load required libraries
library(terra)
library(tidyverse)
library(lubridate)
library(data.table)
library(parallel)
library(ncdf4)
library(ncdf4.helpers)
library(R.matlab) # if needed for mct function

# Load custom functions (assuming mct is defined in your package)
devtools::load_all(".")

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

# Function to process data for a given model and scenario
process_model_cwd <- function(model_name, scenario_type) {
  print(paste("Processing model:", model_name, "Scenario:", scenario_type))

  # Define file paths based on model and scenario
  # Adjust the realization based on the model
  if (model_name %in% c("CNRM-ESM2-1", "UKESM1-0-LL", "CNRM-CM6-1")) {
    realization <- "r1i1p1f2"
  } else if (model_name == "E3SM-1-1") {
    realization <- "r1i1p11f1"
  } else {
    realization <- "r1i1p1f1"
  }

  # Construct file paths using pattern matching
  et_pattern <- paste0("evspsbl_mon_", model_name, "_", scenario_type, "_", realization, "_g025.nc")
  pr_pattern <- paste0("pr_mon_", model_name, "_", scenario_type, "_", realization, "_g025.nc")

  et_files <- list.files("data-raw/cmip6-ng/evspsbl/mon/g025/", pattern = et_pattern, full.names = TRUE)
  pr_files <- list.files("data-raw/cmip6-ng/pr/mon/", pattern = pr_pattern, full.names = TRUE)

  if (length(et_files) == 0 || length(pr_files) == 0) {
    warning(paste("Files not found for model:", model_name, "Scenario:", scenario_type))
    return(NULL)
  }

  # Load the data
  ET <- rast(et_files[1])
  P <- rast(pr_files[1])

  # Special case for CMCC-ESM2 using precipitation from CESM2
  if (model_name == "CMCC-ESM2") {
    P <- rast("data-raw/cmip6-ng/pr/mon/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
  }

  # Calculate water balance
  water_balance_raw <- P - ET

  # Select dates from 1935-01-01 onwards (80 years of data)
  dates <- terra::time(water_balance_raw)
  selection <- which(dates >= as.Date("1935-01-01"))
  water_balance <- subset(water_balance_raw, subset = selection)

  # Rotate to true coordinates (-180;180 instead of 0;360)
  water_balance <- terra::rotate(water_balance)

  # Focus on vegetated land
  land_cover_raw <- rast("data-raw/landcover/landcover_MCD12C1.nc")
  land_cover <- flip(land_cover_raw[[1]], direction = "vertical")
  vegetated_land <- ifel(
    land_cover == 0,
    NA,
    ifel(
      land_cover == 13,  # Urban and Built-up Lands
      NA,
      ifel(
        land_cover > 14,  # Permanent Snow and Ice, Barren
        NA,
        land_cover
      )
    )
  )
  vegetated_land <- resample(vegetated_land, water_balance)
  water_balance <- mask(water_balance, vegetated_land)

  # Transform to dataframe
  df_wb <- terra::as.data.frame(water_balance, xy = TRUE)
  dates <- terra::time(water_balance)

  names(df_wb) <- c("lon", "lat", as.character(dates))

  df_wb_long <- df_wb %>%
    data.table() %>%
    melt(
      measure.vars = as.character(dates),
      variable.name = "date",
      value.name = "water_balance"
    ) %>%
    mutate(date = as.Date(date))

  # Convert water_balance to mm/day
  df_wb_long[, water_balance := water_balance * 86400]  # Convert from kg/m^2/s to mm/day

  # Calculate CWD using mct function
  print(paste("Calculating CWD for model:", model_name, "Scenario:", scenario_type))

  # Group data by lon and lat
  df_grouped_list <- split(df_wb_long, by = c("lon", "lat"))

  # Define function to apply mct to each group
  apply_mct <- function(group_data) {
    result <- tryCatch({
      mct(group_data, "water_balance", "date")
    }, error = function(e) {
      NULL  # Return NULL in case of error
    })
    return(result)
  }

  # Parallel computation using mclapply
  num_cores <- detectCores() - 1
  df_cwd_list <- mclapply(df_grouped_list, apply_mct, mc.cores = num_cores)

  # Remove NULL results
  df_cwd_list <- df_cwd_list[!sapply(df_cwd_list, is.null)]

  if (length(df_cwd_list) == 0) {
    warning(paste("No CWD data calculated for model:", model_name, "Scenario:", scenario_type))
    return(NULL)
  }

  # Combine the results into a single dataframe
  df_cwd_long <- rbindlist(df_cwd_list)

  # Calculate maximum CWD per grid cell
  max_cwd <- df_cwd_long[, .(max_cwd = max(deficit, na.rm = TRUE)), by = .(lon, lat)]
  max_cwd[, `:=`(model = model_name, scenario = scenario_type)]

  # Save the results
  saveRDS(max_cwd, paste0("max_cwd_", model_name, "_", scenario_type, ".rds"), compress = "xz")

  return(max_cwd)
}

# Main script to process all models and scenarios
# Create a list of model-scenario combinations
model_scenario_list <- expand.grid(model_name = model_names, scenario_type = scenario_types, stringsAsFactors = FALSE)

# Use mclapply for parallel processing
num_cores <- detectCores() - 1

results_list <- mclapply(1:nrow(model_scenario_list), function(i) {
  model_name <- model_scenario_list$model_name[i]
  scenario_type <- model_scenario_list$scenario_type[i]

  result <- process_model_cwd(model_name, scenario_type)
  return(result)
}, mc.cores = num_cores)

# Filter out NULL results
results_list <- results_list[!sapply(results_list, is.null)]

# Combine all results into a single data.table
final_results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

# Save the combined results
saveRDS(final_results, "max_cwd_all_models_scenarios.rds", compress = "xz")
