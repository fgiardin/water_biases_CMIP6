# script to calculate critical soil moisture thresholds for CMIP6 models for every pixel using monthly data

# Load libraries
devtools::load_all(".")
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(segmented)
library(data.table)
library(parallel)
library(terra)
library(lubridate)
library(R.matlab) # to read .mat files

# adapt start and end date based on availability of observations
start_date = as.Date("2007-01-01")
end_date = as.Date("2014-12-31")

# variables needed: mrso (total column soil moisture), hfls, rsds, rsus, rlds, rlus
# Rn = (rsds – rsus) + (rlds – rlus) [W/m2]
# EF = hfls [W/m2] / Rn [W/m2]

# load and format model data -----------------------------------------------------

# Use mclapply for parallel processing
num_cores <- detectCores() - 2

# Define the model names
model_names <- c("CESM2",
                 "EC-Earth3-Veg",
                 "IPSL-CM6A-LR",
                 "UKESM1-0-LL",
                 "CNRM-ESM2-1",
                 "E3SM-1-1",
                 "CMCC-ESM2"
                 )

# Define the variables to extract for each model
variables <- c("mrso", "hfls", "rsds", "rsus", "rlds", "rlus")

# List to store the results
results <- list()

results <- mclapply(model_names, function(model_name) {

  print(paste("Processing model", "CESM2"))

  # Load all variables for a model
  data_list <- list()
  for (variable in variables) {

    if (model_name %in% c("CNRM-ESM2-1", "UKESM1-0-LL")) {
      path <- paste0("data-raw/cmip6-ng/", variable, "/mon/g025/", variable, "_mon_", model_name, "_land-hist_r1i1p1f2_g025.nc")
    } else if (model_name == "E3SM-1-1") {
      path <- paste0("data-raw/cmip6-ng/", variable, "/mon/g025/", variable, "_mon_", model_name, "_land-hist_r1i1p11f1_g025.nc")
    } else {
      path <- paste0("data-raw/cmip6-ng/", variable, "/mon/g025/", variable, "_mon_", model_name, "_land-hist_r1i1p1f1_g025.nc")
    }

    data_raw <- rast(path)

    # filter dates
    dates <- terra::time(data_raw)
    selection <- which(dates >= start_date & dates <= end_date) # apply filters on date!
    data_model <- subset(data_raw, subset = selection)

    # Rotate to true coordinates (-180;180 instead of 0;360)
    data_model <- terra::rotate(data_model)

    # Focus on vegetated land
    land_cover_raw <- rast("data-raw/landcover/landcover_MCD12C1.nc")
    land_cover <- flip(land_cover_raw[[1]], direction = "vertical")
    vegetated_land <- ifel(
      land_cover == 0,
      NA,
      ifel(
        land_cover == 13,
        NA,
        ifel(
          land_cover > 14,
          NA,
          land_cover
        )
      )
    )
    vegetated_land <- resample(vegetated_land, data_model)
    data_model <- mask(data_model, vegetated_land)

    data_list[[variable]] <- data_model
  }

  # Calculate Rn and EF using raster operations (faster)
  Rn <- (data_list$rsds - data_list$rsus) + (data_list$rlds - data_list$rlus)
  EF <- data_list$hfls / Rn

  # Filter Rn raster to get only values where Rn > 75 (remove winter months at high latitudes + growing season filter)
  Rn_filtered <- ifel(Rn > 75, Rn, NA)
  EF_filtered <- mask(EF, Rn_filtered)
  mrso_filtered <- mask(data_list$mrso, Rn_filtered)

  # Convert rasters to data frames
  df_Rn <- terra::as.data.frame(Rn_filtered, xy=TRUE)
  df_EF <- terra::as.data.frame(EF_filtered, xy=TRUE)
  df_mrso <- terra::as.data.frame(mrso_filtered, xy=TRUE)

  # extract date and rename columns
  dates <- terra::time(Rn_filtered) # dates should be the same after filtering!
  names(df_Rn) <- c("lon", "lat", as.character(dates))
  names(df_EF) <- c("lon", "lat", as.character(dates))
  names(df_mrso) <- c("lon", "lat", as.character(dates))

  # Pivot df_Rn to long format
  df_Rn_long <- df_Rn %>%
    data.table() %>%
    melt(
      measure.vars = as.character(dates),
      variable.name = "date",
      value.name = "Rn"
    ) %>%
    mutate(date = lubridate::date(as.character(date))) %>% # first convert to date (otherwise "melt" converts it to a factor)
    mutate(date = lubridate::floor_date(date, unit = "month")) %>%  # round down date to first day of each month (days of the month are 01 or 15 or 16 between models, and this would generate duplicate months)
    drop_na()

  # Transform df_EF to long format
  df_EF_long <- df_EF %>%
    data.table() %>%
    melt(
      measure.vars = as.character(dates),
      variable.name = "date",
      value.name = "EF"
    ) %>%
    mutate(date = lubridate::date(as.character(date))) %>% # first convert to date (otherwise "melt" converts it to a factor)
    mutate(date = lubridate::floor_date(date, unit = "month")) %>%
    drop_na()

  # Transform df_mrso to long format
  df_mrso_long <- df_mrso %>%
    data.table() %>%
    melt(
      measure.vars = as.character(dates),
      variable.name = "date",
      value.name = "mrso"
    ) %>%
    mutate(date = lubridate::date(as.character(date))) %>% # first convert to date (otherwise "melt" converts it to a factor)
    mutate(date = lubridate::floor_date(date, unit = "month")) %>%
    drop_na()

  # combine the data with data.table
  result_data <- df_EF_long[
    df_Rn_long,
    on = .(lon, lat, date),
    nomatch = NA
  ][
    df_mrso_long,
    on = .(lon, lat, date),
    nomatch = NA
  ]

  return(result_data)

}, mc.cores = num_cores)

# convert each df in the list() of results to data.table format
results <- lapply(results, as.data.table)
names(results) <- model_names # assign model names for the next step!

# combine results from all models together and automatically create a column corresponding to the id (in this case the model name!)
dt_allresults <- rbindlist(results, idcol = "model") # The idcol argument creates a new column called "model" that contains the list element names (i.e., model names) for each row




# fit bilinear regression per location ------------------------------------

# objective: have a global map of theta_crit per model

dt_final <- dt_allresults

# normalize mrso by model and location
dt_final <- dt_final %>%
  group_by(model, lon, lat) %>%
  mutate(mrso_norm = (mrso - min(mrso)) / (max(mrso) - min(mrso))) %>%
  ungroup() %>%
  setDT()

# split the data into a nested dt with one dataframe per model
split_data <- split(dt_final, dt_final$model)

# Set the number of cores to available - 1
num_cores <- detectCores() - 1

# num_cores <- 3

# Take all possible combinations of model/lon/lat as they appear in data.table
combinations_dt <- unique(dt_final[, .(model, lon, lat)])

# combinations_dt <- combinations_dt[1:4000]

# apply the function in parallel
source("R/fit_bilinear.R")
source("R/fit_bilinear_from_combination.R")

results_list_mrso <- mclapply(1:nrow(combinations_dt), function(i) {

  fit_bilinear_from_combination(combinations_dt[i], split_data, "EF", "mrso_norm")

}, mc.cores = num_cores)

# Filter out list elements that are not data.frames or data.tables (errors)
filtered_list_mrso <- results_list_mrso[sapply(results_list_mrso, function(x) is.data.frame(x) || is.data.table(x))]

# unnest dataframe for plotting
final_results_mrso <- rbindlist(filtered_list_mrso, fill = TRUE)

saveRDS(final_results_mrso, "theta_crit_mrso_allmodels.rds", compress = "xz")


# final_results_mrso <- theta_crit_mrso_allmodels

# count percentage of time under SM limitation ----------------------------
df_count <- dt_final %>%
  left_join(final_results_mrso, by = join_by(lon, lat, model)) %>%
  group_by(lon, lat, model) %>%
  mutate(
    count = ifelse(
      all(is.na(theta_crit)), NA, # Keep NA if all theta_crit values are NA (locations where it's never limited)
      sum(mrso_norm < theta_crit, na.rm = TRUE) / n() # Otherwise, calculate percentage of time under SM limitation
    )
  ) %>%
  ungroup()

saveRDS(df_count, "df_count_mrso.rds", compress = "xz")











