# script to calculate CWD in every grid cell of CMIP6 data and extract maximum CWD per grid cell

# the code was run in paralle on a local machine, but can be also run on a HPC
# given the size of the raw data, we didn't upload it on this repo.
# please refer to the "Data availability statement" in the paper to download the data

# libraries
devtools::load_all(".")
library(terra)
library(tidyverse)
library(lubridate)
library(data.table)
library(parallel)
library(ncdf4)
library(ncdf4.helpers)

# Set model names
model_names <- c("CESM2",
                 "EC-Earth3-Veg",
                 "IPSL-CM6A-LR",
                 "UKESM1-0-LL",
                 "CNRM-ESM2-1",
                 "E3SM-1-1",
                 "CMCC-ESM2")

# Loop through all model names
lapply(model_names, function(model_name) {

  print(model_name)

  # extract data from netCDFs -----------------------------------------------

  # load the data
  if (model_name %in% c("CNRM-ESM2-1", "UKESM1-0-LL")) {
    # "CNRM-ESM2-1" and "UKESM1-0-LL": use realization "r1i1p1f2" ("r1i1p1f1" not available, but they are equivalent)
    ET <- rast(paste0("data-raw/cmip6-ng/evspsbl/mon/g025/evspsbl_mon_", model_name, "_land-hist_r1i1p1f2_g025.nc"))
    P <- rast(paste0("data-raw/cmip6-ng/pr/mon/g025/pr_mon_", model_name, "_land-hist_r1i1p1f2_g025.nc"))

  } else if (model_name == "E3SM-1-1") {
    # "E3SM-1-1": use realization "r1i1p11f1" ("r1i1p11f1" not available, but they are equivalent)
    ET <- rast(paste0("data-raw/cmip6-ng/evspsbl/mon/g025/evspsbl_mon_", model_name, "_land-hist_r1i1p11f1_g025.nc"))
    P <- rast(paste0("data-raw/cmip6-ng/pr/mon/g025/pr_mon_", model_name, "_land-hist_r1i1p11f1_g025.nc"))

  } else if (model_name == "CMCC-ESM2") {
    # "CMCC-ESM2": use precipitation from CESM2 (within "land-hist" it's the same metereological forcing)
    ET <- rast(paste0("data-raw/cmip6-ng/evspsbl/mon/g025/evspsbl_mon_", model_name, "_land-hist_r1i1p1f1_g025.nc"))
    P <- rast("data-raw/cmip6-ng/pr/mon/g025/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")

  } else {
    # for all other models use relization "r1i1p1f1"
    ET <- rast(paste0("data-raw/cmip6-ng/evspsbl/mon/g025/evspsbl_mon_", model_name, "_land-hist_r1i1p1f1_g025.nc"))
    P <- rast(paste0("data-raw/cmip6-ng/pr/mon/g025/pr_mon_", model_name, "_land-hist_r1i1p1f1_g025.nc"))
  }

  # calculate water balance
  water_balance_raw <- P - ET

  # only take 1900-2014
  dates <- terra::time(water_balance_raw)
  selection <- which(dates >= as.POSIXct("1935-01-01")) # create a selection based on a threshold (80 years of data, to be consistent with observational CWD)
  water_balance <- subset(water_balance_raw, subset = selection)

  # rotate to true coordinates (-180;180 instead of 0;360)
  water_balance <- terra::rotate(water_balance)

  # focus on vegetated land
  land_cover_raw <- rast("data-raw/landcover/landcover_MCD12C1.nc") # load land cover
  land_cover <- flip(land_cover_raw[[1]],
                     direction = "vertical") # flip it vertically
  vegetated_land <- ifel(
    land_cover == 0,
    NA,
    ifel(
      land_cover == 13, # 13 = Urban and Built-up Lands
      NA,
      ifel(
        land_cover > 14, # 14 = Cropland/Natural Vegetation Mosaics (keep)
        NA,              # 15 = Permanent Snow and Ice, 16 = Barren
        land_cover
      )
    )
  )

  vegetated_land <- resample(vegetated_land, water_balance) # resample land_cover to match water_balance
  water_balance <- mask(water_balance, vegetated_land) # remove all pixels that are NAs in land_cover
  # water_balance
  # plot(water_balance[[1]])

  # transform to dataframe
  wb_selection <- water_balance

  df_wb <- terra::as.data.frame(wb_selection, xy = TRUE) # xy =TRUE keeps the spatial coordinates
  dates <- terra::time(wb_selection)

  names(df_wb) <- c("lon", "lat", as.character(dates))

  df_wb_long <- df_wb %>% # pivot_longer with data.table
    data.table() %>%
    melt(
      measure.vars = as.character(dates), # name of the columns to be pivoted
      variable.name = "date", # name of the new column where the former column names will be pivoted to (in this case: the dates)
      value.name = "water_balance" # name of the new column where the values will be pivoted to (in this case: the water balance)
    ) %>%
    mutate(date = lubridate::date(as.character(date)))
  saveRDS(df_wb_long, paste0("df_wb_long_", model_name, ".rds"), compress = "xz")


  # calculate CWD -----------------------------------------------------------
  print(paste0("Calculating CWD of ", model_name, " now..."))

  # Convert dataframe to a data.table
  setDT(df_wb_long)

  # convert water_balance to mm
  df_wb_long[, water_balance := water_balance * 86400] # number of seconds in a day (must be further multiplied by 30 -- aka "number of seconds in a month", since we're using monthly data --> multiplied at a later stage of the workflow)

  # Group the data by lon and lat using data.table
  df_grouped <- df_wb_long[, .(grouped_data = list(.SD)), by = .(lon, lat)]

  # Split the data into a list of data.tables based on lon and lat
  df_grouped_list <- split(df_wb_long, by = c("lon", "lat"))

  # Define a function to apply mct_max to each group
  apply_mct <- function(group_data) {
    mct(group_data,
            "water_balance",
            "date")
  }

  # Parallelize the computation using mclapply
  num_cores <- detectCores()-1
  df_cwd_list <- mclapply(df_grouped_list, # list of CWD time series at each lon,lat location
                              apply_mct,
                              mc.cores = num_cores)
  saveRDS(df_cwd_list, paste0("df_cwd_list_", model_name, ".rds"), compress = "xz")


  # Combine the results into a single dataframe
  df_cwd_long <- rbindlist(df_cwd_list)

  # Group the combined data.table by lon and lat and calculate the maximum at each lon,lat location
  max_cwd <- df_cwd_long[, .(max_cwd = max(deficit)), # conver to mm
                            by = .(lon, lat)]

  saveRDS(max_cwd, paste0("max_cwd_", model_name, ".rds"), compress = "xz")

})



# CWDX80 ------------------------------------------------------------------

# process data
CWDX80_raw <- rast("data-raw/CWDX80/cwdx80_halfdeg.nc", "cwdx80")

CWDX80 <- CWDX80_raw

# focus on vegetated land
land_cover_raw <- rast("data-raw/landcover/landcover_MCD12C1.nc") # load land cover
land_cover <- flip(land_cover_raw[[1]], # select correct layer (year: 2018)
                   direction="vertical") # flip it vertically
vegetated_land <- ifel( # only keep vegetated land
  land_cover == 0,
  NA,
  ifel(
    land_cover == 13, # 13 = Urban and Built-up Lands
    NA,
    ifel(
      land_cover > 14, # 14 = Cropland/Natural Vegetation Mosaics (keep)
      NA,              # 15 = Permanent Snow and Ice, 16 = Barren
      land_cover
    )
  )
)
vegetated_land <- terra::resample(vegetated_land, CWDX80) # resample land_cover to match water_balance
CWDX80 <- mask(CWDX80, vegetated_land) # remove all pixels that are NAs in land_cover

# resample CWDX80 to match CMIP6
P <- rast("data-raw/cmip6-ng/pr/mon/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
P <- terra::rotate(P)

CWDX80 <- terra::resample(CWDX80, P[[1]])

# transform to dataframe
df_CWDX80 <- terra::as.data.frame(CWDX80, xy = TRUE) # xy =TRUE keeps the spatial coordinates

names(df_CWDX80) <- c("lon", "lat", "max_cwd") # rename to be consistent with the rest

saveRDS(df_CWDX80, "max_cwd_.ALEXI.rds", compress = "xz")





