# script to calculate deltaSM for every year in every grid cell and
# then take the maximum deltaSM over the entire period

# the code was run in parallel on a local machine, but can be also run on a HPC
# given the size of the raw data, we didn't upload it on this repo.
# please refer to the "Data availability statement" in the paper to download the data

# Load required libraries
library(terra)
library(tidyverse)
library(data.table)
library(lubridate)
library(parallel)
library(ncdf4)
library(ncdf4.helpers)


# CMIP6 -------------------------------------------------------------------

# Use mclapply for parallel processing (locally)
num_cores <- detectCores() - 1 # Get the number of available cores

# Define model names
model_names <- c("CESM2",
                 "EC-Earth3-Veg",
                 "IPSL-CM6A-LR",
                 "UKESM1-0-LL",
                 "CNRM-ESM2-1",
                 "E3SM-1-1",
                 "CMCC-ESM2",
                 "CNRM-CM6-1",
                 "MIROC6"
                 )

# Use lapply to process data for each model
mclapply(model_names, function(model_name) {

  print(model_name)

  # extract data from netCDFs

  # load the data
  if (model_name %in% c("CNRM-CM6-1", "CNRM-ESM2-1", "UKESM1-0-LL")) {
    # CNRM-CM6-1, CNRM-ESM2-1 and UKESM1-0-LL: use realization "r1i1p1f2" ("r1i1p1f1" not available)
    mrso_raw <- rast(paste0("data-raw/cmip6-ng/mrso/mon/g025/mrso_mon_", model_name, "_land-hist_r1i1p1f2_g025.nc"))

  } else if (model_name == "E3SM-1-1") {
    # "E3SM-1-1": use realization "r1i1p11f1" ("r1i1p1f1" not available)
    mrso_raw <- rast(paste0("data-raw/cmip6-ng/mrso/mon/g025/mrso_mon_", model_name, "_land-hist_r1i1p11f1_g025.nc"))

  } else {
    # for all other models use realization "r1i1p1f1"
    mrso_raw <- rast(paste0("data-raw/cmip6-ng/mrso/mon/g025/mrso_mon_", model_name, "_land-hist_r1i1p1f1_g025.nc"))
  }


  # Only take 1900-2014
  dates <- terra::time(mrso_raw)
  selection <- which(dates >= as.POSIXct("2003-01-01")) # match dates with GRACE
  mrso_model <- subset(mrso_raw, subset = selection)

  # Rotate to true coordinates (-180;180 instead of 0;360)
  mrso_model <- terra::rotate(mrso_model)

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
  vegetated_land <- resample(vegetated_land, mrso_model)
  mrso_model <- mask(mrso_model, vegetated_land)

  # Transform to dataframe
  df_model <- terra::as.data.frame(mrso_model, xy = TRUE)
  dates <- terra::time(mrso_model)

  names(df_model) <- c("lon", "lat", as.character(dates))

  df_model_long <- df_model %>%
    data.table() %>%
    melt(
      measure.vars = as.character(dates),
      variable.name = "date",
      value.name = "mrso"
    ) %>%
    mutate(date = lubridate::date(as.character(date)),
           year = lubridate::year(date))

  # Calculate deltaSM per year per pixel
  delta_SM <- df_model_long %>%
    group_by(lon, lat, year) %>%
    summarise(
      deltaSM = max(mrso, na.rm = TRUE) - min(mrso, na.rm = TRUE)
    ) %>%
    ungroup()

  # In every pixel, calculate the max deltaSM over all years
  delta_SMmax <- delta_SM %>%
    group_by(lon, lat) %>%
    summarise(
      deltaSMmax = max(deltaSM, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    setDT() # transform to datatable

  # calculate absolute deltaSM (directly max-min for a specific location across all years, maximum amplitude)
  delta_SMabs <- df_model_long %>%
    group_by(lon, lat) %>%
    summarise(
      deltaSMabs = max(mrso, na.rm = TRUE) - min(mrso, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    setDT() # transform to datatable

  # Merge the two data.tables based on lon and lat
  merged_data_model <- delta_SMmax[delta_SMabs, nomatch=0, on=.(lon, lat)]

  # Save the results
  saveRDS(merged_data_model, paste0("deltaSM_", model_name, ".rds"), compress = "xz")
},
mc.cores = num_cores)

# GRACE -------------------------------------------------------------------
dname = "lwe_thickness"
# see data structure
GRACE_nc <- nc_open("data-raw/GRACE/GRCTellus.JPL.200204_202304.GLO.RL06.1M.MSCNv03CRI.nc")
fillvalue <- ncatt_get(GRACE_nc, dname, "_FillValue")
fillvalue$value # -99999
GRACE_units <- ncatt_get(GRACE_nc, dname,"units")
GRACE_units$value # "cm"
nc_close(GRACE_nc)

# process data
GRACE_raw <- rast("data-raw/GRACE/GRCTellus.JPL.200204_202304.GLO.RL06.1M.MSCNv03CRI.nc", "lwe_thickness") # liquid water equivalent
plot(GRACE_raw[[1]])

# extract dates
dates <- terra::time(GRACE_raw) # monthly data 2002-2023

# # match dates with CMIP
# selection <- which(dates >= as.POSIXct("2003-01-01") &
#                      dates <= as.POSIXct("2014-12-31"))
# GRACE_raw <- subset(GRACE_raw, subset = selection)

# rotate to true coordinates (-180;180 instead of 0;360)
GRACE <- terra::rotate(GRACE_raw)

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
vegetated_land <- resample(vegetated_land, GRACE) # resample land_cover to match water_balance
GRACE <- mask(GRACE, vegetated_land) # remove all pixels that are NAs in land_cover
GRACE
plot(GRACE[[1]])
plot(vegetated_land)

# resample GRACE to match CMIP6
P <- rast("data-raw/cmip6-ng/pr/mon/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
P <- terra::rotate(P)
GRACE <- terra::resample(GRACE, P[[1]])

# transform to dataframe
df_GRACE <- terra::as.data.frame(GRACE, xy = TRUE) # xy =TRUE keeps the spatial coordinates
dates <- terra::time(GRACE)

names(df_GRACE) <- c("lon", "lat", as.character(dates))

df_GRACE_long <- df_GRACE %>% # pivot_longer with data.table (faster)
  data.table() %>%
  melt(
    measure.vars = as.character(dates), # name of the columns to be pivoted
    variable.name = "date", # name of the new column where the former column names will be pivoted to (in this case: the dates)
    value.name = "lwe" # name of the new column where the values will be pivoted to (in this case: the water balance)
  ) %>%
  mutate(date = lubridate::date(as.character(date)),
         year = lubridate::year(date))

# calculate deltaSM per year per pixel
delta_SM_GRACE <- df_GRACE_long %>%
  group_by(lon, lat, year) %>%
  summarise(
    deltaSM = max(lwe, na.rm = TRUE) - min(lwe, na.rm = TRUE)
  ) %>%
  ungroup()

# in every pixel, calculate the max deltaSM over all years
delta_SMmax_GRACE <- delta_SM_GRACE %>%
  group_by(lon, lat) %>%
  summarise(
    deltaSMmax = max(deltaSM, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    deltaSMmax = deltaSMmax * 10 # convert cm to mm
  ) %>%
  setDT() # transform to datatable

# calculate absolute deltaSM (directly max-min for a specific location across all years, maximum amplitude)
delta_SM_GRACE_abs <- df_GRACE_long %>%
  group_by(lon, lat) %>%
  summarise(
    deltaSMabs = max(lwe, na.rm = TRUE) - min(lwe, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    deltaSMabs = deltaSMabs * 10 # convert cm to mm
  ) %>%
  setDT() # transform to datatable

# Merge the two data.tables based on lon and lat
merged_data_GRACE <- delta_SMmax_GRACE[delta_SM_GRACE_abs, nomatch=0, on=.(lon, lat)]

saveRDS(merged_data_GRACE, "deltaSM_.GRACE.rds", compress = "xz")







