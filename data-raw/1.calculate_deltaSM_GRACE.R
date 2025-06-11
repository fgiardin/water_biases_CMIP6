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

# Define model names and scenario types
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

# Use mclapply for parallel processing (locally)
num_cores <- detectCores() - 1 # Get the number of available cores

# source function to process data consistently
source("R/process_model_deltaSM.R")

# Use mclapply to process data for each model in parallel
mclapply(model_names, function(model_name) {
  for (scenario in scenario_types) {
    process_model_deltaSM(model_name, scenario)
  }
}, mc.cores = num_cores)

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

# resample GRACE to match CMIP6
P <- rast("data-raw/cmip6-ng/pr/mon/g025/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
P <- terra::rotate(P)
GRACE <- terra::resample(GRACE, P[[1]])

# apply land mask
vegetated_land <- readRDS("data/land_mask/vegetated_land_mask.rds")
GRACE <- mask(GRACE,
              vegetated_land,
              maskvalues = 0) # indicate to the function that the mask value (indicates which cells should be masked) is 0 (default is NA)

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
         year = lubridate::year(date)) %>%
  dplyr::filter(date >= as.Date("2003-01-01") & date <= as.Date("2014-12-31")) # only take study years (intersection between GRACE and CMIP6)


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











