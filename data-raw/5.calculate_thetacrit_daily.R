# script to extract variables and calculate theta_crit and EFmax for CMIP6 models using DAILY data (only 4 models available)

# Load libraries
devtools::load_all(".")
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(segmented)
library(data.table)
library(parallel)
library(terra)
library(lubridate)
library(R.matlab) # to read .mat files


# find grid cells corresponding to flux locations in CMIP maps ------------------------------------

# load world map of cmip6 as reference
cmip6_world <- terra::rast("data-raw/cmip6-ng/hfls/mon/g025/hfls_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
cmip6_world <- terra::rotate(cmip6_world[[1]]) # rotate to have longitude in -180;180 instead of 0;360 (and only take first layer)

# load list of lon/lat fluxnet
df_flux <- readRDS("data/flux_data/theta_crit_flux_norm.rds") # critical SM thresholds calculated at 128 FLUXNET2015 sites

# Create SpatVector of flux locations
locations <- vect(df_flux[, c("lon", "lat")])
crs(locations) <- crs(cmip6_world) # set projection as in cmip6

# extract locations using terra
cell_coord <- terra::extract(cmip6_world, locations, xy = TRUE, method = "simple") %>%
  rename(lon_cmip6 = x, # rename x and y to lon and lat
         lat_cmip6 = y)

# Bind the cell lon/lat to the original df_flux dataframe
df_key <- cbind(df_flux %>% dplyr::select(sitename, lon, lat), # df_key: df containing exact fluxnet2015 locations with corresponding grid cell
                cell_coord[, c("lon_cmip6", "lat_cmip6")]) %>%
  setDT()

saveRDS(df_key, "df_key.rds", compress = "xz")

# process daily cmip6 data ---------------------------------------------------

# models with daily data and all variables available
models <- c("UKESM1-0-LL", "IPSL-CM6A-LR", "EC-Earth3-Veg", "CNRM-ESM2-1")
lon_mat <- readMat("data-raw/cmip6-ng/mrsol/mon/integrated_1m/lon.mat")
lat_mat <- readMat("data-raw/cmip6-ng/mrsol/mon/integrated_1m/lat.mat")

# source function to extract daily data at flux locations from cmip6 (df_key must be loaded in the environment)
source("R/process_daily_cmip6.R")

# Run the function on each model in parallel
num_cores <- 4 # Use mclapply for parallel processing
model_list <- mclapply(models,
                       process_daily_cmip6, # to run the function df_key (see above) should be loaded in the environment
                       mc.cores = num_cores)

saveRDS(model_list, "daily_cmip6_flux_locations.rds", compress = "xz")

# Combine all the results into a single data.table
final_dt <- rbindlist(model_list)

# only keep EF between 0 and 1 + seasonality filter
final_DT <- final_DT[EF > 0 & EF <= 1.5 & Rn > 100]


# calculate theta_crit per site and model ---------------------------------

filtered_DT <- final_DT %>%
  dplyr::select(-Rn)

# > head(filtered_DT)
# sitename       model       date      lon      lat        EF       SM
# 1:   AR-SLu CNRM-ESM2-1 2000-01-04 -66.4598 -33.4648 0.6837789 194.1172
# 2:   AR-SLu CNRM-ESM2-1 2000-01-06 -66.4598 -33.4648 0.5711509 233.1164
# 3:   AR-SLu CNRM-ESM2-1 2000-01-07 -66.4598 -33.4648 0.7159951 201.5759
# 4:   AR-SLu CNRM-ESM2-1 2000-01-13 -66.4598 -33.4648 0.6808712 209.2566
# 5:   AR-SLu CNRM-ESM2-1 2000-01-15 -66.4598 -33.4648 0.4244146 167.3753
# 6:   AR-SLu CNRM-ESM2-1 2000-01-16 -66.4598 -33.4648 0.7565675 221.0521
saveRDS(data, "filtered_DT.rds", compress = "xz")

# normalize soil moisture by model and location (to compare it with fluxdata!)
filtered_DT <- filtered_DT %>%
  group_by(model, sitename) %>%
  mutate(SM = (SM - min(SM)) / (max(SM) - min(SM))) %>%
  ungroup() %>%
  setDT()

# split the data into a nested dt with one dataframe per model
split_data <- split(filtered_DT, filtered_DT$model)

# Set the number of cores to available - 1
num_cores <- detectCores() - 1

# Take all possible combinations of model/lon/lat as they appear in data.table
combinations_dt <- unique(filtered_DT[, .(model, sitename)])

# focus on less combinations (faster)
# combinations_dt <- combinations_dt[1:4000]

# apply the function in parallel
source("R/fit_bilinear.R")
source("R/fit_bilinear_from_combination.R")

# split_data/filtered_DT contains original lon/lat of fluxnet locations
results_list <- mclapply(1:nrow(combinations_dt), function(i) {

  fit_bilinear_from_combination(combinations_dt[i], split_data, "EF", "SM")

}, mc.cores = num_cores)

# Filter out list elements that are not data.frames or data.tables
filtered_list <- results_list[sapply(results_list, function(x) is.data.frame(x) || is.data.table(x))]

# unnest dataframe for plotting
final_results <- rbindlist(filtered_list)

# attach results to original dataframe + calculate count
dt_count <- filtered_DT %>%
  left_join(final_results, by = join_by(sitename, model)) %>%
  drop_na() %>%
  group_by(model, sitename) %>% # calculate the ratio of days under SM limitation
  mutate(
    count = sum(SM < theta_crit, na.rm = TRUE) / n() # first calculate the number of rows/days where there is SM limitation, then divide by the total number of days!
  ) %>%
  ungroup()

saveRDS(dt_count, "cmip6_daily_theta_crit_count.rds", compress = "xz")


