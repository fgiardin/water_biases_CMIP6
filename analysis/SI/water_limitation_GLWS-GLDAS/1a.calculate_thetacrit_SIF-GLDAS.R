# script to calculate theta_crit in SIF/SIFmax vs GRACE alternative products (GLDAS/GLDAS)
# to be compared to EF vs SM (mrso, total column soil moisture) from CMIP6

# load libraries
library(terra) # using version 1.7.29
library(tidyverse)
library(data.table)
library(R.matlab)
library(segmented)
library(parallel)

# data info
# GRACE: monthly data 2002-2023
# SIF GOME-2: 01-2007 to 09-2016
# focus on 9 years: 01-2007 to 12-2015 (included)


# SIF processing ----------------------------------------------------------

# load world map of cmip6
cmip6_world <- terra::rast("data-raw/cmip6-ng/hfls/mon/g025/hfls_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
cmip6_world <- terra::rotate(cmip6_world[[1]]) # only take sample world map (it's just for resampling) + rotate so longitude is standard

# load SIF
SIF = readMat("data-raw/SIF/GOME2_SIF_Monthly_01_2007_09_2016.mat") # GOME-2 SIF
str(SIF)
# 1x1 degree, lon x lat x months
# take 1-108 to have full years 01/2007 to 12/2015 (9 years)

# reformat data to be read by terra (for consistency with processing of CMIP6)
data_matrix <- SIF$SIF.L3
str(data_matrix) # lon and lat are inverted --> transpose matrix
data_matrix <- aperm(data_matrix, c(2, 1, 3)) # switch second and first dimensions, so that it can easily be read with package terra!
str(data_matrix)

# Define the extent: xmin, xmax, ymin, ymax
extent <- ext(-180, 180, -90, 90)

# Create Date objects for each time slice
dates <- seq(as.Date("2007-01-01"), by = "month", length.out = dim(data_matrix)[3])

# Create an empty SpatRaster with the correct dimensions and extent
r <- rast(nrows=180, ncols=360, nlyr=dim(data_matrix)[3], ext=extent)
r

# Set the resolution to 1 degree
res(r) <- c(1, 1)

# Fill the SpatRaster with data from the matrix
values(r) <- data_matrix

# Set the time attribute
time(r) <- dates
# plot(r[[1]])

SIF_rast <- r
# plot(SIF_rast[[1]])

# only take data from 01/01/2007 to 01/12/2015 (9 full years)
SIF_rast <- subset(SIF_rast, 1:108)

# remove SIF negative values (consistent with previous studies with same dataset)
SIF_rast <- ifel(SIF_rast < 0, NA, SIF_rast)

# focus on growing season
# The growing season is defined here for each pixel as the months in which the climatological mean
# is greater than or equal to half of the climatological maximum.

# calculate (climatological) maximum SIF in every pixel
SIF_max <- max(SIF_rast, na.rm = TRUE)

# take only points that are greater than or equal to 0.5 of the climatological maximum
SIF_season <- ifel(SIF_rast >= 0.5*SIF_max, SIF_rast, NA)

# calculate SIF/SIFmax (dividing by maximum value in time)
SIFmax <- max(SIF_season, na.rm = TRUE) # recalculate maximum (with growing season only)

SIF_norm <- SIF_season / SIFmax # normalize each pixel

# resample to match projection and extent of CMIP6
SIF_norm <- terra::resample(SIF_norm,
                            cmip6_world,
                            method="average")

# apply land mask
vegetated_land <- readRDS("data/land_mask/vegetated_land_mask.rds")
SIF_masked <- mask(SIF_norm,
                   vegetated_land,
                   maskvalues = 0) # indicate to the function that the mask value (indicates which cells should be masked) is 0 (default is NA)

# transform to dataframe for plotting and merging with other data
df_SIF <- terra::as.data.frame(SIF_masked, xy = TRUE)
dates <- terra::time(SIF_masked)
names(df_SIF) <- c("lon", "lat", as.character(dates))

# Pivot to long format
df_SIF_long <- df_SIF %>%
  data.table() %>%
  melt(
    measure.vars = as.character(dates),
    variable.name = "date",
    value.name = "SIF"
  ) %>%
  mutate(date = lubridate::date(as.character(date))) %>% # first convert to date (otherwise "melt" converts it to a factor)
  mutate(date = lubridate::floor_date(date, unit = "month")) %>%  # round down date to first day of each month
  drop_na()

# saveRDS(df_SIF_long, "df_SIF.rds", compress = "xz")


# GLDAS processing --------------------------------------------------------

GLDAS_raw <- rast("data-raw/GLDAS-2_CLSM/SoilMoist_P_tavg_2003_2023.nc", "SoilMoist_P_tavg")
GLDAS <- GLDAS_raw

# extract dates
dates <- terra::time(GLDAS_raw) # monthly data 2002-2023

# normalize TWS by location (to compare it with flux SM)
GLDASmax <- max(GLDAS, na.rm = TRUE)
GLDASmin <- min(GLDAS, na.rm = TRUE)

GLDAS_norm <- (GLDAS - GLDASmin) / (GLDASmax - GLDASmin)
GLDAS_final <- ifel(GLDAS_norm > 0, GLDAS_norm, NA) # remove instances when total soil moisture = 0 (not physically meaningful)

# resample GLDAS to match CMIP6
P <- rast("data-raw/cmip6-ng/pr/mon/g025/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
P <- terra::rotate(P)
GLDAS_final <- terra::resample(GLDAS_final,
                               P[[1]],
                               method = "average")

# apply land mask
vegetated_land <- readRDS("data/land_mask/vegetated_land_mask.rds") # already in the same grid as CMIP6!!
GLDAS_final <- mask(GLDAS_final,
                     vegetated_land,
                     maskvalues = 0)

# transform to dataframe
df_GLDAS <- terra::as.data.frame(GLDAS_final, xy = TRUE) # xy =TRUE keeps the spatial coordinates
dates <- terra::time(GLDAS_final)

names(df_GLDAS) <- c("lon", "lat", as.character(dates))

df_GLDAS_long <- df_GLDAS %>% # pivot_longer with data.table (faster)
  data.table() %>%
  melt(
    measure.vars = as.character(dates), # name of the columns to be pivoted
    variable.name = "date", # name of the new column where the former column names will be pivoted to (in this case: the dates)
    value.name = "TWS" # name of the new column where the values will be pivoted to (in this case: the water balance)
  ) %>%
  mutate(date = lubridate::date(as.character(date)),
         date = floor_date(date, "month") # floor (round) date at the first of the month (so we can merge with SIF)
  )


# fit bilinear regression per location ----------------------------

# merge two dataframes and create a new column (to use the function "fit_bilinear_from_combination")
dt_final <- df_SIF_long %>%
  left_join(df_GLDAS_long,
            by = join_by(lon, lat, date)) %>%
  mutate(model = "OBS") %>%
  drop_na() %>%
  dplyr::filter(date < "2015-01-01") # keep dates until 2014-12-31 (consistent with CMIP6)

# split the data into a nested dt with one dataframe per model
split_data <- split(dt_final, dt_final$model)

# Set the number of cores to available - 1
num_cores <- detectCores() - 1

# Take all possible combinations of model/lon/lat as they appear in data.table
combinations_dt <- unique(dt_final[, .(model, lon, lat)])

# apply the function in parallel
source("R/fit_bilinear.R")
source("R/fit_bilinear_from_combination.R")

results_list_GLDAS <- mclapply(1:nrow(combinations_dt), function(i) {

  fit_bilinear_from_combination(combinations_dt[i], split_data, "SIF", "TWS")

}, mc.cores = num_cores)

# Filter out list elements that are not data.frames or data.tables
filtered_list <- results_list_GLDAS[sapply(results_list_GLDAS, function(x) is.data.frame(x) || is.data.table(x))]

# unnest dataframe for plotting
map_theta_crit_SIF_GLDAS <- rbindlist(filtered_list, fill = TRUE)


saveRDS(map_theta_crit_SIF_GLDAS, "theta_crit_GLDAS.rds", compress = "xz")



# count frequency of SM limitation ----------------------------
df_count <- dt_final %>%
  left_join(map_theta_crit_SIF_GLDAS, by = join_by(lon, lat, model)) %>%
  dplyr::select(-model) %>%
  group_by(lon, lat) %>%
  mutate(
    count = ifelse(
      all(is.na(theta_crit)), NA, # Keep NA if all theta_crit values are NA (locations where it's never limited)
      sum(TWS < theta_crit, na.rm = TRUE) / n() # Otherwise, calculate percentage of time under SM limitation
    )
  ) %>%
  ungroup()

saveRDS(df_count, "df_count_GLDAS.rds", compress = "xz")


