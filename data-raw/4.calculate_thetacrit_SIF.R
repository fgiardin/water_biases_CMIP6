# script to calculate theta_crit in SIF/SIFmax vs GRACE plots
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
cmip6_world <- rotate(cmip6_world[[1]]) # only take sample world map (it's just for resampling) + rotate so longitude is standard

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
vegetated_land <- terra::resample(vegetated_land, SIF_rast)
SIF_rast <- mask(SIF_rast, vegetated_land)

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
SIF_norm <- terra::resample(SIF_norm, cmip6_world, method="bilinear")

# transform to dataframe for plotting and merging with other data
df_SIF <- terra::as.data.frame(SIF_norm, xy = TRUE)
dates <- terra::time(SIF_norm)
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


# GRACE processing --------------------------------------------------------

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
land_cover <- flip(land_cover_raw[[1]], # select correct layer
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
vegetated_land <- terra::resample(vegetated_land, GRACE) # resample land_cover to match water_balance
GRACE <- mask(GRACE, vegetated_land) # remove all pixels that are NAs in land_cover
GRACE
plot(GRACE[[1]])
plot(vegetated_land)

# normalize TWS by location (to compare it with flux SM)
GRACEmax <- max(GRACE, na.rm = TRUE)
GRACEmin <- min(GRACE, na.rm = TRUE)

GRACE_norm <- (GRACE - GRACEmin) / (GRACEmax - GRACEmin)
GRACE_final <- ifel(GRACE_norm > 0, GRACE_norm, NA) # remove instances when total soil moisture = 0 (not physically meaningful)

# resample GRACE to match CMIP6
P <- rast("data-raw/cmip6-ng/pr/mon/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
P <- terra::rotate(P)
GRACE_final <- terra::resample(GRACE_final, P[[1]])

# transform to dataframe
df_GRACE <- terra::as.data.frame(GRACE_final, xy = TRUE) # xy =TRUE keeps the spatial coordinates
dates <- terra::time(GRACE_final)

names(df_GRACE) <- c("lon", "lat", as.character(dates))

df_GRACE_long <- df_GRACE %>% # pivot_longer with data.table (faster)
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
  left_join(df_GRACE_long,
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

results_list_GRACE <- mclapply(1:nrow(combinations_dt), function(i) {

  fit_bilinear_from_combination(combinations_dt[i], split_data, "SIF", "TWS")

}, mc.cores = num_cores)

# Filter out list elements that are not data.frames or data.tables
filtered_list <- results_list_GRACE[sapply(results_list_GRACE, function(x) is.data.frame(x) || is.data.table(x))]

# unnest dataframe for plotting
map_theta_crit_SIF_GRACE <- rbindlist(filtered_list, fill = TRUE)


saveRDS(map_theta_crit_SIF_GRACE, "theta_crit_GRACE.rds", compress = "xz")



# count frequency of SM limitation ----------------------------
df_count <- dt_final %>%
  left_join(map_theta_crit_SIF_GRACE, by = join_by(lon, lat, model)) %>%
  dplyr::select(-model) %>%
  group_by(lon, lat) %>%
  mutate(
    count = ifelse(
      all(is.na(theta_crit)), NA, # Keep NA if all theta_crit values are NA (locations where it's never limited)
      sum(TWS < theta_crit, na.rm = TRUE) / n() # Otherwise, calculate percentage of time under SM limitation
    )
  ) %>%
  ungroup()

saveRDS(df_count, "df_count_GRACE.rds", compress = "xz")


