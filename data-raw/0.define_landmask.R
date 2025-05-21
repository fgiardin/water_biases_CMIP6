# script to define land mask to use across all datasets

# for consistency, land mask will be at 2.5 degrees (same as LMIP data)
# and will be applied after resampling the dataframe at CMIP6 resolution

# load libraries
library(terra) # using version 1.7.29
library(tidyverse)
library(data.table)
library(R.matlab)
library(segmented)
library(parallel)

# load a random LMIP-CMIP6 dataset just for resolution and origin
P <- rast("data-raw/cmip6-ng/pr/mon/g025/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
P <- terra::rotate(P)

# land mask filter: focus on vegetated land
land_cover_raw <- rast("data-raw/landcover/landcover_MCD12C1.nc")
land_cover <- flip(land_cover_raw[[1]], direction = "vertical")
vegetated_land <- ifel(
  land_cover == 0,
  0, # important: use zeros and not NAs here (important for averaging below!)
  ifel(
    land_cover == 13,
    0,
    ifel(
      land_cover > 14,
      0,
      1 # create binary mask (vegetation: yes or no?)
    )
  )
)

# convert to LMIP-CMIP 2.5 degree grid
vegetated_land <- round( # 2) ...and then round to 1 if >0.5 (vegetated land)
  terra::resample( # 1) do average of 0-1 pixels first (while resampling)...
    vegetated_land,
    P,
    method = "average"))

saveRDS(vegetated_land, "vegetated_land_mask.rds", compress = "xz")

# important! it's a 0-1 mask to it should be used accordingly (not NA-1)

