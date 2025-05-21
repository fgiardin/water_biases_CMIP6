# Script to extract and resample the CWDX80 dataframe

# CWDX80 ------------------------------------------------------------------

# process data
CWDX80_raw <- rast("data-raw/CWDX80/cwdx80_halfdeg.nc", "cwdx80")

CWDX80 <- CWDX80_raw

# resample CWDX80 to match CMIP6
P <- rast("data-raw/cmip6-ng/pr/mon/g025/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
P <- terra::rotate(P)

CWDX80 <- terra::resample(CWDX80, P[[1]])

# Focus on vegetated land
vegetated_land <- readRDS("data/land_mask/vegetated_land_mask.rds")
CWDX80 <- mask(
  CWDX80,
  vegetated_land,
  maskvalues = 0)

# transform to dataframe
df_CWDX80 <- terra::as.data.frame(CWDX80, xy = TRUE) # xy =TRUE keeps the spatial coordinates

names(df_CWDX80) <- c("lon", "lat", "max_cwd") # rename to be consistent with the rest

saveRDS(df_CWDX80, "max_cwd_.ALEXI.rds", compress = "xz")
