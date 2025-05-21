# script to put GLWS and GLDAS from preprocessing in right format

### load library and data
library(tidyverse)
library(R.matlab) # to read .mat
library(ggnewscale) # to add multiple scales in same ggplot
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ncdf4)

# load data prepared by Ryan --------------------------------------------------------
GWLS_data <- readMat("data/Comparisons_GRACE/GWLS_2_0_ALWSC.mat")  # deltaSM already calculated

GRACE_data <- readRDS("data/GRACE/summary_merged_deltaSM.rds") %>%
  dplyr::filter(model_name == "GRACE observations") %>%
  dplyr::select(lon, lat, deltaSMmax)

GLDAS_data <- readMat("data/Comparisons_GRACE/GLDAS_CLSM_ALWSC.mat")

# align GRACE and GLDAS dataframes -----------------------------------------
# Convert GLDAS_data to a DataFrame
lon <- GLDAS_data$lon[, 1]
lat <- GLDAS_data$lat[, 1]
dsmldas_max <- GLDAS_data$dsm.max

# Create a dataframe
GLDAS_df <- expand.grid(lon = lon, lat = lat) %>%
  mutate(dsmldas_max = as.vector(dsmldas_max))

# Step 2: Filter GLDAS_data to GRACE_data range
GLDAS_filtered <- GLDAS_df %>%
  filter(lon >= min(GRACE_data$lon), lon <= max(GRACE_data$lon),
         lat >= min(GRACE_data$lat), lat <= max(GRACE_data$lat))

# Step 3: Aggregate GLDAS_data to Match GRACE Resolution
GLDAS_aggregated <- GLDAS_filtered %>%
  mutate(
    lon_bin = floor((lon + 180) / 2.5) * 2.5 - 178.75,
    lat_bin = floor((lat + 90) / 2.5) * 2.5 - 88.75
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(dsmldas_max = mean(dsmldas_max, na.rm = TRUE)) %>%  # average to same resolution as GRACE
  rename(lon = lon_bin, lat = lat_bin) %>%
  ungroup()

# Final GLDAS data
GLDAS_final <- GLDAS_aggregated

range(GLDAS_final$lat)  # Should match [-53.75, 81.25]
range(GLDAS_final$lon)  # Should match [-178.75, 178.75]

# Focus on vegetated land using GRACE pixels
GRACE_pixels <- GRACE_data %>%
  dplyr::select(lon, lat) %>%
  distinct()
GLDAS_matched <- GLDAS_final %>%
  semi_join(GRACE_pixels, by = c("lon", "lat"))

print(dim(GLDAS_matched))
print(dim(GRACE_data))

saveRDS(GLDAS_matched, "GLDAS_2.5deg.rds", compress = "xz")

# align GRACE and GWLS dataframes -----------------------------------------
# Convert GWLS_data to a DataFrame
lon <- GWLS_data$lon[, 1]
lat <- GWLS_data$lat[, 1]
dsmgw_max <- GWLS_data$dsmgw.max

# Create a dataframe
GWLS_df <- expand.grid(lon = lon, lat = lat) %>%
  mutate(dsmgw_max = as.vector(dsmgw_max))

# Step 2: Filter GWLS_data to GRACE_data range
GWLS_filtered <- GWLS_df %>%
  filter(lon >= min(GRACE_data$lon), lon <= max(GRACE_data$lon),
         lat >= min(GRACE_data$lat), lat <= max(GRACE_data$lat))

# Step 3: Aggregate GWLS_data to Match GRACE Resolution
GWLS_aggregated <- GWLS_filtered %>%
  mutate(
    lon_bin = floor((lon + 180) / 2.5) * 2.5 - 178.75,
    lat_bin = floor((lat + 90) / 2.5) * 2.5 - 88.75
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(dsmgw_max = mean(dsmgw_max, na.rm = TRUE)) %>% # average to same resolution of GRACE
  rename(lon = lon_bin, lat = lat_bin) %>%
  ungroup()


# Final GWLS data
GWLS_final <- GWLS_aggregated

range(GWLS_final$lat)  # Should match [-53.75, 81.25]
range(GWLS_final$lon)  # Should match [-178.75, 178.75]

# focus on vegetated land
GRACE_pixels <- GRACE_data %>% # extract right pixels with GRACE
  dplyr::select(lon, lat) %>%
  distinct()

GWLS_matched <- GWLS_final %>%
  semi_join(GRACE_pixels, by = c("lon", "lat"))

print(dim(GWLS_matched))
print(dim(GRACE_data))

saveRDS(GWLS_matched, "GWLS_2.5deg.rds", compress = "xz")
