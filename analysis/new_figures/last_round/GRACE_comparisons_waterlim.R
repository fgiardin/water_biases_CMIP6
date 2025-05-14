# script to compare standardized GRACE to alternative GRACE products

# load libraries
library(terra) # using version 1.7.29
library(tidyverse)
library(data.table)
library(R.matlab)
library(segmented)
library(parallel)
library(ncdf4)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)



# GRACE processing --------------------------------------------------------
# (from script: data-raw/4.calculate_thetacrit_SIF.R)

# process data
GRACE_raw <- rast("data-raw/GRACE/GRCTellus.JPL.200204_202304.GLO.RL06.1M.MSCNv03CRI.nc", "lwe_thickness") # liquid water equivalent
plot(GRACE_raw[[1]])

# extract dates
dates <- terra::time(GRACE_raw) # monthly data 2002-2023

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
vegetated_land <- terra::resample(vegetated_land, GRACE) # resample land_cover to match GRACE
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
P <- rast("data-raw/cmip6-ng/pr/mon/g025/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
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



# processing GLWS ---------------------------------------------------------

GLWS_all <- rast("data-raw/GLWS/v2.0/0.5deg_lat-lon_1m/original/GLWS_2_0_2003-2019.nc")
GLWS_TWS <- rast("data-raw/GLWS/v2.0/0.5deg_lat-lon_1m/original/GLWS_2_0_2003-2019.nc", "Total water storage anomalies")
GLWS_surf <- rast("data-raw/GLWS/v2.0/0.5deg_lat-lon_1m/original/GLWS_2_0_2003-2019.nc", "Surface water")

GLWS_raw <- GLWS_TWS - GLWS_surf
GLWS <- GLWS_raw

# extract dates
nc <- nc_open("data-raw/GLWS/v2.0/0.5deg_lat-lon_1m/original/GLWS_2_0_2003-2019.nc")
tvals      <- ncvar_get(nc, "time") # read the raw time values + units
time_units <- ncatt_get(nc, "time", "units")$value # expressed in years (no comment)
nc_close(nc)
years  <- floor(tvals) # extract year and month
months <- round((tvals - years) * 12 + 0.5)
dates  <- as.Date(sprintf("%04d-%02d-01", years, months)) # build a proper Date vector (effort!)

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
vegetated_land <- terra::resample(vegetated_land, GLWS) # resample land_cover to match GLWS
GLWS <- mask(GLWS, vegetated_land) # remove all pixels that are NAs in land_cover
GLWS
plot(GLWS[[1]])
plot(vegetated_land)

# normalize TWS by location (to compare it with flux SM)
GLWSmax <- max(GLWS, na.rm = TRUE)
GLWSmin <- min(GLWS, na.rm = TRUE)

GLWS_norm <- (GLWS - GLWSmin) / (GLWSmax - GLWSmin)
GLWS_final <- ifel(GLWS_norm > 0, GLWS_norm, NA) # remove instances when total soil moisture = 0 (not physically meaningful)

# resample GLWS to match CMIP6
P <- rast("data-raw/cmip6-ng/pr/mon/g025/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
P <- terra::rotate(P)
GLWS_final <- terra::resample(GLWS_final, P[[1]])

# transform to dataframe
df_GLWS <- terra::as.data.frame(GLWS_final, xy = TRUE) # xy =TRUE keeps the spatial coordinates
names(df_GLWS) <- c("lon", "lat", as.character(dates))

df_GLWS_long <- df_GLWS %>% # pivot_longer with data.table (faster)
  data.table() %>%
  melt(
    measure.vars = as.character(dates), # name of the columns to be pivoted
    variable.name = "date", # name of the new column where the former column names will be pivoted to (in this case: the dates)
    value.name = "TWS" # name of the new column where the values will be pivoted to (in this case: the water balance)
  ) %>%
  mutate(date = lubridate::date(as.character(date)),
         date = floor_date(date, "month") # floor (round) date at the first of the month (so we can merge with SIF)
  )


# processing GLDAS --------------------------------------------------------

GLDAS_raw <- rast("data-raw/GLDAS-2_CLSM/SoilMoist_P_tavg_2003_2023.nc", "SoilMoist_P_tavg")
GLDAS <- GLDAS_raw

# extract dates
dates <- terra::time(GLDAS_raw) # monthly data 2002-2023

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
vegetated_land <- terra::resample(vegetated_land, GLDAS) # resample land_cover to match GLDAS
GLDAS <- mask(GLDAS, vegetated_land) # remove all pixels that are NAs in land_cover
GLDAS
plot(GLDAS[[1]])
plot(vegetated_land)

# normalize TWS by location (to compare it with flux SM)
GLDASmax <- max(GLDAS, na.rm = TRUE)
GLDASmin <- min(GLDAS, na.rm = TRUE)

GLDAS_norm <- (GLDAS - GLDASmin) / (GLDASmax - GLDASmin)
GLDAS_final <- ifel(GLDAS_norm > 0, GLDAS_norm, NA) # remove instances when total soil moisture = 0 (not physically meaningful)

# resample GLDAS to match CMIP6
P <- rast("data-raw/cmip6-ng/pr/mon/g025/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
P <- terra::rotate(P)
GLDAS_final <- terra::resample(GLDAS_final, P[[1]])

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


# put together and plot: GLWS ---------------------------------------------------

# merge
df_merged_GLWS <- inner_join(
  df_GRACE_long %>% rename(TWS_GRACE = TWS),
  df_GLWS_long  %>% rename(TWS_GLWS  = TWS),
  by = c("lon", "lat", "date")
) %>%
  drop_na()

### scatter plot ###
# calculate R2
fit    <- lm(TWS_GLWS ~ TWS_GRACE, data = df_merged_GLWS) # fit linear model
r2     <- summary(fit)$r.squared # extract R2

# raw error metrics
err    <- df_merged_GLWS$TWS_GLWS - df_merged_GLWS$TWS_GRACE
rmse   <- sqrt(mean(err^2))
bias   <- mean(err)

# normalize by the mean of GRACE
mean_grace <- mean(df_merged_GLWS$TWS_GRACE)
rRMSE  <- rmse  / mean_grace
rBias  <- bias  / mean_grace

# scatter plot
scatter_GLWS <- ggplot(df_merged_GLWS, aes(x = TWS_GRACE, y = TWS_GLWS)) +
  geom_point(alpha = 0.4) +
  # 1:1 line
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "grey40") +
  # regression line
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  # labels
  labs(
    x = "GRACE TWS (-)",
    y = "GLWS  TWS (-)",
    title = "GLWS vs GRACE"
  ) +
  annotate("label",
           x     = 0,
           y     = 1,
           label = sprintf(
             "R²    = %.2f\nrRMSE = %.2f\nrBias = %.2f",
             r2, rRMSE, rBias
           ),
           hjust      = 0, vjust = 1,
           fill       = alpha("white", 0.8),  # slightly transparent
           color      = "black",
           fontface   = "bold",
           label.size = 0.2                   # box border thickness
  ) +
  # annotate(
  #   "text",
  #   x    = min(df_merged$TWS_GRACE, na.rm = TRUE),
  #   y    = max(df_merged$TWS_GLWS, na.rm = TRUE),
  #   hjust = 0, vjust = 1,
  #   label = sprintf(
  #     "R²    = %.2f\nrRMSE = %.2f\nrBias = %.2f",
  #     r2, rRMSE, rBias
  #   )
  # ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# save
ggsave("scatter_GRACE-GLWS.png", plot = scatter_GLWS,
       width = 4, height = 4, dpi = 600)


### map of R2 ###
# Compute R2 in every cell cell
df_model <- df_merged_GLWS %>%
  group_by(lon, lat) %>%
  summarize(
    # R2 of TWS_GLWS ~ TWS_GRACE
    count = cor(TWS_GLWS, TWS_GRACE, use = "complete.obs")^2,
    .groups = "drop"
  )

# get coastlines for land lines
coast <- ne_coastline(scale = "medium", returnclass = "sf")

# plot
p <- ggplot() +
  theme_void() +
  geom_tile(data = df_model,
            aes(x = lon, y = lat, fill = count)) +
  geom_sf(data = coast,
          color = "black", fill = NA, size = 0.3) +
  coord_sf(
    xlim   = c(-179.999, 179.999),
    ylim   = c(-60, 88),
    expand = FALSE
  ) +
  scale_fill_gradient(
    name   = expression(R^2),
    low    = "gray80",    # very light
    high   = "#00008B",      # dark
    limits = c(0, 1),     # force legend to span 0–1
    breaks = c(0, 0.25, 0.5, 0.75, 1),     # show only end‐points
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  guides(
    fill = guide_colourbar(
      title.position   = "top",
      title.hjust      = 0.5,
      barwidth         = 15,    # make it long
      barheight        = 0.5,   # thin strip
      frame.colour     = "black",
      ticks            = TRUE,
      ticks.linewidth  = 0.5
    )
  ) +
  theme(
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title.align = 0.5,         # center title above bar
    legend.margin      = margin(t = 5),
    plot.title       = element_text(hjust = 0.5, size = 15),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 15),
    legend.key.width = unit(2.5, "cm"),
    legend.key.height= unit(0.4, "cm"),
    panel.border     = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin      = unit(c(0, -0.3, 0.3, -0.6), "cm")
  ) +
  guides(fill = guide_colourbar(frame.linewidth = 0.5,
                                ticks.linewidth = 0.5,
                                frame.colour = "black",
                                ticks.colour = "black")) +
  labs(
    title = "GLWS vs GRACE", # Global R² per cell:
    fill  = expression(R^2)
  )

# save
ggsave("GRACE-GLWS_R2_Global.png", plot = p, path = "./", width = 6, height = 3.25, dpi = 300)



# put together and plot: GLDAS -------------------------------------------------------------------

df_merged_GLDAS <- inner_join(
  df_GRACE_long %>% rename(TWS_GRACE = TWS),
  df_GLDAS_long %>% rename(TWS_GLDAS = TWS),
  by = c("lon", "lat", "date")
) %>%
  drop_na()


### scatter ###
# fit and stats
fit    <- lm(TWS_GLDAS ~ TWS_GRACE, data = df_merged_GLDAS)
r2     <- summary(fit)$r.squared

err    <- df_merged_GLDAS$TWS_GLDAS - df_merged_GLDAS$TWS_GRACE
n      <- nrow(df_merged_GLDAS)
rmse   <- sqrt(sum(err^2) / n)
bias   <- sum(err)    / n
mean_g <- mean(df_merged_GLDAS$TWS_GRACE)
rRMSE  <- rmse / mean_g
rBias  <- bias / mean_g

# scatter
scatter_GLDAS <- ggplot(df_merged_GLDAS, aes(x = TWS_GRACE, y = TWS_GLDAS)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  labs(
    x     = "GRACE TWS (-)",
    y     = "GLDAS  TWS (-)",
    title = "GLDAS vs GRACE"
  ) +
  annotate("label",
           x     = 0,
           y     = 1,
           label = sprintf(
             "R²    = %.2f\nrRMSE = %.2f\nrBias = %.2f",
             r2, rRMSE, rBias
           ),
           hjust      = 0, vjust = 1,
           fill       = alpha("white", 0.8),  # slightly transparent
           color      = "black",
           fontface   = "bold",
           label.size = 0.2                   # box border thickness
  ) +
  # annotate(
  #   "text",
  #   x    = min(df_merged_GLDAS$TWS_GRACE, na.rm = TRUE),
  #   y    = max(df_merged_GLDAS$TWS_GLDAS, na.rm = TRUE),
  #   hjust = 0, vjust = 1,
  #   label = sprintf(
  #     "R²    = %.2f\nrRMSE = %.2f\nrBias = %.2f",
  #     r2, rRMSE, rBias
  #   )
  # ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# save
ggsave("scatter_GRACE-GLDAS.png", scatter_GLDAS,
       width = 4, height = 4, dpi = 600)


### map ###
# compute per-cell R2
df_model_GLDAS <- df_merged_GLDAS %>%
  group_by(lon, lat) %>%
  summarize(
    count = cor(TWS_GLDAS, TWS_GRACE, use = "complete.obs")^2,
    .groups = "drop"
  )

# coastline for land outlines
coast <- ne_coastline(scale = "medium", returnclass = "sf")

# plot
p_GLDAS <- ggplot() +
  theme_void() +
  geom_tile(data = df_model_GLDAS,
            aes(x = lon, y = lat, fill = count)) +
  geom_sf(data = coast, color = "black", fill = NA, size = 0.3) +
  coord_sf(xlim = c(-179.999, 179.999),
           ylim = c(-60, 88),
           expand = FALSE) +
  scale_fill_gradient(
    name   = expression(R^2),
    low    = "gray80",    # very light
    high   = "#00008B",      # dark
    limits = c(0, 1),     # force legend to span 0–1
    breaks = c(0, 0.25, 0.5, 0.75, 1),     # show only end‐points
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  guides(
    fill = guide_colourbar(
      title.position   = "top",
      title.hjust      = 0.5,
      barwidth         = 15,    # make it long
      barheight        = 0.5,   # thin strip
      frame.colour     = "black",
      ticks            = TRUE,
      ticks.linewidth  = 0.5
    )
  ) +
  theme(
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title.align = 0.5,         # center title above bar
    legend.margin      = margin(t = 5),
    plot.title       = element_text(hjust = 0.5, size = 15),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 15),
    legend.key.width = unit(2.5, "cm"),
    legend.key.height= unit(0.4, "cm"),
    panel.border     = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin      = unit(c(0, -0.3, 0.3, -0.6), "cm")
  ) +
  guides(fill = guide_colourbar(frame.linewidth = 0.5,
                                ticks.linewidth = 0.5,
                                frame.colour = "black",
                                ticks.colour = "black")) +
  labs(
    title = "GLDAS vs GRACE",
    fill  = expression(R^2)
  )

# save
ggsave("GRACE-GLDAS_R2_Global.png", p_GLDAS,
       width = 6, height = 3.25, dpi = 300)


