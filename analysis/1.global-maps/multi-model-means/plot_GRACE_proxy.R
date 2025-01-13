# script to plot alternative to GRACE product for peer review

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
GWLS_data <- readMat("data/GRACE/GWLS_2_0_ALWSC.mat")  # deltaSM already calculated
GRACE_data <- readRDS("data/GRACE/summary_merged_deltaSM.rds") %>%
  dplyr::filter(model_name == "GRACE observations") %>%
  dplyr::select(lon, lat, deltaSMmax)

GWLS_data <- readRDS("data/GLWS/deltaSM_GWLS.rds")



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


# Plot --------------------------------------------------------------------

### Prepare plot
# assign same color to points above upper threshold
upper_threshold <- 1000

# Define the desired breaks for the color scale (around 10 breaks given the upper_threshold)
interval <- upper_threshold / 9
rounded_interval <- round(interval / 100) * 100
breaks <- seq(0, upper_threshold, by = rounded_interval)

# Cap values to 1000 for plotting (outliers)
GWLS_final <- GWLS_data %>%
  mutate(dsmgw_max = pmin(deltaSMmax/10, upper_threshold)) # pmin: replaces all values greater than 1000 with 1000
GRACE_final <- GRACE_data %>%
  mutate(deltaSMmax = pmin(deltaSMmax, upper_threshold))

# download countries
countries <- ne_countries(scale = 50, returnclass = c("sf"))

# download ocean outlines
ocean <- ne_download(
  scale = 50,
  type = "ocean",
  category = "physical",
  returnclass = "sf")

### Plot global map GRACE
a <- ggplot() +
  theme_void() +
  geom_tile(data = GRACE_final,
            aes(x = lon, y = lat,
                color = deltaSMmax, fill = deltaSMmax # this notation extracts the string contained in "plot_variable" from the dataframe
            )) +
  geom_sf(data=ocean, # country borders
          color= "black",
          linetype='solid',
          fill= "white", #cad6df", #D6D6E9
          size=1) +
  scale_color_viridis_c(option = "turbo",
                        breaks = breaks,
                        limits = c(0, upper_threshold)
  ) +
  scale_fill_viridis_c(option = "turbo",
                       breaks = breaks,
                       limits = c(0, upper_threshold)
  ) +
  coord_sf( # cut Antarctica (sorry penguins!)
    xlim = c(-179.999, 179.999),
    ylim = c(-60, 88),
    expand = FALSE) +
  theme(
    plot.title = element_text(hjust = 0.5, size=15), # center title
    # legend.position = "none",
    legend.text = element_text(color = "black", size=12), # family = "Prata"
    legend.title = element_text(color = "black", size=12),
    plot.margin=unit(c(0,-0.3,0.3,-0.6), 'cm'), # unit(c(top, right, bottom, left), units)
    legend.key.width = unit(2.5, "cm"), # control size of colorbar
    legend.key.height = unit(0.4, "cm"),
    panel.border = element_rect(colour = "black", fill= NA, linewidth=0.5)) +
  guides(fill = guide_colourbar(frame.linewidth = 0.5,
                                ticks.linewidth = 0.5,
                                frame.colour = "black",
                                ticks.colour = "black"),
         color = guide_colourbar(frame.linewidth = 0.5,
                                 ticks.linewidth = 0.5,
                                 frame.colour = "black",
                                 ticks.colour = "black")) +
  labs(title = "GRACE", # fill and color same label --> only one colorbar in the legend
       color = expression(Delta*SM[max]~"(mm)"),
       fill = expression(Delta*SM[max]~"(mm)"))


### Plot global map GWLS
b <- ggplot() +
  theme_void() +
  geom_tile(data = GWLS_final,
            aes(x = lon, y = lat,
                color = dsmgw_max, fill = dsmgw_max # this notation extracts the string contained in "plot_variable" from the dataframe
            )) +
  geom_sf(data=ocean, # country borders
          color= "black",
          linetype='solid',
          fill= "white", #cad6df", #D6D6E9
          size=1) +
  scale_color_viridis_c(option = "turbo",
                        breaks = breaks,
                        limits = c(0, upper_threshold)
  ) +
  scale_fill_viridis_c(option = "turbo",
                       breaks = breaks,
                       limits = c(0, upper_threshold)
  ) +
  coord_sf( # cut Antarctica (sorry penguins!)
    xlim = c(-179.999, 179.999),
    ylim = c(-60, 88),
    expand = FALSE) +
  theme(
    plot.title = element_text(hjust = 0.5, size=15), # center title
    # legend.position = "none",
    legend.text = element_text(color = "black", size=12), # family = "Prata"
    legend.title = element_text(color = "black", size=12),
    plot.margin=unit(c(0,-0.3,0.3,-0.6), 'cm'), # unit(c(top, right, bottom, left), units)
    legend.key.width = unit(2.5, "cm"), # control size of colorbar
    legend.key.height = unit(0.4, "cm"),
    panel.border = element_rect(colour = "black", fill= NA, linewidth=0.5)) +
  guides(fill = guide_colourbar(frame.linewidth = 0.5,
                                ticks.linewidth = 0.5,
                                frame.colour = "black",
                                ticks.colour = "black"),
         color = guide_colourbar(frame.linewidth = 0.5,
                                 ticks.linewidth = 0.5,
                                 frame.colour = "black",
                                 ticks.colour = "black")) +
  labs(title = "GWLS: sm + gw", # fill and color same label --> only one colorbar in the legend
       color = expression(Delta*SM[max]~"(mm)"),
       fill = expression(Delta*SM[max]~"(mm)"))




### Put together and plot
plot_list <- list(a, b)

all <- ggarrange(plotlist = plot_list,
                 labels = "auto", # "auto"
                 ncol = 2, nrow = 1,
                 common.legend = TRUE, # have just one common legend
                 legend="bottom")

ggsave("GRACE_GWLS_comparison.png", path = "./", width = 12, height = 3.25, dpi= 600)






