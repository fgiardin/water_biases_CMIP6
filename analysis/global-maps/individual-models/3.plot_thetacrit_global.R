# script to plot global maps of critical soil moisture thresholds per model
# currently not shown in the paper

# load packages
library(tidyverse)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(patchwork)
library(cowplot)
library(ingestr)
library(terra)
library(rgdal)  # for readOGR
library(ggnewscale) # to add multiple scales in same ggplot
library(ggpubr)
sf_use_s2(FALSE) # switch off spherical geometry
library(RColorBrewer)

# load data
df_count_GRACE <- readRDS("data/theta_crit/monthly/df_count_GRACE.rds") %>%
  dplyr::select(-date, -SIF, -TWS) %>%  # remove temporal values
  unique() %>%
  mutate(model = "Observations")

df_count_mrso <- readRDS("data/theta_crit/monthly/df_count_mrso.rds") %>%
  dplyr::select(-date, -EF, -Rn, -mrso, -mrso_norm) %>%
  unique()

summary_thetacrit <- bind_rows(df_count_mrso, df_count_GRACE) %>%
  rename(model_name = model) %>%
  mutate(count = count*100)

plot_variable <- "theta_crit" # either "theta_crit" or "EFmax" (i.e. maximum EF)


# Define breaks for the color scale (around 10 breaks given the upper_threshold)
lower_threshold <- 0.1
upper_threshold <- 1
breaks <- seq(0.1, 1, by = 0.1) # by = rounded_interval

# download countries
countries <- ne_countries(scale = 50, returnclass = c("sf"))

# download ocean outlines
ocean <- ne_download(
  scale = 50,
  type = "ocean",
  category = "physical",
  returnclass = "sf")

background_color = "#cad6df" # set the color corresponding to no data "#ECECED", # "#cad6df", "#D6D6E9"

# create color palette
colors_rdpu <- brewer.pal(9, "PuRd") # Get colors from the RdPu palette
colors_rdpu <- colors_rdpu[2:9] # remove first color(s) as too close to white
color_map <- colorRampPalette(colors_rdpu)(length(breaks)) # Interpolate to get the correct number of colors

# Sort the models alphabetically, keeping "Observations" first
unique_models <- unique(summary_thetacrit$model_name) # Extract unique model names
sorted_models <- c("Observations", sort(unique_models[unique_models != "Observations"]))

# Create separate plots for each model
plot_list <- lapply(sorted_models, function(model) {

  df_model <- filter(summary_thetacrit, model_name == model)

  # Plot global map for each model
  ggplot() +
    theme_void() +
    geom_tile(data = df_model,
              aes(x = lon, y = lat,
                  color = !!sym(plot_variable), fill = !!sym(plot_variable) # this notation extracts the string contained in "plot_variable" from the dataframe
              )) +
    geom_sf(data=ocean, # country borders
            color= "black",
            linetype='solid',
            fill= background_color, #cad6df", #D6D6E9
            size=1) +
    scale_color_gradientn(
      colors = color_map,
      breaks = breaks,
      limits = c(lower_threshold, upper_threshold),
      na.value = "white"
    ) +
    scale_fill_gradientn(
      colors = color_map,
      breaks = breaks,
      limits = c(lower_threshold, upper_threshold),
      na.value = "white"
    ) +
    coord_sf( # cut Antarctica (sorry penguins!)
      xlim = c(-179.999, 179.999),
      ylim = c(-60, 88),
      expand = FALSE) +
    theme(
      plot.title = element_text(hjust = 0.5, size=15), # center title
      legend.text = element_text(color = "black", size=12), # family = "Prata"
      legend.title = element_text(color = "black", size=12),
      plot.margin=unit(c(0,-0.3,0.3,-0.6), 'cm'), # unit(c(top, right, bottom, left), units)
      legend.key.width = unit(2.5, "cm"),
      legend.key.height = unit(0.5, "cm"),
      panel.border = element_rect(colour = "black", fill= NA, linewidth=0.5),
      panel.background = element_rect(fill = background_color, colour = NA) # Setting the background color
    ) +
    guides(fill = guide_colourbar(frame.linewidth = 0.5, ticks.linewidth = 0.5, frame.colour = "black", ticks.colour = "black"),
           color = guide_colourbar(frame.linewidth = 0.5, ticks.linewidth = 0.5, frame.colour = "black", ticks.colour = "black")
    ) +
    labs(title = model,
         color = expression(paste(theta[crit], " (-)   ")),
         fill = expression(paste(theta[crit], " (-)   ")))
         # color = expression(paste(theta[crit], " (m"^3, " m"^-3, ")   ")),
         # fill = expression(paste(theta[crit], " (m"^3, " m"^-3, ")   ")))
})

# Print all plots with one colorbar
all <- ggarrange(plotlist = plot_list,
                 labels = "auto",
                 ncol = 2, nrow = 4,
                 common.legend = TRUE, # have just one common legend
                 legend="bottom")

filename <- paste0("map_", plot_variable, ".png")
ggsave(filename, path = "./", width = 12, height = 11, dpi= 600)
