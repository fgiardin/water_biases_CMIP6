# script to plot frequency of SM limitation (multi-model mean)

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
library(data.table)
library(RColorBrewer)

# load data
df_count_GRACE <- readRDS("data/theta_crit/monthly/df_count_GRACE.rds") %>%
  dplyr::select(-date, -SIF, -TWS) %>%  # remove temporal values
  unique() %>%
  mutate(model = "Observations")

df_count_mrso <- readRDS("data/theta_crit/monthly/df_count_mrso.rds") %>%
  dplyr::select(-date, -EF, -Rn, -mrso, -mrso_norm) %>%
  unique()

dt_count <- bind_rows(df_count_mrso, df_count_GRACE) %>%
  rename(model_name = model) %>%
  mutate(count = count*100) %>%
  mutate(count = ifelse(Intercept < EFmax,
                        count, # keep only rows with meaningful intercept (when the EF vs SM relationship decreases)
                        0))   # this way 0s represent locations that are never water-limited according to our definition!

# make sure all models have same cells that are NA
# NAs = cells where the segmented function couldn't find a breakpoint (not enough data at that cell)
setDT(dt_count) # Convert to data.table if it's not already
na_locations <- unique(dt_count[is.na(theta_crit), .(lon, lat)]) # Identify all unique lon and lat combinations where any of the columns of interest are NA
dt_count[na_locations, # Update the original data.table to set the values to NA for the identified lon and lat, for all models
         on = .(lon, lat),
         `:=` (theta_crit = NA_real_,
               EFmax = NA_real_,
               Slope = NA_real_,
               Intercept = NA_real_,
               count = NA_real_)]

plot_variable <- "count"

# Define breaks for the color scale (around 10 breaks given the upper_threshold)
lower_threshold <- 0
upper_threshold <- 100
# interval <- upper_threshold / 9
# rounded_interval <- round(interval, digits = 2)
breaks <- seq(lower_threshold, upper_threshold, by = 10) # by = rounded_interval

# create color palette
colors_rdpu <- brewer.pal(9, "PuRd") # Get colors from the RdPu palette
colors_rdpu <- colors_rdpu[2:9] # remove first color(s) as too close to white
color_map <- colorRampPalette(colors_rdpu)(length(breaks)) # Interpolate to get the correct number of colors

# download countries
countries <- ne_countries(scale = 50, returnclass = c("sf"))

coastline <- ne_coastline(scale = 50, returnclass = c("sf"))

# download ocean outlines
ocean <- ne_download(
  scale = 50,
  type = "ocean",
  category = "physical",
  returnclass = "sf")


# Sort the models alphabetically, keeping "Observations" first
unique_models <- unique(dt_count$model_name) # Extract unique model names
sorted_models <- c("Observations", sort(unique_models[unique_models != "Observations"]))

background_color = "white" # set the color corresponding to no data "#ECECED", # "#cad6df", "#D6D6E9"


# Create a new (repeated for every model) column with obs data in the original dataframe
grace_data <- dt_count %>%
  filter(model_name == "Observations") %>%
  dplyr::select(lon, lat, count)

# calculate multi-model mean
MMmeans <- dt_count %>%
  filter(model_name != "Observations") %>%
  group_by(lon, lat) %>%
  summarise(
    count = mean(count, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model_name = "Multi-model mean") # add model name to mirror summary_deltaSM

# bind rows to summary_deltaSM
dt_count_MMmeans <- dt_count %>%
  bind_rows(MMmeans)

summary_merged <- full_join(dt_count_MMmeans, # remove GRACE in the original dataframe to avoid repetitions
                            grace_data, by = join_by(lon, lat),
                            suffix = c("", "_GRACE")) # rename new column "medianSMmax" + "_GRACE"

# create the list of models to print
sorted_models <- c("Observations", "Multi-model mean")

# Create separate plots for each model
plot_list <- lapply(sorted_models, function(model) {

  df_model <- filter(summary_merged, model_name == model)

  # calculate stats to compare models-obs
  if(model != "Observations") {

    # calculate R squared
    fit <- lm(count ~ count_GRACE, data = df_model)
    r_squared <- summary(fit)$r.squared
    r_squared_label <- bquote(italic(R)^2 == .(round(r_squared, 2)))

    # calculate abs(mean bias)
    mean_bias_abs <- df_model %>%
      mutate(mean_bias_abs = abs(count - count_GRACE)) %>%
      pull(mean_bias_abs) %>%
      mean(na.rm = TRUE)
    mean_bias_abs_label <- bquote("|Bias|" == .(round(mean_bias_abs, 2)) ~ "(%)")

    # calculate mean bias
    mean_bias <- df_model %>%
      mutate(mean_bias = count - count_GRACE) %>%
      pull(mean_bias) %>%
      mean(na.rm = TRUE)
    mean_bias_label <- bquote("Bias" == .(round(mean_bias, 2)) ~ "(%)")
  }

  # Plot global map for each model
  p <- ggplot() +
    theme_void() +
    geom_tile(data = df_model,
              aes(x = lon, y = lat,
                  color = !!sym(plot_variable), fill = !!sym(plot_variable) # this notation extracts the string contained in "plot_variable" from the dataframe
              )) +
    geom_sf(data=ocean, # country borders
            color= "black",
            linetype='solid',
            fill= background_color,
            size=1) +

    # Orange-Blue color bar
    scale_color_viridis_c(
      option = "turbo",
      breaks = breaks,
      limits = c(lower_threshold, upper_threshold)
    ) +
    scale_fill_viridis_c(
      option = "turbo",
      breaks = breaks,
      limits = c(lower_threshold, upper_threshold)
    ) +
    coord_sf( # cut Antarctica (sorry penguins!)
      xlim = c(-179.999, 179.999),
      ylim = c(-60, 88),
      expand = FALSE) +
    theme(
      plot.title = element_text(hjust = 0.5, size=15), # center title
      # legend.position = "none",
      legend.text = element_text(color = "black", size=12), # family = "Prata"
      legend.title = element_text(color = "black", size=15),
      plot.margin=unit(c(0,-0.3,0.3,-0.6), 'cm'), # unit(c(top, right, bottom, left), units)
      legend.key.width = unit(2.5, "cm"), # control size of colorbar
      legend.key.height = unit(0.4, "cm"),
      panel.border = element_rect(colour = "black", fill= NA, linewidth=0.5),
      panel.background = element_rect(fill = background_color, colour = NA) # Setting the background color
    ) +
    guides(fill = guide_colourbar(frame.linewidth = 0.5, ticks.linewidth = 0.5, frame.colour = "black", ticks.colour = "black"),
           color = guide_colourbar(frame.linewidth = 0.5, ticks.linewidth = 0.5, frame.colour = "black", ticks.colour = "black")
    ) +
    labs(
      title = model,
      color = "Frequency of water limitation (%)",
      fill = "Frequency of water limitation (%)",
      # fill = expression(paste("number of months when SM < ", theta[crit], " (-) "))
    )

  # more specific title rather than just "Observations"
  if(model == "Observations") {
    p <- p + labs(
      title = "GOME-2 SIF and GRACE observations")
  }

  # add stats for models
  if(model != "Observations") {
    p <- p +
      annotate("text", x = -175, y = -25, label = r_squared_label, hjust = 0, vjust = 0, size = 4.3) + # R2
      annotate("text", x = -175, y = -40, label = mean_bias_label, hjust = 0, vjust = 0, size = 4.3) + # mean bias
      annotate("text", x = -175, y = -55, label = mean_bias_abs_label, hjust = 0, vjust = 0, size = 4.3) # abs(mean bias)
  }

  return(p)

})

# Print all plots with one colorbar
all <- ggarrange(plotlist = plot_list,
                 labels = "auto",
                 ncol = 2, nrow = 1,
                 common.legend = TRUE, # have just one common legend
                 legend="bottom")

filename <- paste0("map_", plot_variable, "_MMmean.png")
ggsave(filename, path = "./", width = 12, height = 3.25, dpi= 600)

