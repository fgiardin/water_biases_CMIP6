# script to plot frequency of SM limitation per model (based on script 3.plot_thetacrit_global)

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
  mutate(count = count*100) %>% # convert to %
  mutate(count = ifelse(Intercept < EFmax, # keep only rows with meaningful intercept (when the EF vs SM relationship decreases)
                        count,
                        0))  # this way 0s represent locations that are never water-limited according to our definition

# make sure models and observations have same cells that are NA
# NAs = cells where the segmented function couldn't find a breakpoint (not enough data at that cell)
# i.e. cells where it wasn't possible to find a SM threshold due to data scarcity
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
breaks <- seq(lower_threshold, upper_threshold, by = 10) # by = rounded_interval

# create color palette
colors_rdpu <- brewer.pal(9, "PuRd") # Get colors from the RdPu palette
colors_rdpu <- colors_rdpu[2:9] # remove first color(s) as too close to white
color_map <- colorRampPalette(colors_rdpu)(length(breaks)) # Interpolate to get the correct number of colors

# assign same color to values above/below threshold (to ensure that the structure of the data is highlighted properly)
dt_count[[plot_variable]][dt_count[[plot_variable]] > upper_threshold] <- upper_threshold
dt_count[[plot_variable]][dt_count[[plot_variable]] < lower_threshold] <- lower_threshold

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

summary_merged <- full_join(dt_count, # remove GRACE in the original dataframe to avoid repetitions
                            grace_data, by = join_by(lon, lat),
                            suffix = c("", "_GRACE")) # rename new column "medianSMmax" + "_GRACE"


# calculate WEIGHTS for the bias (to account for spherical coordinates)
res <- 2.5 # resolution in degrees
ref_grid_size <- sin(res * pi / 180) # calculate reference grid size using the sine of resolution converted to radians
summary_merged <- summary_merged %>% # Add weights to the summary_merged dataframe
  mutate(
    lat1 = lat - (res / 2),  # Calculate the southern boundary of the grid cell
    lat2 = lat + (res / 2),  # Calculate the northern boundary of the grid cell
    gridsize = abs(sin(lat1 * pi / 180) - sin(lat2 * pi / 180)),  # Grid cell size based on sine of latitude boundaries
    weights = gridsize / ref_grid_size  # Normalize weights by the reference grid size
  ) %>%
  dplyr::select(-lat1, -lat2, -gridsize)  # Clean up by removing intermediate columns

# Create separate plots for each model
plot_list <- lapply(sorted_models, function(model) {

  df_model <- filter(summary_merged, model_name == model)

  # calculate stats to compare models-obs
  if(model != "Observations") {

    # calculate R squared
    fit <- lm(count ~ count_GRACE, data = df_model)
    r_squared <- summary(fit)$r.squared
    r_squared_label <- bquote(italic(R)^2 == .(round(r_squared, 2)))

    # calculate abs(mean bias) using weights
    mean_bias_abs <- df_model %>%
      filter(!is.na(count) & !is.na(count_GRACE) & !is.na(weights)) %>% # manually remove NAs when using "sum" (num and den should have same number of elements)
      mutate(mean_bias_abs = abs(count - count_GRACE) * weights) %>%
      summarise(weighted_mean_bias_abs = sum(mean_bias_abs, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
      pull(weighted_mean_bias_abs)
    mean_bias_abs_label <- bquote("|Bias|" == .(round(mean_bias_abs, 0)) * "%")

    # calculate mean bias using weights
    mean_bias <- df_model %>%
      filter(!is.na(count) & !is.na(count_GRACE) & !is.na(weights)) %>%
      mutate(mean_bias = (count - count_GRACE) * weights) %>%
      summarise(weighted_mean_bias = sum(mean_bias, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
      pull(weighted_mean_bias)
    mean_bias_label <- bquote("Bias" == .(round(mean_bias, 0)) * "%")
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

    # Orange-Blue color bar (more intuitive)
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
      legend.key.height = unit(0.5, "cm"),
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
                 ncol = 2, nrow = 4,
                 common.legend = TRUE, # have just one common legend
                 legend="bottom")

filename <- paste0("map_", plot_variable, ".png")
ggsave(filename, path = "./", width = 12, height = 11, dpi= 600)

