# script to plot land water storage as quantified by delta_SM using cmip6 data and GRACE
# results per model

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
library(yardstick)
sf_use_s2(FALSE) # switch off spherical geometry

# Define which variable to plot: deltaSMabs or deltaSMmax
# deltaSMabs: absolute deltaSM over 2002-2014
# deltaSMmax: first calculate deltaSM per year and then take maximum among all years (used in main text)
plot_variable <- "deltaSMmax" # (!!! changing to deltaSMabs will still plot the stats of deltaSMmax!)

# Define the directory path where the .rds files are located
directory_path <- "data/deltaSM/"

# get list of files we want to merge (here: output dataframes from script 2)
file_locations <- list.files(directory_path, # define directory where all files are
                             glob2rx("deltaSM_*.rds"), # search for this pattern inside directory
                             full.names = TRUE, # list all paths
                             recursive = FALSE # NOT go through all subdirectories
)

# Initialize an empty list to store dataframes
dataframes_list <- list()

# Loop through the files, read each .rds file, and store the dataframes in the list
for (file_path in file_locations) {
  # Extract the model name from the file name
  model_name <- gsub(".*deltaSM_(.*)\\.rds", "\\1", basename(file_path))

  # Read the dataframe from the .rds file
  df <- readRDS(file_path)

  # Add the model name as a new column in the dataframe
  df$model_name <- model_name

  # Append the dataframe to the list
  dataframes_list[[model_name]] <- df
}

# Combine all dataframes in the list into a single dataframe
summary_deltaSM <- bind_rows(dataframes_list)

# # save it for reproducibility
# saveRDS(summary_deltaSM, "summary_deltaSM_allyears.rds", compress = "xz")

upper_threshold <- 1000

# Define the desired breaks for the color scale (around 10 breaks given the upper_threshold)
interval <- upper_threshold / 9
rounded_interval <- round(interval / 100) * 100
breaks <- seq(0, upper_threshold, by = rounded_interval)

# assign same color to all outliers (to ensure that the structure of the data i highlighted properly)
summary_deltaSM[[plot_variable]][summary_deltaSM[[plot_variable]] > upper_threshold] <- upper_threshold

# download countries
countries <- ne_countries(scale = 50, returnclass = c("sf"))

# download ocean outlines
ocean <- ne_download(
  scale = 50,
  type = "ocean",
  category = "physical",
  returnclass = "sf")

# rename ".GRACE"
summary_deltaSM <- summary_deltaSM %>%
  mutate(model_name = ifelse(model_name == ".GRACE", "GRACE observations", model_name))
model_list <- unique(summary_deltaSM$model_name)

# Create a new (repeated for every model) column with GRACE data in the original dataframe
# only to calculate model-obs stats (not for plotting!)
grace_data <- summary_deltaSM %>%
  filter(model_name == "GRACE observations") %>%
  dplyr::select(-model_name)

summary_merged <- full_join(summary_deltaSM,
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
plot_list <- lapply(model_list, function(model) {

  df_model <- filter(summary_merged, model_name == model)

  # calculate stats to compare models-obs
  if(model != "GRACE Observations") {

    # calculate R squared
    fit <- lm(deltaSMmax ~ deltaSMmax_GRACE, data = df_model)
    r_squared <- summary(fit)$r.squared
    r_squared_label <- bquote(italic(R)^2 == .(round(r_squared, 2)))

    # calculate abs(mean bias) using weights
    mean_bias_abs <- df_model %>%
      filter(!is.na(deltaSMmax) & !is.na(deltaSMmax_GRACE) & !is.na(weights)) %>% # manually remove NAs when using "sum" (num and den should have same number of elements)
      mutate(mean_bias_abs = abs(deltaSMmax - deltaSMmax_GRACE) * weights) %>%
      summarise(weighted_mean_bias_abs = sum(mean_bias_abs, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
      pull(weighted_mean_bias_abs)
    mean_bias_abs_label <- bquote("|Bias|" == .(round(mean_bias_abs, 0)) ~ "mm")

    # calculate mean bias using weights
    mean_bias <- df_model %>%
      filter(!is.na(deltaSMmax) & !is.na(deltaSMmax_GRACE) & !is.na(weights)) %>%
      mutate(mean_bias = (deltaSMmax - deltaSMmax_GRACE) * weights) %>%
      summarise(weighted_mean_bias = sum(mean_bias, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
      pull(weighted_mean_bias)
    mean_bias_label <- bquote("Bias" == .(round(mean_bias, 0)) ~ "mm")
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
      legend.title = element_text(color = "black", size=15),
      plot.margin=unit(c(0,-0.3,0.3,-0.6), 'cm'), # unit(c(top, right, bottom, left), units)
      legend.key.width = unit(2.5, "cm"), # control size of colorbar
      legend.key.height = unit(0.5, "cm"),
      panel.border = element_rect(colour = "black", fill= NA, linewidth=0.5)) +
    guides(fill = guide_colourbar(frame.linewidth = 0.5, ticks.linewidth = 0.5, frame.colour = "black", ticks.colour = "black"),
           color = guide_colourbar(frame.linewidth = 0.5, ticks.linewidth = 0.5, frame.colour = "black", ticks.colour = "black")) +
    labs(title = model,
         color = expression(paste(Delta*SM[max]," (mm)")),
         fill = expression(paste(Delta*SM[max]," (mm)"))) # fill and color same label --> only one colorbar in the legend



  # add stats for model panels
  if(model != "GRACE observations") {

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

