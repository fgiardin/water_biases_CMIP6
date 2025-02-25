# script to plot land water storage as quantified by CWD_max using cmip6 data and ALEXI (S_CWDX)
# multi-model mean
# update to also plot comparison with historical data in a separate plot (goes to SI)

# load packages
devtools::load_all(".")
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
sf_use_s2(FALSE)


# Define the directory path where the .rds files are located
directory_path <- "data/CWD/max_CWD/"

# get list of files we want to merge (here: output dataframes from script 2)
file_locations <- list.files(directory_path, # define directory where all files are
                             glob2rx("max_cwd_*.rds"), # search for this pattern inside directory
                             full.names = TRUE, # list all paths
                             recursive = FALSE # NOT go through all subdirectories
)

# Initialize an empty list to store dataframes
dataframes_list <- list()

# Loop through the files, read each .rds file, and store the dataframes in the list
for (file_path in file_locations) {

  # Extract the model name from the file name
  model_name <- gsub(".*max_cwd_(.*)\\.rds", "\\1", basename(file_path))

  # Read the dataframe from the .rds file
  df <- readRDS(file_path)

  if(model_name == ".ALEXI") {
    # Add a new column to observational data
    df$model <- model_name
  }

  # Append the dataframe to the list
  dataframes_list[[model_name]] <- df
}

# Combine all dataframes in the list into a single dataframe
summary_CWD <- bind_rows(dataframes_list) %>%
  mutate(model = ifelse(model == ".ALEXI", "Observations", model)) %>% # rename observations
  mutate(max_cwd = if_else(model != "Observations", max_cwd * 30, max_cwd)) # convert the units of CMIP6 CWD (number of seconds in a month, see script "calculate_CWD.R")

# Find the overall minimum and maximum values of "max_cwd" across all models (to adjust intervals in plot)
min_max_values <- summary_CWD %>%
  reframe(min_max_cwd = range(max_cwd)) %>%
  pull(min_max_cwd)

# Define the desired breaks for the color scale (0 to 1000 by 100)
breaks <- seq(0, 1000, by = 100)

# download countries
countries <- ne_countries(scale = 50, returnclass = c("sf"))

# download ocean outlines
ocean <- ne_download(
  scale = 50,
  type = "ocean",
  category = "physical",
  returnclass = "sf")

# Create a new (repeated for every model) column with ALEXI data in the original dataframe
# To calculate comparison stats model-ALEXI
alexi_data <- summary_CWD %>%
  filter(model == "Observations") %>%
  dplyr::select(-model,-scenario)

# summary_CWD_hist <- summary_CWD %>%
#   dplyr::filter(scenario == "historical")
# summary_CWD_landhist <- summary_CWD %>%
#   dplyr::filter(scenario == "land-hist")
# quantile(summary_CWD_hist$max_cwd, 0.99) # 34940.76
# quantile(summary_CWD_landhist$max_cwd, 0.99) # 765.6059

# handle outliers
summary_CWD$max_cwd[summary_CWD$max_cwd > 1000] <- 1000
summary_CWD <- summary_CWD %>%
  mutate(max_cwd = if_else(max_cwd == 1000 & scenario == "historical", NA_real_, max_cwd))

# calculate multi-model mean
MMmeans <- summary_CWD %>%
  filter(model != "Observations") %>%
  group_by(lon, lat, scenario) %>%
  summarise(
    max_cwd = mean(max_cwd, na.rm = TRUE)
  ) %>%

  ungroup() %>%
  mutate(model = "Multi-model mean") # add model name to mirror summary_CWD

# bind rows to summary_CWD
summary_CWD_MMmeans <- summary_CWD %>%
  bind_rows(MMmeans)

summary_merged <- summary_CWD_MMmeans %>%
  full_join( # attach ALEXI data as new columns
    alexi_data, by = join_by(lon, lat),
    suffix = c("", "_obs")) # rename new column "medianSMmax" + "_ALEXI"

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



# save dataframe
saveRDS(summary_merged, "summary_maxCWD_MMmean.rds", compress = "xz")



# Get the list of unique scenarios
scenarios <- unique(summary_merged$scenario)
scenarios <- scenarios[!is.na(scenarios)]

model_list <- c("Observations", "Multi-model mean")



# print figure -------------------------------------------------------

# Loop over each scenario
for (scene in scenarios) {

  print(paste("Processing scenario:", scene))

  # Filter data for the current scenario
  first_filter <- dplyr::filter(summary_merged, scenario == scene)

  # Create separate plots for data in model_list
  plot_list <- lapply(model_list, function(current_mdl) {

    if(current_mdl == "Observations") {
      df_model <- dplyr::filter(summary_merged, model == current_mdl)
    } else if (current_mdl != "Observations") {
      # calculate stats to compare models-obs
      df_model <- dplyr::filter(first_filter, model == current_mdl)

      # calculate R squared
      fit <- lm(max_cwd ~ max_cwd_obs, data = df_model)
      r_squared <- summary(fit)$r.squared
      r_squared_label <- bquote(italic(R)^2 == .(round(r_squared, 2)))

      # calculate abs(mean bias) using weights
      mean_bias_abs <- df_model %>%
        dplyr::filter(!is.na(max_cwd) & !is.na(max_cwd_obs) & !is.na(weights)) %>% # manually remove NAs when using "sum" (num and den should have same number of elements)
        mutate(mean_bias_abs = abs(max_cwd - max_cwd_obs) * weights) %>%
        summarise(weighted_mean_bias_abs = sum(mean_bias_abs, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
        pull(weighted_mean_bias_abs)
      mean_bias_abs_label <- bquote("|Bias|" == .(round(mean_bias_abs, 0)) ~ "mm")

      # calculate mean bias using weights
      mean_bias <- df_model %>%
        # dplyr::filter(!(lat >= -23 & lat <= 23)) %>%
        dplyr::filter(!is.na(max_cwd) & !is.na(max_cwd_obs) & !is.na(weights)) %>%
        mutate(mean_bias = (max_cwd - max_cwd_obs) * weights) %>%
        summarise(weighted_mean_bias = sum(mean_bias, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
        pull(weighted_mean_bias)
      # mean_bias_label <- bquote(Bias[extratropics] == .(round(mean_bias, 0)) ~ "mm")
      mean_bias_label <- bquote(Bias == .(round(mean_bias, 0)) ~ "mm")


      # bias focusing on tropics
      mean_bias_tropics <- df_model %>%
        dplyr::filter(lat >= -23 & lat <= 23) %>%
        dplyr::filter(!is.na(max_cwd) & !is.na(max_cwd_obs) & !is.na(weights)) %>%
        mutate(mean_bias = (max_cwd - max_cwd_obs) * weights) %>%
        summarise(weighted_mean_bias = sum(mean_bias, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
        pull(weighted_mean_bias)
      mean_bias_tropics_label <- bquote(Bias[tropics] == .(round(mean_bias_tropics, 0)) ~ "mm")

    }

    # Plot global map for each model
    p <- ggplot() +
      theme_void() +
      geom_tile(data = df_model,
                aes(x = lon, y = lat, color = max_cwd, fill = max_cwd)) +
      geom_sf(data=ocean, # country borders
              color= "black",
              linetype='solid',
              fill= "white", #cad6df", #D6D6E9
              size=1) +
      scale_color_viridis_c(option = "turbo", limits = c(0, 1000), breaks = breaks) + #
      scale_fill_viridis_c(option = "turbo", limits = c(0, 1000), breaks = breaks) + #
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
        legend.key.width = unit(4, "cm"),
        legend.key.height = unit(0.4, "cm"),
        panel.border = element_rect(colour = "black", fill= NA, linewidth=0.5)
      ) +
      guides(fill = guide_colourbar(frame.linewidth = 0.5, ticks.linewidth = 0.5, frame.colour = "black", ticks.colour = "black"),
             color = guide_colourbar(frame.linewidth = 0.5, ticks.linewidth = 0.5, frame.colour = "black", ticks.colour = "black")
      ) +
      labs(title = current_mdl, color = "Max CWD (mm)", fill = "Max CWD (mm)") # fill and color same label --> only one colorbar in the legend

    # add stats for model panels
    if(current_mdl != "Observations") {

      p <- p +
        annotate("text", x = -175, y = -11, label = r_squared_label, hjust = 0, vjust = 0, size = 4.3) + # R2
        annotate("text", x = -175, y = -25, label = mean_bias_label, hjust = 0, vjust = 0, size = 4.3) + # mean bias
        annotate("text", x = -175, y = -40, label = mean_bias_tropics_label, hjust = 0, vjust = 0, size = 4.3) + # tropics
        annotate("text", x = -175, y = -52, label = mean_bias_abs_label, hjust = 0, vjust = 0, size = 4.3) # abs(mean bias)

    }

    return(p)

  })

  # Print all plots with one colorbar
  all <- ggarrange(plotlist = plot_list,
                   labels = "auto",
                   ncol = 2, nrow = 1,
                   common.legend = TRUE, # have just one common legend
                   legend="bottom")

  # Define the filename based on the scenario
  filename <- paste0("map_CWDmax_", scene, "_MMmean.png")
  ggsave(filename, plot = all, path = "./", width = 12, height = 3.25, dpi = 300)
}

# save plot list for combined plot
saveRDS(plot_list, "plot_list_CWD.rds", compress = "xz")
