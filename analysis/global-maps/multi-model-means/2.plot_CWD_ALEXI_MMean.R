# script to plot land water storage as quantified by CWD_max using cmip6 data and ALEXI (S_CWDX)
# multi-model mean

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

  # Add the model name as a new column in the dataframe
  df$model_name <- model_name

  # Append the dataframe to the list
  dataframes_list[[model_name]] <- df
}

# Combine all dataframes in the list into a single dataframe
summary_CWD <- bind_rows(dataframes_list) %>%
  mutate(model_name = ifelse(model_name == ".ALEXI", "Observations", model_name)) %>% # rename observations
  mutate(max_cwd = if_else(model_name != "Observations", max_cwd * 30, max_cwd)) # convert the units of CMIP6 CWD (number of seconds in a month, see script "calculate_CWD.R")

# Find the overall minimum and maximum values of "max_cwd" across all models (to adjust intervals plot)
min_max_values <- summary_CWD %>%
  reframe(min_max_cwd = range(max_cwd)) %>%
  pull(min_max_cwd)

# Define the desired breaks for the color scale (0 to 1000 by 100)
breaks <- seq(0, 1000, by = 100)

# Set values greater than 1000 to 1000 to ensure they are colored using the same color as the maximum break
summary_CWD$max_cwd[summary_CWD$max_cwd > 1000] <- 1000

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
  filter(model_name == "Observations") %>%
  dplyr::select(-model_name)

# calculate multi-model mean
MMmeans <- summary_CWD %>%
  filter(model_name != "Observations") %>%
  group_by(lon, lat) %>%
  summarise(
    max_cwd = mean(max_cwd, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model_name = "Multi-model mean") # add model name to mirror summary_CWD

# bind rows to summary_CWD
summary_CWD_MMmeans <- summary_CWD %>%
  bind_rows(MMmeans)

summary_merged <- summary_CWD_MMmeans %>%
  full_join( # attach ALEXI data as new columns
    alexi_data, by = join_by(lon, lat),
    suffix = c("", "_obs")) # rename new column "medianSMmax" + "_ALEXI"

# save dataframe
saveRDS(summary_merged, "summary_maxCWD_MMmean.rds", compress = "xz")

model_list <- c("Observations", "Multi-model mean")

# Create separate plots for each model
plot_list <- lapply(model_list, function(model) {

  df_model <- filter(summary_merged, model_name == model)

  # calculate stats to compare models-obs
  if(model != "Observations") {

    # calculate R squared
    fit <- lm(max_cwd ~ max_cwd_obs, data = df_model)
    r_squared <- summary(fit)$r.squared
    r_squared_label <- bquote(italic(R)^2 == .(round(r_squared, 2)))

    # calculate abs(mean bias)
    mean_bias_abs <- df_model %>%
      mutate(mean_bias_abs = abs(max_cwd - max_cwd_obs)) %>%
      pull(mean_bias_abs) %>%
      mean(na.rm = TRUE)
    mean_bias_abs_label <- bquote("|Bias|" == .(round(mean_bias_abs, 2)) ~ "mm")

    # calculate mean bias
    mean_bias <- df_model %>%
      mutate(mean_bias = max_cwd - max_cwd_obs) %>%
      pull(mean_bias) %>%
      mean(na.rm = TRUE)
    mean_bias_label <- bquote("Bias" == .(round(mean_bias, 2)) ~ "mm")
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
    labs(title = model, color = "Max CWD (mm)", fill = "Max CWD (mm)") # fill and color same label --> only one colorbar in the legend

  # add stats for model panels
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
                 labels = NULL,
                 ncol = 2, nrow = 1,
                 common.legend = TRUE, # have just one common legend
                 legend="bottom")

ggsave("map_CWDmax_MMmean.png", path = "./", width = 12, height = 3.25, dpi=300)

# save plot list for combined plot
saveRDS(plot_list, "plot_list_CWD.rds", compress = "xz")









