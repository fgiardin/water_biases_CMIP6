# script to plot land water storage as quantified by delta_SM using cmip6 data and GRACE
# multi-model mean
# update to also plot comparison with historical data in a separate plot

# Load required packages
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

# source function to plot and print data
source("R/process_scenario_data.R")

# Define which variable to plot: deltaSMabs or deltaSMmax
# deltaSMabs: absolute deltaSM over 1900-2014
# deltaSMmax: first calculate deltaSM per year and then take maximum among all years
plot_variable <- "deltaSMmax" # (!!! changing to deltaSMabs will still plot the stats of deltaSMmax, need to change manually!)

# Define the upper threshold and breaks (common to both scenarios, depends on GRACE data)
upper_threshold <- 1000

# Define the desired breaks for the color scale (around 10 breaks given the upper_threshold)
interval <- upper_threshold / 9
rounded_interval <- round(interval / 100) * 100
breaks <- seq(0, upper_threshold, by = rounded_interval)

# Download countries and ocean data (do this once)
countries <- ne_countries(scale = 50, returnclass = c("sf"))

# Download ocean outlines
ocean <- ne_download(
  scale = 50,
  type = "ocean",
  category = "physical",
  returnclass = "sf")

# Define the directory path where the .rds files are located
directory_path <- "data/deltaSM/"

# Get list of files we want to merge (here: output dataframes from script 2)
file_locations <- list.files(directory_path, # define directory where all files are
                             pattern = "deltaSM_.*\\.rds", # search for this pattern inside directory
                             full.names = TRUE, # list all paths
                             recursive = FALSE # NOT go through all subdirectories
)

# Initialize empty lists to store dataframes for each scenario
dataframes_land_hist <- list()
dataframes_historical <- list()
grace_data <- NULL

# Loop through the files, read each .rds file, and store the dataframes in the appropriate list
for (file_path in file_locations) {
  # Extract the filename without path
  filename <- basename(file_path)

  # Check if the file is the GRACE data
  if (filename == "deltaSM_.GRACE.rds") {

    # Read the GRACE data
    grace_data <- readRDS(file_path)
    grace_data$model_name <- "GRACE observations"

  } else {

    # Extract the model name and scenario from the file name
    # The filename format is 'deltaSM_<model_name>_<scenario>.rds'
    # So we can extract the parts

    matches <- regmatches(filename, regexec("deltaSM_(.*)_(.*)\\.rds", filename))
    if (length(matches[[1]]) == 3) {
      model_name <- matches[[1]][2]
      scenario <- matches[[1]][3]

      # Read the dataframe from the .rds file
      df <- readRDS(file_path)

      # Add the model name as a new column in the dataframe
      df$model_name <- model_name

      # Store the dataframe in the appropriate list
      if (scenario == "land-hist") {
        dataframes_land_hist[[model_name]] <- df
      } else if (scenario == "historical") {
        dataframes_historical[[model_name]] <- df
      } else {
        warning(paste("Unknown scenario in file:", filename))
      }
    } else {
      warning(paste("Filename does not match expected pattern:", filename))
    }
  }
}

# Process and print 'land-hist' data
if (length(dataframes_land_hist) > 0) {
  # Combine all dataframes into a single dataframe
  summary_deltaSM_land_hist <- bind_rows(dataframes_land_hist)

  # Call the processing function
  process_scenario_data(summary_deltaSM_land_hist, grace_data, "land-hist")
}

# Process and print 'historical' data
if (length(dataframes_historical) > 0) {
  # Combine all dataframes into a single dataframe
  summary_deltaSM_historical <- bind_rows(dataframes_historical)

  # Call the processing function
  process_scenario_data(summary_deltaSM_historical, grace_data, "historical")
}


