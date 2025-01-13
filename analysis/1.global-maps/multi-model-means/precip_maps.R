# script to compare precip between historical and land-hist

# libraries
library(terra)
library(tidyverse)
library(lubridate)
library(data.table)
library(furrr)
library(future)
library(ggplot2)
library(ggpubr)
library(sf)
library(ggnewscale)
sf_use_s2(FALSE)

# Set model names
model_names <- c(
  "CESM2",
  "EC-Earth3-Veg",
  "IPSL-CM6A-LR",
  "UKESM1-0-LL",
  "CNRM-ESM2-1",
  "E3SM-1-1",
  "CMCC-ESM2",
  "CNRM-CM6-1",
  "MIROC6"
)

# Set scenario types
scenario_types <- c("land-hist", "historical")

# Create a list of model-scenario combinations
model_scenario_list <- expand.grid(model_name = model_names, scenario_type = scenario_types, stringsAsFactors = FALSE)

# Function to process precipitation data
process_model_precip <- function(model_name, scenario_type) {
  print(paste("Processing model:", model_name, "Scenario:", scenario_type))

  # Construct file pattern
  pr_pattern <- paste0("pr_mon_", model_name, "_", scenario_type, "_.*\\.nc")

  # List files matching the pattern
  pr_files <- list.files("data-raw/cmip6-ng/pr/mon/g025/", pattern = pr_pattern, full.names = TRUE)

  # Check and load P files
  if (model_name == "CMCC-ESM2") {
    # Special case for CMCC-ESM2: use precipitation from CESM2
    P <- rast("data-raw/cmip6-ng/pr/mon/g025/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
  } else {
    if (length(pr_files) == 0) {
      warning(paste("P files not found for model:", model_name, "Scenario:", scenario_type))
      return(NULL)
    }
    P <- rast(pr_files[1])
  }

  # Select dates from 1935-01-01 onwards
  dates <- terra::time(P)
  selection <- which(dates >= as.Date("2007-01-01")) # for CWD: 1935, for theta_crit: 2007
  P_selected <- subset(P, selection)

  # Rotate to true coordinates (-180 to 180)
  P_selected <- terra::rotate(P_selected)

  # Convert precipitation to mm/day
  P_selected <- P_selected * 86400  # Convert from kg/m^2/s to mm/day

  # Return raster with metadata
  return(P_selected)
}

# Process each model-scenario combination
results <- list()
for (i in seq_len(nrow(model_scenario_list))) {
  model <- model_scenario_list$model_name[i]
  scenario <- model_scenario_list$scenario_type[i]

  result <- process_model_precip(model, scenario)
  if (!is.null(result)) {
    results[[paste(model, scenario, sep = "_")]] <- result
  }
}

# Combine results by scenario to calculate multi-model mean
multi_model_means <- list()

for (scenario in scenario_types) {
  # Filter results by scenario based on the list names
  scenario_rasters <- lapply(names(results), function(key) {
    if (grepl(scenario, key)) {
      return(results[[key]])
    }
    return(NULL)
  })

  # Remove NULL entries
  scenario_rasters <- Filter(Negate(is.null), scenario_rasters)

  # Calculate multi-model mean if there are rasters
  if (length(scenario_rasters) > 0) {
    combined_raster <- rast(scenario_rasters)  # Combine rasters into a single SpatRaster
    multi_model_means[[scenario]] <- terra::mean(combined_raster)  # Calculate multi-model mean
  }
}

# PLOT MAPS ---------------------------------------------------------------

# Load required libraries
library(ggplot2)
library(ggpubr)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggnewscale)
sf_use_s2(FALSE)

# Download countries and ocean outlines for consistent map styling
countries <- ne_countries(scale = 50, returnclass = "sf")
ocean <- ne_download(
  scale = 50,
  type = "ocean",
  category = "physical",
  returnclass = "sf"
)

# Define plot breaks and limits for precipitation
breaks <- seq(0, 10, by = 1)  # Adjust based on expected range
max_limit <- 10  # Set a max value to cap the color scale

# Initialize an empty list for plots
plot_list <- list()

# Iterate over each scenario to create a plot
for (scenario in names(multi_model_means)) {

  # Convert raster to data frame for plotting
  df <- terra::as.data.frame(multi_model_means[[scenario]], xy = TRUE)
  colnames(df) <- c("lon", "lat", "precip_mm_day")

  # Set values greater than the maximum limit to the maximum
  df$precip_mm_day <- pmin(df$precip_mm_day, max_limit)

  # Create the plot
  p <- ggplot() +
    theme_void() +
    geom_tile(data = df, aes(x = lon, y = lat, fill = precip_mm_day)) +
    geom_sf(data = ocean,  # Ocean outlines
            color = "black",
            linetype = 'solid',
            fill = "white",
            size = 1) +
    geom_sf(data = countries,  # Country borders
            color = "black",
            fill = NA,
            size = 0.5) +
    scale_fill_viridis_c(
      option = "turbo",
      limits = c(0, max_limit),
      breaks = breaks,
      name = "Precip (mm/day)",
      na.value = "white"
    ) +
    coord_sf(
      xlim = c(-180, 180),
      ylim = c(-60, 88),
      expand = FALSE
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      plot.margin = unit(c(0, -0.3, 0.3, -0.6), 'cm'),
      legend.key.width = unit(4, "cm"),
      legend.key.height = unit(0.4, "cm")
    ) +
    guides(
      fill = guide_colourbar(
        frame.linewidth = 0.5,
        ticks.linewidth = 0.5,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    labs(
      title = paste("Multi-model mean precipitation (", scenario, ")", sep = ""),
      fill = "Precip (mm/day)"
    )

  # Add to plot list
  plot_list[[scenario]] <- p
}

# Combine the plots into a single figure
combined_plot <- ggarrange(
  plotlist = plot_list,
  ncol = 2,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom"
)

# Save the combined plot
ggsave("multi_model_precip_maps.png", plot = combined_plot, path = "./", width = 12, height = 3.25, dpi=300)


