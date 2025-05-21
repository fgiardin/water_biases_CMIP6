# script to plot land water storage as quantified by CWD_max using cmip6 data and ALEXI (S_CWDX)
# results by model

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
directory_path <- "data/CWD/max_CWD"

# get list of files we want to merge
file_locations <- list.files(directory_path, # define directory where all files are
                             glob2rx("max_cwd_*.rds"), # search for this pattern inside directory
                             full.names = TRUE, # list all paths
                             recursive = FALSE # NOT go through all subdirectories
                             )

# focus only on land-hist scenario and ALEXI observations
file_locations <- file_locations[grepl("land-hist|\\.ALEXI", file_locations)]

# Initialize an empty list to store dataframes
dataframes_list <- list()

# Loop through the files, read each .rds file, and store the dataframes in the list
for (file_path in file_locations) {

  # Read the dataframe from the .rds file
  df <- readRDS(file_path)

  if (grepl("ALEXI", file_path)) { # If it's the ALEXI file, create the column with ".ALEXI"
    df$model_name <- ".ALEXI"

  } else { # Else add the model name as a new column in the dataframe
    df <- df %>%
      dplyr::rename(model_name = model) %>%
      dplyr::select(-scenario)
  }

  # Append the dataframe to the list
  dataframes_list[[df$model_name[1]]] <- df
}

# Combine all dataframes in the list into a single dataframe
summary_CWD <- bind_rows(dataframes_list) %>%
  mutate(max_cwd = if_else(model_name != ".ALEXI", max_cwd * 30, max_cwd)) # adjust the units of CMIP6 CWD (leave ALEXI unchanged)

# Find the overall minimum and maximum values of "max_cwd" across all models
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

# rename ".ALEXI"
summary_CWD <- summary_CWD %>%
  mutate(model_name = ifelse(model_name == ".ALEXI", "ALEXI observations", model_name))
model_list <- unique(summary_CWD$model_name)

# Create a new (repeated for every model) column with ALEXI data in the original dataframe
# only to calculate model-obs stats (not for plotting!)
alexi_data <- summary_CWD %>%
  filter(model_name == "ALEXI observations") %>%
  dplyr::select(-model_name)

summary_merged <- full_join(summary_CWD,
                            alexi_data, by = join_by(lon, lat),
                            suffix = c("", "_ALEXI")) # rename new column "medianSMmax" + "_GRACE"

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
  if(model != "ALEXI observations") {

    # calculate R squared
    fit <- lm(max_cwd ~ max_cwd_ALEXI, data = df_model)
    r_squared <- summary(fit)$r.squared
    r_squared_label <- bquote(italic(R)^2 == .(round(r_squared, 2)))

    # calculate abs(mean bias) using weights
    mean_bias_abs <- df_model %>%
      filter(!is.na(max_cwd) & !is.na(max_cwd_ALEXI) & !is.na(weights)) %>% # manually remove NAs when using "sum" (num and den should have same number of elements)
      mutate(mean_bias_abs = abs(max_cwd - max_cwd_ALEXI) * weights) %>%
      summarise(weighted_mean_bias_abs = sum(mean_bias_abs, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
      pull(weighted_mean_bias_abs)
    mean_bias_abs_label <- bquote("|Bias|" == .(round(mean_bias_abs, 0)) ~ "mm")

    # calculate mean bias using weights
    mean_bias <- df_model %>%
      filter(!is.na(max_cwd) & !is.na(max_cwd_ALEXI) & !is.na(weights)) %>%
      mutate(mean_bias = (max_cwd - max_cwd_ALEXI) * weights) %>%
      summarise(weighted_mean_bias = sum(mean_bias, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
      pull(weighted_mean_bias)
    mean_bias_label <- bquote("Bias" == .(round(mean_bias, 0)) ~ "mm")
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
      legend.key.width = unit(2.5, "cm"), # control size of colorbar
      legend.key.height = unit(0.5, "cm"),
      panel.border = element_rect(colour = "black", fill= NA, linewidth=0.5)
    ) +
    guides(fill = guide_colourbar(frame.linewidth = 0.5, ticks.linewidth = 0.5, frame.colour = "black", ticks.colour = "black"),
           color = guide_colourbar(frame.linewidth = 0.5, ticks.linewidth = 0.5, frame.colour = "black", ticks.colour = "black")
    ) +
    labs(title = ifelse(model == "ALEXI observations",
                        "ALEXI and WATCH-WFDEI observations",
                        model),
         color = "Maximum water storage (mm)",
         fill  = "Maximum water storage (mm)")

  # add stats for model panels
  if(model != "ALEXI observations") {

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
                 ncol = 2, nrow = 5,
                 common.legend = TRUE, # have just one common legend
                 legend="bottom")

ggsave("map_CWDmax_9-models.png", path = "./", width = 12, height = 13.75, dpi=600)


# saveRDS(summary_merged, "summary_max_cwd.rds", compress = "xz")







