# plot GRACE/GLDAS/GLWS vs multi-model mean LMIP

# Load required packages
library(tidyverse)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggpubr)
library(terra)
library(rgdal)
library(ggnewscale)
library(viridis)
sf_use_s2(FALSE)


# Load and prepare data ---------------------------------------------------

# --- User-defined parameters ---
plot_variable <- "deltaSMmax"  # or "deltaSMabs" if desired
upper_threshold <- 1000
interval <- upper_threshold / 9
rounded_interval <- round(interval / 100) * 100
breaks <- seq(0, upper_threshold, by = rounded_interval)

# load GLDAS and GLWS dSMmax
GLDAS <- readRDS("data/Comparisons_GRACE/GLDAS_2.5deg.rds")
GLWS <- readRDS("data/Comparisons_GRACE/GWLS_2.5deg.rds")

# Define the directory where the .rds files are located
directory_path <- "data/deltaSM/"

# List all .rds files in the directory
file_locations <- list.files(directory_path,
                             pattern = "deltaSM_.*\\.rds",
                             full.names = TRUE)

# Initialize an empty list for land-hist data and a variable for GRACE data
dataframes_land_hist <- list()
grace_data <- NULL

# Loop through files and load data
for (file_path in file_locations) {
  filename <- basename(file_path)

  if (filename == "deltaSM_.GRACE.rds") {
    # Read GRACE data and label it
    grace_data <- readRDS(file_path)
    grace_data$model_name <- "GRACE observations"
  } else {
    # Expect filenames of the form "deltaSM_<model>_<scenario>.rds"
    matches <- regmatches(filename, regexec("deltaSM_(.*)_(.*)\\.rds", filename))
    if (length(matches[[1]]) == 3) {
      model_name <- matches[[1]][2]
      scenario <- matches[[1]][3]

      if (scenario == "land-hist") {
        df <- readRDS(file_path)
        df$model_name <- model_name
        dataframes_land_hist[[model_name]] <- df
      }
    } else {
      warning(paste("Filename does not match expected pattern:", filename))
    }
  }
}

# Check if land-hist data exists and combine into one dataframe
if (length(dataframes_land_hist) == 0) {
  stop("No land-hist data found.")
}
summary_deltaSM_land_hist <- bind_rows(dataframes_land_hist)

# Cap values above the upper threshold for the plot variable
summary_deltaSM_land_hist[[plot_variable]][summary_deltaSM_land_hist[[plot_variable]] > upper_threshold] <- upper_threshold
grace_data[[plot_variable]][grace_data[[plot_variable]] > upper_threshold] <- upper_threshold

# Combine model data with GRACE observations
summary_deltaSM_land_hist <- bind_rows(summary_deltaSM_land_hist, grace_data)

# Calculate the multi-model mean (excluding GRACE observations)
MMmeans <- summary_deltaSM_land_hist %>%
  filter(model_name != "GRACE observations") %>%
  group_by(lon, lat) %>%
  summarise(
    deltaSMmax = mean(deltaSMmax, na.rm = TRUE),
    deltaSMabs = mean(deltaSMabs, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model_name = "Multi-model mean")

# Combine the multi-model mean with the rest of the data
summary_deltaSM_MMmeans <- bind_rows(summary_deltaSM_land_hist, MMmeans)

# add GLWS and GLDAS data!
GLDAS_modified <- GLDAS %>% # Process GLDAS: rename, add empty column and model_name
  rename(deltaSMmax = dsmldas_max) %>%
  mutate(deltaSMabs = NA,      # add empty column (NA values)
         model_name = "GLDAS_CLSM025_DA1_D")

GLWS_modified <- GLWS %>% # Process GLWS: rename, add empty column and model_name
  rename(deltaSMmax = dsmgw_max) %>%
  mutate(deltaSMabs = NA,      # add empty column (NA values)
         model_name = "GLWS2.0")

summary_deltaSM_MMmeans <- bind_rows(summary_deltaSM_MMmeans, # Bind the modified GLDAS and GLWS rows to the existing summary_deltaSM_MMmeans
                                     GLDAS_modified,
                                     GLWS_modified)

# Merge with GRACE data (by lon and lat) to allow computation of statistics
summary_merged <- summary_deltaSM_MMmeans %>%
  full_join(grace_data,
            by = join_by(lon, lat),
            suffix = c("", "_GRACE"))

# Calculate grid weights for the bias computation (using a 2.5° resolution)
res <- 2.5
ref_grid_size <- sin(res * pi / 180)
summary_merged <- summary_merged %>%
  mutate(
    lat1 = lat - (res/2),
    lat2 = lat + (res/2),
    gridsize = abs(sin(lat1 * pi / 180) - sin(lat2 * pi / 180)),
    weights = gridsize / ref_grid_size
  ) %>%
  dplyr::select(-lat1, -lat2, -gridsize)


# PLOTS -------------------------------------------------------------------

# --- Download or load spatial boundary data ---
ocean <- ne_download(scale = 50,
                     type = "ocean",
                     category = "physical",
                     returnclass = "sf")

# --- Define left models to compare with the multi-model mean ---
left_obs <- c("GRACE observations", "GLDAS_CLSM025_DA1_D", "GLWS2.0")
reference_model <- "Multi-model mean"

# --- Initialize a list to store the 6 plots ---
plot_list <- list()

# --- Loop through each left model to create paired panels ---
for(model in left_obs) {

  # Extract data for the current model and for the multi-model mean
  df_left <- summary_merged %>% filter(model_name == model)
  df_mm   <- summary_merged %>% filter(model_name == reference_model)

  # For the statistics we join the two datasets by lon abd lat
  # The left dataset is the one to be compared; the multi-model mean gets suffix _mm
  df_stats <- left_join(df_left, df_mm, by = c("lon", "lat"),
                        suffix = c("", "_mm"))

  # --- Calculate statistics ---
  # we want the bias computed as bias = (model - multi-model mean)
    fit <- lm(deltaSMmax_mm ~ deltaSMmax, data = df_stats)
    r_squared <- summary(fit)$r.squared
    r_squared_label <- bquote(italic(R)^2 == .(round(r_squared, 2)))

    df_stats_bias <- df_stats
    # df_stats_bias <- df_stats %>%
    #   filter(!is.na(deltaSMmax), !is.na(deltaSMmax_mm), !is.na(weights))

    mean_bias_abs <- df_stats_bias %>%
      filter(!is.na(deltaSMmax_mm) & !is.na(deltaSMmax) & !is.na(weights)) %>% # manually remove NAs
      mutate(bias_abs = abs(deltaSMmax_mm - deltaSMmax) * weights) %>%
      summarise(weighted_mean_bias_abs = sum(bias_abs, na.rm = TRUE) / sum(weights)) %>%
      pull(weighted_mean_bias_abs)
    mean_bias_abs_label <- bquote("|Bias|" == .(round(mean_bias_abs, 0)) ~ "mm")

    mean_bias <- df_stats_bias %>%
      filter(!is.na(deltaSMmax_mm) & !is.na(deltaSMmax) & !is.na(weights)) %>% # manually remove NAs
      mutate(bias = (deltaSMmax_mm - deltaSMmax) * weights) %>%
      summarise(weighted_mean_bias = sum(bias, na.rm = TRUE) / sum(weights)) %>%
      pull(weighted_mean_bias)
    mean_bias_label <- bquote(Bias == .(round(mean_bias, 0)) ~ "mm")

    mean_bias_tropics <- df_stats_bias %>%
      filter(lat >= -23, lat <= 23) %>%
      filter(!is.na(deltaSMmax_mm) & !is.na(deltaSMmax) & !is.na(weights)) %>% # manually remove NAs
      mutate(bias = (deltaSMmax_mm - deltaSMmax) * weights) %>%
      summarise(weighted_mean_bias = sum(bias, na.rm = TRUE) / sum(weights)) %>%
      pull(weighted_mean_bias)
    mean_bias_tropics_label <- bquote(Bias[tropics] == .(round(mean_bias_tropics, 0)) ~ "mm")

  # --- Create left panel: the individual dataset (without stats annotations) ---
  p_left <- ggplot() +
    theme_void() +
    geom_tile(data = df_left,
              aes(x = lon, y = lat,
                  fill = !!sym(plot_variable),
                  color = !!sym(plot_variable))) +
    geom_sf(data = ocean,
            color = "black",
            linetype = "solid",
            fill = "white",
            size = 1) +
    scale_color_viridis_c(option = "turbo",
                          breaks = breaks,
                          limits = c(0, upper_threshold)) +
    scale_fill_viridis_c(option = "turbo",
                         breaks = breaks,
                         limits = c(0, upper_threshold)) +
    coord_sf(xlim = c(-180, 180), ylim = c(-60, 88), expand = FALSE) +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          legend.text = element_text(color = "black", size = 12),
          legend.title = element_text(color = "black", size = 12),
          plot.margin = unit(c(0, -0.3, 0.3, -0.6), "cm"),
          legend.key.width = unit(2.5, "cm"),
          legend.key.height = unit(0.4, "cm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
    guides(fill = guide_colourbar(frame.linewidth = 0.5,
                                  ticks.linewidth = 0.5,
                                  frame.colour = "black",
                                  ticks.colour = "black"),
           color = guide_colourbar(frame.linewidth = 0.5,
                                   ticks.linewidth = 0.5,
                                   frame.colour = "black",
                                   ticks.colour = "black")) +
    labs(title = model,
         color = expression(Delta*SM[max]~"(mm)"),
         fill = expression(Delta*SM[max]~"(mm)"))

  # --- Create right panel: the multi-model mean reference with stats annotations ---
  p_mm <- ggplot() +
    theme_void() +
    geom_tile(data = df_mm,
              aes(x = lon, y = lat,
                  fill = !!sym(plot_variable),
                  color = !!sym(plot_variable))) +
    geom_sf(data = ocean,
            color = "black",
            linetype = "solid",
            fill = "white",
            size = 1) +
    scale_color_viridis_c(option = "turbo",
                          breaks = breaks,
                          limits = c(0, upper_threshold)) +
    scale_fill_viridis_c(option = "turbo",
                         breaks = breaks,
                         limits = c(0, upper_threshold)) +
    coord_sf(xlim = c(-180, 180), ylim = c(-60, 88), expand = FALSE) +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          legend.text = element_text(color = "black", size = 12),
          legend.title = element_text(color = "black", size = 12),
          plot.margin = unit(c(0, -0.3, 0.3, -0.6), "cm"),
          legend.key.width = unit(2.5, "cm"),
          legend.key.height = unit(0.4, "cm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)) +
    guides(fill = guide_colourbar(frame.linewidth = 0.5,
                                  ticks.linewidth = 0.5,
                                  frame.colour = "black",
                                  ticks.colour = "black"),
           color = guide_colourbar(frame.linewidth = 0.5,
                                   ticks.linewidth = 0.5,
                                   frame.colour = "black",
                                   ticks.colour = "black")) +
    labs(title = reference_model,
         color = expression(Delta*SM[max]~"(mm)"),
         fill = expression(Delta*SM[max]~"(mm)")) +
    annotate("text", x = -175, y = -11, label = r_squared_label,
             hjust = 0, vjust = 0, size = 4.3) +
    annotate("text", x = -175, y = -25, label = mean_bias_label,
             hjust = 0, vjust = 0, size = 4.3) +
    annotate("text", x = -175, y = -40, label = mean_bias_tropics_label,
             hjust = 0, vjust = 0, size = 4.3) +
    annotate("text", x = -175, y = -52, label = mean_bias_abs_label,
             hjust = 0, vjust = 0, size = 4.3)

  # Save both panels in the list. Names include the left model and a suffix for left vs. reference.
  plot_list[[paste0(model, "_left")]] <- p_left
  plot_list[[paste0(model, "_mm")]]   <- p_mm
}

# --- Arrange the six plots (3 rows x 2 columns) with a common legend ---
combined_plot <- ggarrange(plotlist = plot_list,
                           labels = "auto",
                           ncol = 2, nrow = 3,
                           common.legend = TRUE,
                           legend = "bottom")

# --- Save the combined plot ---
ggsave("map_GRACE_comparisons.png", plot = combined_plot,
       width = 12, height = 8.25, dpi = 600)



