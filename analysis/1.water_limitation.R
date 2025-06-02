# script to plot frequency of SM limitation (multi-model mean) [Fig. 1]
# update to also plot comparison with historical data in a separate plot (supplementary figures)

# data processing in following scripts:
# GRACE/SIF data: data-raw/4.calculate_thetacrit_SIF.R
# LMIP-CMIP6 data: data-raw/3.calculate_thetacrit_monthly.R

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
library(matrixStats) # to calculate median of bias
library(RColorBrewer)
setwd("/Users/fgiardina/water_biases_CMIP6")

# *** select scenario ***
# !!! important do one scenario at a time!
scenario_select = "land-hist"
# either "land-hist" or "historical"

# load data
df_count_GRACE <- readRDS("data/theta_crit/monthly/df_count_GRACE.rds") %>%
  dplyr::select(-date, -SIF, -TWS) %>%  # remove temporal values
  unique() %>%
  mutate(model = "Observations")

df_count_mrso <- readRDS("data/theta_crit/monthly/df_count_mrso_allscenarios.rds") %>%
  dplyr::select(-date, -EF, -Rn, -mrso, -mrso_norm) %>%
  dplyr::filter(scenario == scenario_select) %>% # select scenario
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


# plot --------------------------------------------------------------------

plot_variable <- "count"

# Define breaks for the color scale (around 10 breaks given the upper_threshold)
lower_threshold <- 0
upper_threshold <- 100
# interval <- upper_threshold / 9
# rounded_interval <- round(interval, digits = 2)
breaks <- seq(lower_threshold, upper_threshold, by = 10) # by = rounded_interval

# create color palette
colors_rdpu <- brewer.pal(9, "PuRd") # Get colors from the RdPu palette
colors_rdpu <- colors_rdpu[2:9] # remove first color(s) as too close to white and last color as too dark in PDF format
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
  dplyr::filter(model_name == "Observations") %>%
  dplyr::select(lon, lat, count)

# calculate multi-model mean
MMmeans <- dt_count %>%
  dplyr::filter(model_name != "Observations") %>%
  group_by(lon, lat, scenario) %>%
  summarise(
    count = mean(count, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model_name = "Multi-model mean") # add model name to mirror summary_deltaSM

# bind rows to original dt (to add GRACE for plotting panel)
dt_count_MMmeans <- dt_count %>%
  bind_rows(MMmeans)

summary_merged <- full_join(dt_count_MMmeans, # remove GRACE in the original dataframe to avoid repetitions
                            grace_data, by = join_by(lon, lat),
                            suffix = c("", "_GRACE")) # add GRACE column just to calculate comparison stats below

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

saveRDS(summary_merged, "summary_merged_watlim.rds", compress = "xz")


# Get the list of unique scenarios
scenarios <- unique(summary_merged$scenario)
scenarios <- scenarios[!is.na(scenarios)]  # Remove NA if any

# Create the list of panels to print
model_list <- c("Observations", "Multi-model mean")

# function for weighted median
weighted_median <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    ok <- !(is.na(x) | is.na(w))
    x <- x[ok]; w <- w[ok]
  }
  if (!length(x)) return(NA_real_)

  o <- order(x)             # 1. sort the sample
  x <- x[o]; w <- w[o]

  cw <- cumsum(w) / sum(w)  # 2. cumulative weight ∈ [0,1]
  idx <- which(cw >= 0.5)[1]# 3. first cell past 50 %
  x[idx]                     # 4. that value *is* the weighted median
}

# list to collect the fraction of area affected by over/underestimation
area_fraction_stats <- list()   # one element per scenario XXXX

# Loop over each scenario
for (scene in scenarios) {

  print(paste("Processing scenario:", scene))

  # Filter data for the current scenario
  first_filter <- dplyr::filter(summary_merged, scenario == scene)

  # Create separate plots for "Observations" and "Multi-model mean"
  plot_list <- lapply(model_list, function(current_mdl) {

    # Handle "Observations" separately (since it may not have a scenario)
    if (current_mdl == "Observations") {
      df_model <- dplyr::filter(summary_merged, model_name == current_mdl)
    } else {
      df_model <- dplyr::filter(first_filter, model_name == current_mdl)
    }

    # Calculate statistics to compare models with observations
    if (current_mdl != "Observations") {

      # Calculate R-squared
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

      # calculate mean bias (raw bias) using weights
      mean_bias_int <- df_model %>% # save intermediary raw bias to calculate percentages below
        filter(!is.na(count) & !is.na(count_GRACE) & !is.na(weights)) %>%
        mutate(mean_bias = (count - count_GRACE) * weights)

      mean_bias <- mean_bias_int %>%
        summarise(weighted_mean_bias = sum(mean_bias, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
        pull(weighted_mean_bias)
      mean_bias_label <- bquote("Bias" == .(round(mean_bias, 0)) * "%")

      # calculate %land affected by over and underestimation
      raw_bias <- mean_bias_int$mean_bias
      n_tot      <- length(raw_bias)
      n_over     <- sum(raw_bias >  0, na.rm = TRUE)
      n_under    <- sum(raw_bias <  0, na.rm = TRUE)

      frac_df <- tibble(
        scenario   = scene,
        frac_over  = paste0(round(n_over  / n_tot * 100), "%"),
        frac_under = paste0(round(n_under / n_tot * 100), "%")
      )
      print(frac_df)


      # bias focusing on tropics
      mean_bias_tropics_int <- df_model %>%
        filter(lat >= -23 & lat <= 23) %>%
        filter(!is.na(count) & !is.na(count_GRACE) & !is.na(weights)) %>%
        mutate(mean_bias = (count - count_GRACE) * weights)

      mean_bias_tropics <- mean_bias_tropics_int %>%
        summarise(weighted_mean_bias = sum(mean_bias, na.rm = TRUE) / sum(weights, na.rm = TRUE)) %>%
        pull(weighted_mean_bias)
      mean_bias_tropics_label <- bquote(Bias[tropics] == .(round(mean_bias_tropics, 0)) * "%")

      # calculate %land affected by over and underestimation **IN TROPICS**
      raw_bias_tr <- mean_bias_tropics_int$mean_bias
      n_tot_tr      <- length(raw_bias_tr)
      n_over_tr    <- sum(raw_bias_tr >  0, na.rm = TRUE)
      n_under_tr    <- sum(raw_bias_tr <  0, na.rm = TRUE)

      frac_df_tropics <- tibble(
        scenario   = scene,
        frac_over_tropics  = paste0(round(n_over_tr  / n_tot_tr * 100), "%"),
        frac_under_tropics = paste0(round(n_under_tr / n_tot_tr * 100), "%")
      )
      print(frac_df_tropics)

     # median using weights
      median_bias <- df_model %>%
        filter(!is.na(count) & !is.na(count_GRACE) & !is.na(weights)) %>%
        transmute(bias = count - count_GRACE, w = weights) %>%
        summarise(weighted_median_bias = weighted_median(bias, w)) %>%
        # mutate(bias = count - count_GRACE) %>%
        # summarise(weighted_median_bias = weightedMedian(bias, w = weights, na.rm = TRUE)) %>%
        pull(weighted_median_bias)
      median_bias_label <- bquote("Median bias" == .(round(median_bias, 0)) * "%")
      print(median_bias_label) # just print (not to overpopulate figure)

      # median using weights: tropics
      median_bias_tropics <- df_model %>%
        filter(lat >= -23 & lat <= 23) %>%
        filter(!is.na(count) & !is.na(count_GRACE) & !is.na(weights)) %>%
        transmute(bias = count - count_GRACE, w = weights) %>%
        summarise(weighted_median_bias = weighted_median(bias, w)) %>%
        # mutate(bias = count - count_GRACE) %>%
        # summarise(weighted_median_bias = weightedMedian(bias, w = weights, na.rm = TRUE)) %>%
        pull(weighted_median_bias)
      median_bias_tropics_label <- bquote("Median bias"[tropics] == .(round(median_bias_tropics, 0)) * "%")
      print(median_bias_tropics_label)
    }

    # Plot global map for each model
    p <- ggplot() +
      theme_void() +
      geom_tile(data = df_model,
                aes(x = lon, y = lat, color = count, fill = count)) +
      geom_sf(data = ocean,
              color = "black",
              linetype = 'solid',
              fill = "white",
              size = 1) +

      # Orange-Blue color bar
      scale_color_viridis_c(
        option = "turbo",
        breaks = breaks,
        limits = c(lower_threshold, upper_threshold),
        na.value = "#8C8C8C"  # NA color #969696 or #A0A0A0
      ) +
      scale_fill_viridis_c(
        option = "turbo",
        breaks = breaks,
        limits = c(lower_threshold, upper_threshold),
        na.value = "#8C8C8C"  # NA color
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
        title = current_mdl,
        color = "Frequency of water limitation (%)",
        fill = "Frequency of water limitation (%)",
        # fill = expression(paste("number of months when SM < ", theta[crit], " (-) "))
      )

    # more specific title rather than just "Observations"
    if(current_mdl == "Observations") {
      p <- p + labs(
        title = "GOME-2 SIF and GRACE observations")
    }

    # add stats for models
    if(current_mdl != "Observations") {
      p <- p +
        # annotate("text", x = -175, y = 11, label = r_squared_label, hjust = 0, vjust = 0, size = 3.8) + # R2
        # annotate("text", x = -175, y = 0, label = median_bias_label, hjust = 0, vjust = 0, size = 3.8) + # median of bias
        # annotate("text", x = -175, y = -13,  label = median_bias_tropics_label, hjust = 0, vjust = 0, size = 3.8) + # median of bias in the tropics
        annotate("text", x = -175, y = -11, label = r_squared_label, hjust = 0, vjust = 0, size = 3.8) + # R2
        annotate("text", x = -175, y = -25, label = mean_bias_label, hjust = 0, vjust = 0, size = 3.8) + # mean bias
        annotate("text", x = -175, y = -40, label = mean_bias_tropics_label, hjust = 0, vjust = 0, size = 3.8) + # tropics
        annotate("text", x = -175, y = -52, label = mean_bias_abs_label, hjust = 0, vjust = 0, size = 3.8) # abs(mean bias)
    }

    return(p)
  })

  # Arrange and save the plots for the current scenario
  all <- ggarrange(
    plotlist = plot_list,
    labels = "auto",
    ncol = 2, nrow = 1,
    common.legend = TRUE,
    legend = "bottom"
  )

  # print in png
  filename <- paste0("map_watlim_", scene, "_MMmean.png")
  ggsave(filename,
         plot = all,
         path = "./",
         width = 12,
         height = 3.25,
         dpi = 300)

  # print in PDF
  filename <- paste0("map_watlim_", scene, "_MMmean.pdf")
  ggsave(filename,
         plot = all,
         device = cairo_pdf, # save in PDF vectographic format (for publishing)
         path = "./",
         width = 12,
         height = 3.25)
}


