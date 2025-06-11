# Function to process and plot data for a given scenario (land-hist and historical)

process_scenario_data <- function(summary_deltaSM_scenario, grace_data, scenario_name) {
  # Assign same color to points above upper threshold
  summary_deltaSM_scenario[[plot_variable]][summary_deltaSM_scenario[[plot_variable]] > upper_threshold] <- upper_threshold

  # include GRACE data in summary_deltaSM_scenario
  grace_data_capped <- grace_data
  grace_data_capped[[plot_variable]][grace_data_capped[[plot_variable]] > upper_threshold] <- upper_threshold

  # Combine model data and GRACE data
  summary_deltaSM_scenario <- bind_rows(summary_deltaSM_scenario, grace_data_capped)

  # Calculate multi-model mean (excluding GRACE observations)
  MMmeans <- summary_deltaSM_scenario %>%
    filter(model_name != "GRACE observations") %>%
    group_by(lon, lat) %>%
    summarise(
      deltaSMmax = mean(deltaSMmax, na.rm = TRUE),
      deltaSMabs = mean(deltaSMabs, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(model_name = "Multi-model mean") # add model name to mirror summary_deltaSM

  # Bind rows to summary_deltaSM
  summary_deltaSM_MMmeans <- summary_deltaSM_scenario %>%
    bind_rows(MMmeans)

  # Merge with GRACE data for calculating statistics
  summary_merged <- summary_deltaSM_MMmeans %>%
    # filter(model_name != "GRACE observations") %>% # Exclude GRACE observations to avoid self-comparison
    full_join( # attach GRACE data as new columns
      grace_data, by = join_by(lon, lat),
      suffix = c("", "_GRACE")) # rename new column "deltaSMmax" + "_GRACE"

  # Calculate WEIGHTS for the bias (to account for spherical coordinates)
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

  # Create the list of models to plot
  model_list <- c("GRACE observations", "Multi-model mean")

  # Create separate plots for each model
  plot_list <- lapply(model_list, function(model) {

    df_model <- filter(summary_deltaSM_MMmeans, model_name == model)

    # Calculate stats to compare models with GRACE
    if(model != "GRACE observations") {

      df_stats <- filter(summary_merged, model_name == model)

      # Calculate R squared
      fit <- lm(deltaSMmax ~ deltaSMmax_GRACE, data = df_stats)
      r_squared <- summary(fit)$r.squared
      r_squared_label <- bquote(italic(R)^2 == .(round(r_squared, 2)))

      # Calculate absolute mean bias using weights
      mean_bias_abs <- df_stats %>%
        filter(!is.na(deltaSMmax) & !is.na(deltaSMmax_GRACE) & !is.na(weights)) %>% # manually remove NAs
        mutate(mean_bias_abs = abs(deltaSMmax - deltaSMmax_GRACE) * weights) %>%
        summarise(weighted_mean_bias_abs = sum(mean_bias_abs, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
        pull(weighted_mean_bias_abs)
      mean_bias_abs_label <- bquote("|Bias|" == .(round(mean_bias_abs, 0)) ~ "mm")

      # Calculate mean bias using weights
      mean_bias <- df_stats %>%
        # dplyr::filter(!(lat >= -23 & lat <= 23)) %>%
        dplyr::filter(!is.na(deltaSMmax) & !is.na(deltaSMmax_GRACE) & !is.na(weights)) %>%
        mutate(mean_bias = (deltaSMmax - deltaSMmax_GRACE) * weights) %>%
        summarise(weighted_mean_bias = sum(mean_bias, na.rm = TRUE) / sum(weights), na.rm = TRUE) %>%
        pull(weighted_mean_bias)
      mean_bias_label <- bquote(Bias == .(round(mean_bias, 0)) ~ "mm")
      # mean_bias_label <- bquote(Bias[extratropics] == .(round(mean_bias, 0)) ~ "mm")

      # bias focusing on tropics
      mean_bias_tropics <- df_stats %>%
        dplyr::filter(lat >= -23 & lat <= 23) %>%
        dplyr::filter(!is.na(deltaSMmax) & !is.na(deltaSMmax_GRACE) & !is.na(weights)) %>%
        mutate(mean_bias = (deltaSMmax - deltaSMmax_GRACE) * weights) %>%
        summarise(weighted_mean_bias = sum(mean_bias, na.rm = TRUE) / sum(weights, na.rm = TRUE)) %>%
        pull(weighted_mean_bias)
      mean_bias_tropics_label <- bquote(Bias[tropics] == .(round(mean_bias_tropics, 0)) ~ "mm")
    }

    # Plot global map for each model
    p <- ggplot() +
      theme_void() +
      geom_tile(data = df_model,
                aes(x = lon, y = lat,
                    color = !!sym(plot_variable), fill = !!sym(plot_variable)
                )) +
      geom_sf(data = ocean, # country borders
              color = "black",
              linetype = 'solid',
              fill = "white",
              size = 1) +
      scale_color_viridis_c(option = "turbo",
                            breaks = breaks,
                            limits = c(0, upper_threshold)
      ) +
      scale_fill_viridis_c(option = "turbo",
                           breaks = breaks,
                           limits = c(0, upper_threshold)
      ) +
      coord_sf( # cut Antarctica
        xlim = c(-179.999, 179.999),
        ylim = c(-60, 88),
        expand = FALSE) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 15), # center title
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_text(color = "black", size = 12),
        plot.margin = unit(c(0, -0.3, 0.3, -0.6), 'cm'), # unit(c(top, right, bottom, left), units)
        legend.key.width = unit(2.5, "cm"), # control size of colorbar
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
      # geom_hline(yintercept = 10, color = "red", linetype = "dashed", size = 1) +
      # geom_hline(yintercept = -10, color = "blue", linetype = "dashed", size = 1) +
      # geom_hline(yintercept = 15, color = "red", linetype = "dashed", size = 1) +
      # geom_hline(yintercept = -15, color = "blue", linetype = "dashed", size = 1) +
      # geom_hline(yintercept = 23, color = "red", linetype = "dashed", size = 1) +
      # geom_hline(yintercept = -23, color = "blue", linetype = "dashed", size = 1)

    # Add stats for model panels
    if(model != "GRACE observations") {
      p <- p +
        annotate("text", x = -175, y = -11, label = r_squared_label, hjust = 0, vjust = 0, size = 4.3) + # R2
        annotate("text", x = -175, y = -25, label = mean_bias_label, hjust = 0, vjust = 0, size = 4.3) + # mean bias
        annotate("text", x = -175, y = -40, label = mean_bias_tropics_label, hjust = 0, vjust = 0, size = 4.3) + # tropics
        annotate("text", x = -175, y = -52, label = mean_bias_abs_label, hjust = 0, vjust = 0, size = 4.3) # abs(mean bias)
    }

    return(p)
  })

  # Arrange and save the plots
  all <- ggarrange(plotlist = plot_list,
                   labels = "auto",
                   ncol = 2, nrow = 1,
                   common.legend = TRUE,
                   legend = "bottom")

  # Save the figure
  filename <- paste0("map_", plot_variable, "_MMmean_", scenario_name, ".png")
  ggsave(filename, plot = all, path = "./", width = 12, height = 3.25, dpi = 600)

  # Optionally, save the plot list
  # saveRDS(plot_list, paste0("plot_list_SMmax_", scenario_name, ".rds"), compress = "xz")
}
