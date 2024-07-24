
# script to plot EF vs SM at three eddy-covariance sites, focusing on DAILY fluxnet vs models

# using flux data and model data to show discrepancies and illustrate our methods

# load libraries
library(tidyverse)
library(ggpubr)
library(data.table)
library(rnaturalearth)
library(sf)
library(patchwork)


# load daily data
df_daily_cmip6 <- readRDS("data/theta_crit/cmip6_daily_theta_crit_count.rds") %>%
  # EF must drop of at least 30%, otherwise we consider no water limitation (see Methods)
  # drop is defined as the difference between EFmax and the intercept with the y-axis
  dplyr::mutate(count = ifelse(Intercept < EFmax - 0.3, count, NA))

df_daily_flux <- readRDS("data/flux_data/df_flux_allsites.rds") %>%
  rename(SM = soil_moisture) %>%
  # EF must drop of at least 30%, otherwise we consider no water limitation
  # drop is defined as the difference between EFmax and the intercept with the y-axis
  dplyr::mutate(
    # theta_crit = ifelse(Intercept < EFmax - 0.3, theta_crit, NA),
    count = ifelse(Intercept < EFmax - 0.3, count, NA))

# specify sites to plot
site1 = "GF-Guy"
site2 = "BR-Sa3"
site3 = "US-Ton"
plot_model = c("UKESM1-0-LL", "IPSL-CM6A-LR", "EC-Earth3-Veg", "CNRM-ESM2-1")



# plots by-model -------------------------------------------------------------------


# function to create plots at specific locations
source("R/plot_flux.R")

plot_list <- list()

# flux
plot_list[['site1_flux']] <- plot_flux(df_daily_flux, site1, show_x = FALSE, show_y = TRUE)
plot_list[['site2_flux']] <- plot_flux(df_daily_flux, site2, show_x = FALSE, show_y = TRUE)
plot_list[['site3_flux']] <- plot_flux(df_daily_flux, site3, show_x = TRUE, show_y = TRUE)

# Loop over models
for (pmodel in plot_model) {
  # Filter data for the current model
  filtered_data <- df_daily_cmip6 %>% filter(model == pmodel)

  # Generate and save plot for site1
  plot_name_site1 <- paste("site1", pmodel, sep = "_")
  plot_list[[plot_name_site1]] <- plot_flux(filtered_data, site1, show_x = FALSE, show_y = FALSE)

  # Generate and save plot for site2
  plot_name_site2 <- paste("site2", pmodel, sep = "_")
  plot_list[[plot_name_site2]] <- plot_flux(filtered_data, site2, show_x = FALSE, show_y = FALSE)

  # Generate and save plot for site3
  plot_name_site3 <- paste("site3", pmodel, sep = "_")
  plot_list[[plot_name_site3]] <- plot_flux(filtered_data, site3, show_x = TRUE, show_y = FALSE)
}

# Extract plots in the desired order
plots_in_order <- list(
  plot_list[["site1_flux"]], plot_list[["site1_UKESM1-0-LL"]], plot_list[["site1_IPSL-CM6A-LR"]], plot_list[["site1_EC-Earth3-Veg"]], plot_list[["site1_CNRM-ESM2-1"]],
  plot_list[["site2_flux"]], plot_list[["site2_UKESM1-0-LL"]], plot_list[["site2_IPSL-CM6A-LR"]], plot_list[["site2_EC-Earth3-Veg"]], plot_list[["site2_CNRM-ESM2-1"]],
  plot_list[["site3_flux"]], plot_list[["site3_UKESM1-0-LL"]], plot_list[["site3_IPSL-CM6A-LR"]], plot_list[["site3_EC-Earth3-Veg"]], plot_list[["site3_CNRM-ESM2-1"]]
)



# Create one plot
all <- wrap_plots(plots_in_order, ncol = 5, nrow = 3) + plot_layout(guides = "collect")

# Create column titles
column_titles <- c("FLUXNET2015", "UKESM1-0-LL", "IPSL-CM6A-LR", "EC-Earth3-Veg", "CNRM-ESM2-1")
title_plots <- lapply(column_titles, function(title) {
  ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = title, size = 6.5, fontface = "bold")
})

# Combine the titles and the plot grid
title_row <- wrap_plots(title_plots, ncol = 5)
final_plot <- wrap_plots(title_row, all, ncol = 1, heights = c(0.08, 1))

ggsave("EFvsSM_sites.png", path = "./", width = 14, height = 8.7, dpi= 300) # width = 15, height = 9.2

# ggsave("EFvsSM_sites.png", path = "./", width = 18, height = 10.8, dpi= 300)

