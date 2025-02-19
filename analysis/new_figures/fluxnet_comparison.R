# script to plot EF vs SM at three eddy-covariance sites, focusing on DAILY fluxnet vs models
# using flux data vs model data to show discrepancies + illustrate our methods

# load libraries
library(tidyverse)
library(ggpubr)
library(data.table)
library(rnaturalearth)
library(sf)
library(patchwork)
library(ragg) # to save plots with a lot of points faster
library(cowplot)


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
plot_list[['site1_flux']] <- plot_flux(df_daily_flux, site1, show_x = FALSE, show_y = TRUE) # show_x = FALSE, show_y = TRUE
plot_list[['site2_flux']] <- plot_flux(df_daily_flux, site2, show_x = FALSE, show_y = TRUE) # show_x = FALSE, show_y = TRUE
plot_list[['site3_flux']] <- plot_flux(df_daily_flux, site3, show_x = TRUE, show_y = TRUE) # show_x = TRUE, show_y = TRUE

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
  plot_list[[plot_name_site3]] <- plot_flux(filtered_data, site3, show_x = TRUE, show_y = FALSE) # show_x = TRUE, show_y = FALSE
}

# Extract plots in the desired order
plots_in_order <- list(
  plot_list[["site1_flux"]], plot_list[["site1_UKESM1-0-LL"]], plot_list[["site1_IPSL-CM6A-LR"]], plot_list[["site1_EC-Earth3-Veg"]], plot_list[["site1_CNRM-ESM2-1"]],
  plot_list[["site2_flux"]], plot_list[["site2_UKESM1-0-LL"]], plot_list[["site2_IPSL-CM6A-LR"]], plot_list[["site2_EC-Earth3-Veg"]], plot_list[["site2_CNRM-ESM2-1"]],
  plot_list[["site3_flux"]], plot_list[["site3_UKESM1-0-LL"]], plot_list[["site3_IPSL-CM6A-LR"]], plot_list[["site3_EC-Earth3-Veg"]], plot_list[["site3_CNRM-ESM2-1"]]
)



# Plot --------------------------------------------------------------------

### Create plot of all panels
all <- (wrap_plots(plots_in_order, ncol = 5, nrow = 3) +
          plot_layout(guides = "collect") +
          plot_annotation(tag_levels = "a")) &
  theme(plot.tag = element_text(face = "bold"))

# ggsave("EFvsSM_sites_justplot.png", plot = all, device = agg_png, width = 14, height = 8.7, dpi = 300) # width = 15, height = 9.2


### Create plot of column titles
obs_plot <- ggplot() + # Create the left title: "Observations"
  theme_void() +
  annotate("text", x = 0.5, y = 0.5,
           label = "Observations",
           size = 7.5, fontface = "bold", color = "black") +
  labs(tag = "")

models_plot <- ggplot() + # Create the right title: "CMIP6 models"
  theme_void() +
  annotate("text", x = 0.5, y = 0.5,
           label = "CMIP6 models",
           size = 7.5, fontface = "bold", color = "#4D4D4D") +
  labs(tag = "")

# Combine the two plots into one row
title_row <- obs_plot + models_plot +
  plot_layout(ncol = 2, widths = c(1.3, 4)) # ensure proportions with the rest of the plot


### Combine the title row and the main grid using cowplot
# final_plot <- title_row / all +
#   plot_layout(heights = c(0.08, 1))
final_plot <- plot_grid(title_row, all, ncol = 1, rel_heights = c(0.08, 1))

# add common x and y axes
final_plot_no_axis <- final_plot &  # Remove all axes title themes from final_plot
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0, 5.5, 15, 15)  # top, right, bottom, left (in points or use unit("cm"))
  )

final_with_axis <- ggdraw(final_plot_no_axis, clip = "off") +  # Use ggdraw() to overlay overall labels
  draw_label("Normalized land water storage (-)", x = 0.5, y = 0.005, hjust = 0.5, vjust = 0, size = 18) +
  draw_label("Evaporative Fraction (-)", x = 0.01, y = 0.5, hjust = 0.5, angle = 90, size = 18)


ggsave("EFvsSM_sites.png", plot = final_with_axis, device = agg_png, width = 14, height = 8.7, dpi = 300) # width = 15, height = 9.2

# ggsave("EFvsSM_sites.png", path = "./", width = 18, height = 10.8, dpi= 300)







