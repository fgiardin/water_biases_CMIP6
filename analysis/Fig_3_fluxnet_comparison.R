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

# models
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



# Plot --------------------------------------------------------------------

### Create plot of all panels
all <- (wrap_plots(plots_in_order, ncol = 5, nrow = 3) +
          plot_layout(guides = "collect") +
          plot_annotation(tag_levels = "a")) &
  theme(plot.tag = element_text(face = "bold"))



### Create plot of column titles
obs_plot <- ggplot() +
  # Right boundary line
  annotate("segment", x = 1, xend = 1, y = 0, yend = 1, # y = 0.25, yend = 0.75,
           color = "black", size = 1) +
  # Text in the middle
  annotate("text", x = 0.5, y = 0.5,
           label = "Observations",
           size = 7.5, fontface = "bold", color = "black") +
  # Make sure ggplot shows the entire [0..1] range in both x and y
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  labs(tag = "")
  # # Optionally disable clipping so nothing outside xlim/ylim is cut
  # coord_cartesian(clip = "off") +

  # theme(
  #   plot.margin = margin(0, 0, 0, 0)  # remove outer margins if desired
  # ) +



models_plot <- ggplot() +
  # annotate("rect", xmin = 0, xmax = 1, ymin = 0.25, ymax = 0.75,
  #          color = "black", fill = NA, size = 1) +
  annotate("text", x = 0.5, y = 0.5,
           label = "CMIP6 models",
           size = 7.5, fontface = "bold", color = "#4D4D4D") +
  theme_void() +
  labs(tag = "")


# Combine the two plots into one row
title_row <- obs_plot + models_plot +
  plot_layout(ncol = 2, widths = c(1.385, 4)) # ensure proportions with the rest of the plot
title_row

### Combine the title row and the main grid using cowplot
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


# save in PDF
ggsave("EFvsSM_sites.pdf",
       plot = final_with_axis,
       device = cairo_pdf, # save in PDF vectographic format (for publishing)
       path = "./",
       width = 14,
       height = 8.7) # width = 14, height = 8.7

# save in png
ggsave("EFvsSM_sites.png",
       plot = final_with_axis,
       dpi = 300,
       path = "./",
       width = 14,
       height = 8.7) # width = 14, height = 8.7









