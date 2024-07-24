# function to plot EF vs SM at a specific location (FLUXNET2015 site locations)

library(ggplot2)
library(dplyr)

plot_flux <- function(df, site_name, show_x = TRUE, show_y = TRUE) {

  df_site <- df %>%
    filter(sitename == site_name)

  p <- ggplot(df_site, aes(x = SM, y = EF)) +
    geom_point(alpha = 0.2, size = 1) + # alpha = 0.23, size = 1.5
    theme_minimal() +
    theme(
      axis.text.x = if(show_x) element_text(size = 14) else element_blank(),
      axis.title.x = if(show_x) element_text(size = 16) else element_blank(),
      axis.ticks.x = if(show_x) element_line(color = "black") else element_blank(), # add axes ticks
      axis.text.y = if(show_y) element_text(size = 14) else element_blank(),
      axis.title.y = if(show_y) element_text(size = 16) else element_blank(),
      axis.ticks.y = if(show_y) element_line(color = "black") else element_blank(), # add axes ticks
      legend.text = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16),
      # plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.background = element_blank(),  # Set background to white
      # panel.background = element_rect(fill = "white"),  # Set background to white
      # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) # Add a border around the plot
      panel.border = element_blank(), # Remove the border around the plot
      axis.line = element_line(color = "black") # Add axis lines
    ) +
    labs(
      title = site_name,
      x = "ALW (-)",
      y = "EF (-)"
    ) +
    scale_y_continuous(breaks = seq(0, 1.4, 0.4), limits = c(0, 1.5), expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(0, 1, 0.3), limits = c(0, 1), expand = c(0, 0))

  # plot the segments (results of bilinear regression) only if provided in original df
  if (all(c("theta_crit", "Intercept", "Slope") %in% names(df_site))) {
    p <- p +
      geom_segment(data = data.frame(x = min(df_site$SM), xend = df_site$theta_crit,
                                     y = df_site$Intercept + df_site$Slope * min(df_site$SM), yend = df_site$Slope * df_site$theta_crit + df_site$Intercept),
                   aes(x = x, y = y, xend = xend, yend = yend),
                   color = "red",
                   size = 1) +
      geom_segment(data = data.frame(x = df_site$theta_crit, xend = max(df_site$SM),
                                     y = df_site$EFmax, yend = df_site$EFmax),
                   aes(x = x, y = y, xend = xend, yend = yend),
                   color = "red",
                   size = 1)
  }

  if ("count" %in% names(df_site)) {
    count_value <- round(df_site$count * 100, 2)  # Round and transform to percentage

    # Use ifelse to check for NA and adjust the label accordingly
    label_text <- ifelse(is.na(count_value),
                         "Water limitation: NA",
                         paste0("Water limitation: ", count_value, "%"))

    p <- p +
      annotate("label", x = 0.5, y = 1.35, label = label_text,
               hjust = 0.5, # Center text
               vjust = 0.5,
               size = 4.6, color = "black",
               fill = "white",
               label.size = 0.5,
               label.r = unit(0.15, "lines"))
  }

  return(p)
}
