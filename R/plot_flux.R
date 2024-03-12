
# function to plot EF vs SM at a specific location (FLUXNET2015 site locations)

library(ggplot2)
library(dplyr)

plot_flux <- function(df, site_name) {

  df_site <- df %>%
    filter(sitename == site_name)

  p <- ggplot(df_site, aes(x = SM, y = EF)) +
    geom_point(alpha = 0.2, size = 1) + # alpha = 0.23, size = 1.5
    theme_minimal() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
      axis.ticks = element_line(color = "black"), # add axes ticks
      # panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.background = element_rect(fill = "white"),  # Set background to white
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # Add a border around the plot
    ) +
    labs(
      title = site_name,
      x = "Total water storage (-)",
      y = "Evaporative Fraction (-)"
    ) +
    scale_y_continuous(breaks = seq(0, 1.4, 0.2), limits = c(0, 1.5), expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0))

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
                         paste("Water limitation:", count_value, "%"))

    p <- p +
      annotate("text", x = 0.3, y = 0.3, label = label_text,
               hjust = 0, # Align text to the left
               vjust = 1.1,
               size = 4.4, color = "black")
  }

  return(p)
}
