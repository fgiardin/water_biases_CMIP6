
simple_map <- function(data_to_plot, variable_to_plot, print_map = TRUE, save_map = FALSE) {

  # load packages
  library(tidyverse)
  library(raster)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(sf)

  data_to_plot <- theta_crit_mrsol_CESM2
  variable_to_plot <- "theta_crit"

  # download countries
  countries <- ne_countries(scale = 50, returnclass = c("sf"))

  p <- ggplot() +
    theme_void() +
    geom_tile(data = data_to_plot,
              aes(x = lon, y = lat,
                         color = variable_to_plot,
                         fill = variable_to_plot)) +
    geom_sf(data=countries, # country borders
            color= "grey23",
            linetype='solid',
            fill= NA,
            size=0.3)

  # Either print or save the plot
  if (print_map) {
    print(p)
  }

  if (save) {
    ggsave("simple_map.png", plot = p, path = "./", width = 11)
  }
}

# Example usage:
# simple_map(theta_crit_mrsol_CESM2, "theta_crit")


