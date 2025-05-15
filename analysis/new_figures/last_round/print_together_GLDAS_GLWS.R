# script to simply print together in one figure results of GLDAS and GLWS for water limitation

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
library(RColorBrewer)
setwd("/Users/fgiardina/water_biases_CMIP6")


# load data
plot_list_GLDAS <- readRDS("data/Comparisons_GRACE/water_limitation/plot_list_GLDAS.rds")
plot_list_GLWS <- readRDS("data/Comparisons_GRACE/water_limitation/plot_list_GLWS.rds")

all_plots <- c(plot_list_GLDAS, plot_list_GLWS)

# Arrange plots
all <- ggarrange(
  plotlist = all_plots,
  labels = "auto",
  ncol = 2, nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)

# save
ggsave("Fig_SI_watlim_GLDAS_GLWS.png", plot = all, path = "./", width = 12, height = 6, dpi = 300)
