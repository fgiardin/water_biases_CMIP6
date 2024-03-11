# summary figure of GRACE and ALEXI methods vs CMIP6 multi-model mean grouped by IPCC region!

# load packages
devtools::load_all(".")
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
sf_use_s2(FALSE)

# load data
plot_list_CWD <- readRDS("data/summary-figures/plot_list_maxCWD_IPCCregions.rds")
plot_list_SMmax <- readRDS("data/summary-figures/plot_list_SMmax_IPCCregions.rds")

# combine lists
combined_list <- c(plot_list_SMmax, plot_list_CWD)

# adjust panels titles
combined_list[[1]] <- combined_list[[1]] +
  labs(title = expression(paste(Delta*SM[max]," and ", Delta*TWS[max]))) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

combined_list[[2]] <- combined_list[[2]] +
  labs(title = expression(paste(CWD[max])),
       x = "ALEXI and WATCH-WFDEI CWD (mm)") +
  theme(axis.title.x = element_text(size = 12),  # Smaller x axis label
        axis.title.y = element_text(size = 12))

# Print all plots with one colorbar
all <- ggarrange(plotlist = combined_list,
                 labels = "auto", # "auto"
                 ncol = 2, nrow = 2,
                 common.legend = TRUE, # have just one common legend
                 legend="bottom")

ggsave("Fig_2.png", path = "./", width = 8.8, height = 9.6, dpi=300)

