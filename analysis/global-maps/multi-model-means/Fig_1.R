# summary figure of GRACE and ALEXI methods vs CMIP6 multi-model mean (4 panels)

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
plot_list_CWD <- readRDS("data/summary-figures/plot_list_CWD.rds")
plot_list_SMmax <- readRDS("data/summary-figures/plot_list_SMmax.rds")

# combine lists
combined_list <- c(plot_list_SMmax, plot_list_CWD)

# adjust panels titles
combined_list[[1]] <- combined_list[[1]] +
  labs(title = expression(paste(Delta*TWS[max],": GRACE observations"))) +
  theme(legend.title = element_text(color = "black", size=15))

combined_list[[2]] <- combined_list[[2]] +
  labs(title = expression(paste(Delta*SM[max],": CMIP6 multi-model mean")))+
  theme(legend.title = element_text(color = "black", size=15))

combined_list[[3]] <- combined_list[[3]] +
  labs(title = expression(paste(CWD[max],": ALEXI and WATCH-WFDEI observations"))) +
  theme(legend.title = element_text(color = "black", size=15))

combined_list[[4]] <- combined_list[[4]] +
  labs(title = expression(paste(CWD[max],": CMIP6 multi-model mean"))) +
  theme(legend.title = element_text(color = "black", size=15))

# adjust legend titles
for(i in seq_along(combined_list)) {
  combined_list[[i]] <- combined_list[[i]] +
    labs(fill = "Maximum water storage (mm)", color = "Maximum water storage (mm)")
}

# Print all plots with one colorbar
all <- ggarrange(plotlist = combined_list,
                 labels = "auto", # "auto"
                 ncol = 2, nrow = 2,
                 common.legend = TRUE, # have just one common legend
                 legend="bottom")

ggsave("Fig_1.png", path = "./", width = 12, height = 6, dpi=300)

