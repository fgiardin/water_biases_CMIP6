# Script compare daily GRACE values with soil moisture simulated at Fluxnet locations of Fig. 3
# plots comparison scatter plot and broken stick regression like fig. 3

# load libraries
library(tidyverse)
library(ggpubr)
library(data.table)
library(rnaturalearth)
library(sf)
library(patchwork)
library(ragg) # to save plots with a lot of points faster
library(cowplot)
library(segmented)

# load daily GRACE
df_daily_GRACE <- read_csv("data/GRACE_daily/Grace_daily.csv",     # automatically guesses header, delim, types
                           na = c("", "NA")) %>%
  pivot_longer(  # conver to same format as fluxnet
    cols      = c("GF-Guy", "BR-Sa3", "US-Ton"),  # the columns to stack
    names_to  = "sitename",                       # new column for the old names
    values_to = "SM_GRACE"                           # new column for the values
  ) %>%
  rename(date = Date) %>%
  mutate(
    date = dmy(date) # convert to proper date format
  ) %>%
  group_by(sitename) %>% # normalize SM (like for fluxnet)
  mutate(SM_GRACE = (SM_GRACE - min(SM_GRACE)) / (max(SM_GRACE) - min(SM_GRACE))) %>%
  ungroup()

# load flux data
df_daily_flux <- readRDS("data/flux_data/df_flux_allsites.rds") %>%
  rename(SM = soil_moisture) %>%
  dplyr::select(date, sitename, SM) %>%
  dplyr::filter(sitename %in% c("GF-Guy", "BR-Sa3", "US-Ton"))

# merge two datasets
df_merged <- df_daily_flux %>%
  inner_join(df_daily_GRACE,
             by = c("date", "sitename"))

# scatter plot -----------------------------------------------------------

# Compute statistics
stats <- df_merged %>%
  summarise(
    R2     = cor(SM, SM_GRACE)^2,
    Bias   = mean(SM - SM_GRACE),
    RMSE   = sqrt(mean((SM - SM_GRACE)^2)),
    Mean   = mean(c(SM, SM_GRACE)),
    rBias  = Bias / Mean,
    rRMSE  = RMSE / Mean
  )

xpos <- 0.01 # Choose annotation positions (10% in from the top/right)
ypos <- 0.99


# repeat with subset of data
df_sub <- df_merged %>% # Define subset
  filter(SM >= 0.25, SM <= 0.75)

stats_sub <- df_sub %>% # Compute stats on the subset
  summarise(
    R2    = cor(SM, SM_GRACE)^2,
    Bias  = mean(SM - SM_GRACE),
    RMSE  = sqrt(mean((SM - SM_GRACE)^2)),
    Mean  = mean(c(SM, SM_GRACE)),
    rBias = Bias / Mean,
    rRMSE = RMSE / Mean
  )

xpos_sub <- 0.25 # Choose position
ypos_sub <- 0.99


# Plot
a_plot <- ggplot(df_merged, aes(x = SM, y = SM_GRACE)) +
  geom_point(alpha = 0.5) +                         # all points
  # geom_point(data = df_sub, aes(x = SM, y = SM_GRACE),   # highlight subset in green
  #            color = "forestgreen", size = 2, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +           # 1:1 line
  geom_smooth(method = "lm", se = FALSE, color = "blue") +          # regression line
  annotate("text",
           x = xpos, y = ypos,
           hjust = 0, vjust = 1,
           # color = "red",
           label = paste0(
             "R² = ",  round(stats$R2,   2), "\n",
             "rBias = ", round(stats$rBias, 2), "\n",
             "rRMSE = ",round(stats$rRMSE, 2)
           )
  ) +
  # annotate("text",
  #          x = xpos_sub, y = ypos_sub,
  #          hjust = 0, vjust = 1,
  #          color = "forestgreen",
  #          label = paste0(
  #            "R² = ", round(stats_sub$R2,   2), "\n",
  #            "rBias = ", round(stats_sub$rBias, 3), "\n",
  #            "rRMSE = ",round(stats_sub$rRMSE,3)
  #          )
  # ) +
  labs(
    x = "SPLASH soil moisture (-)",
    y = "GRACE daily (-)"
  ) +
  theme_minimal()


ggsave("GRACE_daily_flux_comparison.png", plot = a_plot,
       width = 4, height = 4, dpi = 600)



# EF vs SM plots ----------------------------------------------------------

# load daily GRACE
df_daily_GRACE <- read_csv("data/GRACE_daily/Grace_daily.csv",     # automatically guesses header, delim, types
                           na = c("", "NA")) %>%
  pivot_longer(  # conver to same format as fluxnet
    cols      = c("GF-Guy", "BR-Sa3", "US-Ton"),  # the columns to stack
    names_to  = "sitename",                       # new column for the old names
    values_to = "SM_GRACE"                           # new column for the values
  ) %>%
  rename(date = Date) %>%
  mutate(
    date = dmy(date) # convert to proper date format
  ) %>%
  group_by(sitename) %>% # normalize SM (like for fluxnet)
  mutate(SM_GRACE = (SM_GRACE - min(SM_GRACE)) / (max(SM_GRACE) - min(SM_GRACE))) %>%
  ungroup()

# load EF data from FLUXNET (processing consistent with script 3.fluxnet_comparison)
df_daily_flux_EF <- readRDS("data/flux_data/df_flux_allsites.rds") %>%
  rename(SM = soil_moisture) %>%
  # EF must drop of at least 30%, otherwise we consider no water limitation
  # drop is defined as the difference between EFmax and the intercept with the y-axis
  dplyr::mutate(
    # theta_crit = ifelse(Intercept < EFmax - 0.3, theta_crit, NA),
    count = ifelse(Intercept < EFmax - 0.3, count, NA)) %>%
  dplyr::select(date, sitename, EF)

# merge two datasets
df_merged_EF <- df_daily_flux_EF %>%
  inner_join(df_daily_GRACE,
             by = c("date", "sitename"))


# fit bilinear for all three sites simultaneously
source("R/fit_bilinear.R")

df_out_all <- df_merged_EF %>%
  # 1) group & pack each site's rows into a list‐column called `data`
  nest_by(sitename, .keep = FALSE) %>%

  # 2) for each site, run fit_bilinear() on its data
  mutate(params = list(fit_bilinear(data, "EF", "SM_GRACE"))) %>%

  # 3) unpack back out (still per site!!!):
  unnest(data) %>%            # brings back date, EF, SM_GRACE
  unnest_wider(params) %>%    # spreads theta_crit, EFmax, Slope, Intercept into columns

  ungroup() %>%
  rename(SM = SM_GRACE) %>% # rename to use plot_flux function
  group_by(sitename) %>%
  mutate(
    count = ifelse(
      all(is.na(theta_crit)), NA, # Keep NA if all theta_crit values are NA (locations where it's never limited)
      sum(SM < theta_crit, na.rm = TRUE) / n() # Otherwise, calculate percentage of time under SM limitation
    )
  ) %>%
  ungroup() %>%
  dplyr::mutate(count = ifelse(Intercept < EFmax - 0.3, count, NA)) # set water_lim to zero if less than 30% drop in EF (see Methods)


# function to create plots at specific locations
source("R/plot_flux.R")

# (only these threes sites here)
site1 = "GF-Guy"
site2 = "BR-Sa3"
site3 = "US-Ton"

plot_list <- list()

# flux
plot_list[['site1_flux']] <- plot_flux(df_out_all, site1, show_y = TRUE, show_x = TRUE) + labs(x = "ALW GRACE (-)", y = "EF (-)")
plot_list[['site2_flux']] <- plot_flux(df_out_all, site2, show_y = TRUE, show_x = TRUE) + labs(x = "ALW GRACE (-)", y = "")
plot_list[['site3_flux']] <- plot_flux(df_out_all, site3, show_y = TRUE, show_x = TRUE) + labs(x = "ALW GRACE (-)", y = "")

# plot
all <- ggarrange(
  plotlist = plot_list,
  labels = "auto",
  ncol = 3, nrow = 1,
  common.legend = TRUE,
  legend = "bottom"
)

# save
ggsave("EFvsSM_sites_GRACE.png", plot = all, device = agg_png, width = 8.7, height = 3, dpi = 300)

