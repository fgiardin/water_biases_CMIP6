# Script to plot daily GRACE values and compared them with soil moisture simulated at Fluxnet locations of Fig. 3

# load libraries
library(tidyverse)
library(ggpubr)
library(data.table)
library(rnaturalearth)
library(sf)
library(patchwork)
library(ragg) # to save plots with a lot of points faster
library(cowplot)

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

# filter extreme values
df_merged <- df_merged %>%
  dplyr::filter()


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



