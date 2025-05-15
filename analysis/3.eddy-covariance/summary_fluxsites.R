# summary plot showing results for all 128 fluxnet2015 sites (Supplementary Information)

# load libraries
library(terra)
library(tidyverse)
library(ggpubr)
library(scales)

# load count data for flux and cmip6 daily
df_count_raw <- readRDS("data/theta_crit/cmip6_daily_theta_crit_count.rds")
df_count_flux_raw <- readRDS("data/flux_data/theta_crit_flux_norm.rds")
df_key <- readRDS("data/flux_data/df_key.rds")

df_count_flux <- df_count_flux_raw %>%
  # keep only rows with meaningful intercept (when the relationship decreases)
  dplyr::mutate(count = ifelse(Intercept < EFmax - 0.3, count, 0)) %>%
  dplyr::select(sitename, count) %>%
  mutate(model = "FLUXNET2015") %>%
  unique()

df_count_models <- df_count_raw %>%
  # keep only rows with meaningful intercept (when the EF vs SM relationship decreases)
  dplyr::mutate(count = ifelse(Intercept < EFmax - 0.3, count, 0)) %>%
  dplyr::select(sitename, model, count) %>%
  unique() %>%
  group_by(sitename) %>%
  summarize(count = mean(count)) %>%  # calculate multi-model mean
  ungroup() %>%
  mutate(model = "CMIP6 multi-model mean")

df <- rbind(df_count_models, df_count_flux)

# create 2 panels to plot
df <- df[order(df$sitename), ] # sort the dataframe by 'sitename'
index_split <- ceiling(nrow(df) / 2) # create an index for splitting the dataframe in half alphabetically
df$panel_group <- ifelse(seq_along(df$sitename) <= index_split, "Group 1", "Group 2") # Assign 'panel_group' based on the position in the sorted dataframe

# plots -------------------------------------------------------------------

ggplot(df, aes(x = sitename, y = count * 100, color = model)) +
  geom_point() +
  geom_line(aes(group = model)) +
  scale_color_manual(values = c("CMIP6 multi-model mean" = "red", "FLUXNET2015" = "black")) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "lightgrey"),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        strip.text = element_blank()) +  # Remove facet strip labels
  labs(y = "Frequency of water limitation (%)", x = "") +
  facet_wrap(~panel_group, scales = "free_x", ncol = 1)

ggsave("summary_flux.png", path = "./", width = 10, height = 6, dpi= 600)













