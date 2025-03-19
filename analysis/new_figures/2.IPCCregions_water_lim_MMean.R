# script to plot water limitation results grouped by IPCC region
# multi-model mean

# load packages
library(tidyverse)
library(terra)
library(ggpubr)
library(data.table)

# load plot data
df_count_GRACE <- readRDS("data/theta_crit/monthly/df_count_GRACE.rds") %>%
  dplyr::select(-date, -SIF, -TWS) %>%  # remove temporal values
  unique() %>%
  mutate(model = "Observations")

df_count_mrso <- readRDS("data/theta_crit/monthly/df_count_mrso_allscenarios.rds") %>%
  dplyr::filter(scenario == "land-hist") %>%
  dplyr::select(-date, -EF, -Rn, -mrso, -mrso_norm, -scenario) %>%
  unique()

dt_count <- bind_rows(df_count_mrso, df_count_GRACE) %>%   # NAs = cells where the segmented function couldn't find a breakpoint (not enough data at that cell)
  rename(model_name = model) %>%
  mutate(count = count*100) %>%
  mutate(count = ifelse(Intercept < EFmax,
                        count, # keep only rows with meaningful intercept (when the EF vs SM relationship decreases)
                        0))   # this way 0s represent locations that are never water-limited according to our definition!

# make sure all models have same cells that are NA
setDT(dt_count) # Convert to data.table if it's not already
na_locations <- unique(dt_count[is.na(theta_crit), .(lon, lat)]) # Identify all unique lon and lat combinations where any of the columns of interest are NA
dt_count[na_locations, # Update the original data.table to set the values to NA for the identified lon and lat, for all models
         on = .(lon, lat),
         `:=` (theta_crit = NA_real_,
               EFmax = NA_real_,
               Slope = NA_real_,
               Intercept = NA_real_,
               count = NA_real_)]

# calculate multi-model mean
MMmeans <- dt_count %>%
  dplyr::filter(model_name != "Observations") %>%
  group_by(lon, lat) %>%
  summarise(
    count = mean(count, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model_name = "Multi-model mean") # add model name to mirror summary_deltaSM

dt_count <- dt_count %>% # bind rows to original dt
  bind_rows(MMmeans) %>%
  dplyr::select(-theta_crit, -EFmax, -Slope, -Intercept)

# load IPCC data
IPCC_mask <- rast("data-raw/IPCC_regions/ar6_land_regions_g025_2D.nc")
IPCC_mask <- terra::rotate(IPCC_mask)


# process and merge data --------------------------------------------------

df_IPCC <- terra::as.data.frame(IPCC_mask, xy = TRUE) %>%
  rename(lon = x,
         lat = y,
         ID = mask)

# Create a dataframe with the region codes (including all land IPCC regions)
regions <- c('GIC', 'NWN', 'NEN', 'WNA', 'CNA', 'ENA', 'NCA', 'SCA', 'CAR',
             'NWS', 'NSA', 'NES', 'SAM', 'SWS', 'SES', 'SSA', 'NEU', 'WCE',
             'EEU', 'MED', 'SAH', 'WAF', 'CAF', 'NEAF', 'SEAF', 'WSAF', 'ESAF',
             'MDG', 'RAR', 'WSB', 'ESB', 'RFE', 'WCA', 'ECA', 'TIB', 'EAS',
             'ARP', 'SAS', 'SEA', 'NAU', 'CAU', 'EAU', 'SAU', 'NZ', 'EAN', 'WAN')

full_region_names <- c('GIC'  = 'Greenland/Iceland', 'NWN'  = 'N.W.North-America', 'NEN'  = 'N.E.North-America',
                       'WNA'  = 'W.North-America', 'CNA'  = 'C.North-America', 'ENA'  = 'E.North-America',
                       'NCA'  = 'N.Central-America', 'SCA'  = 'S.Central-America', 'CAR'  = 'Caribbean',
                       'NWS'  = 'N.W.South-America','NSA'  = 'N.South-America', 'NES'  = 'N.E.South-America',
                       'SAM'  = 'South-American-Monsoon', 'SWS'  = 'S.W.South-America', 'SES'  = 'S.E.South-America',
                       'SSA'  = 'S.South-America', 'NEU'  = 'N.Europe', 'WCE'  = 'West&Central-Europe',
                       'EEU'  = 'E.Europe', 'MED'  = 'Mediterranean', 'SAH'  = 'Sahara',
                       'WAF'  = 'Western-Africa', 'CAF'  = 'Central-Africa', 'NEAF' = 'N.Eastern-Africa',
                       'SEAF' = 'S.Eastern-Africa', 'WSAF' = 'W.Southern-Africa', 'ESAF' = 'E.Southern-Africa',
                       'MDG'  = 'Madagascar', 'RAR'  = 'Russian-Arctic', 'WSB'  = 'W.Siberia', 'ESB'  = 'E.Siberia',
                       'RFE'  = 'Russian-Far-East', 'WCA'  = 'W.C.Asia', 'ECA'  = 'E.C.Asia', 'TIB'  = 'Tibetan-Plateau',
                       'EAS'  = 'E.Asia', 'ARP'  = 'Arabian-Peninsula', 'SAS'  = 'S.Asia', 'SEA'  = 'S.E.Asia',
                       'NAU'  = 'N.Australia', 'CAU'  = 'C.Australia', 'EAU'  = 'E.Australia', 'SAU'  = 'S.Australia',
                       'NZ'   = 'New-Zealand', 'EAN'  = 'E.Antarctica', 'WAN'  = 'W.Antarctica')
numbers <- 0:45
df <- data.frame(ID = numbers, Region = regions)

# merge to original IPCC df
df_IPCC <- df_IPCC %>%
  left_join(df, by = join_by(ID)) %>%
  dplyr::filter(ID > 0)

# merge to count df
df_plot <- dt_count %>%
  left_join(df_IPCC, by = join_by(lon, lat)) %>%
  # drop_na() %>%
  group_by(Region, model_name) %>% # in every IPCC cell and per every model, calculate mean/median count
  # dplyr::filter(n() >= 10) %>% # only retain Regions with values from at least 10 grid cells
  summarise(median_count = mean(count, na.rm = TRUE)) %>%
  ungroup()

# Filter out the data for the observations (count based on SIF vs GRACE)
grace_data <- df_plot %>%
  filter(model_name == "Observations") %>%
  dplyr::select(-model_name) # remove model name column (as GRACE data will have its own column in df_merged)

# Create a new (repeated for every model) column with GRACE data in the original dataframe
df_merged <- full_join(df_plot %>% dplyr::filter(model_name != "Observations"), # remove GRACE in the original dataframe to avoid repetitions
                       grace_data, by = "Region",
                       suffix = c("", "_obs")) # rename new column "count" + "_obs"

# Get the unique model names (excluding .GRACE) and create scatter plots.
unique_models <- unique(df_merged$model_name)

# Extract unique regions from df_merged
unique_regions_in_data <- unique(df_merged$Region)

# manually remove Sahara and Arabian Peninsula because not vegetated
# they're not removed automatically in previous step because there still are very few pixels in those areas
unique_regions_in_data <- unique_regions_in_data[!unique_regions_in_data %in% c("SAH", "ARP")]

# Filter the original region vector to include only those present in df_merged
filtered_regions <- regions[regions %in% unique_regions_in_data]

# Reorder the Region factor in df_merged according to the filtered_regions
df_merged$Region <- factor(df_merged$Region, levels = filtered_regions)


# manage colors -----------------------------------------------------------

# Base colors for major regions (one color per region)
base_colors <- list(
  "North and Central America" = "#00008B", # DarkBlue
  "South America"             = "#4CAF4C", # DarkGreen
  "Africa"                    = "#CD3333", # DarkRed
  "Europe"                    = "#FFD700", # Gold
  "Russia/Asia"               = "#BF3EFF", # Purple
  "Australia, New Zealand and South East Asia" = "#FF7F00" # Orange
)

# sub-regions in order of appearance within each major region
macro_regions <- list(
  "North and Central America" = c('NWN', 'NEN', 'WNA', 'CNA', 'ENA', 'NCA', 'SCA', 'CAR'),
  "South America" = c('NWS', 'NSA', 'NES', 'SAM', 'SWS', 'SES', 'SSA'),
  "Europe" = c('NEU', 'WCE', 'EEU', 'MED'),
  "Africa" = c('SAH', 'WAF', 'CAF', 'NEAF', 'SEAF', 'WSAF', 'ESAF', 'MDG'),
  "Russia/Asia" = c('RAR', 'WSB', 'ESB', 'RFE', 'WCA', 'ECA', 'TIB', 'EAS', 'ARP', 'SAS'),
  "Australia, New Zealand and South East Asia" = c('SEA', 'NAU', 'CAU', 'EAU', 'SAU', 'NZ')
)
macro_regions <- lapply(macro_regions, intersect, filtered_regions) # only keep regions that appear in the final dataset

# Assigning the base color to each sub-region within the macro regions
region_colors <- unlist(lapply(names(macro_regions), function(region_name) {
  setNames(rep(base_colors[[region_name]], length(macro_regions[[region_name]])), macro_regions[[region_name]])
}))

# Assign back the names of the regions
names(region_colors) <- filtered_regions


# manage shape of points --------------------------------------------------

# Initialize a vector to hold the shapes for each region
shapes <- numeric(length(filtered_regions))

# vector to be replicated (9 values max)
# numbers correspond to the shape numbers of ggplot
# x = 7:18
x = c(15:17, 4, 8, 9, 10, 14, 18)

# For each macro-region, assign shapes and reset for each macro-region
for (macro_region in names(macro_regions)) {
  sub_regions <- macro_regions[[macro_region]]
  shapes[which(filtered_regions %in% sub_regions)] <- rep(x, length.out = length(sub_regions))
}

df_merged <- df_merged %>%
  drop_na()



# scatter plots -----------------------------------------------------------

plots <- list()
unique_models <- "Multi-model mean"

for(model in unique_models) {

  # filter model info only of the current model
  df_model <- df_merged %>% filter(model_name == model)

  # calculate R2 within each model
  fit <- lm(median_count ~ median_count_obs, data = df_model)
  r_squared <- summary(fit)$r.squared
  r_squared_label <- bquote(italic(R)^2 == .(round(r_squared, 2)))

  plots[[model]] <- ggplot(data = df_model,
                           aes(x = median_count_obs, # same x-axis for all sub-plots
                               y = median_count,
                               color = Region,
                               shape = Region # as.factor(unique_id)
                           )) +

    # add R2
    annotate("text", x = 0.7, y = 93, label = r_squared_label, hjust = 0, vjust = 0, size = 5) +

    # Add the scatter plot points
    geom_point(size = 2) +  # Add the scatter plot points
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add the y=x dashed line
    labs(x = "Observations (%)",
         y = paste(model, "(%)"), # Adding (%) to the y-axis title
         color = "IPCC WGI reference regions", # legend title
         shape = "IPCC WGI reference regions",
         title = "Frequency of water limitation"
    ) +
    theme_minimal(base_size = 14) +  # consistently increase base font
    theme(axis.ticks = element_line(color = "black"), # add axes ticks
          panel.grid.major = element_blank(),  # Remove major gridlines
          panel.grid.minor = element_blank(),  # Remove minor gridlines
          panel.background = element_rect(fill = "white"),  # Set background to white
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # Add a border around the plot
          plot.title = element_text(size = 14, hjust = 0.5)
    ) +

    scale_color_manual(
      name = "IPCC WGI reference regions", # to share same legend, the two scale_*_manual need to have same name, labels and limits
      labels = full_region_names[filtered_regions], # use full-length region names
      limits = filtered_regions, # display regions in right order
      values = region_colors, # Use the custom color vector
    ) +

    scale_shape_manual( # assign different shapes to different regions in each region group
      name = "IPCC WGI reference regions",
      labels = full_region_names[filtered_regions],
      limits = filtered_regions, # display regions in right order
      values = shapes
    ) +
    guides(
      color = guide_legend(
        # title = "IPCC WGI reference regions",
        label.position = "right",
        title.position = "top",
        title.hjust = 0.5),
      shape = guide_legend(
        # title = "IPCC WGI reference regions",
        label.position = "right",
        title.position = "top",
        title.hjust = 0.5)
    ) +
    ylim(0, 100) +
    xlim(0, 100)

}

# calculate distance from 1:1 line

df_merged <- df_merged %>%
  filter(model_name == "Multi-model mean") %>%
  mutate(distance_from_1to1 = (median_count - median_count_obs) / sqrt(2))

# save plot
empty_plot <- ggplot() + # add letter to empty space in plot
  theme_void() +
  ggtitle("")  # remove any default title
plots_with_empty <- c(plots, list(empty_plot))

all <- ggarrange(plotlist = plots_with_empty,
                 labels = "auto",
                 ncol = 3, nrow = 1,
                 common.legend = TRUE, # have just one common legend
                 legend= "bottom")

ggsave("IPCC_regions_water-lim.png", path = "./", width = 11, height = 6, dpi= 600) # width = 8, height = 15,


# original phrasing in text
# The most affected include the Congo rainforest (CAF) and Southern Africa (WSAF, ESAF),
# Australia (SAU, EAU, CAU), the Amazon (NWS, NSA, SAM) and the Southern Cone (SSA, SWS),
# and also Eastern Central Asia (ECA) (Fig. 2).









