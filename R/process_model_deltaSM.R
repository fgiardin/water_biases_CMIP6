# Function to process data for a given model and scenario
# It calculates deltaSM based on mrso (total column soil moisture)

process_model_deltaSM <- function(model_name, scenario) {

  # Construct file path using pattern matching
  file_pattern <- paste0("mrso_mon_", model_name, "_", scenario, "_.*_g025.nc")
  file_list <- list.files("data-raw/cmip6-ng/mrso/mon/g025/", pattern = file_pattern, full.names = TRUE)

  if (length(file_list) == 0) {
    warning(paste("No file found for model", model_name, "and scenario", scenario))
    return(NULL)
  }

  file_path <- file_list[1]  # Assuming there's only one file per model and scenario

  print(paste("Processing:", file_path))

  # Load the data
  mrso_raw <- rast(file_path)

  # Only take data from 2003 onwards to match GRACE data
  dates <- terra::time(mrso_raw)
  selection <- which(dates >= as.POSIXct("2003-01-01"))
  mrso_model <- subset(mrso_raw, subset = selection)

  # Rotate to true coordinates (-180;180 instead of 0;360)
  mrso_model <- terra::rotate(mrso_model)

  # apply land mask
  vegetated_land <- readRDS("data/land_mask/vegetated_land_mask.rds")
  mrso_model <- mask(mrso_model,
                vegetated_land,
                maskvalues = 0)

  # Transform to dataframe
  df_model <- terra::as.data.frame(mrso_model, xy = TRUE)
  dates <- terra::time(mrso_model)

  names(df_model) <- c("lon", "lat", as.character(dates))

  df_model_long <- df_model %>%
    data.table() %>%
    melt(
      measure.vars = as.character(dates),
      variable.name = "date",
      value.name = "mrso"
    ) %>%
    mutate(
      date = lubridate::date(as.character(date)),
      year = lubridate::year(date)
    )

  # Calculate deltaSM per year per pixel
  delta_SM <- df_model_long %>%
    group_by(lon, lat, year) %>%
    summarise(
      deltaSM = max(mrso, na.rm = TRUE) - min(mrso, na.rm = TRUE)
    ) %>%
    ungroup()

  # In every pixel, calculate the max deltaSM over all years
  delta_SMmax <- delta_SM %>%
    group_by(lon, lat) %>%
    summarise(
      deltaSMmax = max(deltaSM, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    setDT()

  # Calculate absolute deltaSM (maximum amplitude across all years)
  delta_SMabs <- df_model_long %>%
    group_by(lon, lat) %>%
    summarise(
      deltaSMabs = max(mrso, na.rm = TRUE) - min(mrso, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    setDT()

  # Merge the two data.tables based on lon and lat
  merged_data_model <- delta_SMmax[delta_SMabs, nomatch = 0, on = .(lon, lat)]

  # Save the results with model name and scenario in the filename
  saveRDS(merged_data_model, paste0("deltaSM_", model_name, "_", scenario, ".rds"), compress = "xz")

  return(NULL)
}
