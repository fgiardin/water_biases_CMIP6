# Function to extract and process CMIP6 data for a given model and scenario
# to prepare them for the calculation of moisture limitation thresholds

process_moisture_thresholds <- function(model_name, scenario) {
  print(paste("Processing model:", model_name, "Scenario:", scenario))

  # Initialize a list to store data for each variable
  data_list <- list()

  # Loop over each variable and load the data
  for (variable in variables) {

    # Construct file pattern to match the correct file
    file_pattern <- paste0(variable, "_mon_", model_name, "_", scenario, "_.*_g025.nc")
    file_list <- list.files(
      path = paste0("data-raw/cmip6-ng/", variable, "/mon/g025/"),
      pattern = file_pattern,
      full.names = TRUE
    )

    if (length(file_list) == 0) {
      warning(paste("No file found for variable", variable, "model", model_name, "and scenario", scenario))
      return(NULL)
    }

    file_path <- file_list[1]  # Assuming there's only one file per variable, model, and scenario

    # Load the data
    data_raw <- rast(file_path)

    # Filter dates
    dates <- terra::time(data_raw)
    selection <- which(dates >= start_date & dates <= end_date)
    data_model <- subset(data_raw, subset = selection)

    # Rotate to true coordinates (-180;180 instead of 0;360)
    data_model <- terra::rotate(data_model)

    # Focus on vegetated land
    vegetated_land <- readRDS("data/land_mask/vegetated_land_mask.rds")
    data_model <- mask(
        data_model,
        vegetated_land,
        maskvalues = 0)

    # Store the processed data
    data_list[[variable]] <- data_model
  }

  # Check if all variables are loaded
  if (length(data_list) != length(variables)) {
    warning(paste("Missing variables for model", model_name, "and scenario", scenario))
    return(NULL)
  }

  # Calculate Rn and EF using raster operations
  Rn <- (data_list$rsds - data_list$rsus) + (data_list$rlds - data_list$rlus)
  EF <- data_list$hfls / Rn

  # Filter Rn raster to get only values where Rn > 75
  Rn_filtered <- ifel(Rn > 75, Rn, NA)
  EF_filtered <- mask(EF, Rn_filtered)
  mrso_filtered <- mask(data_list$mrso, Rn_filtered)

  # Convert rasters to data frames
  df_Rn <- terra::as.data.frame(Rn_filtered, xy = TRUE)
  df_EF <- terra::as.data.frame(EF_filtered, xy = TRUE)
  df_mrso <- terra::as.data.frame(mrso_filtered, xy = TRUE)

  # Extract dates and rename columns
  dates <- terra::time(Rn_filtered)
  names(df_Rn) <- c("lon", "lat", as.character(dates))
  names(df_EF) <- c("lon", "lat", as.character(dates))
  names(df_mrso) <- c("lon", "lat", as.character(dates))

  # Pivot data frames to long format
  df_Rn_long <- df_Rn %>%
    data.table() %>%
    melt(
      measure.vars = as.character(dates),
      variable.name = "date",
      value.name = "Rn"
    ) %>%
    mutate(
      date = lubridate::date(as.character(date)),
      date = lubridate::floor_date(date, unit = "month")
    ) %>%
    drop_na()

  df_EF_long <- df_EF %>%
    data.table() %>%
    melt(
      measure.vars = as.character(dates),
      variable.name = "date",
      value.name = "EF"
    ) %>%
    mutate(
      date = lubridate::date(as.character(date)),
      date = lubridate::floor_date(date, unit = "month")
    ) %>%
    drop_na()

  df_mrso_long <- df_mrso %>%
    data.table() %>%
    melt(
      measure.vars = as.character(dates),
      variable.name = "date",
      value.name = "mrso"
    ) %>%
    mutate(
      date = lubridate::date(as.character(date)),
      date = lubridate::floor_date(date, unit = "month")
    ) %>%
    drop_na()

  # Combine the data using data.table
  result_data <- df_EF_long[
    df_Rn_long,
    on = .(lon, lat, date),
    nomatch = NA
  ][
    df_mrso_long,
    on = .(lon, lat, date),
    nomatch = NA
  ]

  # Add model and scenario information
  result_data[, `:=`(model = model_name, scenario = scenario)]

  return(result_data)
}
