# Function to process data for a given model and scenario

# make sure you only have the realizations you need in the data folder
# the script only looks at the model name and scenario

process_model_cwd <- function(model_name, scenario_type) {

  print(paste("Processing model:", model_name, "Scenario:", scenario_type))

  ### LOAD DATA ###
  # Construct simplified file patterns (only model name and scenario type are required)
  et_pattern <- paste0("evspsbl_mon_", model_name, "_", scenario_type, "_.*\\.nc")
  pr_pattern <- paste0("pr_mon_", model_name, "_", scenario_type, "_.*\\.nc")

  # List files matching the pattern
  et_files <- list.files("data-raw/cmip6-ng/evspsbl/mon/g025/", pattern = et_pattern, full.names = TRUE)
  pr_files <- list.files("data-raw/cmip6-ng/pr/mon/g025/", pattern = pr_pattern, full.names = TRUE)

  # Check for ET files
  if (length(et_files) == 0) {
    warning(paste("ET files not found for model:", model_name, "Scenario:", scenario_type))
    return(NULL)
  }

  # Load ET data
  ET <- rast(et_files[1])

  # Check and load P files
  if (model_name == "CMCC-ESM2") {
    # Special case for CMCC-ESM2: use precipitation from CESM2
    P <- rast("data-raw/cmip6-ng/pr/mon/g025/pr_mon_CESM2_land-hist_r1i1p1f1_g025.nc")
  } else {
    # For other models, check for P files
    if (length(pr_files) == 0) {
      warning(paste("P files not found for model:", model_name, "Scenario:", scenario_type))
      return(NULL)
    }
    # Load P data
    P <- rast(pr_files[1])
  }

  ### CWD CALCULATION ###
  # Calculate water balance
  water_balance_raw <- P - ET

  # Select dates from 1935-01-01 onwards (80 years of data)
  dates <- terra::time(water_balance_raw)
  selection <- which(dates >= as.Date("1935-01-01"))
  water_balance <- subset(water_balance_raw, subset = selection)

  # Rotate to true coordinates (-180;180 instead of 0;360)
  water_balance <- terra::rotate(water_balance)

  # Focus on vegetated land
  vegetated_land <- readRDS("data/land_mask/vegetated_land_mask.rds")
  water_balance <- mask(
    water_balance,
    vegetated_land,
    maskvalues = 0)

  # Transform to dataframe
  df_wb <- terra::as.data.frame(water_balance, xy = TRUE)
  dates <- terra::time(water_balance)

  names(df_wb) <- c("lon", "lat", as.character(dates))

  df_wb_long <- df_wb %>%
    data.table() %>%
    melt(
      measure.vars = as.character(dates),
      variable.name = "date",
      value.name = "water_balance"
    ) %>%
    mutate(date = as.Date(date))

  # Convert water_balance to mm/day
  df_wb_long[, water_balance := water_balance * 86400]  # Convert from kg/m^2/s to mm/day

  # Calculate CWD using mct function
  print(paste("Calculating CWD for model:", model_name, "Scenario:", scenario_type))

  # Group data by lon and lat
  df_grouped_list <- split(df_wb_long, by = c("lon", "lat"))

  # Define function to apply mct to each group
  apply_mct <- function(group_data) {
    result <- tryCatch({
      mct(group_data, "water_balance", "date")
    }, error = function(e) {
      NULL  # Return NULL in case of error
    })
    return(result)
  }

  # Set up parallel processing plan for grid cells
  plan(multisession, workers = future::availableCores() - 1)

  # Parallel computation using future_map
  df_cwd_list <- future_map(df_grouped_list, apply_mct)

  # Revert to sequential plan
  plan(sequential)

  # Remove NULL results
  df_cwd_list <- df_cwd_list[!sapply(df_cwd_list, is.null)]

  if (length(df_cwd_list) == 0) {
    warning(paste("No CWD data calculated for model:", model_name, "Scenario:", scenario_type))
    return(NULL)
  }

  # Combine the results into a single dataframe
  df_cwd_long <- rbindlist(df_cwd_list)

  # Calculate maximum CWD per grid cell
  max_cwd <- df_cwd_long[, .(max_cwd = max(deficit, na.rm = TRUE)), by = .(lon, lat)]
  max_cwd[, `:=`(model = model_name, scenario = scenario_type)]

  # Save the results
  saveRDS(max_cwd, paste0("max_cwd_", model_name, "_", scenario_type, ".rds"), compress = "xz")

  return(max_cwd)
}
