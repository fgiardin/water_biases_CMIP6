# function based on fit_bilinear designed to iterate over all combinations of models, longitude and latitude

fit_bilinear_from_combination <- function(combination, split_data, y_var, x_var) {

  # Get the column names of the "combination" dataframe
  col_names <- colnames(combination)

  # Create the key for split_data using model and scenario
  key <- paste(combination$model[[1]], combination$scenario[[1]], sep = ".")

  # Check if the key exists in split_data
  if (!key %in% names(split_data)) {
    cat("No data found for key:", key, "\n")
    return(NULL)
  }

  dt <- split_data[[key]] # Extract the data table for the model and scenario

  # Focus on grid cell / flux site
  if ("lon" %in% col_names & "lat" %in% col_names) { # If lon and lat are provided, use those to filter our dt

    subset_data <- dt[lon == combination$lon[[1]] & lat == combination$lat[[1]]] # Take all available points at a given location

  } else if ("sitename" %in% col_names) { # Instead if a sitename is provided, filter using that
    subset_data <- dt[sitename == combination$sitename[[1]]]

  } else {
    # Generate an error message if none of the specified columns exist
    stop("Error: Neither 'lon' and 'lat' columns nor 'sitename' column exists in combinations_dt")
  }

  # Don't fit segmented regression if there's not enough data available at that location
  if (nrow(subset_data) < 25) {
    cat("Not enough data for model:", combination$model[[1]],
        "scenario:", combination$scenario[[1]],
        "lat:", combination$lat[[1]],
        "lon:", combination$lon[[1]], "\n")
    return(data.table(
      model = combination$model[[1]],
      scenario = combination$scenario[[1]],
      lon = combination$lon[[1]],
      lat = combination$lat[[1]],
      theta_crit = NA,
      EFmax = NA
    ))
  }

  # Run the bilinear algorithm
  results <- fit_bilinear(subset_data, y_var, x_var)

  # Check if results contain an error
  if ("error" %in% names(results)) {
    cat("Error fitting bilinear model for model:", combination$model[[1]], "scenario:", combination$scenario[[1]], "lat:", combination$lat[[1]], "lon:", combination$lon[[1]], "\n", results$error, "\n")
    return(data.table(
      model = combination$model[[1]],
      scenario = combination$scenario[[1]],
      lon = combination$lon[[1]],
      lat = combination$lat[[1]],
      theta_crit = NA,
      EFmax = NA
    ))
  }

  # Return the results
  return(data.table(
    model = combination$model[[1]],
    scenario = combination$scenario[[1]],
    lon = combination$lon[[1]],
    lat = combination$lat[[1]],
    theta_crit = results$theta_crit,
    EFmax = results$EFmax,
    Slope = results$Slope,
    Intercept = results$Intercept
  ))
}
