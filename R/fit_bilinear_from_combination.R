
# function based on fit_bilinear designed to iterate over all combinations of models, longitude and latitude

fit_bilinear_from_combination <- function(combination, split_data, y_var, x_var) {

  # to run manually
  # combination <- combinations_dt[3]
  # y_var <- "SIF_over_SIFmax"
  # x_var <- "TWS_norm"

  # Get the column names of the "combination" dataframe
  col_names <- colnames(combination)

  dt <- split_data[[combination$model[[1]]]] # Extracting value from data.table (focus on one model at a time)

  # focus on grid cell / fluxsite
  if ("lon" %in% col_names & "lat" %in% col_names) { # if lon and lat are provided, use those to filter our dt

    subset_data <- dt[lon == combination$lon[[1]] & lat == combination$lat[[1]]] # take all available points at a given location

  } else if ("sitename" %in% col_names) { # instead if a sitename is provided, filter using that
    subset_data <- dt[sitename == combination$sitename[[1]]]

  } else {
    # Generate an error message if none of the specified columns exist
    stop("Error: Neither 'lon' and 'lat' columns nor 'sitename' column exists in combinations_dt")
  }


  # don't fit segmented regression if there's not enough data available at that location
  # trade off between having enough points so that the function can find a meaningful break point and not excluding too many cells on global grid
  if (nrow(subset_data) < 25) { # FYI - different models may have different points with not enough data, because the Rn > 75 filter depends on the Rn of the specific model

    cat("Not enough data for model:", combination$model[[1]], "lat:", combination$lat[[1]], "lon:", combination$lon[[1]], "\n")

    return(data.table(model = combination$model[[1]], lon = combination$lon[[1]], lat = combination$lat[[1]], # return statement will terminate the execution of the function
                      theta_crit = NA, EFmax = NA))
  }



  # run the bilinear algorithm
  results <- fit_bilinear(subset_data, y_var, x_var)

  # Check if results contain an error (propagate error message from fit_bilinear)
  if ("error" %in% names(results)) {
    cat("Error fitting bilinear model for model:", combination$model[[1]], "lat:", combination$lat[[1]], "lon:", combination$lon[[1]], "\n", results$error, "\n")
    return(data.table(model = combination$model[[1]], lon = combination$lon[[1]], lat = combination$lat[[1]],
                      theta_crit = NA, EFmax = NA))
  }


  # adjust results based on input
  if ("lon" %in% col_names & "lat" %in% col_names) { # if lon and lat are provided, use those to filter our dt

    return(data.table(model = combination$model[[1]], lon = combination$lon[[1]], lat = combination$lat[[1]],
                      theta_crit = results$theta_crit, EFmax = results$EFmax, Slope = results$Slope, Intercept = results$Intercept))

  } else if ("sitename" %in% col_names) { # instead if a sitename is provided, filter using that

    return(data.table(model = combination$model[[1]], sitename = combination$sitename[[1]],
                      theta_crit = results$theta_crit, EFmax = results$EFmax, Slope = results$Slope, Intercept = results$Intercept))

  } else {
    # Generate an error message if none of the specified columns exist
    stop("Error: Neither 'lon' and 'lat' columns nor 'sitename' column exists in combinations_dt")
  }

}


