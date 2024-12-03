# function based on the segmented package to extract theta_crit and EFmax (Zheng Fu et al.)

fit_bilinear <- function(data, y_var, x_var) {

  y <- data[[y_var]]
  x <- -(data[[x_var]]) # Negative as per the definition of the segmented function

  # Check if y or x contains only NA values
  if (all(is.na(y)) || all(is.na(x))) {
    return(list(error = "y or x contains only NA values"))
  }

  out.lm <- lm(y ~ 1)

  # Try to fit the segmented model
  out.segmented <- tryCatch({
    segmented(out.lm, seg.Z = ~x)
  }, error = function(e) {

    # Return a list with the error if any
    return(list(error = paste("ERROR:", conditionMessage(e))))
  })

  # Check if an error occurred
  if ("error" %in% names(out.segmented)) {
    return(out.segmented) # Return the error
  }

  # Extract results
  theta_crit <- -out.segmented$psi[1, 2] # Setting positive again
  EFmax <- out.segmented$coefficients["(Intercept)"] # Intercept of the flat line
  Slope <- -slope(out.segmented)$x["slope2", "Est."] # Setting positive again
  Intercept <- EFmax - Slope * theta_crit # Intercept with y-axis of the non-flat line

  return(data.table(
    theta_crit = theta_crit,
    EFmax = EFmax,
    Slope = Slope,
    Intercept = Intercept
  ))
}
