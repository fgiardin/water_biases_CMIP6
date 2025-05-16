# function based on the segmented package to extract theta_crit and EFmax (Zheng Fu et al.)

fit_bilinear <- function(data, y_var, x_var) {

  # 1) ensure 'segmented' is installed
  if (!requireNamespace("segmented", quietly = TRUE)) {
    stop(
      "The 'segmented' package is not installed.\n",
      "Please install it with:\n",
      "  install.packages('segmented')"
    )
  }
  # 2) ensure it's loaded
  if (!"segmented" %in% .packages()) {
    stop(
      "The 'segmented' package is installed but not loaded.\n",
      "Please load it before calling fit_bilinear():\n",
      "  library(segmented)"
    )
  }

  y <- data[[y_var]]
  x <- -(data[[x_var]]) # negative as per the definition of the segmented function (see documentation)

  out.lm <- lm(y ~ 1)
  out.segmented <- tryCatch({
    segmented(out.lm, seg.Z=~x)

  }, error=function(e) {
    # return a list with the error if any
    return(list(error = paste("ERROR:", conditionMessage(e))))
  })

  # # detect if segmented() failed and exit (uncomment if you have an error and the function can't run, this lets you see the error)
  # if (!is.null(out.segmented$error)) {
  #   return(data.table(
  #     theta_crit = NA_real_,
  #     EFmax      = NA_real_,
  #     Slope      = NA_real_,
  #     Intercept  = NA_real_,
  #     error      = out.segmented$error
  #   ))
  # }

  theta_crit <- -out.segmented$psi[1,2] # setting positive again
  EFmax <- out.segmented[["coefficients"]][["(Intercept)"]] # intercept with y axis of flat line
  Slope <- -slope(out.segmented)[["x"]]["slope2", "Est."] # setting positive again
  Intercept = EFmax-Slope*theta_crit # intercept with y axis of the non-flat line

  # ggplot(data, aes(x = SM, y = EF)) +
  #   geom_density_2d_filled(contour_var = "ndensity", alpha = 0.8) + # Adjust alpha for transparency
  #   geom_hline(yintercept = EFmax, color = "red", size = 2) +       # Horizontal line for EFmax
  #   geom_abline(intercept = EFmax - Slope * theta_crit,
  #               slope = Slope, color = "blue", size = 2) +          # Line with Slope
  #   labs(x = "Soil Moisture (SM)", y = "Evapotranspiration Fraction (EF)") +
  #   ylim(0, 1.2) +                                                  # Set y-axis limits
  #   theme_minimal()                                                 # A minimal theme

  return(data.table(theta_crit = theta_crit, EFmax = EFmax, Slope = Slope, Intercept = Intercept))
}
