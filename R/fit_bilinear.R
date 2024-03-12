
# function based on the segmented package to extract theta_crit and EFmax (Zheng Fu et al.)

fit_bilinear <- function(data, y_var, x_var) {

  y <- data[[y_var]]
  x <- -(data[[x_var]]) # negative as per the definition of the segmented function (see documentation)

  out.lm <- lm(y ~ 1)
  out.segmented <- tryCatch({
    segmented(out.lm, seg.Z=~x)

  }, error=function(e) {
    # return a list with the error if any
    return(list(error = paste("ERROR:", conditionMessage(e))))
  })

  theta_crit <- -out.segmented$psi[1,2] # setting positive again
  EFmax <- out.segmented[["coefficients"]][["(Intercept)"]] # intercept with y axis of flat line
  Slope <- -slope(out.segmented)[["x"]]["slope2", "Est."] # setting positive again
  Intercept = EFmax-Slope*theta_crit # intercept with y axis of the non-flat line

  return(data.table(theta_crit = theta_crit, EFmax = EFmax, Slope = Slope, Intercept = Intercept))
}
