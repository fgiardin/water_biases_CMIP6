mct_max <- function(df, varname_wbal, varname_date) {
  # make sure that the data is sorted by date in ascending order
  df <- df[order(df[[varname_date]]), ]

  # Calculate the cumulative water deficit (CWD) for the entire timeseries
  df$deficit <- cumsum(df[[varname_wbal]])

  # Find the maximum CWD value over the entire period
  max_cwd <- max(df$deficit, na.rm = TRUE)

  # Return a data.table with the maximum CWD value for the lon-lat combination
  result <- data.table(max_cwd = max_cwd)
  return(result)
}
