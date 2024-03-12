
# function to extract Rn, EF and SM from .mat files and create a long datatable for plotting

# df_key must be loaded in the R environment

library(tidyverse)
library(terra)
library(R.matlab)
library(data.table)

process_daily_cmip6 <- function(model_name) {

  # Parent directory where all the daily model data are
  path_prefix <- "data-raw/cmip6-ng/DAILY/"

  # Load model-specific data
  # size is 144x72x5475 (lon x lat x days), 15 years of data (2000-2014)
  Rn_model <- readMat(paste0(path_prefix, "Rn_daily_", model_name, ".mat"))
  EF_model <- readMat(paste0(path_prefix, "EF_daily_", model_name, ".mat"))
  SM_model <- readMat(paste0(path_prefix, "SM1m_daily_", model_name, ".mat"))

  # extract values
  EF_array <- array(EF_model[[1]], # only extract values (ignore attributes of the list, as they will increase the size of the array like crazy)
                    dim = c(144, 72, 5475))
  Rn_array <- array(Rn_model[[1]],
                    dim = c(144, 72, 5475))
  SM_array <- array(SM_model[[1]],
                    dim = c(144, 72, 5475))

  # free up memory
  rm(Rn_model)
  rm(EF_model)
  rm(SM_model)

  # convert to SpatRasters
  str(EF_array) # lon and lat are inverted --> transpose matrix
  EF_array <- aperm(EF_array, c(2, 1, 3)) # switch second and first dimensions, so that it can easily be read with package terra!
  SM_array <- aperm(SM_array, c(2, 1, 3))
  Rn_array <- aperm(Rn_array, c(2, 1, 3))

  # reverse order in the first dimension (the map is flipped on the latitude otherwise)
  EF_array <- EF_array[ # subset based on the reversed order
    rev( # reverse the order
      seq_len( # create a sequence from 1 till the argument (i.e. total number of elements)
        dim(EF_array)[1])), , ] # take the total number of elements in the first dimension
  SM_array <- SM_array[rev(seq_len(dim(SM_array)[1])), , ]
  Rn_array <- Rn_array[rev(seq_len(dim(Rn_array)[1])), , ]

  # Define the extent: xmin, xmax, ymin, ymax
  extent <- ext(-180, 180, -90, 90)

  # create a sequence of dates
  date_seq <- seq(from = as.Date("2000-01-01"), to = as.Date("2014-12-31"), by = "day")
  dates <- date_seq[!format(date_seq, "%m-%d") == "02-29"] # remove 29/02 from leap years to match CMIP6 dates

  # Create an empty SpatRaster with the correct dimensions and extent
  EF_rast <- rast(nrows=180, ncols=360, nlyr=dim(EF_array)[3], ext=extent)
  SM_rast <- rast(nrows=180, ncols=360, nlyr=dim(SM_array)[3], ext=extent)
  Rn_rast <- rast(nrows=180, ncols=360, nlyr=dim(SM_array)[3], ext=extent)

  # Set the resolution to 2.5 degrees
  res(EF_rast) <- c(2.5, 2.5)
  res(SM_rast) <- c(2.5, 2.5)
  res(Rn_rast) <- c(2.5, 2.5)

  # Fill the SpatRaster with data from the matrix
  values(EF_rast) <- EF_array
  values(SM_rast) <- SM_array
  values(Rn_rast) <- Rn_array

  # set dates
  time(EF_rast) <- dates
  time(SM_rast) <- dates
  time(Rn_rast) <- dates

  # Create a SpatVector with all points from df_key
  points <- vect(df_key, geom = c("lon_cmip6", "lat_cmip6"), crs = "EPSG:4326")

  # Extract data for all points at once
  EF_data <- extract(EF_rast, points)
  SM_data <- extract(SM_rast, points)
  Rn_data <- extract(Rn_rast, points)

  # Convert the extracted data to dataframes
  names(EF_data) <- c("ID", as.character(dates)) # rename the columns to the dates
  EF_data <- EF_data %>% # attach sitename list
    bind_cols(df_key %>%
                dplyr::select(sitename)) %>%
    dplyr::select(-ID)
  dt_EF <- melt(EF_data %>% setDT(), # pivot longer
                  id.vars = "sitename",
                  variable.name = "date",
                  value.name = "EF")

  names(SM_data) <- c("ID", as.character(dates)) # rename the columns to the dates

  SM_data <- SM_data %>% # attach sitename list
    bind_cols(df_key %>%
                dplyr::select(sitename)) %>%
    dplyr::select(-ID)
  dt_SM <- melt(SM_data %>% setDT(), # pivot longer
                id.vars = "sitename",
                variable.name = "date",
                value.name = "SM")

  names(Rn_data) <- c("ID", as.character(dates)) # rename the columns to the dates

  Rn_data <- Rn_data %>% # attach sitename list
    bind_cols(df_key %>%
                dplyr::select(sitename)) %>%
    dplyr::select(-ID)
  dt_Rn <- melt(Rn_data %>% setDT(), # pivot longer
                id.vars = "sitename",
                variable.name = "date",
                value.name = "Rn")

  # join the dataframes
  df_combined <- dt_EF[dt_SM, on=.(date, sitename), nomatch=0]  # First left join with dt_SM
  df_combined <- df_combined[dt_Rn, on=.(date, sitename), nomatch=0]  # Second left join with dt_Rn

  # converting date from character to date object and adding model column
  df_combined[, `:=`(date = as.Date(date),
                     model = model_name)]

  # saveRDS(df_combined, "df_combined_UKESM1.rds", compress = "xz")

  return(df_combined)
}
