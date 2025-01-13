# script read GLDAS data, average it from daily to monthly and put in the same format as GRACE and GWLS

# Load libraries
library(ncdf4)
library(tidyverse)
library(lubridate)

# List all .nc4 files recursively in directory
file_list <- list.files(
  path = "data-raw/GLDAS-2_CLSM",
  pattern = "\\.nc4$",
  full.names = TRUE,
  recursive = TRUE
)

# Open the first file to explore
nc_data <- nc_open(file_list[1])

# Print the file summary
print(nc_data)

# Check variable names
variables <- names(nc_data$var)
print(variables)
