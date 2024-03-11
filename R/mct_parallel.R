# mct function adapted to take as input a global dataframe and return a maximum values of CWD around the globe

# libraries
library(data.table)
library(furrr)
library(future)
library(future.apply)

mct_parallel <- function(df, varname_wbal, varname_date, thresh_terminate = 0.0, thresh_drop = 0.9) {

  if (thresh_terminate > thresh_drop)
    rlang::abort("Aborting. thresh_terminate must be smaller or equal thresh_drop.")

  setDT(df)  # Convert dataframe to data.table

  # Initialize columns
  df[, iinst := NA_integer_]
  df[, dday := NA_integer_]
  df[, deficit := 0]

  # Define function for within-group processing
  mct_group <- function(data) {
    idx <- 0
    iinst <- 1

    # Search all dates
    while (idx <= (nrow(data) - 1)) {
      idx <- idx + 1

      # If the water balance (deficit = prec - et) is negative, start accumulating deficit
      if (data[[varname_wbal]][idx] < 0) {
        dday <- 0
        deficit <- 0
        max_deficit <- 0
        iidx <- idx
        done_finding_dropday <- FALSE

        # Continue accumulating deficit as long as the deficit has not fallen below (thresh_terminate) times the maximum deficit attained in this event
        while (iidx <= (nrow(data) - 1) && (deficit - data[[varname_wbal]][iidx] > thresh_terminate * max_deficit)) {
          dday <- dday + 1
          deficit <- deficit - data[[varname_wbal]][iidx]

          # Record the maximum deficit attained in this event
          if (deficit > max_deficit) {
            max_deficit <- deficit
            done_finding_dropday <- FALSE
          }

          # Record the day when deficit falls below (thresh_drop) times the current maximum deficit
          if (deficit < (max_deficit * thresh_drop) && !done_finding_dropday) {
            iidx_drop <- iidx
            done_finding_dropday <- TRUE
          }

          # Once deficit has fallen below threshold, all subsequent dates are dropped (dday set to NA)
          if (done_finding_dropday) {
            data$iinst[iidx] <- NA_integer_
            data$dday[iidx] <- NA_integer_
          } else {
            data$iinst[iidx] <- iinst
            data$dday[iidx] <- dday
            iidx_drop <- iidx
          }

          data$deficit[iidx] <- deficit

          iidx <- iidx + 1
        }

        # Record instance
        this_inst <- data.table(
          idx_start = idx,
          len = iidx_drop - idx,
          iinst = iinst,
          date_start = data[[varname_date]][idx],
          date_end = data[[varname_date]][iidx_drop - 1],
          deficit = max_deficit
        )

        inst <<- rbindlist(list(inst, this_inst))  # Append to inst (using <<- to assign to the global variable)

        # Update
        iinst <- iinst + 1
        dday <- 0
        idx <- iidx
      }
    }

    return(data)
  }

  # Set up future backend for parallel execution
  plan(future::multisession, workers = availableCores()-1)

  # Apply the group-wise function using data.table's fast grouping and parallel processing with furrr
  df <- df[, furrr::future_map(.SD, mct_group), by = c("lon", "lat")]

  return(df)
}
