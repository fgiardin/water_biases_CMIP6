# MCT 2
# supply P - ET

mct <- function(df, varname_wbal, varname_date, thresh_terminate = 0.0, thresh_drop = 0.9){

  if (thresh_terminate > thresh_drop) rlang::abort("Aborting. thresh_terminate must be smaller or equal thresh_drop.")

  inst <- tibble()
  idx <- 0
  iinst <- 1 #CWD event

  df <- df %>%
    ungroup() %>%
    # dplyr::select(date, !!varname_wbal) %>%
    mutate(iinst = NA, dday = NA, deficit = 0)

  ## search all dates
  while (idx <= (nrow(df)-1)){
    idx <- idx + 1

    ## if the water balance (deficit = prec - et) is negative, start accumulating deficit
    ## cumulate negative water balances (deficits)
    if (df[[ varname_wbal ]][idx]<0){

      dday <- 0
      deficit <- 0
      max_deficit <- 0
      iidx <- idx
      done_finding_dropday <- FALSE

      ## continue accumulating deficit as long as the deficit has not fallen below (thresh_terminate) times the maximum deficit attained in this event
      while (iidx <= (nrow(df)-1) && (deficit - df[[ varname_wbal ]][iidx] > thresh_terminate * max_deficit)){
        dday <- dday + 1
        deficit <- deficit - df[[ varname_wbal ]][iidx]

        ## record the maximum deficit attained in this event
        if (deficit > max_deficit){
          max_deficit <- deficit
          done_finding_dropday <- FALSE
        }

        ## record the day when deficit falls below (thresh_drop) times the current maximum deficit
        if (deficit < (max_deficit * thresh_drop) && !done_finding_dropday){
          iidx_drop <- iidx
          done_finding_dropday <- TRUE
        }

        ## once deficit has fallen below threshold, all subsequent dates are dropped (dday set to NA)
        if (done_finding_dropday){
          df$iinst[iidx] <- NA
          df$dday[iidx]  <- NA
        } else {
          df$iinst[iidx] <- iinst
          df$dday[iidx]  <- dday
          iidx_drop <- iidx
        }

        df$deficit[iidx] <- deficit

        iidx <- iidx + 1

      }

      ## record instance
      this_inst <- tibble( idx_start = idx, len = iidx_drop-idx, iinst = iinst, date_start=df[[varname_date]][idx], date_end = df[[varname_date]][iidx_drop-1], deficit = max_deficit )
      inst <- inst %>% bind_rows(this_inst)

      ## update
      iinst <- iinst + 1
      dday <- 0
      idx <- iidx
    }

  }
  return(df=df)
}
