#' Process Uniaxial Minute-to-Minute Accelerometer Data
#' 
#' Calculates a variety of physical activity variables based on uniaxial 
#' minute-to-minute accelerometer count values for individual participants. A 
#' data dictionary for the variables returned is available here: 
#' \url{https://sites.google.com/site/danevandomelen/r-package-accelerometry/data-dictionary}.
#' 
#' 
#' @param counts Integer vector with accelerometer count values.
#' 
#' @param steps Integer vector with steps.
#' 
#' @param nci_methods Logical value for whether to set all arguments such as to 
#' replicate the data processing methods used in the NCI's SAS programs. More 
#' specifically: 
#' 
#' \code{valid_days = 4}
#' 
#' \code{valid_week_days = 0}
#' 
#' \code{valid_weekend_days = 0}
#' 
#' \code{int_cuts = c(100, 760, 2020, 5999)}
#' 
#' \code{cpm_nci = TRUE}
#' 
#' \code{days_distinct = TRUE}
#' 
#' \code{nonwear_window = 60}
#' 
#' \code{nonwear_tol = 2}
#' 
#' \code{nonwear_tol.upper = 100}
#' 
#' \code{nonwear_nci = TRUE}
#' 
#' \code{weartime_minimum = 600}
#' 
#' \code{weartime_maximum = 1440}
#' 
#' \code{partialday_minimum = 1440}
#' 
#' \code{active_bout_length = 10}
#' 
#' \code{active_bout_tol = 2}
#' 
#' \code{mvpa_bout_tol_lower = 0}
#' 
#' \code{vig_bout_tol_lower = 0}
#' 
#' \code{active_bout_nci = TRUE}
#' 
#' \code{sed_bout_tol = 0}
#' 
#' \code{sed_bout_tol_maximum = 759}
#' 
#' \code{artifact_thresh = 32767}
#' 
#' \code{artifact_action = 3}
#' 
#' If \code{TRUE}, you can still specify non-default values for \code{brevity}, 
#' \code{weekday_weekend}, and \code{return_form}.
#' 
#' @param start_date Date of first day of monitoring (must be of class 'date'). 
#' Only used to extract day of week, so if day of week is known but date is not, 
#' user can enter any date that corresponds to that day of the week. The default 
#' date corresponds to the first Sunday in January 2014. 
#' 
#' @param start_time Character vector indicating start time for monitoring. If 
#' not specified it is assumed to be 00:00:00, i.e. the very beginning of the 
#' first day.
#' 
#' @param id Numeric value specifying ID number of participant. 
#' 
#' @param brevity Integer value controlling the number of physical activity 
#' variables that are returned. Choices are \code{1} for basic indicators of 
#' physical activity volume, \code{2} for addditional indicators of activity 
#' intensities, activity bouts, sedentary behavior, and peak activity, and 
#' \code{3} for additional hourly count averages.
#' 
#' @param valid_days Integer value specifying minimum number of valid days to 
#' be considered valid for analysis.
#' 
#' @param valid_week_days Integer value specifying minimum number of valid 
#' weekdays to be considered valid for analysis.
#' 
#' @param valid_weekend_days Integer value specifying minimum number of valid 
#' weekend days to be considered valid for analysis.
#' 
#' @param int_cuts Numeric vector with four cutpoints from which five intensity 
#' ranges are derived. For example, \code{int_cuts = c(100, 760, 2020, 5999)} 
#' creates: 0-99 = intensity 1; 100-759 = intensity level 2; 760-2019 = 
#' intensity 3; 2020-5998 = intensity 4; >= 5999 = intensity 5. Intensities 1-5 
#' are typically viewed as sedentary, light, lifestyle, moderate, and vigorous.
#' 
#' @param cpm_nci Logical value for whether to calculate average counts per 
#' minute by dividing average daily counts by average daily weartime, as opposed 
#' to taking the average of each day's counts per minute value. Strongly 
#' recommend leave as \code{FALSE} unless you wish to replicate the NCI's SAS 
#' programs.
#' 
#' @param days_distinct Logical value for whether to treat each day of data as 
#' distinct, as opposed to analyzing the entire monitoring period as one 
#' continuous segment.
#' 
#' @param nonwear_window Integer value specifying minimum length of a non-wear 
#' period.
#' 
#' @param nonwear_tol Integer value specifying tolerance for non-wear algorithm, 
#' i.e. number of minutes with non-zero counts allowed during a non-wear 
#' interval.
#' 
#' @param nonwear_tol_upper Integer value specifying maximum count value for a 
#' minute with non-zero counts during a non-wear interval.
#' 
#' @param nonwear_nci Logical value for whether to use non-wear algorithm from 
#' NCI's SAS programs.
#' 
#' @param weartime_minimum Integer value specifying minimum number of wear time 
#' minutes for a day to be considered valid.
#' 
#' @param weartime.maximum Integer value specifying maximum number of wear time 
#' minutes for a day to be considered valid. The default is 1440, but you may 
#' want to use a lower value (e.g. 1200) if participants were instructed to 
#' remove devices for sleeping, but often did not.
#' 
#' @param partialday_minimum Integer value specifying minimum number of minutes 
#' for a partial day of monitoring to be processed and potentially considered 
#' valid for analysis. Generally applies only to the first and last days of 
#' monitoring, which may not cover full 24-hour periods. Included because some 
#' researchers may prefer to exclude a day that only has data from, say, 1 pm to 
#' midnight, even if sufficient wear time was achieved during that 11-hour time 
#' frame.
#' 
#' @param active_bout_length Integer value specifying minimum length of an 
#' active bout.
#' 
#' @param active_bout_tol Integer value specifying number of minutes with counts 
#' outside the required range to allow during an active bout. If non-zero and 
#' \code{active_bout_nci = FALSE}, specifying non-zero values for 
#' \code{mvpa.bout.tol.lower} and \code{vig.bout.tol.lower} is highly 
#' recommended. Otherwise minutes immediately before and after an active bout 
#' will tend to be classified as part of the bout.
#' 
#' @param mvpa_bout_tol_lower Integer value specifying lower cut-off for count 
#' values outside of required intensity range for an MVPA bout.
#' 
#' @param vig_bout_tol_lower Integer value specifying lower cut-off for count 
#' values outside of required intensity range for a VPA bout.
#' 
#' @param active_bout_nci Logical value for whether to use algorithm from the 
#' NCI's SAS programs for classifying active bouts.
#' 
#' @param sed_bout_tol Integer value specifying number of minutes with counts 
#' outside sedentary range to allow during a sedentary bout.
#' 
#' @param sed_bout_tol_maximum Integer value specifying upper cut-off for count 
#' values outside sedentary range during a sedentary bout.
#' 
#' @param artifact_thresh Integer value specifying the smallest count value that 
#' should be considered an artifact.
#' 
#' @param artifact_action Integer value controlling method of correcting 
#' artifacts. Choices are \code{1} to exclude days with one or more artifacts, 
#' \code{2} to lump artifacts into non-wear time, \code{3} to replace artifacts 
#' with the average of neighboring count values, and \code{4} to take no action
#' (see \bold{Note}).
#' 
#' @param weekday_weekend Logical value for whether to calculate averages for 
#' weekdays and weekend days separately (in addition to all valid days). 
#' 
#' @param return_form Integer value controlling how variables are returned. 
#' Choices are \code{1} for per-person basis (i.e. averages for each 
#' participant across valid days), \code{2} for per-day basis, and \code{3} for 
#' both. 
#' 
#' 
#' @return
#' Numeric matrix or list of two numeric matrices, depending on 
#' \code{return_form}.
#'
#' 
#' @examples
#' # Load accelerometer data for first 5 participants in NHANES 2003-2004
#' data(unidata)
#' 
#' # Get data from ID number 21005
#' id.part1 <- unidata[unidata[, "seqn"] == 21005, "seqn"]
#' counts.part1 <- unidata[unidata[, "seqn"] == 21005, "paxinten"]
#' 
#' # Create vector of all 10-minute moving averages
#' all.movingaves <- movingaves(x = counts.part1, window = 10, integer = TRUE)
#' 
#' # Calculate maximum 10-minute moving average
#' max.movingave <- movingaves(x = counts.part1, window = 10, integer = TRUE, 
#'                             max = TRUE)
#' 
#' 
#' @export
process_uni <-
  function(counts, steps = NULL, 
           nci_methods = FALSE, 
           start_date = as.Date("2014-01-05"), start_time = "00:00:00", 
           id = NULL, 
           brevity = 1, 
           valid_days = 1, valid_week_days = 0, valid_weekend_days = 0, 
           int_cuts = c(100, 760, 2020, 5999), 
           cpm_nci = FALSE, 
           days_distinct = FALSE, 
           nonwear_window = 60, nonwear_tol = 0, nonwear_tol_upper = 99, 
           nonwear_nci = FALSE, 
           weartime_minimum = 600, weartime_maximum = 1440, 
           partialday_minimum = 1440, 
           active_bout_length = 10, active_bout_tol = 0, 
           mvpa_bout_tol_lower = 0, vig_bout_tol_lower = 0, 
           active_bout_nci = FALSE, sed_bout_tol = 0, 
           sed_bout_tol_maximum = int.cuts[2] - 1, 
           artifact_thresh = 25000, artifact_action = 1, 
           weekday_weekend = FALSE, return_form = 2) {
  
  # Get number of minutes of data
  datalength <- length(counts)
  
  # If nci.methods is TRUE, set inputs to replicate data processing done by NCI's SAS programs
  if (nci.methods == TRUE) {
    
    # Set certain inputs to match NCI methods
    valid.days <- 4
    valid.week.days <- 0
    valid.weekend.days <- 0
    int.cuts <- c(100, 760, 2020, 5999)
    cpm.nci <- TRUE
    days.distinct <- TRUE
    nonwear.window <- 60
    nonwear.tol <- 2
    nonwear.tol.upper <- 100
    nonwear.nci <- TRUE
    weartime.minimum <- 600
    weartime.maximum <- 1440
    partialday.minimum <- 1440
    active.bout.length <- 10
    active.bout.tol <- 2
    mvpa.bout.tol.lower <- 0
    vig.bout.tol.lower <- 0
    active.bout.nci <- TRUE
    sed.bout.tol <- 0
    sed.bout.tol.maximum <- 759
    artifact.thresh <- 32767
    artifact.action <- 3
    
  }
  
  # Get start/stop minutes for each day of monitoring
  extratime <- max(1, round(as.numeric(difftime(as.POSIXct(paste(start.date, "24:00:00")), as.POSIXct(paste(start.date, start.time)), units = "mins"))))
  startmins <- 1
  stopmins <- min(extratime, datalength)
  if (stopmins < datalength) {
    startmins <- c(startmins, seq(stopmins+1, datalength, 1440))
    stopmins <- c(stopmins, startmins[2:length(startmins)]+1439)
    stopmins[length(stopmins)] <- datalength
  }
  
  # If id value or vector is provided, get first value
  if (is.null(id)) {
    id <- 1
  } else {
    id <- id[1]
  }
  
  # Calculate number of full days of data
  numdays <- length(startmins)
  
  # Initialize matrix to save daily physical activity variables
  dayvars <- matrix(NA, ncol = 68, nrow = numdays)
  
  # If artifact.action = 3, replace minutes with counts > artifact.thresh with average of surrounding minutes
  if (artifact.action == 3) {
    counts <- accel.artifacts(counts = counts, thresh = artifact.thresh)
  }
  
  # Call weartime.flag function to flag minutes valid for analysis
  wearflag <- accel.weartime(counts = counts,
                             window = nonwear.window,
                             tol = nonwear.tol,
                             tol.upper = nonwear.tol.upper,
                             nci = nonwear.nci,
                             days.distinct = days.distinct)
  
  # If artifact.action = 2, consider minutes with counts >= artifact.thresh as non-weartime
  if (artifact.action == 2) {
    artifact.locs <- which(counts >= artifact.thresh)
    wearflag[artifact.locs] <- 0
    counts[artifact.locs] <- 0
  }
  
  # Identify bouts of MVPA, VPA, and sedentary time
  if (brevity == 2 | brevity == 3) {
    boutedMVPA <- accel.bouts(counts = counts,
                              weartime = wearflag,
                              bout.length = active.bout.length,
                              thresh.lower = int.cuts[3],
                              tol = active.bout.tol,
                              tol.lower = mvpa.bout.tol.lower,
                              nci = active.bout.nci,
                              days.distinct = days.distinct,
                              skipchecks = TRUE)
    boutedvig <- accel.bouts(counts = counts,
                             weartime = wearflag,
                             bout.length = active.bout.length,
                             thresh.lower = int.cuts[4],
                             tol = active.bout.tol,
                             tol.lower = vig.bout.tol.lower,
                             nci = active.bout.nci,
                             days.distinct = days.distinct,
                             skipchecks = TRUE)
    boutedsed10 <- accel.bouts(counts = counts,
                               weartime = wearflag,
                               bout.length = 10,
                               thresh.upper = int.cuts[1]-1,
                               tol = sed.bout.tol,
                               tol.upper = sed.bout.tol.maximum,
                               days.distinct = days.distinct,
                               skipchecks = TRUE)
    boutedsed30 <- accel.bouts(counts = counts,
                               weartime = wearflag,
                               bout.length = 30,
                               thresh.upper = int.cuts[1]-1,
                               tol = sed.bout.tol,
                               tol.upper = sed.bout.tol.maximum,
                               days.distinct = days.distinct,
                               skipchecks = TRUE)
    boutedsed60 <- accel.bouts(counts = counts,
                               weartime = wearflag,
                               bout.length = 60,
                               thresh.upper = int.cuts[1]-1,
                               tol = sed.bout.tol,
                               tol.upper = sed.bout.tol.maximum,
                               days.distinct = days.distinct,
                               skipchecks = TRUE)
  }
  
  # Get day of week
  currentday <- weekdays(start.date-1)
  if (currentday == "Sunday") {
    currentday <- 1
  } else if (currentday == "Monday") {
    currentday <- 2
  } else if (currentday == "Tuesday") {
    currentday <- 3
  } else if (currentday == "Wednesday") {
    currentday <- 4
  } else if (currentday == "Thursday") {
    currentday <- 5
  } else if (currentday == "Friday") {
    currentday <- 6
  } else if (currentday == "Saturday") {
    currentday <- 7
  }
  
  # Loop through accelerometer data for i days
  for (i in 1:numdays) { 
    
    # Update day of week
    currentday <- currentday + 1
    if (currentday == 8) {
      currentday <- 1
    }
    
    # Load accelerometer data from day i    
    day.counts <- counts[startmins[i]:stopmins[i]]
    day.wearflag <- wearflag[startmins[i]:stopmins[i]]
    if (brevity == 2 | brevity == 3) {
      day.boutedMVPA <- boutedMVPA[startmins[i]:stopmins[i]]
      day.boutedvig <- boutedvig[startmins[i]:stopmins[i]]
      day.boutedsed10 <- boutedsed10[startmins[i]:stopmins[i]]
      day.boutedsed30 <- boutedsed30[startmins[i]:stopmins[i]]
      day.boutedsed60 <- boutedsed60[startmins[i]:stopmins[i]]
    }
    if (!is.null(steps)) {
      day.steps <- steps[startmins[i]:stopmins[i]]
    }
    
    # Calculate constants that are used more than once
    daywear <- sum(day.wearflag)
    maxcount <- max(day.counts)
    daylength <- length(day.counts)
    
    # ID number
    dayvars[i,1] <- id
    
    # Day of week
    dayvars[i,2] <- currentday
    
    # Check whether day is valid for analysis; if not, mark as invalid
    if (daywear < weartime.minimum | daywear > weartime.maximum | (artifact.action == 1 & maxcount >= artifact.thresh) |
        daylength < partialday.minimum) {
      dayvars[i,3] <- 0
    } else {
      dayvars[i,3] <- 1
    }
    
    # Minutes of valid wear time
    dayvars[i,4] <- daywear
    
    # Store day.counts[day.wearflag == 1] into its own vector
    day.counts.valid <- day.counts[day.wearflag == 1]
    
    # Total counts during wear time
    dayvars[i,5] <- sum(day.counts.valid)
    
    # Counts per minute - calculated as total counts during wear time divided by wear time
    dayvars[i,6] <- dayvars[i,5]/dayvars[i,4]
    
    # Steps
    if (!is.null(steps)) {
      dayvars[i,7] <- sum(day.steps[day.wearflag == 1])
    }
    
    if (brevity == 2 | brevity == 3) {
      
      # Minutes in various intensity levels
      intensities <- accel.intensities(counts = day.counts.valid, thresh = int.cuts)
      dayvars[i,8:15] <- intensities[1:8]
      
      # Proportions of daily wear time in each intensity level
      dayvars[i,16:23] <- dayvars[i,8:15]/daywear
      
      # Counts accumulated during wear time in each intensity level  
      dayvars[i,24:31] <- intensities[9:16]
      
      # Bouted sedentary time
      dayvars[i,32] <- sum(day.boutedsed10)
      dayvars[i,33] <- sum(day.boutedsed30)
      dayvars[i,34] <- sum(day.boutedsed60)
      
      # Sedentary breaks
      dayvars[i,35] <- accel.sedbreaks(counts = day.counts, weartime = day.wearflag, thresh = int.cuts[1], skipchecks = TRUE)
      
      # Maximum 1-min, 5-min, 10-min, and 30-min count averages
      dayvars[i,36] <- maxcount
      dayvars[i,37] <- movingaves(x = day.counts, window = 5, return.max = TRUE, skipchecks = TRUE)
      dayvars[i,38] <- movingaves(x = day.counts, window = 10, return.max = TRUE, skipchecks = TRUE)
      dayvars[i,39] <- movingaves(x = day.counts, window = 30, return.max = TRUE, skipchecks = TRUE)
      
      # MVPA and vigorous physical activity in >= 10-min bouts
      dayvars[i,42] <- sum(day.boutedMVPA)
      dayvars[i,43] <- sum(day.boutedvig)
      dayvars[i,44] <- sum(dayvars[i,42:43])
      
      if (dayvars[i,42] > 0) {
        dayvars[i,40] <- sum(rle2(day.boutedMVPA)[,1] == 1)
      } else {
        dayvars[i,40] <- 0
      }
      if (dayvars[i,43] > 0) {
        dayvars[i,41] <- sum(rle2(day.boutedMVPA)[,1] == 1)
      } else {
        dayvars[i,41] <- 0
      }
      
      if (brevity == 3) {
        
        # Hourly counts/min averages
        if (daylength == 1440) {
          dayvars[i,45:68] <- blockaves(x = day.counts, window = 60, skipchecks = TRUE)
        }
        
      }
    }
    
  }
  
  # Format matrix of daily physical activity variables
  colnames(dayvars) <- c("id", "day", "valid_day", "valid_min", "counts", "cpm", "steps", "sed_min", "light_min", "life_min", 
                         "mod_min", "vig_min", "lightlife_min", "mvpa_min", "active_min", "sed_percent", "light_percent", 
                         "life_percent", "mod_percent", "vig_percent", "lightlife_percent", "mvpa_percent", "active_percent", 
                         "sed_counts", "light_counts", "life_counts", "mod_counts", "vig_counts", "lightlife_counts", 
                         "mvpa_counts", "active_counts", "sed_bouted_10min", "sed_bouted_30min", "sed_bouted_60min", 
                         "sed_breaks", "max_1min_counts", "max_5min_counts", "max_10min_counts", "max_30min_counts", 
                         "num_mvpa_bouts", "num_vig_bouts", "mvpa_bouted", "vig_bouted", "guideline_min", "cpm_hour1", 
                         "cpm_hour2", "cpm_hour3", "cpm_hour4", "cpm_hour5", "cpm_hour6", "cpm_hour7", "cpm_hour8", "cpm_hour9", 
                         "cpm_hour10", "cpm_hour11", "cpm_hour12", "cpm_hour13", "cpm_hour14", "cpm_hour15", "cpm_hour16", 
                         "cpm_hour17", "cpm_hour18", "cpm_hour19", "cpm_hour20", "cpm_hour21", "cpm_hour22", "cpm_hour23", 
                         "cpm_hour24")
  
  # Drop variables according to brevity setting
  if (brevity == 1) {
    dayvars <- dayvars[,1:7]
  } else if (brevity == 2) {
    dayvars <- dayvars[,1:44]
  }
  
  # Drop steps if NULL
  if (is.null(steps)) {
    dayvars <- dayvars[, -7, drop = FALSE]
  }
  
  # Calculate daily averages
  locs.valid <- which(dayvars[,3] == 1)
  locs.valid.wk <- which(dayvars[,3] == 1 & dayvars[,2] %in% 2:6)
  locs.valid.we <- which(dayvars[,3] == 1 & dayvars[,2] %in% c(1, 7))
  
  # If weekday.weekend is FALSE, just calculate overall averages
  if (weekday.weekend == FALSE) {
    averages <- c(id = id, valid_days = length(locs.valid), include = 0, colMeans(x = dayvars[locs.valid,4:ncol(dayvars), drop = FALSE]))
    if (length(locs.valid) >= valid.days & length(locs.valid.wk) >= valid.week.days & length(locs.valid.we) >= valid.weekend.days) {
      averages[3] <- 1
    }
  } else {
    
    # Otherwise calculate averages by all valid days and by valid weekdays and valid weekend days
    averages <- c(id = id, valid_days = length(locs.valid), valid_week_days = length(locs.valid.wk),
                  valid_weekend_days = length(locs.valid.we), include = 0,
                  colMeans(x = dayvars[locs.valid,4:ncol(dayvars)]),
                  colMeans(x = dayvars[locs.valid.wk,4:ncol(dayvars)]),
                  colMeans(x = dayvars[locs.valid.we,4:ncol(dayvars)]))
    if (length(locs.valid) >= valid.days & length(locs.valid.wk) >= valid.week.days & length(locs.valid.we) >= valid.weekend.days) {
      averages[5] <- 1
    }
    
    # Modify variable names for weekdays and weekends
    numvars <- (length(averages)-5)/3
    names(averages)[(6+numvars):(6+2*numvars-1)] <- paste("wk_", names(averages)[(6+numvars):(6+2*numvars-1)], sep = "")
    names(averages)[(6+2*numvars):(6+3*numvars-1)] <- paste("we_", names(averages)[(6+2*numvars):(6+3*numvars-1)], sep = "")
    
  }
  
  # If cpm.nci is TRUE, re-calculate averages for cpm variables
  if (cpm.nci == TRUE) {
    averages["cpm"] <- averages["counts"]/averages["valid_min"]
  }
  
  # Convert averages to matrix for the hell of it
  averages <- matrix(averages, nrow = 1, dimnames = list(NULL, names(averages)))
  
  
  # Return data frame(s)
  if (return.form == 1) {
    return(averages)
  } else if (return.form == 2) {
    return (dayvars)
  } else if (return.form == 3) {
    retlist <- list(averages = averages, dayvars = dayvars)
    return(retlist)
  }
  }