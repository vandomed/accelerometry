#' Process Uniaxial Minute-to-Minute Accelerometer Data
#' 
#' Calculates a variety of physical activity variables based on uniaxial 
#' minute-to-minute accelerometer count values for individual participants. 
#' Assumes first 1440 minutes are day 1, next 1440 are day 2, and so on. If 
#' final day has less than 1440 minutes, it is excluded. A data dictionary for 
#' the variables created is available here: 
#' \url{https://vandomed.github.io/process_uni_dictionary}.
#' 
#' 
#' @param counts Integer vector with accelerometer count values.
#' 
#' @param steps Integer vector with steps.
#' 
#' @param nci_methods Logical value for whether to set all arguments so as to 
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
#' \code{nonwear_tolupper = 100}
#' 
#' \code{nonwear_nci = TRUE}
#' 
#' \code{weartime_minimum = 600}
#' 
#' \code{weartime_maximum = 1440}
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
#' If \code{TRUE}, you can still specify non-default values for \code{brevity} 
#' and \code{weekday_weekend}.
#' 
#' @param start_day Integer value specifying day of week for first day of 
#' monitoring, with 1 = Sunday, ..., 7 = Satuday.
#' 
#' @param start_date Date for first day of monitoring, which function can use to 
#' figure out \code{start_day}. 
#' 
#' @param id Numeric value specifying ID number of participant. 
#' 
#' @param brevity Integer value controlling the number of physical activity 
#' variables generated. Choices are 1 for basic indicators of physical activity 
#' volume, 2 for addditional indicators of activity intensities, activity bouts, 
#' sedentary behavior, and peak activity, and 3 for additional hourly count 
#' averages.
#' 
#' @param hourly_var Character string specifying what hourly activity variable 
#' to record, if \code{brevity = 3}. Choices are "counts", "cpm", "sed_min", 
#' "sed_bouted_10min", and "sed_breaks". 
#' 
#' @param hourly_wearmin Integer value specifying minimum number of wear time 
#' minutes needed during a given hour to record a value for the hourly activity 
#' variable.
#' 
#' @param hourly_normalize Logical value for whether to normalize hourly 
#' activity by number of wear time minutes.
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
#' minute by dividing average daily counts by average daily wear time, as 
#' opposed to taking the average of each day's counts per minute value. Strongly 
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
#' @param weartime_maximum Integer value specifying maximum number of wear time 
#' minutes for a day to be considered valid. The default is 1440, but you may 
#' want to use a lower value (e.g. 1200) if participants were instructed to 
#' remove devices for sleeping, but often did not.
#' 
#' @param active_bout_length Integer value specifying minimum length of an 
#' active bout.
#' 
#' @param active_bout_tol Integer value specifying number of minutes with counts 
#' outside the required range to allow during an active bout. If non-zero and 
#' \code{active_bout_nci = FALSE}, specifying non-zero values for 
#' \code{mvpa_bout_tol_lower} and \code{vig_bout_tol_lower} is highly 
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
#' artifacts. Choices are 1 to exclude days with one or more artifacts, 2 to 
#' lump artifacts into non-wear time, 3 to replace artifacts with the average of 
#' neighboring count values, and 4 to take no action.
#' 
#' @param weekday_weekend Logical value for whether to calculate averages for 
#' weekdays and weekend days separately (in addition to all valid days). 
#' 
#' @param return_form Integer value controlling how variables are returned. 
#' Choices are 1 for per-person basis (i.e. averages for each participant across 
#' valid days), 2 for per-day basis, and 3 for both. 
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
#' # Process data from ID 21005 and request per-day variables
#' accel.days <- process_uni(counts = counts.part1, id = id.part1)
#' 
#' # Process data from ID 21005 and request daily averages
#' accel.averages <- process_uni(counts = counts.part1, id = id.part1, 
#'                               return_form = 1)
#'                               
#' # Process data from ID 21005 and request per-day variables and daily averages
#' accel.list <- process_uni(counts = counts.part1, id = id.part1, 
#'                           return_form = 3)
#'                           
#' # Process data according to methods used in NCI's SAS programs
#' accel.nci1 <- process_uni(counts = counts.part1, id = id.part1, brevity = 2, 
#'                           valid_days = 4, cpm_nci = TRUE, days_distinct = TRUE, 
#'                           nonwear_tol = 2, nonwear_tol_upper = 100, 
#'                           nonwear_nci = TRUE, weartime_maximum = 1440, 
#'                           active_bout_tol = 2, active_bout_nci = TRUE, 
#'                           artifact_thresh = 32767, artifact_action = 3, 
#'                           return_form = 3)
#'                           
#' # Repeat, but use nci_methods input for convenience
#' accel.nci2 <- process_uni(counts = counts.part1, id = id.part1, 
#'                           nci_methods = TRUE, brevity = 2, return_form = 3)
#'                           
#' # Results are identical
#' all(accel.nci1[[1]] == accel.nci2[[1]])
#' all(accel.nci1[[2]] == accel.nci2[[2]])
#' 
#' 
#' @export
process_uni <- function(counts, steps = NULL, 
                        nci_methods = FALSE, 
                        start_day = 1, start_date = NULL, 
                        id = NULL, 
                        brevity = 1, 
                        hourly_var = "cpm", hourly_wearmin = 0, 
                        hourly_normalize = FALSE, 
                        valid_days = 1, valid_week_days = 0, valid_weekend_days = 0, 
                        int_cuts = c(100, 760, 2020, 5999), 
                        cpm_nci = FALSE, 
                        days_distinct = FALSE, 
                        nonwear_window = 60, nonwear_tol = 0, nonwear_tol_upper = 99, 
                        nonwear_nci = FALSE, 
                        weartime_minimum = 600, weartime_maximum = 1440, 
                        active_bout_length = 10, active_bout_tol = 0, 
                        mvpa_bout_tol_lower = 0, vig_bout_tol_lower = 0, 
                        active_bout_nci = FALSE, sed_bout_tol = 0, 
                        sed_bout_tol_maximum = int_cuts[2] - 1, 
                        artifact_thresh = 25000, artifact_action = 1, 
                        weekday_weekend = FALSE, return_form = 2) {
  
  # If requested, set inputs to mimic NCI's SAS programs
  if (nci_methods) {
    
    valid_days <- 4
    valid_week_days <- 0
    valid_weekend_days <- 0
    int_cuts <- c(100, 760, 2020, 5999)
    cpm_nci <- TRUE
    days_distinct <- TRUE
    nonwear_window <- 60
    nonwear_tol <- 2
    nonwear_tol_upper <- 100
    nonwear_nci <- TRUE
    weartime_minimum <- 600
    weartime_maximum <- 1440
    active_bout_length <- 10
    active_bout_tol <- 2
    mvpa_bout_tol_lower <- 0
    vig_bout_tol_lower <- 0
    active_bout_nci <- TRUE
    sed_bout_tol <- 0
    sed_bout_tol_maximum <- 759
    artifact_thresh <- 32767
    artifact_action <- 3
    
  }
    
  # Get start/end indices for each day of monitoring
  n.minutes <- length(counts)
  end.indices <- seq(1440, n.minutes, 1440)
  start.indices <- end.indices - 1439
  n.days <- length(start.indices)
  
  # Get single value for ID
  id <- ifelse(is.null(id), 1, id[1])
  
  # Calculate acceptable range for wear time
  weartime.range <- c(weartime_minimum, weartime_maximum)
  
  # Create variable for whether there steps is specified
  have.steps <- ! is.null(steps)
  
  # Initialize matrix to save daily physical activity variables
  day.vars <- matrix(NA, ncol = 68, nrow = n.days)
  
  # If artifact_action = 3, replace minutes with counts >= artifact_thresh with 
  # average of surrounding minutes
  if (artifact_action == 3) {
    counts <- artifacts(counts = counts, thresh = artifact_thresh)
  }
  
  # Call weartime function to flag minutes valid for analysis
  wearflag <- weartime(counts = counts,
                       window = nonwear_window,
                       tol = nonwear_tol,
                       tol_upper = nonwear_tol_upper,
                       nci = nonwear_nci,
                       days_distinct = days_distinct)
  
  # If artifact_action = 2, consider minutes with counts >= artifact_thresh  
  # non-wear time
  if (artifact_action == 2) {
    artifact.locs <- which(counts >= artifact_thresh)
    wearflag[artifact.locs] <- 0
    counts[artifact.locs] <- 0
  }
  
  # Identify bouts of MVPA, VPA, and sedentary time
  if (brevity %in% c(2, 3)) {
    
    bouted.MVPA <- 
      bouts(counts = counts, weartime = wearflag, 
            bout_length = active_bout_length, thresh_lower = int_cuts[3],
            tol = active_bout_tol, tol_lower = mvpa_bout_tol_lower,
            nci = active_bout_nci, days_distinct = days_distinct)
    
    bouted.VPA <- 
      bouts(counts = counts, weartime = wearflag, 
            bout_length = active_bout_length, thresh_lower = int_cuts[4],
            tol = active_bout_tol, tol_lower = vig_bout_tol_lower,
            nci = active_bout_nci, days_distinct = days_distinct)
    
    bouted.sed10 <- 
      bouts(counts = counts, weartime = wearflag, bout_length = 10,
            thresh_upper = int_cuts[1] - 1, tol = sed_bout_tol, 
            tol_upper = sed_bout_tol_maximum, days_distinct = days_distinct)
    
    bouted.sed30 <- 
      bouts(counts = counts, weartime = wearflag, bout_length = 30,
            thresh_upper = int_cuts[1] - 1, tol = sed_bout_tol, 
            tol_upper = sed_bout_tol_maximum, days_distinct = days_distinct)
    
    bouted.sed60 <- 
      bouts(counts = counts, weartime = wearflag, bout_length = 60,
            thresh_upper = int_cuts[1] - 1, tol = sed_bout_tol, 
            tol_upper = sed_bout_tol_maximum, days_distinct = days_distinct)

  }
  
  # Get days of week (1 = Sunday, ..., 7 = Saturday)
  if (! is.null(start_date)) {
    days <- as.POSIXlt(seq(start_date, start_date + n.days - 1, 1))$wday + 1
  } else {
    days <- start_day: (start_day + n.days - 1)
    days <- ifelse(days > 7, days - 7, days)
  }
  
  # Loop through each day of data
  for (ii in 1: n.days) { 
    
    # Day of week
    day.ii <- days[ii]
    
    # Load data for iith day
    start.ii <- start.indices[ii]
    end.ii <- end.indices[ii]
    counts.ii <- counts[start.ii: end.ii]
    wearflag.ii <- wearflag[start.ii: end.ii]
    if (brevity %in% c(2, 3)) {
      bouted.MVPA.ii <- bouted.MVPA[start.ii: end.ii]
      bouted.VPA.ii <- bouted.VPA[start.ii: end.ii]
      bouted.sed10.ii <- bouted.sed10[start.ii: end.ii]
      bouted.sed30.ii <- bouted.sed30[start.ii: end.ii]
      bouted.sed60.ii <- bouted.sed60[start.ii: end.ii]
    }
    if (! is.null(steps)) {
      steps.ii <- steps[start.ii: end.ii]
    }
    
    # Calculate constants that are used more than once
    sum_wearflag.ii <- sum(wearflag.ii)
    max_counts.ii <- max(counts.ii)
    
    # ID number and day of week
    day.vars[ii, 1: 2] <- c(id, day.ii)
    
    # Valid day
    day.vars[ii, 3] <- 
      ifelse(! inside(sum_wearflag.ii, weartime.range) || 
               (artifact_action == 1 & max_counts.ii >= artifact_thresh), 0, 1)
    
    # Wear time minutes
    day.vars[ii, 4] <- sum_wearflag.ii
    
    # Create vector of counts during wear time
    wearcounts.ii <- counts.ii[wearflag.ii == 1]
    
    # Total counts during wear time
    sum_wearcounts.ii <- sum(wearcounts.ii)
    day.vars[ii, 5] <- sum_wearcounts.ii
    
    # Counts per minute
    day.vars[ii, 6] <- sum_wearcounts.ii / sum_wearflag.ii
    
    # Steps
    if (have.steps) {
      day.vars[ii, 7] <- sum(steps.ii[wearflag.ii == 1])
    }
    
    if (brevity %in% c(2, 3)) {
      
      # Minutes in each intensity
      intensities.ii <- intensities(counts = wearcounts.ii, int_cuts = int_cuts)
      intensities.min.ii <- intensities.ii[1: 8]
      day.vars[ii, 8: 15] <- intensities.min.ii
      
      # Percentage of wear time in each intensity
      day.vars[ii, 16: 23] <- intensities.min.ii / sum_wearflag.ii
      
      # Counts accumulated from each intensity
      day.vars[ii, 24: 31] <- intensities.ii[9: 16]
      
      # Bouted sedentary time
      day.vars[ii, 32: 34] <- c(sum(bouted.sed10.ii), 
                                sum(bouted.sed30.ii), 
                                sum(bouted.sed60.ii))
      
      # Sedentary breaks
      day.vars[ii, 35] <- sedbreaks(counts = counts.ii, weartime = wearflag.ii, 
                                    thresh = int_cuts[1])
      
      # Maximum 1-min, 5-min, 10-min, and 30-min count average
      day.vars[ii, 36] <- max_counts.ii
      day.vars[ii, 37] <- movingaves(x = counts.ii, window = 5, max = TRUE)
      day.vars[ii, 38] <- movingaves(x = counts.ii, window = 10, max = TRUE)
      day.vars[ii, 39] <- movingaves(x = counts.ii, window = 30, max = TRUE)
      
      # Minutes of bouted MVPA and bouted VPA
      sum_bouted.MVPA.ii <- sum(bouted.MVPA.ii)
      sum_bouted.VPA.ii <- sum(bouted.VPA.ii)
      day.vars[ii, 42] <- sum_bouted.MVPA.ii
      day.vars[ii, 43] <- sum_bouted.VPA.ii
      
      # Guideline minutes
      day.vars[ii, 44] <- sum_bouted.MVPA.ii + sum_bouted.VPA.ii
      
      # Number of MVPA and VPA bouts
      if (sum_bouted.MVPA.ii > 0) {
        day.vars[ii, 40] <- sum(rle2(bouted.MVPA.ii)[, 1] == 1)
      } else {
        day.vars[ii, 40] <- 0
      }
      if (sum_bouted.VPA.ii > 0) {
        day.vars[ii, 41] <- sum(rle2(bouted.VPA.ii)[, 1] == 1)
      } else {
        day.vars[ii, 41] <- 0
      }
      
      if (brevity == 3) {
        
        # Hourly variable - first analyze all minutes
        if (hourly_var == "counts") {
          day.vars[ii, 45: 68] <- blocksums(x = counts.ii * wearflag.ii, 
                                            window = 60)
        } else if (hourly_var == "cpm") {
          day.vars[ii, 45: 68] <- blockaves(x = counts.ii * wearflag.ii, 
                                            window = 60)
        } else if (hourly_var == "sed_min") {
          day.vars[ii, 45: 68] <- blocksums(x = counts.ii < int_cuts[1], 
                                            window = 60)
        } else if (hourly_var == "sed_bouted_10min") {
          day.vars[ii, 45: 68] <- blocksums(x = bouted.sed10.ii, window = 60)
        } else if (hourly_var == "sed_breaks") {
          sedbreaks.ii <- sedbreaks(counts = counts.ii, weartime = wearflag.ii, 
                                    flags = TRUE)
          day.vars[ii, 45: 68] <- blocksums(x = sedbreaks.ii, window = 60)
        }
        
        # Calculate hourly wear time if necessary
        if (hourly_normalize || hourly_wearmin > 0) {
          hourly.weartime <- blocksums(x = wearflag.ii, window = 60)
        }
        
        # Normalize for wear time if requested
        if (hourly_normalize) {
          day.vars[ii, 45: 68] <- day.vars[ii, 45: 68] / hourly.weartime
        }
        
        # Replace with NAs where wear time is insufficient
        if (hourly_wearmin > 0) {
          day.vars[ii, 45: 68] <- day.vars[ii, 45: 68] * 
            ifelse(hourly.weartime >= hourly_wearmin, 1, NA)
        }
        
      }
    }
    
  }
  
  # Format day.vars
  if (brevity == 1) {
    
    day.vars <- day.vars[, 1: 6]
    colnames(day.vars) <- 
      c("id", "day", "valid_day", "valid_min", "counts", "cpm")
  
  } else if (brevity == 2) {
    
    day.vars <- day.vars[, 1: 44]
    colnames(day.vars) <- 
      c("id", "day", "valid_day", "valid_min", "counts", "cpm", "steps", 
        "sed_min", "light_min", "life_min", "mod_min", "vig_min", 
        "lightlife_min", "mvpa_min", "active_min", "sed_percent", 
        "light_percent", "life_percent", "mod_percent", "vig_percent", 
        "lightlife_percent", "mvpa_percent", "active_percent", "sed_counts", 
        "light_counts", "life_counts", "mod_counts", "vig_counts", 
        "lightlife_counts", "mvpa_counts", "active_counts", "sed_bouted_10min", 
        "sed_bouted_30min", "sed_bouted_60min", "sed_breaks", "max_1min_counts", 
        "max_5min_counts", "max_10min_counts", "max_30min_counts", 
        "num_mvpa_bouts", "num_vig_bouts", "mvpa_bouted", "vig_bouted", 
        "guideline_min")
    
  } else if (brevity == 3) {
    
    colnames(day.vars) <- 
      c("id", "day", "valid_day", "valid_min", "counts", "cpm", "steps", 
        "sed_min", "light_min", "life_min", "mod_min", "vig_min", 
        "lightlife_min", "mvpa_min", "active_min", "sed_percent", 
        "light_percent", "life_percent", "mod_percent", "vig_percent", 
        "lightlife_percent", "mvpa_percent", "active_percent", "sed_counts", 
        "light_counts", "life_counts", "mod_counts", "vig_counts", 
        "lightlife_counts", "mvpa_counts", "active_counts", "sed_bouted_10min", 
        "sed_bouted_30min", "sed_bouted_60min", "sed_breaks", "max_1min_counts", 
        "max_5min_counts", "max_10min_counts", "max_30min_counts", 
        "num_mvpa_bouts", "num_vig_bouts", "mvpa_bouted", "vig_bouted", 
        "guideline_min", paste(hourly_var, "_hour", 1: 24, sep = ""))
    
  }
  
  # Drop steps if NULL
  if (! have.steps) {
    day.vars <- day.vars[, -7, drop = FALSE]
  }
  
  # Calculate daily averages
  locs.valid <- which(day.vars[, 3] == 1)
  locs.valid.wk <- which(day.vars[, 3] == 1 & day.vars[, 2] %in% 2: 6)
  locs.valid.we <- which(day.vars[, 3] == 1 & day.vars[, 2] %in% c(1, 7))
  n.valid <- length(locs.valid)
  n.valid.wk <- length(locs.valid.wk)
  n.valid.we <- length(locs.valid.we)
  
  if (! weekday_weekend) {
    
    # Just calculate averages across all valid days
    averages <- 
      c(id = id, valid_days = n.valid, include = 0, 
        colMeans(x = day.vars[locs.valid, 4: ncol(day.vars), drop = FALSE], 
                 na.rm = TRUE))
    
    if (n.valid >= valid_days &&  n.valid.wk >= valid_week_days && 
        n.valid.we >= valid_weekend_days) {
      averages[3] <- 1
    }
    
    # If cpm_nci is TRUE, re-calculate daily average for counts per minute
    if (cpm_nci) {
      averages["cpm"] <- averages["counts"] / averages["valid_min"]
    }
    
  } else {
    
    # Calculate averages for all valid days and for weekdays/weekends
    ncol_day.vars <- ncol(day.vars)
    averages <- 
      c(id = id, valid_days = n.valid, valid_week_days = n.valid.wk, 
        valid_weekend_days = n.valid.we, include = 0, 
        colMeans(x = day.vars[locs.valid, 4: ncol_day.vars, drop = FALSE], na.rm = TRUE),
        colMeans(x = day.vars[locs.valid.wk, 4: ncol_day.vars, drop = FALSE], na.rm = TRUE),
        colMeans(x = day.vars[locs.valid.we, 4: ncol_day.vars, drop = FALSE], na.rm = TRUE))
    
    if (n.valid >= valid_days && n.valid.wk >= valid_week_days && 
        n.valid.we >= valid_weekend_days) {
      averages[5] <- 1
      
    }
    
    # Modify variable names for weekdays and weekends
    n.vars <- (length(averages) - 5) / 3
    
    indices <- (6 + n.vars): (6 + 2 * n.vars - 1)
    var.names <- names(averages)[indices]
    names(averages)[indices] <- paste("wk_", var.names, sep = "")
    
    indices2 <- (6 + 2 * n.vars): (6 + 3 * n.vars - 1)
    names(averages)[indices2] <- paste("we_", var.names, sep = "")
    
    # If cpm_nci is TRUE, re-calculate daily average for counts per minute
    if (cpm_nci) {
      averages["cpm"] <- averages["counts"] / averages["valid_min"]
      averages["wk_cpm"] <- averages["wk_counts"] / averages["wk_valid_min"]
      averages["we_cpm"] <- averages["we_counts"] / averages["we_valid_min"]
    }
    
  }
  
  # Convert averages to matrix for the hell of it
  averages <- matrix(averages, nrow = 1, dimnames = list(NULL, names(averages)))
  
  
  # Return data frame(s)
  if (return_form == 1) {
    return(averages)
  } else if (return_form == 2) {
    return (day.vars)
  } else if (return_form == 3) {
    return(list(averages = averages, day.vars = day.vars))
  }
}