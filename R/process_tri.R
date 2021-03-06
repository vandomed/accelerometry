#' Process Triaxial Minute-to-Minute Accelerometer Data
#' 
#' Calculates a variety of physical activity variables based on triaxial 
#' minute-to-minute accelerometer count values for individual participants. 
#' Assumes first 1440 minutes are day 1, next 1440 are day 2, and so on. If 
#' final day has less than 1440 minutes, it is excluded. A data dictionary for 
#' the variables created is available here: 
#' \url{https://github.com/vandomed/accelerometry/blob/master/process_tri_dictionary.csv}.
#' 
#' 
#' @param counts Integer matrix with three columns of count values, e.g. 
#' vertical-axis counts, anteroposterior (AP)-axis counts, and 
#' mediolateral (ML)-axis counts. 
#' 
#' @param steps Integer vector with steps.
#' 
#' @param nci_methods Logical value for whether to set all arguments so as to 
#' replicate the data processing methods used in the NCI's SAS programs. More 
#' specifically: 
#' 
#' \code{valid_days = 4}
#' 
#' \code{valid_wk_days = 0}
#' 
#' \code{valid_we_days = 0}
#' 
#' \code{int_axis = "vert"}
#' 
#' \code{int_cuts = c(100, 760, 2020, 5999)}
#' 
#' \code{cpm_nci = TRUE}
#' 
#' \code{days_distinct = TRUE}
#' 
#' \code{nonwear_axis = "vert"}
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
#' to record, if \code{brevity = 3}. Choices are "counts_vert", "counts_ap", 
#' "counts_ml", "counts_sum", "counts_vm", "cpm_vert", "cpm_ap", "cpm_ml", 
#' "sed_min", "sed_bouted_10min", and "sed_breaks". 
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
#' @param valid_wk_days Integer value specifying minimum number of valid 
#' weekdays to be considered valid for analysis.
#' 
#' @param valid_we_days Integer value specifying minimum number of valid weekend 
#' days to be considered valid for analysis.
#' 
#' @param int_axis Character string specifying which axis should be used to 
#' classify intensities. Choices are "vert", "ap", "ml", "sum" (for triaxial 
#' sum), and "vm (for triaxial vector magnitude).
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
#' @param nonwear_axis Character string specifying which axis should be used to 
#' classify non-wear time. Choices are "vert", "ap", "ml", "sum" (for triaxial 
#' sum), and "vm" (for triaxial vector magnitude).
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
#' values outside of required intensity range for a vigorous bout.
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
#' @param artifact_axis Character string specifying which axis should be used to 
#' identify artifacts (impossibly high count values). Choices are "vert", "ap", 
#' "ml", "sum" (for triaxial sum), and "vm" (for triaxial vector magnitude).
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
#' @param return_form Character string controlling how variables are returned. 
#' Choices are "daily" for per-day summaries, "averages" for averages across 
#' all valid days, and "both" for a list containing both.
#' 
#' 
#' @return
#' Numeric matrix or list of two numeric matrices, depending on 
#' \code{return_form}.
#' 
#' 
#' @references
#' National Cancer Institute. Risk factor monitoring and methods: SAS programs 
#' for analyzing NHANES 2003-2004 accelerometer data. Available at: 
#' \url{http://riskfactor.cancer.gov/tools/nhanes_pam}. Accessed Aug. 19, 2018.
#' 
#' 
#' @examples
#' # Note that the 'tridata' dataset contains 7 days of fake triaxial 
#' # accelerometer data
#' 
#' # Process data using default parameters and request per-day variables
#' accel.days <- process_tri(
#'   counts = tridata, 
#'   return_form = "daily"
#' )
#' 
#' # Repeat, but request averages across all valid days
#' accel.averages <- process_tri(
#'   counts = tridata,
#'   return_form = "averages"
#' )
#' 
#' # Create per-day summary again, but with many more variables
#' accel.days2 <- process_tri(
#'   counts = tridata, 
#'   brevity = 2,
#'   return_form = "daily"
#' )
#' names(accel.days2)
#' 
#' 
#' @export
process_tri <- function(counts, 
                        steps = NULL, 
                        nci_methods = FALSE, 
                        start_day = 1, 
                        start_date = NULL, 
                        id = NULL, 
                        brevity = 1, 
                        hourly_var = "cpm", 
                        hourly_wearmin = 0, 
                        hourly_normalize = FALSE, 
                        valid_days = 1, 
                        valid_wk_days = 0, 
                        valid_we_days = 0, 
                        int_axis = "vert", 
                        int_cuts = c(100, 760, 2020, 5999), 
                        cpm_nci = FALSE, 
                        days_distinct = FALSE, 
                        nonwear_axis = "vert", 
                        nonwear_window = 60, 
                        nonwear_tol = 0, 
                        nonwear_tol_upper = 99, 
                        nonwear_nci = FALSE, 
                        weartime_minimum = 600, 
                        weartime_maximum = 1440, 
                        active_bout_length = 10, 
                        active_bout_tol = 0, 
                        mvpa_bout_tol_lower = 0, 
                        vig_bout_tol_lower = 0, 
                        active_bout_nci = FALSE, 
                        sed_bout_tol = 0, 
                        sed_bout_tol_maximum = int_cuts[2] - 1, 
                        artifact_axis = "vert", 
                        artifact_thresh = 25000, 
                        artifact_action = 1, 
                        weekday_weekend = FALSE, 
                        return_form = "daily") {
  
  # If requested, set inputs to mimic NCI's SAS programs
  if (nci_methods) {
    
    valid_days <- 4
    valid_wk_days <- 0
    valid_we_days <- 0
    int_axis <- "vert"
    int_cuts <- c(100, 760, 2020, 5999)
    cpm_nci <- TRUE
    days_distinct <- TRUE
    nonwear_axis <- "vert"
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
    artifact_axis <- "vert"
    artifact_thresh <- 32767
    artifact_action <- 3
    
  }
  
  # Get start/end indices for each day of monitoring
  n.minutes <- nrow(counts)
  end.indices <- seq(1440, n.minutes, 1440)
  start.indices <- end.indices - 1439
  n.days <- length(start.indices)
  
  # Get single value for ID
  id <- ifelse(is.null(id), 1, id[1])
  
  # Calculate acceptable range for wear time
  weartime.range <- c(weartime_minimum, weartime_maximum)
  
  # Create variable for whether there steps is specified
  have.steps <- ! is.null(steps)
  
  # Separate out counts for each axis, calculate sum and vector magnitude, and 
  # create 5-column matrix with all 5 measures
  counts.vert <- counts[, 1]
  counts.ap <- counts[, 2]
  counts.ml <- counts[, 3]
  counts.sum <- counts.vert + counts.ap + counts.ml
  counts.vm <- sqrt(counts.vert^2 + counts.ap^2 + counts.ml^2)
  
  # Determine which counts vector to use for artifact detection
  if (artifact_axis == "vert") {
    counts.artifacts <- counts.vert
  } else if (artifact_axis == "ap") {
    counts.artifacts <- counts.ap
  } else if (artifact_axis == "ml") {
    counts.artifacts <- counts.ml
  } else if (artifact_axis == "sum") {
    counts.artifacts <- counts.sum
  } else if (artifact_axis == "vm") {
    counts.artifacts <- counts.vm
  }
  
  # If artifact_action = 3, replace minutes with counts >= artifact_thresh with 
  # average of surrounding minutes
  max_artifact <- max(counts.artifacts)
  if (artifact_action == 3 && max_artifact >= artifact_thresh) {
    counts.vert <- artifacts(counts = counts.vert, 
                             thresh = artifact_thresh, 
                             counts_classify = counts.artifacts)
    counts.ap <- artifacts(counts = counts.ap, 
                           thresh = artifact_thresh, 
                           counts_classify = counts.artifacts)
    counts.ml <- artifacts(counts = counts.ml, 
                           thresh = artifact_thresh, 
                           counts_classify = counts.artifacts)
    counts.sum <- artifacts(counts = counts.sum, 
                            thresh = artifact_thresh, 
                            counts_classify = counts.artifacts)
    counts.vm <- artifacts(counts = counts.vm, 
                           thresh = artifact_thresh, 
                           counts_classify = counts.artifacts)
  }
  
  # Determine which counts vector to use for wear-time detection
  if (nonwear_axis == "vert") {
    counts.nonwear <- counts.vert
  } else if (nonwear_axis == "ap") {
    counts.nonwear <- counts.ap
  } else if (nonwear_axis == "ml") {
    counts.nonwear <- counts.ml
  } else if (nonwear_axis == "sum") {
    counts.nonwear <- counts.sum
  } else if (nonwear_axis == "vm") {
    counts.nonwear <- counts.vm
  }
  
  # Call weartime function to flag minutes valid for analysis
  wearflag <- weartime(
    counts = counts.nonwear,
    window = nonwear_window,
    tol = nonwear_tol,
    tol_upper = nonwear_tol_upper,
    nci = nonwear_nci,
    days_distinct = days_distinct
  )
  
  # If artifact_action = 2, consider minutes with counts >= artifact_thresh  
  # non-wear time
  if (artifact_action == 2 && max_artifact >= artifact_thresh) {
    artifact.locs <- which(counts.artifacts >= artifact_thresh)
    wearflag[artifact.locs] <- 0
    counts.vert[artifact.locs, ] <- 0
    counts.ap[artifact.locs, ] <- 0
    counts.ml[artifact.locs, ] <- 0
    counts.sum[artifact.locs, ] <- 0
    counts.vm[artifact.locs, ] <- 0
  }
  
  # Determine which counts vector to use for intensities and activity bouts
  if (int_axis == "vert") {
    counts.int <- counts.vert
  } else if (int_axis == "ap") {
    counts.int <- counts.ap
  } else if (int_axis == "ml") {
    counts.int <- counts.ml
  } else if (int_axis == "sum") {
    counts.int <- counts.sum
  } else if (int_axis == "vm") {
    counts.int <- counts.vm
  }
  
  # Identify MVPA, vigorous, and sedentary bouts
  if (brevity %in% c(2, 3)) {
    
    bouted.MVPA <- bouts(
      counts = counts.int, 
      weartime = wearflag, 
      bout_length = active_bout_length, 
      thresh_lower = int_cuts[3],
      tol = active_bout_tol, 
      tol_lower = mvpa_bout_tol_lower,
      nci = active_bout_nci, 
      days_distinct = days_distinct
    )
    
    bouted.vig <- bouts(
      counts = counts.int, 
      weartime = wearflag, 
      bout_length = active_bout_length, 
      thresh_lower = int_cuts[4],
      tol = active_bout_tol, 
      tol_lower = vig_bout_tol_lower,
      nci = active_bout_nci, 
      days_distinct = days_distinct
    )
    
    # Note - may want to add a sed_axis option in future. May want to use VM
    # to classifiy sedentary time, for example.
    bouted.sed10 <- bouts(
      counts = counts.int, 
      weartime = wearflag, 
      bout_length = 10,
      thresh_upper = int_cuts[1] - 1, 
      tol = sed_bout_tol, 
      tol_upper = sed_bout_tol_maximum, 
      days_distinct = days_distinct
    )
    
    bouted.sed30 <- bouts(
      counts = counts.int, 
      weartime = wearflag, 
      bout_length = 30,
      thresh_upper = int_cuts[1] - 1, 
      tol = sed_bout_tol, 
      tol_upper = sed_bout_tol_maximum, 
      days_distinct = days_distinct
    )
    
    bouted.sed60 <- bouts(
      counts = counts.int, 
      weartime = wearflag, 
      bout_length = 60,
      thresh_upper = int_cuts[1] - 1, 
      tol = sed_bout_tol, 
      tol_upper = sed_bout_tol_maximum, 
      days_distinct = days_distinct
    )
    
  }
  
  # Get days of week (1 = Sunday, ..., 7 = Saturday)
  if (! is.null(start_date)) {
    days <- as.POSIXlt(seq(start_date, start_date + n.days - 1, 1))$wday + 1
  } else {
    days <- start_day: (start_day + n.days - 1)
    days <- ifelse(days > 7, days - 7, days)
  }
  
  # Initialize matrix to save daily physical activity variables
  if (brevity == 1) {
    day.vars <- matrix(NA, nrow = n.days, ncol = 15)
  } else if (brevity == 2) {
    day.vars <- matrix(NA, nrow = n.days, ncol = 100)
  } else {
    day.vars <- matrix(NA, nrow = n.days, ncol = 124)
  }
  
  # Re-assemble counts matrix as vert, AP, ML, sum, vm
  counts <- cbind(counts.vert, counts.ap, counts.ml, counts.sum, counts.vm)
  
  # Loop through each day of data
  for (ii in 1: n.days) { 
    
    # Load data for iith day
    daylocs.ii <- start.indices[ii]: end.indices[ii]
    counts.ii <- counts[daylocs.ii, ]
    
    wearflag.ii <- wearflag[daylocs.ii]
    wearlocs.ii <- which(wearflag.ii == 1)
    wcounts.ii <- counts.ii[wearlocs.ii, ]
    
    counts.artifacts.ii <- counts.artifacts[daylocs.ii]
    
    if (brevity %in% c(2, 3)) {
      bouted.MVPA.ii <- bouted.MVPA[daylocs.ii]
      bouted.vig.ii <- bouted.vig[daylocs.ii]
      bouted.sed10.ii <- bouted.sed10[daylocs.ii]
      bouted.sed30.ii <- bouted.sed30[daylocs.ii]
      bouted.sed60.ii <- bouted.sed60[daylocs.ii]
    }
    if (have.steps) {
      steps.ii <- steps[daylocs.ii]
    }
    
    # Calculate constants that are used more than once
    sum_wearflag.ii <- sum(wearflag.ii)
    max_artifact.ii <- max(counts.artifacts.ii)
    
    # ID number and day of week
    day.vars[ii, 1: 2] <- c(id, days[ii])
    
    # Valid day
    day.vars[ii, 3] <- ifelse(inside(sum_wearflag.ii, weartime.range), 1, 0)
    if (artifact_action == 1) {
      if (max(counts.artifacts[daylocs.ii]) >= artifact_thresh) {
        day.vars[ii, 3] <- 0
      }
    }
    
    # Wear time minutes
    day.vars[ii, 4] <- sum_wearflag.ii
    
    # Total counts during wear time
    sum_wcounts.ii <- apply(wcounts.ii, 2, sum)
    day.vars[ii, 5: 9] <- sum_wcounts.ii
    
    # Counts per minute
    day.vars[ii, 10: 14] <- sum_wcounts.ii / sum_wearflag.ii
    
    # Steps
    if (have.steps) {
      day.vars[ii, 15] <- sum(steps.ii[wearlocs.ii])
    }
    
    if (brevity %in% c(2, 3)) {
      
      # Minutes in each intensity
      counts.int.ii <- counts.int[daylocs.ii]
      wcounts.int.ii <- counts.int.ii[wearlocs.ii]
      intensities.ii <- intensities(counts = wcounts.int.ii, int_cuts = int_cuts)
      intensities.min.ii <- intensities.ii[1: 8]
      day.vars[ii, 16: 23] <- intensities.min.ii
      
      # Percentage of wear time in each intensity
      day.vars[ii, 24: 31] <- intensities.min.ii / sum_wearflag.ii
      
      # Counts accumulated from each intensity
      sedlocs.ii <- which(wcounts.int.ii < int_cuts[1])
      day.vars[ii, 32: 36] <- apply(wcounts.ii[sedlocs.ii, , drop = FALSE], 2, sum)
      
      lightlocs.ii <- which(wcounts.int.ii >= int_cuts[1] & wcounts.int.ii < int_cuts[2])
      day.vars[ii, 37: 41] <- apply(wcounts.ii[lightlocs.ii, , drop = FALSE], 2, sum)
      
      lifelocs.ii <- which(wcounts.int.ii >= int_cuts[2] & wcounts.int.ii < int_cuts[3])
      day.vars[ii, 42: 46] <- apply(wcounts.ii[lifelocs.ii, , drop = FALSE], 2, sum)
      
      modlocs.ii <- which(wcounts.int.ii >= int_cuts[3] & wcounts.int.ii < int_cuts[4])
      day.vars[ii, 47: 51] <- apply(wcounts.ii[modlocs.ii, , drop = FALSE], 2, sum)
      
      viglocs.ii <- which(wcounts.int.ii >= int_cuts[4])
      day.vars[ii, 52: 56] <- apply(wcounts.ii[viglocs.ii, , drop = FALSE], 2, sum)
      
      lightlifelocs.ii <- c(lightlocs.ii, lifelocs.ii)
      day.vars[ii, 57: 61] <- apply(wcounts.ii[lightlifelocs.ii, , drop = FALSE], 2, sum)
      
      mvpalocs.ii <- c(modlocs.ii, viglocs.ii)
      day.vars[ii, 62: 66] <- apply(wcounts.ii[mvpalocs.ii, , drop = FALSE], 2, sum)
      
      activelocs.ii <- c(lightlifelocs.ii, mvpalocs.ii)
      day.vars[ii, 67: 71] <- apply(wcounts.ii[activelocs.ii, , drop = FALSE], 2, sum)
      
      # Bouted sedentary time
      day.vars[ii, 72: 74] <- c(sum(bouted.sed10.ii), 
                                sum(bouted.sed30.ii), 
                                sum(bouted.sed60.ii))
      
      # Sedentary breaks
      day.vars[ii, 75] <- sedbreaks(counts = counts.int.ii, 
                                    weartime = wearflag.ii, 
                                    thresh = int_cuts[1])
      
      # Maximum 1-min, 5-min, 10-min, and 30-min count average in each axis
      day.vars[ii, 76: 80] <- apply(counts.ii, 2, max)
      day.vars[ii, 81: 85] <- apply(counts.ii, 2, function(x) movingaves(x = x, window = 5, max = TRUE))
      day.vars[ii, 86: 90] <- apply(counts.ii, 2, function(x) movingaves(x = x, window = 10, max = TRUE))
      day.vars[ii, 91: 95] <- apply(counts.ii, 2, function(x) movingaves(x = x, window = 30, max = TRUE))
      
      # Minutes in MVPA and vigorous bouts
      sum_bouted.MVPA.ii <- sum(bouted.MVPA.ii)
      sum_bouted.vig.ii <- sum(bouted.vig.ii)
      day.vars[ii, 96] <- sum_bouted.MVPA.ii
      day.vars[ii, 97] <- sum_bouted.vig.ii
      
      # Guideline minutes
      day.vars[ii, 98] <- sum_bouted.MVPA.ii + sum_bouted.vig.ii
      
      # Number of MVPA and vigorous bouts
      if (sum_bouted.MVPA.ii > 0) {
        day.vars[ii, 99] <- sum(rle2(bouted.MVPA.ii)[, 1] == 1)
      } else {
        day.vars[ii, 99] <- 0
      }
      if (sum_bouted.vig.ii > 0) {
        day.vars[ii, 100] <- sum(rle2(bouted.vig.ii)[, 1] == 1)
      } else {
        day.vars[ii, 100] <- 0
      }
      
      if (brevity == 3) {
        
        # Hourly variable - first analyze all minutes
        if (hourly_var == "counts_vert") {
          hourly.ii <- blocksums(x = counts.ii[, 1] * wearflag.ii, window = 60)
        } else if (hourly_var == "counts_ap") {
          hourly.ii <- blocksums(x = counts.ii[, 2] * wearflag.ii, window = 60)
        } else if (hourly_var == "counts_ml") {
          hourly.ii <- blocksums(x = counts.ii[, 3] * wearflag.ii, window = 60)
        } else if (hourly_var == "counts_sum") {
          hourly.ii <- blocksums(x = counts.ii[, 4] * wearflag.ii, window = 60)
        } else if (hourly_var == "counts_vm") {
          hourly.ii <- blocksums(x = counts.ii[, 5] * wearflag.ii, window = 60)
        } else if (hourly_var == "cpm_vert") {
          hourly.ii <- blockaves(x = counts.ii[, 1] * wearflag.ii, window = 60)
        } else if (hourly_var == "cpm_ap") {
          hourly.ii <- blockaves(x = counts.ii[, 2] * wearflag.ii, window = 60)
        } else if (hourly_var == "cpm_ml") {
          hourly.ii <- blockaves(x = counts.ii[, 3] * wearflag.ii, window = 60)
        } else if (hourly_var == "cpm_sum") {
          hourly.ii <- blockaves(x = counts.ii[, 4] * wearflag.ii, window = 60)
        } else if (hourly_var == "cpm_vm") {
          hourly.ii <- blockaves(x = counts.ii[, 5] * wearflag.ii, window = 60)
        } else if (hourly_var == "sed_min") {
          hourly.ii <- blocksums(x = counts.int.ii < int_cuts[1], window = 60)
        } else if (hourly_var == "sed_bouted_10min") {
          hourly.ii <- blocksums(x = bouted.sed10.ii, window = 60)
        } else if (hourly_var == "sed_breaks") {
          sedbreaks.ii <- sedbreaks(counts = counts.int.ii, 
                                    weartime = wearflag.ii, 
                                    flags = TRUE)
          hourly.ii <- blocksums(x = sedbreaks.ii, window = 60)
        }
        
        # Calculate hourly wear time if necessary
        if (hourly_normalize || hourly_wearmin > 0) {
          hourly.weartime <- blocksums(x = wearflag.ii, window = 60)
        }
        
        # Normalize by wear time if requested
        if (hourly_normalize) {
          hourly.ii <- hourly.ii / hourly.weartime
        }
        
        # Replace with NAs where wear time is insufficient
        if (hourly_wearmin > 0) {
          hourly.ii <- hourly.ii * ifelse(hourly.weartime >= hourly_wearmin, 1, 
                                          NA)
        }
        day.vars[ii, 101: 124] <- hourly.ii
        
      }
    }
    
  }
  
  # Format day.vars
  axis.labels <- c("vert", "ap", "ml", "sum", "vm")
  if (brevity == 1) {
  
    colnames(day.vars) <- c("id", "day", "valid_day", "valid_min", 
                            paste("counts", axis.labels, sep = "_"), 
                            paste("cpm", axis.labels, sep = "_"), 
                            "steps")
    
  } else if (brevity == 2) {
    
    int.labels <- c("sed", "light", "life", "mod", "vig", "lightlife", "mvpa", "active")
    colnames(day.vars) <- 
      c("id", "day", "valid_day", "valid_min", 
        paste("counts", axis.labels, sep = "_"), 
        paste("cpm", axis.labels, sep = "_"), 
        "steps", 
        paste(int.labels, "min", sep = "_"), 
        paste(int.labels, "percent", sep = "_"), 
        paste(rep(int.labels, each = 5), "counts", axis.labels, sep = "_"), 
        paste("sed_bouted_", c(10, 30, 60), "min", sep = ""),
        "sed_breaks", 
        paste("max_", rep(c(1, 5, 10, 30), each = 5), "min_", axis.labels, sep = ""), 
        "num_mvpa_bouts", "num_vig_bouts", "mvpa_bouted", "vig_bouted", 
        "guideline_min")
    
  } else if (brevity == 3) {
    
    int.labels <- c("sed", "light", "life", "mod", "vig", "lightlife", "mvpa", "active")
    colnames(day.vars) <- 
      c("id", "day", "valid_day", "valid_min", 
        paste("counts", axis.labels, sep = "_"), 
        paste("cpm", axis.labels, sep = "_"), 
        "steps", 
        paste(int.labels, "min", sep = "_"), 
        paste(int.labels, "percent", sep = "_"), 
        paste(rep(int.labels, each = 5), "counts", axis.labels, sep = "_"), 
        paste("sed_bouted_", c(10, 30, 60), "min", sep = ""),
        "sed_breaks", 
        paste("max_", rep(c(1, 5, 10, 30), each = 5), "min_", axis.labels, sep = ""), 
        "num_mvpa_bouts", "num_vig_bouts", "mvpa_bouted", "vig_bouted", 
        "guideline_min", 
        paste(hourly_var, "_hour", 1: 24, sep = ""))
    
  }
  
  # # Drop steps if NULL (Removed because it makes process_nhanes code easier)
  # if (! have.steps) {
  #   day.vars <- day.vars[, -7, drop = FALSE]
  # }
  
  # Calculate daily averages if requested
  if (return_form %in% c("averages", "both")) {
    
    locs.valid <- which(day.vars[, 3] == 1)
    locs.valid.wk <- which(day.vars[, 3] == 1 & day.vars[, 2] %in% 2: 6)
    locs.valid.we <- which(day.vars[, 3] == 1 & day.vars[, 2] %in% c(1, 7))
    n.valid <- length(locs.valid)
    n.valid.wk <- length(locs.valid.wk)
    n.valid.we <- length(locs.valid.we)
    
    if (! weekday_weekend) {
      averages <- c(id = id, valid_days = n.valid, valid_wk_days = n.valid.wk, 
                    valid_we_days = n.valid.we, 
                    include = ifelse(n.valid >= valid_days && 
                                       n.valid.wk >= valid_wk_days && 
                                       n.valid.we >= valid_we_days, 1, 0), 
                    colMeans(x = day.vars[locs.valid, -c(1: 3), drop = FALSE], na.rm = TRUE))
    } else {
      averages <- c(id = id, valid_days = n.valid, valid_wk_days = n.valid.wk, 
                    valid_we_days = n.valid.we, 
                    include = ifelse(n.valid >= valid_days && 
                                       n.valid.wk >= valid_wk_days && 
                                       n.valid.we >= valid_we_days, 1, 0), 
                    colMeans(x = day.vars[locs.valid, -c(1: 3), drop = FALSE], na.rm = TRUE), 
                    colMeans(x = day.vars[locs.valid.wk, -c(1: 3), drop = FALSE], na.rm = TRUE), 
                    colMeans(x = day.vars[locs.valid.we, -c(1: 3), drop = FALSE], na.rm = TRUE))
      cols <- colnames(day.vars)[-c(1: 3)]
      names(averages) <- c(names(averages)[1: 5], cols, 
                           paste("wk_", cols, sep = ""), 
                           paste("we_", cols, sep = ""))
    }
    
    # If cpm_nci is TRUE, re-calculate daily averages for counts per minute
    if (cpm_nci) {
      averages["cpm_vert"] <- averages["counts_vert"] / averages["valid_min"]
      averages["cpm_ap"] <- averages["counts_ap"] / averages["valid_min"]
      averages["cpm_ml"] <- averages["counts_ml"] / averages["valid_min"]
      averages["cpm_sum"] <- averages["counts_sum"] / averages["valid_min"]
      averages["cpm_vm"] <- averages["counts_vm"] / averages["valid_min"]
      
      if (weekday_weekend) {
        averages["wk_cpm_vert"] <- averages["wk_counts_vert"] / averages["wk_valid_min"]
        averages["wk_cpm_ap"] <- averages["wk_counts_ap"] / averages["wk_valid_min"]
        averages["wk_cpm_ml"] <- averages["wk_counts_ml"] / averages["wk_valid_min"]
        averages["wk_cpm_sum"] <- averages["wk_counts_sum"] / averages["wk_valid_min"]
        averages["wk_cpm_vm"] <- averages["wk_counts_vm"] / averages["wk_valid_min"]
        
        averages["we_cpm_vert"] <- averages["we_counts_vert"] / averages["we_valid_min"]
        averages["we_cpm_ap"] <- averages["we_counts_ap"] / averages["we_valid_min"] 
        averages["we_cpm_ml"] <- averages["we_counts_ml"] / averages["we_valid_min"] 
        averages["we_cpm_sum"] <- averages["we_counts_sum"] / averages["we_valid_min"] 
        averages["we_cpm_vm"] <- averages["we_counts_vm"] / averages["we_valid_min"] 
      }
    }
    
    # Convert averages to matrix for the hell of it
    averages <- matrix(averages, nrow = 1, dimnames = list(NULL, names(averages)))
    
  }
  
  # Return data frame(s)
  if (return_form == "averages") {
    return(averages)
  } else if (return_form == "daily") {
    return (day.vars)
  } else {
    return(list(averages = averages, day.vars = day.vars))
  }
}
