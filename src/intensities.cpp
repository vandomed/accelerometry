#include <Rcpp.h>
using namespace Rcpp;

//' Physical Activity Intensities
//' 
//' Given a vector of accelerometer count values, calculates time spent in 5 
//' mutually exclusive user-defined intensity levels (typically representing
//' sedentary, light, lifestyle, moderate, and vigorous) as well as the total 
//' counts accumulated in various intensities. Non-wear time should be removed 
//' from \code{counts} before calling \code{intensities} to avoid overestimating 
//' sedentary time.
//' 
//' 
//' @inheritParams artifacts
//' 
//' @param thresh Numeric vector with four cutpoints from which five intensity 
//' ranges are derived. For example, \code{thresh = c(100, 760, 2020, 5999)} 
//' creates: 0-99 = intensity 1; 100-759 = intensity level 2; 760-2019 = 
//' intensity 3; 2020-5998 = intensity 4; >= 5999 = intensity 5.
//' 
//' 
//' @return Integer vector of length 16 in which the first eight values are 
//' minutes in intensities 1, 2, 3, 4, 5, 2-3, 4-5, and 2-5, and the next eight 
//' are counts accumulated during time spent in each of those intensities.
//' 
//' 
//' @examples
//' # Load accelerometer data for first 5 participants in NHANES 2003-2004
//' data(unidata)
//' 
//' # Get data from ID number 21005
//' counts.part1 <- unidata[unidata[, "seqn"] == 21005, "paxinten"]
//' 
//' # Create vector of counts during valid wear time only
//' # counts.part1.weartime <- counts.part1[weartime(counts = counts.part1) == 1]
//' counts.part1.weartime <- counts.part1
//' 
//' # Calculate physical activity intensity variables
//' intensity.variables <- intensities(counts = counts.part1.weartime)
//' 
//' 
//' @export
// [[Rcpp::export]]
IntegerVector intensities(IntegerVector counts, 
                          IntegerVector thresh = IntegerVector::create(100, 760, 2020, 5999)) {
  
  // Get length(counts) and initialize output vector
  int n = counts.size();
  IntegerVector out(16);
  
  // Get intensity cutpoints
  int thresh_1 = thresh[0]; 
  int thresh_2 = thresh[1];
  int thresh_3 = thresh[2];
  int thresh_4 = thresh[3];
  
  // Loop through each count value
  for (int a = 0; a < n; ++a) {
    int counts_a = counts[a];
    if (counts_a < thresh_1) {
      out[0] += 1;
      out[8] += counts_a;
    }
    else if (counts_a >= thresh_1 && counts_a < thresh_1) {
      out[1] += 1;
      out[5] += 1;
      out[7] += 1;
      out[9] += counts_a;
      out[13] += counts_a;
      out[15] += counts_a;
    }
    else if (counts_a >= thresh_2 && counts_a < thresh_3) {
      out[2] += 1;
      out[5] += 1;
      out[7] += 1;
      out[10] += counts_a;
      out[13] += counts_a;
      out[15] += counts_a;
    }
    else if (counts_a >= thresh_3 && counts_a < thresh_4) {
      out[3] += 1;
      out[6] += 1;
      out[7] += 1;
      out[11] += counts_a;
      out[14] += counts_a;
      out[15] += counts_a;
    }
    else if (counts_a >= thresh_4) {
      out[4] += 1;
      out[6] += 1;
      out[7] += 1;
      out[12] += counts_a;
      out[14] += counts_a;
      out[15] += counts_a;
    }
  }
  
  // Return output vector
  return(out);
}