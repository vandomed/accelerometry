#' Estimate a Participant's Probability of Adhering to "N Days Per Week" Type of 
#' Physical Activity Guideline (Dodd's Method)
#' 
#' Implements the Bayesian approach developed by Dodd and used in the landmark 
#' Troiano et al. paper (\emph{MSSE} 2008). 
#' 
#' The approach aims to estimate a participant's probability of meeting 
#' guidelines of the form "at least x minutes per day for at least y days per 
#' week" based on observing X active days out of n monitoring days. We 
#' illustrate here with the "5+ active days per week" guideline that motivated 
#' the approach.
#' 
#' The prior assumption for the participant's \emph{daily} adherence probability 
#' is:
#' 
#' p_d ~ Uni(0, 1)
#' 
#' Given p_d, the number of active days out of n monitoring days is distributed:
#' 
#' X|p_d ~ Bin(n, p_d)
#' 
#' It can be shown that the posterior for p_d is:
#' 
#' p_d|X ~ Beta(X + 1, n - X + 1) 
#' 
#' Under a somewhat questionable independence assumption, the weekly adherence 
#' probability is p_w = P(Y >= 5) with Y ~ Bin(7, p_d). Dodd estimates p_w as:
#' 
#' p_w.hat = P(p_d >= 5/7 | X)
#' 
#' which can be calculated using \code{\link[stats:Beta]{pbeta}}.
#' 
#' In my view, the quantity P(p_d >= 5/7 | X) is not a good estimator for p_w. 
#' Consider what would happen in a really long protocol. The Beta posterior for 
#' p_d would be very tightly centered around the true p_d, and 
#' p_w.hat = P(p_d >= 5/7 | X) would be very close to either 0 or 1 -- not very 
#' close to what we're trying to estimate, p_w.
#' 
#' A solution is to define p_d.hat as the posterior mean, median, or mode, and
#' map that estimate to p_w, i.e. p_w.hat = P(Y >= 5) with Y ~ Bin(7, p_d.hat). 
#' So there is an option for that.
#' 
#' 
#' @param n Number of monitoring days.
#' @param x Number of active days.
#' @param posterior Can be \code{NULL} for original Dodd method or 
#' \code{"mean"} or \code{"median"} for modified version described above.
#' @param n.rec Denominator for recommendation.
#' @param x.rec Numerator for recommendation.
#' @param goal Recommended number of active days out of \code{n}.
#' 
#' 
#' @examples
#' # Generate data from hypothetical study with 1000 subjects, valid days 
#' # randomly sampled from 1-7, and p_d's drawn from Beta(1, 2).
#' set.seed(1)
#' n <- sample(1: 7, size = 1000, replace = TRUE)
#' p_d <- rbeta(n = 1000, shape1 = 1, shape2 = 2)
#' x <- rbinom(n = 1000, size = n, prob = p_d)
#' 
#' # Estimate p_w's using Dodd's method
#' p_w.hat <- adherence_dodd(n = n, x = )
#' 
#' # First step: Estimate (alpha, beta) via maximum likelihood. Have to change 
#' # 0's to 0.01 and 1's to 0.99 to avoid Inf's
#' p_d.hat <- x / n
#' p_d.hat[p_d.hat == 0] <- 0.01
#' p_d.hat[p_d.hat == 1] <- 0.99
#' mles <- mles_beta(x = p_d.hat)
#' 
#' # Estimate each subject's weekly adherence probability
#' p_w.hat <- adherence_garriguet(n = n, x = x, alpha = mles$par[1], beta = mles$par[2])
#' 
#' # Note that the mean p_w.hat differs considerably from the true mean p_w, 
#' # reflecting bias in the estimator.
#' mean(p_w.hat)
#' mean(pbinom(q = 4, size = 7, prob = p_d, lower.tail = FALSE))
#' 
#' 
#' @references
#' Dodd, K. (2008). Estimation of the population prevalence of adherence to 
#' physical activity recommendations based on NHANES accelerometry measurements. 
#' Technical Report. Available at: 
#' \url{https://epi.grants.cancer.gov/nhanes_pam/bayesian_adherence_estimation.pdf}. 
#' Accessed Nov. 13, 2018.
#' 
#' Troiano, R.P., Berrigan, D., Dodd, K.W., Masse, L.C. and McDowell, M. (2008). 
#' Physical activity in the United States measured by accelerometer. Medicine \& 
#' Science in Sports \& Exercise 40(1): 181--188.
#' 
#' 
#' @export
adherence_dodd <- function(n, x, n.rec = 7, x.rec = 5, posterior = NULL) {
  
  # Original Dodd method
  if (is.null(posterior)) {
    return(pbeta(q = goal / 7, shape1 = x + 1, shape2 = n - x + 1, lower.tail = FALSE))
  }
  
  # Modified Dodd method
  alpha.hat <- x + 1
  beta.hat <- n - x + 1
  if (posterior == "mean") {
    return(pbinom(q = x.rec - 1, 
                  size = n.rec, 
                  prob = alpha.hat / (alpha.hat + beta.hat), 
                  lower.tail = FALSE))
  }
  if (posterior == "median") {
    return(pbinom(q = x.rec - 1, 
                  size = n.rec, 
                  prob = qbeta(p = 0.5, shape1 = alpha.hat, shape2 = beta.hat), 
                  lower.tail = FALSE))
  }
  
}
