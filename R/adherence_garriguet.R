#' Estimate a Participant's Probability of Adhering to "N Days Per Week" Type of 
#' Physical Activity Guideline (Garriguet's Method)
#' 
#' Implements the Bayesian approach described by Garriguet 
#' (\emph{Statistics Canada} 2016).
#' 
#' The approach aims to estimate a participant's probability of meeting 
#' guidelines of the form "at least x minutes per day for at least y days per 
#' week" based on observing X active days out of n monitoring days.
#' 
#' The prior assumption for the participant's \emph{daily} adherence probability 
#' is: 
#' 
#' p_d ~ Beta(alpha, beta)
#' 
#' where alpha and beta are estimated via maximum likelihood using the observed 
#' sample proportions if active days for all study participants. This can be 
#' done separately via \code{\link{mles_beta}}. 
#' 
#' Given p_d, the number of active days out of n monitoring days is distributed:
#' 
#' X|p_d ~ Bin(n, p_d)
#' 
#' It can be shown that the posterior for p_d is:
#' 
#' p_d|X ~ Beta(alpha2 = alpha + X, beta2 = beta + n - X)
#' 
#' Garriguet then uses the Beta-binomial distribution, which describes binomial 
#' data with success probability randomly drawn from Beta(alpha, beta). The 
#' weekly adherence estimator is defined as:
#' 
#' p_w.hat <- P(Y >= 5) with Y ~ Betabin(7, alpha2, beta2)
#' 
#' which can be calculated using \code{\link{mles_beta}}.
#' 
#' 
#' @param n Number of monitoring days.
#' @param x Number of exercise days.
#' @param alpha Parameter in p_d ~ Beta(alpha, beta). Corresponds to 
#' \code{shape1} in \code{\link[stats]{Beta}} functions.
#' @param beta Parameter in p_d ~ Beta(alpha, beta). Corresponds to 
#' \code{shape2} in \code{\link[stats]{Beta}} functions.
#' @param n.rec Denominator for recommendation.
#' @param x.rec Numerator for recommendation.
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
#' Garriguet, D. (2016). Using a betabinomial distribution to estimate the 
#' prevalence of adherence to physical activity guidelines among children and 
#' youth. Statistics Canada, Catalogue no. 82-003-X. Health Reports 27(4): 3-9. 
#' Available at: \url{https://www150.statcan.gc.ca/n1/pub/82-003-x/2016004/article/14489-eng.pdf}.
#' 
#' 
#' @export
adherence_garriguet <- function(n, x, alpha, beta, n.rec = 7, x.rec = 5) {
  pbbinom(q = x.rec - 1, size = n.rec, alpha = alpha + x, beta = beta + n - x, lower.tail = FALSE)
}