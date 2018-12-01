#' Get Maximum Likelihood Estimates for Beta Distribution
#' 
#' Same idea as \code{fitdistr} function in \strong{MASS}, but has default 
#' starting values and uses \code{\link[stats]{nlminb}} rather than 
#' \code{\link[stats]{optim}}. 
#' 
#' @param x Observations assumed to be iid Beta(alpha, beta).
#' @param start Starting values for alpha and beta.
#' 
#' @examples
#' # Generate data from Beta(1, 2) and get MLE's
#' set.seed(1)
#' x <- rbeta(n = 1000, shape1 = 1, shape2 = 2)
#' mles <- mles_beta(x)
#' mles$par
#' 
#' @export
mles_beta <- function(x, start = c(0.5, 0.5)) {
  ll.f <- function(theta) {
    -sum(dbeta(x = x, shape1 = theta[1], shape2 = theta[2], log = TRUE))
  }
  nlminb(start = start, objective = ll.f, lower = c(0.001, 0.001))
}
