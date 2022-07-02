#'@useDynLib numbat
#'@import Rcpp
NULL

## refer to https://github.com/evanbiederstedt/poilogcpp

#' Returns the density for the Poisson lognormal distribution with parameters mu and sig
#' 
#' @param x vector of integers, the observations
#' @param mu mean of lognormal distribution
#' @param sig standard deviation of lognormal distribution
#' @param log boolean Return the log density if TRUE (default=FALSE)
#' @return NULL
#' @keywords internal
dpoilog <- function(x, mu, sig, log=FALSE){
  if (!(length(x) == length(mu) & length(x) == length(sig))) stop('All parameters must be same length') 
  if (any((x[x!=0]/trunc(x[x!=0]))!=1)) stop('all x must be integers')
  if (any(x<0)) stop('one or several values of x are negative')
  if (!all(is.finite(c(mu,sig)))) {
    stop('all parameters should be finite')
  }
  if (any(is.na(c(x,mu,sig)))) stop('Parameters cannot be NA')
  if (any(sig<=0)) {
      stop(c('sig is not larger than 0', unique(sig)))
  }

  p = poilog1(x=as.integer(x), my=as.double(mu), sig=as.double(sig^2))

  p[p == 0] = 1e-15

  if (log) {
    return(log(p))
  } else {
    return(p)
  }
}