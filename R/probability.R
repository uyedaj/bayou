#' Conditional Poisson distribution
#' 
#' \code{cdpois} calculates the probability density of a value \code{k} from a Poisson distribution with a maximum \code{kmax}. \code{rdpois} draws random numbers from a conditional Poisson distribution.
#' 
#' @rdname cdpois
#' @param k random variable value
#' @param n number of samples to draw
#' @param kmax maximum value of the conditional Poisson distribution
#' @param log log transformed density
#' @param lambda rate parameter of the Poisson distribution
#' @param ... additional parameters passed to \code{dpois} or \code{rpois}
#' @export
#' @examples
#' cdpois(10,1,10)
#' cdpois(11,1,10)
#' rdpois(5,10,10)
cdpois <- function(k,lambda,kmax,log=TRUE){
  kmax <- ceiling(kmax)
  i <- 0:kmax
  R <- sum(dpois(i,lambda))
  if(k<=kmax){
    num <- dpois(k,lambda)
  } else {num <- 0}
  if(log){
    log(num/R)
  } else {num/R }
}
#' @rdname cdpois
rdpois <- function(n,lambda,kmax,...){
  kmax <- ceiling(kmax)
  i=rep(kmax+1,n)
  j=0
  while(any(i>kmax)){
    i[i>kmax] <- rpois(sum(i>kmax),lambda)
    j <- j+1
    if(j>100){stop ("Lambda too high relative to kmax")}
  }
  return(i)
}

