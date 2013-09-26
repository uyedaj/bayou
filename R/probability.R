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
  if(kmax < lambda) stop("lambda is too high relative to kmax")
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

dsb <- function(sb, ntips=ntips, bmax=1, prob=1, log=TRUE){
  if(any(!(bmax %in% c(0,1,Inf)))) stop("Number of shifts allowed per branch must be 0, 1, or Inf") 
  if(length(bmax)==1) bmax <- rep(bmax, 2*ntips-2)
  if(length(bmax)!=(2*ntips-2)) stop ("bmax not a multiple of the number of branches")
  sbt <- table(sb)
  if(any(sbt > bmax[as.numeric(names(sbt))])){
    dens <- 0
    if(log) return(log(0)) else 0
  } else {
    if(max(bmax)==1){
      if(length(prob)>1) warning("cannot sample unequal probabilities without replacement, assuming equal probabilities for each branch")
      dens <- 1/choose(sum(bmax),sum(sbt))
      if(log) return(log(dens)) else return(dens)
    } else {
      if(any(!(bmax %in% c(0,Inf)))) stop("Cannot sample unequal probabilities without replacement")
      if(length(prob)==1) prob <- rep(1,2*ntips-2)
      if(length(prob)!=2*ntips-2) stop("Number of probabilities provided must equal number of branches")  
      prob[bmax==0] <- 0
      sbp.all <- prob/sum(prob)
      sbp <- c(sbp.all[as.numeric(names(sbt))],1-sum(sbp.all[as.numeric(names(sbt))]))
      if(log) return(dmultinom(c(sbt,0),prob=sbp,log=TRUE)) else return(dmultinom(c(sbt,0),prob=sbp))
    }
  }    
}

rsb <- function(k, ntips=ntips, bmax=1, prob=1, log=TRUE){
  if(any(!(bmax %in% c(0,1,Inf)))) stop("Number of shifts allowed per branch must be 0, 1, or Inf") 
  if(length(bmax)==1) bmax <- rep(bmax, 2*ntips-2)
  if(length(bmax)!=(2*ntips-2)) stop ("bmax not a multiple of the number of branches")
  if(max(bmax)==1){
    if(length(prob)>1) warning("cannot sample unequal probabilities without replacement, assuming equal probabilities for each branch")
    sb <- .sample((1:(2*ntips-2))[bmax==1], k, replace=FALSE)
    return(sb)
  } else {
    if(any(!(bmax %in% c(0,Inf)))) stop("Cannot sample unequal probabilities without replacement")
    if(length(prob)==1) prob <- rep(1,2*ntips-2)
    if(length(prob)!=2*ntips-2) stop("Number of probabilities provided must equal number of branches")  
    prob[bmax==0] <- 0
    sbp.all <- prob/sum(prob)
    sb <- suppressWarnings(.sample((1:(2*ntips-2)), k, prob=sbp.all, replace=TRUE))
    return(sb)
  }
}    
  


dsb.equal <- function(sb,ntips=ntips,log=TRUE){
  if(log){
    return(log(1/choose(2*ntips-2,length(sb))))
  } else {return(1/choose(2*ntips-2,length(sb)))}
}
dsb.equalnotips <- function(sb,ntips=ntips,log=TRUE){
  if(log){
    return(log(1/choose(ntips-2,length(sb))))
  } else {return(1/choose(ntips-2,length(sb)))}
}
dsb.free <- function(sb,ntips=ntips,edge.length=edge.length,log=TRUE){
  dd <- tapply(sb,sb,length)
  nd <- unique(sb)
  pd <- edge.length[nd]/sum(edge.length)
  dmultinom(c(dd,0),prob=c(pd,1-sum(pd)),log=log)
} 

rsb.equal <- function(k, ntips=ntips, exclude.branches=NULL){
  if(!is.null(exclude.branches)){
    br <- (1:(2*ntips-2))[-c(exclude.branches)]
  }
  return(.sample(1:(2*ntips-2),k,replace=FALSE))
}
rsb.equalnotips <- function(k, ntips=ntips, tree=tree, exclude.branches=NULL){
  if(!is.null(exclude.branches)){
    exclude.branches <- c(exclude.branches,tree$edge[,2]<=ntips)
  } else {
    exclude.branches <- which(tree$edge[,2]<=ntips)
  }
  br <- (1:(2*ntips-2))[-c(exclude.branches)]
  return(.sample(br,k,replace=FALSE))
}
rsb.free <- function(k, ntips=ntips, edge.length=edge.length, exclude.branches=NULL){
  if(!is.null(exclude.branches)){
    br <- (1:(2*ntips-2))[-c(exclude.branches)]
    edge.length <- edge.length[-c(exclude.branches)]
  } else {br <- (1:(2*ntips-2))}
  prb <- edge.length/sum(edge.length)
  sb <- suppressWarnings(.sample(br,k,prob=prb, replace=TRUE))
  return(sb)
}

dloc <- function(loc,min=0,max=1,log=TRUE) if(log) return (rep(0,length(loc))) else return(rep(1,length(loc)))

rloc <- function(k,min=0,max=1){
  return(runif(k))
}