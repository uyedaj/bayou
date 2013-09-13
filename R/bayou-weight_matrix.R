#' Calculate the weight matrix of a set of regimes on a phylogeny
#' 
#' These functions calculate weight matrices from regimes specified in phytools' simmap format.
#' \code{simmap.W} calculates the weight matrix for a set of regimes from a phylogeny
#' with a stored regime history. \code{.simmap.W} calculates the same matrix, but without checks and is 
#' generally run internally. 
#' 
#' @rdname simmap.W
#' @param tree either a tree of class "phylo" or a cache object produced by bayOU's internal 
#' functions. Must include list element 'maps' which is a simmap reconstruction of regime history.
#' @param pars a list of the parameters used to calculate the weight matrix. Only pars$alpha is
#' necessary to calculate the matrix, but others can be present.
#' @param TotExp a parameter dependent on the total tree height and the alpha parameter that determines
#' the influence of the root state. 
#' 
#' @details \code{.simmap.W} is more computationally efficient within a mcmc and is used internally. The value
#' of \code{TotExp} is supplied to speed computation and reduce redundancy, and cache objects must be supplied as
#' the phylogeny, and the parameter \code{ntheta} must be present in the list \code{pars}.
simmap.W <- function(tree,pars){
  if(class(tree)=="phylo"){
    X <- rep(NA,length(tree$tip.label))
    names(X) <- tree$tip.label
    cache <- .prepare.ou.univariate(tree,X)
  } else {cache <- tree}
  if(is.null(pars$ntheta)){
    pars$ntheta <- length(unique(names(unlist(cache$maps))))
  }
  TotExp <- exp(-cache$height*pars$alpha)
  nbranch <- length(cache$edge.length)
  maps <- cache$maps
  shifts <- unlist(lapply(maps,length),F,F)-1
  irow <- rep(1:nbranch,shifts+1)
  csbase <- cache$nH[irow]
  csmaps <- csbase+unlist(lapply(maps,cumsum),FALSE,TRUE)
  multips <- which(irow[2:length(irow)]==irow[1:(length(irow)-1)])
  csbase[multips+1] <- csmaps[multips]
  tmp <- (exp(csmaps*pars$alpha)-exp(csbase*pars$alpha))*TotExp
  bW <- matrix(0,nrow=nbranch,ncol=pars$ntheta)
  index <- irow + (as.integer(names(tmp))-1)*nbranch
  if(any(shifts>1)){
    tmp <- tapply(tmp,index,sum)
    bW[as.numeric(names(tmp))] <- tmp
  } else {bW[index] <- tmp}
  W=cache$branchtrace%*%bW
  W[,1] <- W[,1]+TotExp
  return(W)
}
#' @rdname simmap.W
.simmap.W <- function(cache,pars,TotExp){
  nbranch <- length(cache$edge.length)
  maps <- cache$maps
  shifts <- unlist(lapply(maps,length),F,F)-1
  irow <- rep(1:nbranch,shifts+1)
  csbase <- cache$nH[irow]
  csmaps <- csbase+unlist(lapply(maps,cumsum),FALSE,TRUE)
  multips <- which(irow[2:length(irow)]==irow[1:(length(irow)-1)])
  csbase[multips+1] <- csmaps[multips]
  tmp <- (exp(csmaps*pars$alpha)-exp(csbase*pars$alpha))*TotExp
  bW <- matrix(0,nrow=nbranch,ncol=pars$ntheta)
  index <- irow + (as.integer(names(tmp))-1)*nbranch
  if(any(shifts>1)){
    tmp <- tapply(tmp,index,sum)
    bW[as.numeric(names(tmp))] <- tmp
  } else {bW[index] <- tmp}
  W=cache$branchtrace%*%bW
  W[,1] <- W[,1]+TotExp
  return(W)
}

#' Calculate the weight matrix of a set of regimes on a phylogeny
#' 
#' @rdname edgemap.W
edgemap.W <- function(tree, pars,emap,alpha=NULL){
  if(class(tree)=="phylo"){
    ntips <- length(tree$tip.label)
    X <- rep(1,ntips)
    names(X) <- tree$tip.label
    cache <- .prepare.ou.univariate(tree,X)
  } else {cache <- tree}
  if(is.null(pars)){
    pars$alpha <- alpha
    pars$ntheta <- sum(emap$sh)+1
  }
  TotExp <- exp(-cache$height*pars$alpha)
  if(pars$ntheta>1){
    tmp <- exp(cbind(cache$nH,cache$nH+emap$r1,cache$nH+emap$r1+emap$r2)*pars$alpha)
    tmp2 <- (tmp[,2:3]-tmp[,1:2])*TotExp
    bW <- matrix(0,nrow=dim(emap)[1],ncol=pars$ntheta)
    index <- 1:dim(emap)[1]+(emap$t2-1)*dim(emap)[1]
    bW[index] <- tmp2[,2]
    index <- 1:dim(emap)[1]+(emap$t1-1)*dim(emap)[1]
    bW[index] <- tmp2[,1]                    
    W=cache$branchtrace%*%bW#t(sapply(cache$branchtrace,function(x) .colSums(bW[x,],sum(x),pars$ntheta)))
    W[,1] <- W[,1]+TotExp
    W[W==Inf] <- NA
    if(mean(cache$nH/(log(2)/pars$alpha))>600 & sum(is.na(W))>0){
      sub <- subset(emap,!is.na(emap$tip))
      t2.tmp <- sub$t2
      names(t2.tmp) <- sub$tip
      t2.tmp <- t2.tmp[cache$tip.label]
      for(k in 1:length(t2.tmp)){
        if(is.na(W[k,t2.tmp[k]])){
          W[k,t2.tmp[k]] <- 1
        }
      }
      W[is.na(W)] <- 0
    }
  } else {W <-  rep(1,cache$ntips)}
  W
}
#' @rdname edgemap.W
.edgemap.W <- function(cache, pars,emap,TotExp,alpha=NULL){
  if(pars$ntheta>1){
    tmp <- exp(cbind(cache$nH,cache$nH+emap$r1,cache$nH+emap$r1+emap$r2)*pars$alpha)
    tmp2 <- (tmp[,2:3]-tmp[,1:2])*TotExp
    bW <- matrix(0,nrow=dim(emap)[1],ncol=pars$ntheta)
    index <- 1:dim(emap)[1]+(emap$t2-1)*dim(emap)[1]
    bW[index] <- tmp2[,2]
    index <- 1:dim(emap)[1]+(emap$t1-1)*dim(emap)[1]
    bW[index] <- tmp2[,1]                    
    W=cache$branchtrace%*%bW#t(sapply(cache$branchtrace,function(x) .colSums(bW[x,],sum(x),pars$ntheta)))
    W[,1] <- W[,1]+TotExp
    W[W==Inf] <- NA
    if(mean(cache$nH/(log(2)/pars$alpha))>600 & sum(is.na(W))>0){
      sub <- subset(emap,!is.na(emap$tip))
      t2.tmp <- sub$t2
      names(t2.tmp) <- sub$tip
      t2.tmp <- t2.tmp[cache$tip.label]
      for(k in 1:length(t2.tmp)){
        if(is.na(W[k,t2.tmp[k]])){
          W[k,t2.tmp[k]] <- 1
        }
      }
      W[is.na(W)] <- 0
    }
  } else {W <-  rep(1,cache$ntips)}
  W
}
