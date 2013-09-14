emap2simmap <- function(emap,tree){
  foo <- function(x){
    tmp <- unlist(emap[x,c('r1','r2')])
    names(tmp) <- c(emap$t1[x],emap$t2[x])
    tmp[tmp>0]
  }
  nb <- sum(emap$sh)
  if(nb>0){
    col <- c("#000000",rainbow(nb))
  } else {col <- 1}
  names(col) <- 1:(nb+1)
  tree$maps <- lapply(1:dim(emap)[1],foo)
  tree$col <- col
  tree
}

QG.alpha <- function(pars){
  pars$h2*pars$P/(pars$P+pars$w2*pars$P)
}
QG.sig2 <- function(pars){
  (pars$h2*pars$P)/pars$Ne
}

OU.repar <- function(pars){
  alpha <- log(2)/pars$halflife
  sig2 <- (2*log(2)/(pars$halflife))*pars$Vy
  return(list(alpha=alpha,sig2=sig2))
}
