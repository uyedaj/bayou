#' Simulates parameters from bayou models
#' 
#' \code{priorSim} Simulates parameters from the prior distribution specified by \code{make.prior}
#' 
#' @param prior A prior function created by \code{bayou::make.prior}
#' @param tree A tree of class 'phylo'
#' @param plot A logical indicating whether the simulated parameters should be plotted
#' @param nsim The number of parameter sets to be simulated
#' @param shiftpars A vector of parameters that split upon a shift, default is "theta"
#' @param ... Parameters passed on to \code{plotSimmap(...)}
#' 
#' @return A list of bayou parameter lists
#' 
#' @export
priorSim <- function(prior, tree, plot=TRUE, nsim=1, shiftpars="theta", ...){
  tree <- reorder(tree,'postorder')
  model <- attributes(prior)$model
  dists <- attributes(prior)$dist
  fixed <- which(attributes(prior)$dist=="fixed")
  allnames <- gsub('^[a-zA-Z]',"", names(attributes(prior)$dist))
  notfixed <- which(attributes(prior)$dist!="fixed")
  dists <- dists[notfixed]
  prior.params <- attributes(prior)$param
  rdists <- lapply(dists,function(x) gsub('^[a-zA-Z]',"r",x))
  prior.params <- lapply(prior.params,function(x) x[-which(names(x)=="log")])
  rdists.fx <- lapply(rdists,get)
  rdists.fx <- lapply(1:length(rdists.fx),function(x) .set.defaults(rdists.fx[[x]],defaults=prior.params[[x]]))
  names(rdists.fx) <- gsub('^[a-zA-Z]',"r",names(rdists))
  N <- setNames(rep(1, length(allnames)), allnames)
  N[fixed] <- "fixed"
  N[shiftpars] <- "ntheta"
  N["loc"] <- "loc"
  
  simpar <- lapply(1:nsim, function(x){ll <- lapply(1:length(allnames), function(y) numeric(0)); names(ll) <- allnames; ll})
  simpar <- lapply(simpar, function(x){names(x) <- allnames; x})
  
  if(N["k"]=="fixed"){
    simpar <- lapply(1:nsim, function(x) {simpar[[x]][["k"]] <- attributes(prior)$fixed$k; simpar[[x]]})
  } else {
    simpar <- lapply(1:nsim, function(x) {simpar[[x]][["k"]] <- rdists.fx[["rk"]](1); simpar[[x]]})
  }
  if(N["sb"]=="fixed"){
    simpar <- lapply(1:nsim, function(x) {simpar[[x]][["sb"]] <- attributes(prior)$fixed$sb; simpar[[x]]})
  } else {
    simpar <- lapply(1:nsim, function(x) {simpar[[x]][["sb"]] <- rdists.fx[["rsb"]](simpar[[x]]$k); simpar[[x]]})
  }
  if(!"ntheta" %in% names(attributes(prior)$fixed)){
    simpar <- lapply(1:nsim, function(x) {simpar[[x]][["ntheta"]] <- simpar[[x]]$k + 1; simpar[[x]]})
  } else {
    simpar <- lapply(1:nsim, function(x) {simpar[[x]][["ntheta"]] <- attributes(prior)$fixed$ntheta; simpar[[x]]})
  }
  if(!"t2" %in% names(attributes(prior)$fixed)){
      simpar <- lapply(1:nsim, function(x) {
        if(simpar[[x]]$k == 0){simpar[[x]][["t2"]] <- numeric(0)} else {simpar[[x]][["t2"]] <-  2:(simpar[[x]]$k+1)};
        simpar[[x]]})
  } else {
    simpar <- lapply(1:nsim, function(x) {simpar[[x]][["t2"]] <- attributes(prior)$fixed$t2; simpar[[x]]})
  }

  for(i in 1:length(allnames)){
    if(!allnames[i] %in% c("k", "sb")){
      if(N[i]=="fixed"){
        simpar <- lapply(1:nsim, function(x) {simpar[[x]][[allnames[i]]] <- attributes(prior)$fixed[[allnames[i]]]; simpar[[x]]})
      } else {
        rdistname <- paste("r", allnames[i], sep="")
        if(N[i]=="loc"){
          n <- sapply(simpar, function(x) x$k)
          simpar <- lapply(1:nsim, function(x) {simpar[[x]][[allnames[i]]] <- rdists.fx[[rdistname]](n[x])*tree$edge.length[simpar[[x]]$sb]; simpar[[x]]})
        } else {
          if(N[i]=="1"){
            n <- rep(1, nsim)
          }
          if(N[i]=="ntheta"){
            n <- sapply(simpar, function(x) x$ntheta)
          }
          simpar <- lapply(1:nsim, function(x) {simpar[[x]][[allnames[i]]] <- rdists.fx[[rdistname]](n[x]); simpar[[x]]})
        }
      }
    }
  }
  if(plot){
    for(i in 1:nsim){
      if(simpar[[i]]$k > 0){
      maps <- pars2simmap(simpar[[i]],tree)
      plotRegimes(maps$tree, ...)
      } else {
        plot(tree, ...)
      }
    }
  }
  return(list(pars=simpar,tree=tree))
}

#' Simulates data from bayou models
#' 
#' \code{dataSim} Simulates data for a given bayou model and parameter set
#' 
#' @param pars A bayou formated parameter list
#' @param model The type of model specified by the parameter list (either "OU", "OUrepar" or "QG").
#' @param tree A tree of class 'phylo'
#' @param map.type Either "pars" if the regimes are taken from the parameter list, or "simmap" if taken from the stored simmap in the tree
#' @param SE A single value or vector equal to the number of tips specifying the measurement error that should be simulated at the tips
#' @param phenogram A logical indicating whether or not the simulated data should be plotted as a phenogram
#' @param ... Optional parameters passed to \code{phenogram(...)}.
#'
#' @description This function simulates data for a given set of parameter values.
#' 
#' @export
dataSim <- function(pars, model, tree, map.type="pars", SE=0, 
                    phenogram=TRUE, ...){
  if(model %in% c("QG")){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model %in% c("OUrepar")){
    p <- OU.repar(pars)
    pars$alpha <- p$alpha
    pars$sig2 <- p$sig2
  }
  if(map.type=="simmap"){
    print("Using mapped regimes from ape tree file")
    maps <- tree$maps
  }
  if(map.type=="pars"){
    print("Using mapped regimes from parameter list")
    maps <- pars2simmap(pars, tree)$tree$maps
  }
  dummy <- rep(0, length(tree$tip.label))
  names(dummy) <- tree$tip.label
  tree$maps <- maps
  cache <- .prepare.ou.univariate(tree,dummy)
  cache$maps <- maps
  W <- .simmap.W(cache,pars)
  if(pars$k > 0){
    E.th <- W%*%pars$theta
  } else {E.th <- W*pars$theta}
  Sigma <- .ouMatrix(vcv.phylo(tree),pars$alpha)*pars$sig2
  diag(Sigma) <- diag(Sigma)+SE
  X <- mvrnorm(1,E.th,Sigma)
  if(phenogram){
    col <- c(1,rainbow(pars$k))
    names(col) <- 1:pars$ntheta
    phenogram(cache$phy,X,colors=col, spread.labels=FALSE, ...)
  }
  return(list(W=W, E.th=E.th,dat=X))
}

