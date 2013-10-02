#' Make a prior function for bayOU
#' 
#' This function generates a prior function to be used for bayOU according to user specifications.
#'
#' @param tree A tree object of class "phylo"
#' @param dists A list providing the function names of the distribution functions describing the prior distributions of parameters (see details). If no
#' distributions are provided for a parameter, default values are given. Note that the names are provided as text strings, not the functions themselves.
#' @param param A list providing the parameter values of the prior distributions (see details). 
#' @param plot.prior A logical indicating whether the prior distributions should be plotted.
#' @param model One of three specifications of the OU parameterization used. 
#' Takes values \code{"OU"} (alpha & sig2), \code{"QG"} (h2, P, w2, Ne), or \code{"OUrepar"} (halflife,Vy)
#' @param type Specifies the type of input for regime mapping used. Takes values 
#' \code{"pars"} (regime shits specified in the parameter list by \code{pars$sb} and \code{pars$loc}),
#' \code{"emap"} (regime shifts specified by a provided edge map), or \code{"simmap"} (regime shifts obtained from \code{tree$maps})
#' 
#' @details Default distributions and parameter values are given as follows:
#' OU: \code{list(dists=list("dalpha"="dlnorm","dsig2"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="equal","dloc"="dunif"),
#'    param=list("dalpha"=list(),"dsig2"=list(),"dtheta"=list(),"dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list()))}
#' QG: \code{list(dists=list("dh2"="dbeta","dP"="dlnorm","dw2"="dlnorm","dNe"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="equal","dloc"="dunif"),
#'    param=list("dh2"=list(shape1=1,shape2=1),"dP"=list(),"dw2"=list(),"dNe"=list(),"dtheta"=list(),"dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list()))}
#' OUrepar: \code{list(dists=list("dhalflife"="dlnorm","dVy"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="equal","dloc"="dunif"),
#'    param=list("dhalflife"=list("meanlog"=0.25,"sdlog"=1.5),"dVy"=list("meanlog"=1,"sdlog"=2),"dk"=list(lambda=1,kmax=2*ntips-2),"dtheta"=list(),"dloc"=list(min=0,max=1)),"dsb"=list())}
#' 
#' \code{dalpha, dsig2, dh2, dP, dw2, dNe, dhalflife},  and \code{dVy} must be positive continuous distributions and provide the parameters used to calculate alpha and sigma^2 of the OU model. 
#' \code{dtheta} must be continuous and describes the prior distribution of the optima. dk is the prior distribution for the number of shifts. For Poisson and conditional Poisson (cdpois) are provided
#' the parameter \code{lambda}, which provides the total number of shifts expected on the tree (not the rate per unit branch length). Otherwise, \code{dk} can take any positive, discrete distribution.  
#' dsb indicates the prior probability of a given set of branches having shifts. Three options are possible: \code{"equal"},\code{"equalnotips"} or \code{"free"}. \code{"equal"} assigns equal probability to each branch
#' sampled without replacement. \code{"equalnotips"} but restricts shifts to only internal branches. \code{"free"} allows multiple shifts per branch, and determines the prior probability of a set of shift locations
#' using a multinomial distribution with probabilities proportional to branch lengths taken from the tree. \code{"dloc"} indicates the prior probability of a shift on a given branch. Currently, all locations are
#' given equal density of 1. All distributions are set to return log-transformed probability densities. 
#' 
#' @return returns a prior function that calculates the log prior density for a set of parameter values provided in a list with correctly named values.
#' 
#' @examples 
#' ## Load data
#' data(chelonia.simmap)
#' tree <- chelonia.simmap$tree
#' dat <- chelonia.simmap$dat
#' 
#' ##Make prior function
#' prior <- make.prior(tree,dists=list(dalpha="dunif",dsig2="dunif"),param=list(dalpha=list(min=0,max=10),dsig2=list(min=0,max=10)),type="pars")
#' 
#' ##Define a set of parameter values
#' emap <- chelonia.simmap$emap
#' pars <- list(alpha=0.1, sig2=1, k=16, theta=c(3,4,5,6), sb=which(emap$sh==1), loc=emap$r1[which(emap$sh==1)])
#' 
#' ##Calculate prior probability
#' prior(pars)
#' 
#' ##Return a list of functions used to calculate prior
#' attributes(prior)$functions
#' 
#' ##Return parameter values used in prior distribution
#' attributes(prior)$parameters
#' 
#' ##Alternative parameterization using phylogenetic half-life and stationary variance
#' pars <- list(halflife=20, Vy=1, k=16, theta=c(3,4,5,6), sb=which(emap$sh==1), loc=emap$r1[which(emap$sh==1)])
#' prior2 <- make.prior(tree,dists=list(dhalflife="dlnorm",dVy="dlnorm"),param=list(dhalflife=list(meanlog=3,sdlog=1),dVy=list(meanlog=0,sdlog=2)),model="OUrepar",type="pars")
#' prior2(pars)
#' 
#' ##Specify different types of inputs
#' prior.pars <- make.prior(tree, model="OUrepar", type="pars")
#' prior.emap <- make.prior(tree, model="OUrepar", type="emap")
#' prior.simmap <- make.prior(tree, model="OUrepar", type="simmap")
#' 
#' prior.pars(pars)
#' prior.emap(pars,emap=emap)
#' prior.simmap(pars,tree)

make.prior <- function(tree,dists=list(),param=list(),plot.prior=TRUE,model="OU",type="pars"){
  tree <- reorder.phylo(tree, "postorder")
  nH <- max(nodeHeights(tree))
  ntips <- length(tree$tip.label)
  TH <- sum(tree$edge.length)
  default.OU <- list(dists=list("dalpha"="dlnorm","dsig2"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="dsb","dloc"="dunif"),param=list("dalpha"=list(),"dsig2"=list(),"dtheta"=list(),"dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list(ntips=ntips, bmax=1, prob=1)))
  default.QG <- list(dists=list("dh2"="dbeta","dP"="dlnorm","dw2"="dlnorm","dNe"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="dsb","dloc"="dunif"),param=list("dh2"=list(shape1=1,shape2=1),"dP"=list(),"dw2"=list(),"dNe"=list(),"dtheta"=list(),"dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list(ntips=ntips, bmax=1, prob=1)))
  default.OUrepar <- list(dists=list("dhalflife"="dlnorm","dVy"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="dsb","dloc"="dunif"),param=list("dhalflife"=list("meanlog"=0.25,"sdlog"=1.5),"dVy"=list("meanlog"=1,"sdlog"=2),"dk"=list(lambda=1,kmax=2*ntips-2),"dtheta"=list(),"dloc"=list(min=0,max=1),"dsb"=list(ntips=ntips, bmax=1, prob=1)))
  #default.OUcpp <- list(dists=list("dalpha"="dlnorm","dsig2"="dlnorm","dsig2jump"="dlnorm","dk"="dpois","dtheta"="dnorm","dloc"="dunif"),param=list("dalpha"=NULL,"dsig2"=list(),"dsig2jump"=list(),"dtheta"=list(),"dk"=list(lambda=1),"dloc"=list(min=0,max=TH)))
  #default.QGcpp <- list(dists=list("dh2"="dbeta","dP"="dlnorm","dw2"="dlnorm","dNe"="dlnorm","dk"="dpois","dtheta"="dnorm","dloc"="dunif"),param=list("dh2"=list(shape1=1,shape2=1),"dP"=list(),"dw2"=list(),"dNe"=list(),"dsig2jump"=list(),"dtheta"=list(),"dk"=list(lambda=1),"dloc"=list(min=0,max=TH)))
  #default.OUreparcpp <- list(dists=list("dhalflife"="dlnorm","dVy"="dlnorm","dsig2jump"="dlnorm","dk"="dpois","dtheta"="dnorm","dloc"="dunif"),param=list("dhalflife"=list("meanlog"=0.25,"sdlog"=1.5),"dVy"=list("meanlog"=1,"sdlog"=2),"dk"=list(lambda=1),"dsig2jump"=list(),"dtheta"=list(),"dloc"=list(min=0,max=TH)))
  default <- switch(model,"OU"=default.OU,"QG"=default.QG,"OUrepar"=default.OUrepar)#,"OUcpp"=default.OUcpp,"QGcpp"=default.QGcpp,"OUreparcpp"=default.OUreparcpp)
  notprovided <- setdiff(names(default$dist),names(dists))
  pars.notprovided <- setdiff(names(default$param),names(param))
  dists[notprovided] <- default$dists[notprovided]
  param[pars.notprovided] <- default$param[pars.notprovided]
  if(length(setdiff(names(default$param),names(param)))>0)
    stop("Provided parameters are not in the model")
  if(length(setdiff(names(default$dists),names(dists)))>0)
    stop("Provided parameters are not in the model")
  if(dists$dsb=="dsb") param$dsb$ntips <- ntips
  if(dists$dloc=="dunif" & grep("dsb",dists$dsb)==1){
        dists$dloc <- "dloc"
    }
  prior.fx <- lapply(dists,get)
  param <- suppressWarnings(lapply(param,function(x){ x$log = TRUE; x}))
  prior.param <- param[match(names(prior.fx),names(param))]
  prior.fx <- lapply(1:length(prior.param),function(x) .set.defaults(prior.fx[[x]],defaults=prior.param[[x]]))
  names(prior.fx) <- names(prior.param)
  #if(model %in% c("OUcpp","QGcpp","OUreparcpp")){
    #droot <- prior.fx$dtheta
    #prior.fx$dtheta <- function(x){
      #droot(x[1])
    #}
  #}
  par.names <- gsub('^[a-zA-Z]',"",names(dists))
  
  if(plot.prior){
    par(mfrow=rep(ceiling(sqrt(length(dists))),2))
    rfx <- lapply(gsub('^[a-zA-Z]',"r",dists),function(x) try(get(x),silent=TRUE))
    rprior.param <- prior.param[1:(length(prior.param))]
    rprior.param <- lapply(rprior.param, function(x) x[-length(x)])
    if(dists$dsb=="dsb" & any(rprior.param$dsb$bmax==1)){rprior.param$dsb$bmax[rprior.param$dsb$bmax==1] <- Inf; rprior.param$dsb$prob <- 1}
    rfx <- lapply(1:length(rprior.param),function(x) try(.set.defaults(rfx[[x]],defaults=rprior.param[[x]]),silent=TRUE))
    plot.names<-par.names[sapply(rfx,class)=="function"]
    rfx <- rfx[sapply(rfx,class)=="function"]
    names(rfx) <- names(rprior.param)
    nsim <-500000
    for(i in 1:length(rfx)){
      if(names(rfx)[i]=="dsb"){
        curve(sapply(x,function(y) prior.fx$dsb(y,log=FALSE)),xlim=c(1,(2*ntips-2)),ylab="Density",main="branches")
      } else {
        x <- rfx[[i]](nsim)
        qq <- quantile(x,c(0.001,0.999))
        plot(density(x),xlim=qq, main=plot.names[i],lwd=2)
      }
    }
  }    
  
  if(type=="pars"){
    priorFUN <- function(pars,cache){
      if(any(!(par.names %in% names(pars)))) stop(paste("Missing parameters: ", paste(par.names[!(par.names %in% names(pars))],collapse=" ")))
      pars.o <- pars[match(par.names,names(pars))]
      pars.o <- pars.o[!is.na(names(pars.o))]
      densities <- sapply(1:length(pars.o),function(x) prior.fx[[x]](pars.o[[x]]))
      names(densities) <- par.names
      lnprior <- sum(unlist(densities,F,F))
      return(lnprior)
    }
  }
  if(type=="emap"){
    priorFUN <- function(pars,cache,emap){
      pars$sb <- which(emap$sh==1)
      pars$loc <- emap$r1[pars$sb]
      if(any(!(par.names %in% names(pars)))) stop(paste("Missing parameters: ", paste(par.names[!(par.names %in% names(pars))],collapse=" ")))
      pars.o <- pars[match(par.names,names(pars))]
      pars.o <- pars.o[!is.na(names(pars.o))]
      densities <- sapply(1:length(pars.o),function(x) prior.fx[[x]](pars.o[[x]]))
      names(densities) <- par.names
      lnprior <- sum(unlist(densities,F,F))
      return(lnprior)
    }
  }
  if(type=="simmap"){
    priorFUN <- function(pars,cache){
      pars$sb <- rep(1:length(cache$edge.length),sapply(cache$maps,length)-1)
      pars$loc <- sapply(pars$sb,function(x) cache$maps[[x]][-1])
      pars$t2 <- names(pars$loc)
      pars$loc <- unname(pars$loc)
      if(any(!(par.names %in% names(pars)))) stop(paste("Missing parameters: ", paste(par.names[!(par.names %in% names(pars))],collapse=" ")))
      pars.o <- pars[match(par.names,names(pars))]
      pars.o <- pars.o[!is.na(names(pars.o))]
      densities <- sapply(1:length(pars.o),function(x) prior.fx[[x]](pars.o[[x]]))
      names(densities) <- par.names
      lnprior <- sum(unlist(densities,F,F))
      return(lnprior)
    }
  }
  attributes(priorFUN) <- list("model"=model,"parnames"=par.names,"distributions"=dists,"parameters"=prior.param,"functions"=prior.fx)
  class(priorFUN) <- c("priorFn","function")
  return(priorFUN)
}
