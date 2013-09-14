#' Make a prior function for bayOU
#' 
#' This function generates a prior function to be used for bayOU according to user specifications.

make.prior <- function(tree,dists=list(),param=list(),plot.prior=FALSE,model="OU"){
  nH <- max(nodeHeights(tree))
  ntips <- length(tree$tip.label)
  TH <- sum(tree$edge.length)
  default.OU <- list(dists=list("dalpha"="dlnorm","dsig2"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="equal","dloc"="dunif"),param=list("dalpha"=list(),"dsig2"=list(),"dtheta"=list(),"dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list()))
  default.QG <- list(dists=list("dh2"="dbeta","dP"="dlnorm","dw2"="dlnorm","dNe"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="equal","dloc"="dunif"),param=list("dh2"=list(shape1=1,shape2=1),"dP"=list(),"dw2"=list(),"dNe"=list(),"dtheta"=list(),"dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list()))
  default.OUrepar <- list(dists=list("dhalflife"="dlnorm","dVy"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="equal","dloc"="dunif"),param=list("dhalflife"=list("meanlog"=0.25,"sdlog"=1.5),"dVy"=list("meanlog"=1,"sdlog"=2),"dk"=list(lambda=1,kmax=2*ntips-2),"dtheta"=list(),"dloc"=list(min=0,max=1)))
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
  if(model %in% c("OU","QG","OUrepar")){
  if(dists$dsb %in% c("equal","equalnotips","free")){
      if(dists$dsb=="equal"){
        dsb <- function(sb,ntips=ntips,log=TRUE){
          if(log){
            return(log(1/choose(2*ntips-2,length(sb))))
          } else {return(1/choose(2*ntips-2,length(sb)))}
        }
        maxL <- 1
      }
      if(dists$dsb=="equalnotips"){
        dsb <- function(sb,ntips=ntips,log=TRUE){
          if(log){
          return(log(1/choose(ntips-2,length(sb))))
          } else {return(1/choose(ntips-2,length(sb)))}
        }
        maxL <- 1
      }
      if(dists$dsb=="free"){
        dsb <- function(sb,ntips=ntips,log=TRUE) if(log) return(0) else return(1)
        maxL <- TH
      }
      dists$dsb <- "dsb"
      param$dsb <- list("ntips"=ntips)
    }
    if(dists$dloc=="dunif" & dists$dsb=="dsb"){
      dloc <- function(x,min=0,max=maxL,log=TRUE) if(log) return(0) else return(1)
      param$dloc <- list(min=0,max=maxL)
      dists$dloc <- "dloc"
    }
  }
  prior.fx <- lapply(dists,get)
  param <- suppressWarnings(lapply(param,function(x){ x$log = TRUE; x}))
  prior.param <- param[match(names(prior.fx),names(param))]
  if(dists$dk %in% c("dpois","cdpois")){
    prior.param$dk$lambda <- prior.param$dk$lambda*TH
  }
  prior.fx <- lapply(1:length(prior.param),function(x) set.defaults(prior.fx[[x]],defaults=prior.param[[x]]))
  names(prior.fx) <- names(prior.param)
  #if(model %in% c("OUcpp","QGcpp","OUreparcpp")){
    #droot <- prior.fx$dtheta
    #prior.fx$dtheta <- function(x){
      #droot(x[1])
    #}
  #}
  
  if(plot.prior){
  }
  par.names <- gsub('^[a-zA-Z]',"",names(dists))
  priorFUN <- function(pars,cache){
    if(any(!(par.names %in% names(pars)))) stop(paste("Missing parameters: ", paste(par.names[!(par.names %in% names(pars))],collapse=" ")))
    pars.o <- pars[match(par.names,names(pars))]
    pars.o <- pars.o[!is.na(names(pars.o))]
    densities <- sapply(1:length(pars.o),function(x) prior.fx[[x]](pars.o[[x]]))
    names(densities) <- par.names
    lnprior <- sum(unlist(densities,F,F))
    return(lnprior)
  }
  attributes(priorFUN) <- list("model"=model,"parnames"=par.names,"distributions"=dists,"parameters"=prior.param,"functions"=prior.fx)
  return(priorFUN)
}
