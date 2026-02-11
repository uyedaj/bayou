#' Make a reference function in bayou
#'
#' This function generates a reference function from a mcmc chain for use in marginal likelihood
#' estimation.
#'
#' @param chain An mcmc chain produced by \code{bayou.mcmc()} and loaded with \code{load.bayou()}
#' @param model A string specifying the model ("OU", "QG", "OUrepar") or a model parameter list
#' @param priorFn The prior function used to generate the mcmc chain
#' @param burnin The proportion of the mcmc chain to be discarded when generating the reference function
#' @param plot Logical indicating whether or not a plot should be created
#'
#' @details Distributions are fit to each mcmc chain and the best-fitting distribution is chosen as
#' the reference distribution for that parameter using the method of Fan et al. (2011). For positive
#' continuous parameters \code{alpha, sigma^2, halflife, Vy, w2, Ne}, Log-normal, exponential, gamma and weibull
#' distributions are fit. For continuous distributions \code{theta}, Normal, Cauchy and Logistic distributions
#' are fit. For discrete distributions, \code{k}, negative binomial, poisson and geometric distributions are fit.
#' Best-fitting distributions are determined by AIC.
#'
#' @export
#' @return Returns a reference function of class "refFn" that takes a parameter list and returns the log density
#' given the reference distribution. If \code{plot=TRUE}, a plot is produced showing the density of variable parameters
#' and the fitted distribution from the reference function (in red).
make.refFn <- function(chain, model, priorFn, burnin=0.3, plot=TRUE){
  if(is.character(model)){
    model.pars <- switch(model, "OU"=model.OU, "QG"=model.QG, "OUrepar"=model.OUrepar)#, "bd"=model.bd)
  } else {
    model.pars <- model
    model <- "Custom"
  }
  contdists <- c("norm", "cauchy", "logis")
  poscontdists <- c("lnorm", "exp", "gamma", "weibull")
  discdists <- c("nbinom", "pois", "geom")
  bounddists <- c("beta")
  parorder <- model.pars$parorder
  postburn <- max(c(1,round(burnin*length(chain[[1]]),0))):length(chain[[1]])
  dists <- attributes(priorFn)$functions
  dists <- dists[!(names(dists) %in% paste("d", model.pars$shiftpars, sep=""))]
  distnames <- gsub('^[d]', "", names(dists))

  ## This code tests each prior for whether it's continuous, positive continuous, bounded or discrete...not perfect; only tests bounded
  ## distributions between 0 and 1, etc....
  test.dists <- suppressWarnings(apply(is.finite(sapply(dists, function(x) c(x(-0.5), x(1.5), x(0.5), x(1)))), 2, function(x) paste(as.numeric(x), collapse="")))
  dist.types <- rep(NA, length(dists))
  names(dist.types) <- names(dists)
  dist.types[which(test.dists=="0111")] <- "pcdist"
  dist.types[which(test.dists=="1111")]  <- "cdist"
  dist.types[which(test.dists=="0010" |test.dists== "0011")] <- "bdist"
  dist.types[which(test.dists=="0001")] <- "ddist"
  refFx <- list()
  refNames <- list()
  dists <- list()
  parameters <- list()
  x <- NULL
  for(i in 1:length(dist.types)){
    parname <- gsub('^[a-zA-Z]',"",names(dist.types)[i])
    xx <- unlist(chain[[parname]][postburn])
    fitdists <- switch(dist.types[i], "ddist"=discdists, "pcdist"=poscontdists, "bdist"=bounddists, "cdist"=contdists)
    {
      tmpFits <- lapply(fitdists, function(x) suppressWarnings(try(fitdistrplus::fitdist(xx, x), silent=TRUE)))
      tmpFits[sapply(tmpFits, function(x) (class(x)=="try-error"))] <- lapply(which(sapply(tmpFits, function(x) (class(x)=="try-error"))), function(j) suppressWarnings(try(fitdistrplus::fitdist(xx, fitdists[j], method="mme"), silent=TRUE)))
      tmpFits <- tmpFits[sapply(tmpFits, function(x) !(class(x)=="try-error"))]
      aic <- sapply(tmpFits, function(x) x$aic)
      fit <- tmpFits[[which(aic==min(aic,na.rm=TRUE))]]
      ## Fix for problem with negative binomial distribution
      if(fit$distname == "nbinom"){fit$estimate <- c(fit$estimate, prob=unname(fit$estimate['size']/(fit$estimate['size']+fit$estimate['mu'])))}
      fitPars <- as.list(fit$estimate)
      fitPars$log <- TRUE
      fitName <- fit$distname
      fitfx <- get(paste("d",fitName, sep=""))
      refFx[[i]] <- .set.defaults(fitfx, defaults=fitPars)
    }

  dists[[i]] <- fitName
  parameters[[i]] <- fitPars
  }
names(refFx) <- gsub('^[a-zA-Z]',"", names(dist.types))
par.names <- names(refFx)
names(dists) <- par.names
names(parameters) <- par.names
if(attributes(priorFn)$distributions$dsb!="fixed"){
  dists$dsb <- "dsb"
  parameters$dsb <- attributes(priorFn)$parameters$dsb
  refFx$dsb <- attributes(priorFn)$functions$dsb
  par.names <- c(par.names, "sb")
}
if(attributes(priorFn)$distributions$dloc!="fixed"){
  dists$dloc = "dloc"
  parameters$dloc <- attributes(priorFn)$parameters$dloc
  refFx$dloc <- attributes(priorFn)$functions$dloc
  par.names <- c(par.names, "loc")
}

if(plot){
  pars2plot <- par.names[!(par.names %in% c("sb", "loc"))]
  par(mfrow=c(ceiling(length(pars2plot)/2),2))
  for(i in 1:length(pars2plot)){
    plot(density(unlist(chain[[pars2plot[i]]][postburn])), main=pars2plot[i])
    if(pars2plot[i]=="k"){
      points(seq(ceiling(par('usr')[1]),floor(par('usr')[2]),1), refFx[[pars2plot[i]]](seq(ceiling(par('usr')[1]),floor(par('usr')[2]),1),log=FALSE),pch=21,bg="red")
    } else {x <- NULL; curve(refFx[[pars2plot[i]]](x,log=FALSE), add=TRUE, col="red")}
  }
}
refFUN <- function(pars,cache){
  if(any(!(par.names %in% names(pars)))) stop(paste("Missing parameters: ", paste(par.names[!(par.names %in% names(pars))],collapse=" ")))
  pars.o <- pars[match(par.names,names(pars))]
  pars.o <- pars.o[!is.na(names(pars.o))]
  densities <- sapply(1:length(pars.o),function(x) refFx[[x]](pars.o[[x]]))
  names(densities) <- par.names
  lnprior <- sum(unlist(densities,F,F))
  return(lnprior)
}
attributes(refFUN) <- list("model"=model,"parnames"=par.names,"distributions"=dists,"parameters"=parameters,"functions"=refFx)
class(refFUN) <- c("refFn","function")
return(refFUN)
}


#' Makes a power posterior function in bayou
#'
#' This function generates a power posterior function for estimation of marginal likelihood using the stepping stone method
#'
#' @param Bk The sequence of steps to be taken from the reference function to the posterior
#' @param priorFn The prior function to be used in marginal likelihood estimation
#' @param refFn The reference function generated using \code{make.refFn()} from a preexisting mcmc chain
#' @param model A string specifying the model type ("OU", "OUrepar", "QG") or a model parameter list
#'
#' @details For use in stepping stone estimation of the marginal likelihood using the method of Fan et al. (2011).
#' @export
#' @return A function of class "powerposteriorFn" that returns a list of four values: \code{result} (the log density of the power posterior),
#' \code{lik} (the log likelihood), \code{prior} (the log prior), \code{ref} the log reference density.
make.powerposteriorFn <- function(Bk, priorFn, refFn, model){
  #Turn these off for now, need to add back in checks
  #model <- attributes(priorFn)$model
  #if(model != attributes(refFn)$model) stop("Error: prior and reference function are not of same type")
  if(is.character(model)){
    model.pars <- switch(model, "OU"=model.OU, "QG"=model.QG, "OUrepar"=model.OUrepar)#, "bd"=model.bd)
  } else {
    model.pars <- model
    model <- "Custom"
  }
  powerposteriorFn <- function(k, Bk, pars, cache, dat, model=model.pars){
    lik <- model$lik.fn(pars, cache, dat)$loglik
    prior <- priorFn(pars, cache)
    ref <- refFn(pars, cache)
    coeff <- c(Bk[k],Bk[k],(1-Bk[k]))
    result <- c(lik, prior, ref)
    result[coeff==0] <- 0
    result <- result*coeff
    result <- sum(result)
    return(list(result=result, lik=lik, prior=prior, ref=ref))
  }
  class(powerposteriorFn) <- c("powerposteriorFn", "function")
  return(powerposteriorFn)
}

powerPosteriorFn <- function(k, Bk, lik, prior, ref){
  coeff <- c(Bk[k],Bk[k],(1-Bk[k]))
  result <- c(lik, prior, ref)
  result[coeff==0] <- 0
  result <- result*coeff
  result <- sum(result)
  return(result)
}


#' S3 method for printing ssMCMC objects
#'
#' @param x An ssMCMC object
#' @param ... Optional arguments passed to print
#'
#' @return **No return value**, called for **side effects**.
#' This function prints a summary of an `ssMCMC` object, including the estimated marginal likelihood,
#' the number of power posteriors run, and the log marginal likelihood contributions across steps.
#'
#' @export
#' @method print ssMCMC
print.ssMCMC <- function(x, ...){
  cat("Stepping stone estimation of marginal likelihood\n")
  cat("Marginal Likelihood:\n")
  print(x$lnr, ...)
  cat(paste("A total of ", length(x$Bk), " power posteriors were run along the sequence: ",paste(round(x$Bk,5), collapse="\t\t"), "\n", sep=""))
  cat("lnr_k", round(unlist(x$lnrk),2))

}
#' S3 method for plotting ssMCMC objects
#'
#' @param x An 'ssMCMC' object
#' @param ... Additional arguments passed to \code{plot}
#'
#' @details Produces 4 plots. The first 3 plot the prior, reference function and likelihood. Different colors
#' indicate different power posteriors for each. These chains should appear to be well mixed. The final plot
#' shows the sum of the marginal likelihood across each of the steps in the stepping stone algorithm.
#' @return No return value, called for side effects. This function generates
#' diagnostic plots for an `ssMCMC` object, including log-likelihood, log-prior,
#' log-reference function, and cumulative marginal likelihood across power posteriors.

#' @export
#' @method plot ssMCMC
plot.ssMCMC <- function(x, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  par(mfrow=c(2,2))
  if(is.null(attributes(x)$burnin)){
    start <- 1
  } else {
    start <- round(attributes(x)$burnin*length(x$chains[[1]][[1]]),0)
  }
  postburn <- start:length(x$chains[[1]][[1]])
  lnL <- lapply(x$chains, function(x) x$lnL[postburn])
  rangelnL <- c(min(unlist(lnL))-2, max(unlist(lnL))+2)
  plot(0,0,type="n", xlim=c(0,length(unlist(lnL))), ylim=rangelnL,xaxt="n",xlab="",ylab="lnL", main="lnL",...)
  xindex <- lapply(1:length(lnL), function(x) (x-1)*length(lnL[[1]]) + 1:length(lnL[[1]]))
  sapply(1:length(lnL), function(x) lines(xindex[[x]], lnL[[x]], col=x))
  abline(v=seq(0,length(unlist(lnL)), length.out=length(lnL)+1),lty=2)

  pr <- lapply(x$chains, function(x) x$prior[postburn])
  rangepr <- c(min(unlist(pr))-2, max(unlist(pr))+2)
  plot(0,0,type="n", xlim=c(0,length(unlist(pr))), ylim=rangepr,xaxt="n",xlab="",ylab="Ln prior", main="ln prior",...)
  xindex <- lapply(1:length(pr), function(x) (x-1)*length(pr[[1]]) + 1:length(pr[[1]]))
  sapply(1:length(pr), function(x) lines(xindex[[x]], pr[[x]], col=x))
  abline(v=seq(0,length(unlist(pr)), length.out=length(pr)+1),lty=2)

  ref <- lapply(x$chains, function(x) x$ref[postburn])
  rangeref <- c(min(unlist(ref))-2, max(unlist(ref))+2)
  plot(0,0,type="n", xlim=c(0,length(unlist(ref))), ylim=rangeref,xaxt="n",xlab="",ylab="Ln ref", main="ln ref",...)
  xindex <- lapply(1:length(ref), function(x) (x-1)*length(ref[[1]]) + 1:length(ref[[1]]))
  sapply(1:length(ref), function(x) lines(xindex[[x]], ref[[x]], col=x))
  abline(v=seq(0,length(unlist(ref)), length.out=length(ref)+1),lty=2)

  plot(x$Bk, c(0, cumsum(x$lnrk)), ylab="ln r", xlab="power posterior",pch=21, bg=1:length(ref),cex=1.5, ...)
  lines(x$Bk, c(0,cumsum(x$lnrk)))
}

.pull.rsample <- function(samp, chain){
  #pars.list <- lapply(samp,function(y) pull.pars(y,chain,model=model))
  #emap.list <- lapply(samp,function(y) read.emap(chain$branch.shift[[y]],chain$location[[y]],chain$t2[[y]],cache$phy)$emap)
  L <- chain$lnL[samp]+chain$prior[samp]-chain$ref[samp]
  Lmax <- max(L)
  Lfactored <- L-Lmax
  return(list(Lmax=Lmax,Lfactored=Lfactored))
}

## Compute marginal likelihood
##
## \code{computelnr} computes the marginal likelihood of a set of chains estimated via stepping stone
## sampling and produced by the function \code{steppingstone}
.computelnr <- function(Kchains,Bk,samp){
  lnr <- list()
  for(i in 1:(length(Bk)-1)){
    Lk <- .pull.rsample(samp, Kchains[[i]])
    lnr[[i]] <- (Bk[i+1]-Bk[i])*Lk$Lmax+log(1/length(Lk$Lfactored)*sum(exp(Lk$Lfactored)^(Bk[i+1]-Bk[i])))
  }
  return(list("lnr"=sum(unlist(lnr)),"lnrk"=lnr))
}




