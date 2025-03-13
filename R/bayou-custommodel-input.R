getPreValues <- function(cache, col){
  V <- phytools::vcvPhylo(cache$phy, anc.nodes=FALSE)
  X <- cache$pred[,col]
  unknown <- is.na(X)
  known <- !unknown
  Vkk <- V[known, known]
  Vuu <- V[unknown, unknown]
  Vku <- V[known, unknown]
  Vuk <- V[unknown, known]
  iVkk <- solve(Vkk)
  sigmabar <- as.matrix(Matrix::forceSymmetric(Vuu - Vuk%*%iVkk%*%Vku))
  cholSigmabar <- chol(sigmabar)
  mubarmat <- Vuk%*%iVkk
  return(list(V=V, X=X, unknown=unknown, known=known, Vkk=Vkk, Vuu=Vuu, Vku=Vku, Vuk=Vuk, iVkk=iVkk, sigmabar=sigmabar, mubarmat=mubarmat, cholSigmabar=cholSigmabar))
}

# Proposal function to simulate conditional draws from a multivariate normal distribution
.imputePredBM <- function(cache, pars, d, move,ct=NULL, prevalues=pv, prior=prior){
  #(tree, dat, sig2, plot=TRUE, ...){
  X <- prevalues$X
  Vuk <- pars$pred.sig2*prevalues$Vuk
  iVkk <- (1/pars$pred.sig2)*prevalues$iVkk
  Vku <- pars$pred.sig2*prevalues$Vku
  Vuu <- pars$pred.sig2*prevalues$Vuu
  known <- prevalues$known
  unknown <- prevalues$unknown
  mu <- rep(pars$pred.root, cache$n)
  muk <- mu[known]
  muu <- mu[unknown]
  mubar <- t(muu + Vuk%*%iVkk%*%(X[known]-muk))
  sigmabar <- Vuu - Vuk%*%iVkk%*%Vku
  res <- MASS::mvrnorm(1, mubar, sigmabar)
  pars.new <- pars
  pars.new$missing.pred <- res
  hr=Inf
  type="impute"
  return(list(pars=pars.new, hr=hr, decision = type))
}

.make.monitorFn <- function(model, noMonitor=c("missing.pred", "ntheta"), integers=c("gen","k"), verbose=TRUE){
  parorder <- model$parorder
  rjpars <- model$rjpars
  exclude <- which(parorder %in% noMonitor)
  if(length(exclude) > 0){
    pars2monitor <- parorder[-exclude]
  } else {pars2monitor <- parorder}
  if(length(rjpars) > 0){
    rjp <- which(pars2monitor %in% rjpars)
    pars2monitor[rjp] <- paste("r", pars2monitor[rjp], sep="")
  }
  pars2monitor <- c("gen", "lnL", "prior", pars2monitor)
  type <- rep(".2f", length(pars2monitor))
  type[which(pars2monitor %in% integers)] <- "i"
  string <- paste(paste("%-8", type, sep=""), collapse="")
  monitor.fn = function(i, lik, pr, pars, accept, accept.type, j){
    names <- pars2monitor
    #names <- c("gen", "lnL", "prior", "alpha" , "sig2", "rbeta1", "endo", "k")
    #string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
    acceptratios <- unlist(accept/accept.type) #tapply(accept, accept.type, mean)
    names <- c(names, names(acceptratios))
    if (verbose){
      if(j==0){
        cat(sprintf("%-7.7s", names), "\n", sep=" ")
      }
      cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$beta1[1], pars$endo, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
    }
  }
}

.getTipMap <- function(pars, cache){
  map <- .pars2map(pars,cache)
  tipreg <- rev(map$theta)
  ntipreg <- rev(map$branch)
  #ntipreg <- names(map$theta)
  dups <- !duplicated(ntipreg) & ntipreg %in% (1:nrow(cache$edge))[cache$externalEdge]
  tipreg <- tipreg[which(dups)]
  ntipreg <- ntipreg[which(dups)]
  o <- order(cache$edge[as.numeric(ntipreg), 2])
  betaID <- tipreg[o]
}


.colorRamp <- function(trait, pal, nn){
  strait <- (trait-min(trait))/max(trait-min(trait))
  itrait <- round(strait*nn, 0)+1
  return(pal(nn+1)[itrait])
}

.addColorBar <- function(x, y, height, width, pal, trait, ticks, adjx=0, n=100,cex.lab=1,pos=2, text.col="black"){
  legend_image <- as.raster(matrix(rev(pal(n)),ncol = 1))
  #text(x = 1.5, y = round(seq(range(ave.Div)[1], range(ave.Div)[2], l = 5), 2), labels = seq(range(ave.Div)[1], range(ave.Div)[2], l = 5))
  seqtrait <- seq(min(trait), max(trait), length.out=nrow(legend_image))
  mincut <- n-which(abs(seqtrait - min(ticks))==min(abs(seqtrait-min(ticks))))
  maxcut <- n-which(abs(seqtrait - max(ticks))==min(abs(seqtrait-max(ticks))))
  legend_cut <- legend_image[maxcut:mincut,]
  legend_cut <- rbind(matrix(rep(legend_image[1,1],round(0.05*n,0)),ncol=1), legend_cut)
  rasterImage(legend_cut, x, y, x+width, y+height)
  ticklab <- format(ticks, digits=2, trim=TRUE)
  ticklab[length(ticklab)] <- paste(">", ticklab[length(ticklab)], sep="")
  text(x+adjx, y=seq(y, y+height, length.out=length(ticks)), labels=ticklab, pos=pos,cex=cex.lab, col=text.col)
}


#' This function makes a bayou model object that can be used for customized allometric regression models.
#'
#' @param f A formula describing the relationship between the data and one or more predictors (use 'dat'
#' for the dependent variable)
#' @param rjpars A character vector of parameters to split at the mapped shifts on the tree
#' @param tree A phylogenetic tree
#' @param dat A named vector of trait data (dependent variable)
#' @param pred A matrix or data frame with named columns with predictor data represented in the specified
#' formula
#' @param prior A prior function made by the 'make.prior' function
#' @param SE A single value or vector of measurement error estimates
#' @param impute The name of a single predictor for which missing values will be imputed using BM (see details).
#' Default is NULL.
#' @param startpar An optional list of starting parameters for the model. If not provided, the model will simulate
#' starting values from the prior function.
#' @param moves An optional list of moves to be passed on to bayou.makeMCMC.
#' @param control.weights An optional list of control weights to be passed on to bayou.makeMCMC.
#' @param D A vector of tuning parameters to be passed on to bayou.makeMCMC.
#' @param shiftpars The names of the parameters defining the map of shifts (for now, always c("sb", "loc", "t2")).
#' @param model The parameterization of the OU model, either "OU", "OUrepar" or "QG".
#' @param slopechange "immediate", "alphaWeighted" or "fullPGLS"
#'
#' @details This function generates a list with the '$model', which provides the specifications of the regression
#' model and '$startpar', which provides starting values to input into bayou.makeMCMC. Note that this model assumes
#' that predictors immediately affect trait values at a shift. In other words, regardless of the past history of the
#' predictor, only the current value affects the current expected trait value. This is only reasonable for allometric
#' models, although it may be appropriate for other models if phylogenetic inertia is very low (short half-lives).
#'
#' One predictor variable may include missing data (coded as "NA"). The model will assume the maximum-likelihood
#' best-fit BM model and simulate the missing predictor values throughout the course of the MCMC. These values will
#' then be used to calculate the likelihood given the parameters for each MCMC step.
#' @return
#' A list with two elements:
#' \describe{
#'   \item{model}{A list containing MCMC settings, likelihood functions, monitoring functions, and other parameters for the Bayesian OU model.}
#'   \item{startpar}{A list of starting parameter values for the model. If missing data is imputed, it includes `missing.pred`.}
#' }
#' @export
makeBayouModel <- function(f, rjpars, tree, dat, pred, prior, SE=0, slopechange="immediate", impute=NULL, startpar=NULL, moves=NULL, control.weights=NULL, D=NULL, shiftpars=c("sb", "loc", "t2"), model="OU"){
  cache <- .prepare.ou.univariate(tree, dat, SE=SE, pred=pred)
  vars <- terms(f)
  cache$pred <- as.data.frame(cache$pred)
  dep <-  rownames(attr(vars, "factors"))[attr(vars, "response")]
  mf <- cbind(cache$dat, cache$pred)
  colnames(mf)[1] <- dep
  MF <- model.frame(f, data=mf, na.action=na.pass)
  MM <- model.matrix(f, MF)
  colnames(MM) <- gsub(":", "x", colnames(MM))
  parnames <- paste("beta", colnames(MM)[-1], sep="_")
  if(length(rjpars) > 0){
    rjpars2 <- c(rjpars, paste("beta", rjpars, sep="_"))
    rj <- which(colnames(MM) %in% rjpars2)-1
    if(slopechange=="alphaWeighted"){
      expFn <- function(pars, cache){
        W <- C_weightmatrix(cache, pars)$W
        if(length(impute)>0){
          MF[is.na(MF[,impute]),impute] <- pars$missing.pred #$impute
          MM <- model.matrix(f, MF)
        }
        parframe <- lapply(pars[parnames], function(x) return(x))
        parframe[rj] <- lapply(parframe[rj], function(x) W%*%x)
        ExpV <- apply(sapply(1:length(parframe), function(x) parframe[[x]]*MM[,x+1]), 1, sum)
        return(ExpV)
      }
    } else {
      expFn <- function(pars, cache){
        betaID <- .getTipMap(pars, cache)
        if(length(impute)>0){
          MF[is.na(MF[,impute]),impute] <- pars$missing.pred #$impute
          MM <- model.matrix(f, MF)
        }
        parframe <- lapply(pars[parnames], function(x) return(x))
        parframe[rj] <- lapply(parframe[rj], function(x) x[betaID])
        ExpV <- apply(sapply(1:length(parframe), function(x) parframe[[x]]*MM[,x+1]), 1, sum)
        return(ExpV)
      }
    }
  } else {
    rjpars2 <- numeric(0)
    expFn <- function(pars, cache){
      #betaID <- getTipMap(pars, cache)
      if(length(impute)>0){
        MF[is.na(MF[,impute]),impute] <- pars$missing.pred #$impute
        MM <- model.matrix(f, MF)
      }
      parframe <- lapply(pars[parnames], function(x) return(x))
      #parframe[rjpars] <- lapply(parframe[rjpars], function(x) x[betaID])
      ExpV <- apply(sapply(1:length(parframe), function(x) parframe[[x]]*MM[,x+1]), 1, sum)
      return(ExpV)
    }
  }
  varnames <- switch(model, "OU"=c("alpha","sig2"), "OUrepar"=c("halflife", "Vy"), "QG"=c("h2", "P", "Ne", "w2"))
  if(model=="OU"){
    varTransform <- function(pars) return(pars)
  }
  if(model=="OUrepar"){
    varTransform <- function(pars){
      repar <- OU.repar(pars)
      pars$alpha <- repar$alpha
      pars$sig2 <- repar$sig2
      return(pars)
    }
  }
  if(model=="QG"){
    varTransform <- function(pars){
      pars$alpha <- QG.alpha(pars)
      pars$sig2 <- QG.sig2(pars)
      return(pars)
    }
  }
  likFn <- function(pars, cache, X, model="Custom"){
    n <- cache$n
    X <- cache$dat
    pred <- cache$pred
    ## Permit alternative OU parameterizations
    pars <- varTransform(pars)
    ## Specify the model here
    X = X - expFn(pars, cache)
    cache$dat <- X
    ### The part below mostly does not change
    pars2 <- pars
    if(slopechange=="fullPGLS"){pars2$alpha = 1e10}
    X.c <- C_weightmatrix(cache, pars2)$resid
    transf.phy <- C_transf_branch_lengths(cache, 1, X.c, pars$alpha)
    transf.phy$edge.length[cache$externalEdge] <- transf.phy$edge[cache$externalEdge] + cache$SE[cache$phy$edge[cache$externalEdge, 2]]^2*(2*pars$alpha)/pars$sig2
    comp <- C_threepoint(list(n=n, N=cache$N, anc=cache$phy$edge[, 1], des=cache$phy$edge[, 2], diagMatrix=transf.phy$diagMatrix, P=X.c, root=transf.phy$root.edge, len=transf.phy$edge.length))
    if(pars$alpha==0){
      inv.yVy <- comp$PP/pars$sig2
      detV <- comp$logd + n*log(pars$sig2)
    } else {
      inv.yVy <- comp$PP*(2*pars$alpha)/(pars$sig2)
      detV <- comp$logd+n*log(pars$sig2/(2*pars$alpha))
    }
    llh <- -0.5*(n*log(2*pi)+detV+inv.yVy)
    return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
  }
  monitorFn <- function(i, lik, pr, pars, accept, accept.type, j, verbose=TRUE){
    names <- c("gen", "lnL", "prior", varnames, parnames, "rtheta", "k")
    format <- c("%-8i",rep("%-8.2f",4), rep("%-8.2f", length(parnames)), "%-8.2f","%-8i")
    acceptratios <- unlist(accept/accept.type) #tapply(accept, accept.type, mean)
    names <- c(names, names(acceptratios))
    if (verbose) {
      if(j==0){
        cat(sprintf("%-7.7s", names), "\n", sep=" ")
      }
      item <- c(i, lik, pr, pars[[varnames[1]]], pars[[varnames[2]]], sapply(pars[parnames], function(x) x[1]), pars$theta[1], pars$k)
      cat(sapply(1:length(item), function(x) sprintf(format[x], item[x])), sprintf("%-8.2f", acceptratios),"\n", sep="")
    }
  }
  rdists <- .getSimDists(prior)
  ## Set default moves if not specified.
  if(length(rjpars) > 0){
    if(is.null(moves)){
      moves =  c(switch(model, "OU" = list(alpha=".multiplierProposal", sig2=".multiplierProposal"),
                        "OUrepar" = list(halflife=".jointHalflifeVyProposal", Vy=".jointHalflifeVyProposal")),
                 as.list(setNames(rep(".vectorSlidingWindowSplit", length(parnames)), parnames)),
                 c(theta=".vectorSlidingWindowSplit", k=".splitmergePrior", slide=".slide2"))
    }
    if(is.null(control.weights)){
      control.weights <- setNames(rep(1, length(parnames)+5), c(varnames, parnames, "k", "theta", "slide"))
      control.weights[c(varnames[1], parnames)] <- 2
      control.weights[c("theta", parnames[rj])] <- 10
      control.weights["k"] <- 5
      control.weights <- as.list(control.weights)
    }

    if(is.null(D)){
      D <- lapply(rdists[!(names(rdists) %in% shiftpars)], function(x) sd(x(1000))/50)
      D$k <- rep(1, length(rjpars))
      D$slide <- 1
    }
    {
    parorder <- c(varnames, parnames,"theta", "k","ntheta")
    rjord <- which(parorder %in% rjpars2)
    fixed <- names(attributes(prior)$fixed)
    if(length(rjord > 0)){
      parorder <- c(parorder[-rjord], fixed , parorder[rjord])
    } else {
      parorder <- c(parorder, fixed)
    }
    parorder <- parorder[!duplicated(parorder) & !(parorder %in% shiftpars)]

    }
    if(is.null(startpar)){
      startpar <- priorSim(prior, cache$phy, shiftpars=rjpars2)$pars[[1]]
      startpar <- startpar[c(parorder, shiftpars)]
      #simdists <- rdists[parorder[!(parorder %in% c(rjpars2,shiftpars, "ntheta"))]]
      #if(length(attributes(prior)$fixed)>0){
      #  simdists[names(attributes(prior)$fixed)] <- lapply(1:length(attributes(prior)$fixed), function(x) function(n) attributes(prior)$fixed[[x]])
      #  fixed.pars <- attributes(prior)$fixed
      #  fixed <- TRUE
      #} else {fixed <- FALSE}
      #simdists <- simdists[!is.na(names(simdists))]
      #startpar <- lapply(simdists, function(x) x(1))
      #startpar$ntheta <- startpar$k+1
      #startpar[parorder[(parorder %in% c(rjpars2))]] <- lapply(rdists[parorder[(parorder %in% c(rjpars2))]], function(x) x(startpar$ntheta))
      #startpar <- c(startpar, list(sb=sample(1:length(cache$bdesc), startpar$k, replace=FALSE, prob = sapply(cache$bdesc, length)), loc=rep(0, startpar$k), t2=2:startpar$ntheta))
      #startpar <- startpar[c(parorder, shiftpars)]
    }
  } else {
    rj <- numeric(0)
    if(is.null(moves)){
      moves =  c(switch(model, "OU" = list(alpha=".multiplierProposal", sig2=".multiplierProposal"),
                        "OUrepar" = list(halflife=".jointHalflifeVyProposal", Vy=".jointHalflifeVyProposal")),
                 as.list(setNames(rep(".vectorSlidingWindowSplit", length(parnames)), parnames)),
                 c(theta=".vectorSlidingWindowSplit"))
    }
    if(is.null(control.weights)){
      control.weights <- setNames(rep(1, length(parnames)+5), c(varnames, parnames, "k", "theta", "slide"))
      control.weights[c(varnames[1], parnames)] <- 2
      control.weights[c("theta", parnames[rj])] <- 6
      control.weights[c("k","slide")] <- 0
      control.weights <- as.list(control.weights)
    }
    if(is.null(D)){
      D <- lapply(rdists[!(names(rdists) %in% shiftpars)], function(x) sd(x(1000))/50)
      D$k <- 1
      D$slide <- 1
    }
    {
    parorder <- c(varnames, parnames,"theta", "k","ntheta")
    rjord <- which(parorder %in% rjpars2)
    fixed <- names(attributes(prior)$fixed)
    if(length(rjord > 0)){
      parorder <- c(parorder[-rjord], fixed , parorder[rjord])
    } else {
      parorder <- c(parorder, fixed)
    }
    parorder <- parorder[!duplicated(parorder) & !(parorder %in% shiftpars)]
    }
    if(is.null(startpar)){
      startpar <- priorSim(prior, cache$phy, shiftpars = rjpars2)$pars[[1]]
      startpar <- startpar[c(parorder, shiftpars)]
      #simdists <- rdists[parorder[!(parorder %in% c(rjpars2,shiftpars, "ntheta"))]]
      #if(length(attributes(prior)$fixed)>0){
      #  simdists[names(attributes(prior)$fixed)] <- lapply(1:length(attributes(prior)$fixed), function(x) function(n) attributes(prior)$fixed[[x]])
      #  fixed.pars <- attributes(prior)$fixed
      #}
      #simdists <- simdists[!is.na(names(simdists))]
      #startpar <- lapply(simdists, function(x) x(1))
      #if(!("k" %in% fixed)){
      #  startpar$k <- 0
      #  startpar$ntheta <- startpar$k+1
      #  if(startpar$k==0) startpar$t2 <- numeric(0) else startpar$t2 <- 2:(startpar$ntheta)
      #} else {
      #  startpar$ntheta <- startpar$k+1
      #  if(startpar$k==0) startpar$t2 <- numeric(0) else startpar$t2 <- 2:(startpar$ntheta)
      #}
      #startpar <- startpar[c(parorder, shiftpars)]
    }
  }
  rjpars[!(rjpars %in% "theta")] <- paste("beta",rjpars[!(rjpars %in% "theta")], sep="_")
  model <- list(moves=moves, control.weights=control.weights, D=D, rjpars=rjpars, parorder=parorder, shiftpars=shiftpars, monitor.fn=monitorFn, call=f, expFn=expFn, lik.fn=likFn)
  if(length(impute)>0){
    missing <- which(is.na(cache$pred[,impute])) #$impute
    pv <- getPreValues(cache, impute) #$impute
    model$moves$missing.pred <- ".imputePredBM"
    model$control.weights$missing.pred <- 1
    model$D$missing.pred <- 1
    startpar <- .imputePredBM(cache, startpar, d=1, NULL, ct=NULL, prevalues=pv)$pars#$impute
    bp <- which(names(startpar)=="pred.root")
    model$parorder <- c(parorder[1:bp], "missing.pred", if(length(parorder)>bp)parorder[(bp+1):length(parorder)] else NULL)
    startpar <- startpar[c(parorder, names(startpar)[!names(startpar) %in% parorder])]
    model$prevalues <- pv
  }

  #try(prior(startpar))
  #try(likFn(startpar, cache, cache$dat))
  return(list(model=model, startpar=startpar))
}


.getSimDists <- function(prior){
  dists <- attributes(prior)$dist
  fixed <- which(attributes(prior)$dist=="fixed")
  notfixed <- which(attributes(prior)$dist!="fixed")
  dists <- dists[notfixed]
  prior.params <- attributes(prior)$param
  rdists <- lapply(dists,function(x) gsub('^[a-zA-Z]',"r",x))
  prior.params <- lapply(prior.params,function(x) x[-which(names(x)=="log")])
  rdists.fx <- lapply(rdists,get)
  rdists.fx <- lapply(1:length(rdists.fx),function(x) .set.defaults(rdists.fx[[x]],defaults=prior.params[[x]]))
  names(rdists.fx) <- gsub('^[a-zA-Z]',"",names(rdists))
  return(rdists.fx)
}
