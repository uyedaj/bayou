#SE=0; ngen=1000; samp=10; chunk=100; control=NULL; tuning=NULL; new.dir=TRUE; plot.freq=NULL; outname="bayou"; ticker.freq=1000; tuning.int=NULL; moves=NULL; control.weights=NULL; lik.fn=NULL; plot.fn=NULL
#model <- model.Impute; plot.fn <- NULL
#startpar=list(alpha=0.1, sig2=3, beta1=1, k=1, ntheta=2, theta=c(4,4), sb=200, loc=0, t2=2)
# Bayesian sampling of multi-optima OU models
#
# @description Runs a reversible-jump Markov chain Monte Carlo on continuous phenotypic data on a phylogeny,
# sampling possible shift locations and shift magnitudes, and shift numbers.
#
# @param tree a phylogenetic tree of class 'phylo'
# @param dat a named vector of continuous trait values matching the tips in tree
# @param SE The standard error of the data. Either a single value applied to all
# the data, or a vector of length(dat).
# @param model The parameterization of the OU model used. Either "OU" for standard parameterization with
# alpha and sigma^2; "OUrepar" for phylogenetic half-life and stationary variance (Vy), or "QG" for the
# Lande model, with parameters h^2 (heritability), P (phenotypic variance), omega^2 (width of adaptive
# landscape), and Ne (effective population size)
# @param prior A prior function of class 'priorFn' that gives the prior distribution of all parameters
# @param ngen The number of generations to run the Markov Chain
# @param samp The frequency at which Markov samples are retained
# @param chunk The number of samples retained in memory before being written to a file
# @param control A list providing a control object governing how often and which proposals are used
# @param tuning A named vector that governs how liberal or conservative proposals are that equals the
# number of proposal mechanisms.
# @param new.dir If TRUE, then results are stored in a new temporary directory. If FALSE, results are
# written to the current working directory. If a character string,
# then results are written to that working directory.
# @param plot.freq How often plots should be made during the mcmc. If NULL, then plots are not produced
# @param outname The prefix given to files created by the mcmc
# @param plot.fn Function used in plotting, defaults to phytools::phenogram
# @param ticker.freq How often a summary log should be printed to the screen
# @param tuning.int How often the tuning parameters should be adjusted as a fraction of the total
# number of generations (currently ignored)
# @param startpar A list with the starting parameters for the mcmc. If NULL, starting parameters are
# simulated from the prior distribution
# @param moves A named list providing the proposal functions to be used in the mcmc. Names correspond to
# the parameters to be modified in the parameter list. See 'details' for default values.
# @param control.weights A named vector providing the relative frequency each proposal mechanism is to
# be used during the mcmc
# @param lik.fn Likelihood function to be evaluated. Defaults to \code{bayou.lik}.
# @param verbose Determines whether information is outputted to the console for the user to view
#
# @return
# For `bayou.mcmc`, a list of class `"bayouFit"` containing:
# \describe{
#   \item{model}{The model parameterization used.}
#   \item{dir.name}{Directory where results are stored.}
#   \item{dir}{Full directory path for stored results.}
#   \item{outname}{Filename prefix for output files.}
#   \item{accept}{Vector of acceptance rates for different proposals.}
#   \item{accept.type}{Vector indicating types of accepted proposals.}
#   \item{tree}{The phylogenetic tree used.}
#   \item{dat}{Continuous trait data used in the MCMC simulation.}
#   \item{tmpdir}{Logical indicating if a temporary directory was used.}
#   \item{startpar}{List of starting parameter values.}
# }
#
# For `bayou.makeMCMC`, a list of class `"bayouMCMCFn"` containing:
# \describe{
#   \item{model}{The model parameterization used.}
#   \item{model.pars}{List of model parameters used in the MCMC simulation.}
#   \item{dir.name}{Directory where results are stored.}
#   \item{dir}{Full directory path for stored results.}
#   \item{outname}{Filename prefix for output files.}
#   \item{tree}{The phylogenetic tree used.}
#   \item{dat}{Continuous trait data used in the MCMC simulation.}
#   \item{pred}{Predictor data matrix or `NULL` if not used.}
#   \item{SE}{Measurement error values provided.}
#   \item{tmpdir}{Logical indicating if a temporary directory was used.}
#   \item{startpar}{List of starting parameter values.}
#   \item{run}{Function to execute the MCMC simulation.}
#   \item{steppingstone}{Function for performing stepping-stone sampling for marginal likelihood estimation.}
#   \item{load}{Function to load the MCMC chain after the simulation.}
# }
#
# @name bayou-deprecated
# @section \code{bayou.mcmc}: This function is deprecated, please use \code{\link{bayou.makeMCMC}}.
# @useDynLib bayou
# @details
# By default, the alpha, sig2 (and various reparameterizations of these parameters) are adjusted with
# multiplier proposals, theta are adjusted with sliding window proposals,
# and the number of shifts is adjusted by splitting and merging, as well as sliding the shifts both
# within and between branches. Allowed shift locations are specified by the
# prior function (see make.prior()).

#model="bd"; tree <- phy; SE=0; ngen=1000; samp=10; chunk=100; control=NULL;tuning=NULL; new.dir=TRUE;plot.freq=100; outname="bayou";ticker.freq=1000; tuning.int=c(0.1,0.2,0.3); startpar=pars; moves=NULL; control.weights=NULL; lik.fn <- bdSplit.lik; plot.fn <- plotSimmap

bayou.mcmc <- function(tree, dat, SE=0, model="OU", prior, ngen=10000, samp=10, chunk=100,
                       control=NULL, tuning=NULL, new.dir=FALSE, plot.freq=500, outname="bayou",
                       plot.fn=phenogram, ticker.freq=1000, tuning.int=c(0.1,0.2,0.3), startpar=NULL,
                       moves=NULL, control.weights=NULL, lik.fn=NULL, verbose=TRUE){
  .Deprecated("bayou.makeMCMC")
  if(FALSE){
  fixed <- gsub('^[a-zA-Z]',"",names(attributes(prior)$distributions)[which(attributes(prior)$distributions=="fixed")])
  if("loc" %in% fixed){
    fixed <- c(fixed,"slide")
  }
  if(is.null(moves)){
    moves <- switch(model,"QG"=list(h2=".multiplierProposal",P=".multiplierProposal",w2=".multiplierProposal",Ne=".multiplierProposal",k=".splitmerge",theta=".adjustTheta",slide=".slide"),
                          "OU"=list(alpha=".multiplierProposal",sig2=".multiplierProposal",k=".splitmerge",theta=".adjustTheta",slide=".slide"),
                          "OUrepar"=list(halflife=".multiplierProposal",Vy=".multiplierProposal",k=".splitmerge",theta=".adjustTheta",slide=".slide"),
                          "bd"=list(r=".vectorMultiplier", eps=".vectorMultiplier", k=".splitmergebd")) #,"OUcpp"=list(alpha=".multiplierProposal",sig2=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"), "QGcpp"=list(h2=".multiplierProposal",P=".multiplierProposal",w2=".multiplierProposal",Ne=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"),"OUreparcpp"=list(halflife=".multiplierProposal",Vy=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"))
    moves <- moves[which(!(names(moves) %in% fixed))]
  }

  cache <- .prepare.ou.univariate(tree,dat, SE=SE)
  dat <- cache$dat
  if(is.null(startpar)){
    if(any(fixed %in% c("h2", "P", "w2", "Ne", "halflife", "Vy"))){
      stop(paste("Parameters '", paste(fixed[fixed %in% c("h2", "P", "w2", "Ne", "halflife", "Vy")], collapse=" "), "' are set to be fixed but no starting values are supplied.
                   Please specify starting parameter values",sep=""))
    }
    startpar <- priorSim(prior,cache$phy,model,nsim=1,plot=FALSE, exclude.branches=NULL)$pars[[1]]
    if(length(fixed)>0){
      assumed <- sapply(fixed, function(x) switch(x, "slide"="", "sb"="sb=numeric(0)", "k"= "k=0", "alpha"="alpha=0", "sig2"="sig2=0", "loc"="0.5*edge.length"))
      if(verbose){
        print(paste("Warning: Fixed parameters '", paste(fixed,collapse=", "), "' not specified, assuming values: ", paste(assumed,collapse=", "),sep="" ))
      }
    }
  }
  if(length(fixed)==0 & is.null(control.weights)){
      control.weights <- switch(model,"OU"=list("alpha"=4,"sig2"=2,"theta"=4,"slide"=2,"k"=10),"QG"=list("h2"=5,"P"=2,"w2"=5,"Ne"=5,"theta"=5,"slide"=3,"k"=20),"OUrepar"=list("halflife"=5,"Vy"=3,"theta"=5,"slide"=3,"k"=20), "bd"=list("r"=2, "eps"=1, "k"=5, slide=0))
      ct <- .buildControl(startpar, prior, control.weights)
    } else {
      if(is.null(control.weights)){
        control.weights <- switch(model,"OU"=list("alpha"=4,"sig2"=2,"theta"=4,"slide"=2,"k"=10),"QG"=list("h2"=5,"P"=2,"w2"=5,"Ne"=5,"theta"=5,"slide"=3,"k"=20),"OUrepar"=list("halflife"=5,"Vy"=3,"theta"=5,"slide"=3,"k"=20), "bd"=list("r"=2, "eps"=1, "k"=5, slide=0))
        #"OUcpp"=list("alpha"=3,"sig2"=3,"sig2jump"=3,"theta"=3,"slide"=5,"k"=10),"QGcpp"=list("h2"=1,"P"=1,"w2"=2,"Ne"=2,"sig2jump"=3,"theta"=3,"slide"=5,"k"=10),"OUreparcpp"=list("halflife"=3,"Vy"=3,"sig2jump"=3,"theta"=3,"slide"=5,"k"=10)
        control.weights[fixed[fixed %in% names(control.weights)]] <- 0
      } else {control.weights <- control.weights}
      ct <- .buildControl(startpar, prior, move.weights=control.weights)
    }
  if(is.null(tuning)){
    D <- switch(model, "OU"=list(alpha=1, sig2= 1, k = 4,theta=2,slide=1), "QG"=list(h2=1, P=1, w2=1, Ne=1, k = 4, theta=2, slide=1), "OUrepar"=list(halflife=1, Vy=1, k=4, theta=2, slide=1), "bd"=list(r=1, eps=1, k=4))#,"OUcpp"=list(alpha=1, sig2= 1,sig2jump=2, k = 4,theta=2,slide=1),"QGcpp"=list(h2=1,P=1,w2=1,Ne=1,sig2jump=2,k=4,theta=2,slide=1),"OUreparcpp"=list(halflife=1,Vy=1,sig2jump=2,k=4,theta=2,slide=1))
  } else {D <- tuning}
  if(is.logical(new.dir)){
    if(new.dir){
      dir.name <- paste(sample(LETTERS,10),collapse="")
      dir <- paste(tempdir(),"/",dir.name,"/",sep="")
      dir.create(dir)
      } else {
      dir <- paste(getwd(),"/",sep="")
      dir.name <- sapply(strsplit(dir,'/'), function(x) x[length(x)])
    }
  } else {
    if(is.character(new.dir)){
      dir.name <- paste(sample(LETTERS,10),collapse="")
      dir <- paste(new.dir,"/",dir.name,"/",sep="")
      dir.create(dir)
    }
  }


  #mapsb <<- file(paste(dir, outname,".sb",sep=""),open="w")
  #mapsloc <<- file(paste(dir, outname,".loc",sep=""),open="w")
  #mapst2 <<- file(paste(dir, outname,".t2",sep=""),open="w")
  #pars.output <<- file(paste(dir, outname,".pars",sep=""),open="w")
  files <- list(mapsb=file(paste(dir, outname,".sb",sep=""),open="a"),
                mapsloc=file(paste(dir, outname,".loc",sep=""),open="a"),
                mapst2=file(paste(dir, outname,".t2",sep=""),open="a"),
                pars.output=file(paste(dir, outname,".pars",sep=""),open="a"))

  oldpar <- startpar
  store <- list("out"=list(), "sb"=list(), "loc"=list(), "t2"=list())

  if(is.null(lik.fn)) lik.fn <- bayou.lik #.OU.lik
  oll  <- lik.fn(oldpar, cache, dat, model=model)$loglik
  pr1 <- prior(oldpar,cache)
  parorder <- switch(model,"QG"=c("h2","P","w2","Ne","k","ntheta","theta"), "OU"=c("alpha","sig2","k","ntheta","theta"),"OUrepar"=c("halflife","Vy","k","ntheta","theta"),"OUcpp"=c("alpha","sig2","sig2jump","k","ntheta","theta"), "bd"=c("r", "eps", "k", "ntheta"))#,"QGcpp"=c("h2","P","w2","Ne","sig2jump","k","ntheta","theta"),"OUreparcpp"=c("halflife","Vy","sig2jump","k","ntheta","theta"))

  accept.type <- NULL
  accept <- NULL
  if(!is.null(plot.freq)){
    tr <- .toSimmap(.pars2map(oldpar, cache),cache)
    tcols <- makeTransparent(rainbow(oldpar$ntheta),alpha=200)
    names(tcols)<- 1:oldpar$ntheta
    plot.fn(tr,dat,colors=tcols,ftype="off")
    plot.dim <- list(par('usr')[1:2],par('usr')[3:4])
  }
  #tuning.int <- round(tuning.int*ngen,0)
  mcmc.loop <- function(verbose=TRUE){
    for (i in 1:ngen){
      ct <- .updateControl(ct, oldpar, fixed)
      u <- runif(1)
      prop <- .proposalFn(u,ct,D,moves,cache,oldpar)
      new.pars <- prop$pars
      #new.cache <- prop$prop$cache
      accept.type <- c(accept.type,paste(prop$decision,prop$move,sep="."))
      pr2 <- prior(new.pars,cache)
      hr <- prop$hr
      new <- lik.fn(new.pars, cache, dat, model=model)
      nll <- new$loglik
      nll <- ifelse(is.na(nll), -Inf, nll)
      if (runif(1) < exp(nll-oll+pr2-pr1+hr)){
        oldpar <- new.pars
        pr1 <- pr2
        oll <- nll
        accept <- c(accept,1)
      } else {
        accept <- c(accept,0)
      }
      store <- .store.bayou(i, oldpar, oll, pr1, store, samp, chunk, parorder,files)
      if(!is.null(plot.freq)){
        if(i %% plot.freq==0){
          tr <- .toSimmap(.pars2map(oldpar, cache),cache)
          tcols <- makeTransparent(rainbow(oldpar$ntheta),alpha=200)
          names(tcols)<- 1:oldpar$ntheta
          plot(plot.dim[[1]],plot.dim[[2]],type="n",xlab="time",ylab="phenotype")
          mtext(paste("gens = ",i," lnL = ",round(oll,2)),3)
          #regime.plot probably doesn't work for simmaps
          #try(regime.plot(oldpar,tr$tree,tcols,type="density",model=model),silent=TRUE)
          plot.fn(tr,dat,colors=tcols,ftype="off",add=TRUE)
        }
      }
      #if(i %in% tuning.int){
      #   D <- tune.D(D,accept,accept.type)$D
      #}
      if(i%%ticker.freq==0){
        alpha <- switch(model,"QG"=QG.alpha(oldpar),"OU"=oldpar$alpha,"OUrepar"=OU.repar(oldpar)$alpha,"OUcpp"=oldpar$alpha,"QGcpp"=QG.alpha(oldpar),"OUrepar"=OU.repar(oldpar)$alpha, "bd"=oldpar$r)
        sig2 <- switch(model,"QG"=QG.sig2(oldpar),"OU"=oldpar$sig2,"OUrepar"=OU.repar(oldpar)$sig2,"OUcpp"=oldpar$sig2,"QGcpp"=QG.sig2(oldpar),"OUrepar"=OU.repar(oldpar)$sig2, "bd"=oldpar$eps)
        if(model=="bd"){
          tick <- c(i,oll,pr1,median(alpha),median(sig2),oldpar$k,tapply(accept,accept.type,mean))
          tick[-1] <- round(tick[-1],2)
          names(tick)[1:6] <- c('gen','lnL','prior','r','eps','K')
        } else {
          tick <- c(i,oll,pr1,log(2)/alpha,sig2/(2*alpha),oldpar$k,tapply(accept,accept.type,mean))
          tick[-1] <- round(tick[-1],2)
          names(tick)[1:6] <- c('gen','lnL','prior','half.life','Vy','K')
        }
        if(verbose) {
          if(i==ticker.freq ){
            cat(c(names(tick),'\n'),sep='\t\t\t')
          }
          cat(c(tick,'\n'),sep='\t\t\t')
        }
      }
    }
  }
  mcmc.loop()
  lapply(files, close)
  out <- list('model'=model, 'dir.name'=dir.name,'dir'=dir, 'outname'=outname, 'accept'=accept,'accept.type'=accept.type, 'tree'=tree, 'dat'=dat, 'tmpdir'=ifelse(new.dir==TRUE, TRUE, FALSE), 'startpar'=startpar)
  class(out) <- c("bayouFit", "list")
  return(out)
  }
}

#' Builds a bayouMCMCFn object that can run MCMC and stepping stone analyses.
#'
#' @description Runs a reversible-jump Markov chain Monte Carlo on continuous phenotypic data
#' on a phylogeny, sampling possible shift locations and
#' shift magnitudes, and shift numbers.
#'
#' @param tree a phylogenetic tree of class 'phylo'
#' @param dat a named vector of continuous trait values matching the tips in tree
#' @param pred A matrix or data frame with named columns with predictor data represented in the
#' specified
#' formula
#' @param SE The standard error of the data. Either a single value applied to all the data, or a
#' vector of length(dat).
#' @param model The parameterization of the OU model used. Either "OU" for standard parameterization
#' with alpha and sigma^2; "OUrepar" for phylogenetic half-life and stationary variance (Vy), or
#' "QG" for the Lande model, with parameters h^2 (heritability), P (phenotypic variance), omega^2
#' (width of adaptive landscape), and Ne (effective population size)
#' @param prior A prior function of class 'priorFn' that gives the prior distribution of all
#' parameters
#' @param samp The frequency at which Markov samples are retained
#' @param chunk The number of samples retained in memory before being written to a file
#' @param control A list providing a control object governing how often and which proposals are used
#' @param tuning A named vector that governs how liberal or conservative proposals are that equals
#' the number of proposal mechanisms.
#' @param file.dir If a character string, then results are written to that working directory. If NULL,
#' then results are not saved to files, but instead held in memory. Default is `tempdir()`, which
#' writes to an R temporary directory.
#' @param plot.freq How often plots should be made during the mcmc. If NULL, then plots are not
#' produced
#' @param outname The prefix given to files created by the mcmc
#' @param ticker.freq How often a summary log should be printed to the screen
#' @param plot.fn Function used in plotting, defaults to phytools::phenogram
#' @param startpar A list with the starting parameters for the mcmc. If NULL, starting parameters
#' are simulated from the prior distribution
#' @param moves A named list providing the proposal functions to be used in the mcmc. Names correspond
#'  to the parameters to be modified in the parameter list. See 'details' for default values.
#' @param control.weights A named vector providing the relative frequency each proposal mechanism is
#' to be used during the mcmc
#' @param lik.fn Likelihood function to be evaluated. Defaults to \code{bayou.lik}.
#' @param perform.checks A logical indicating whether to use bayou.checkModel to validate model inputs.
#' @param verbose A logical that determines whether information is outputted to the console for the user to view
#'
#' @useDynLib bayou
#' @export
#' @details
#' By default, the alpha, sig2 (and various reparameterizations of these parameters) are adjusted
#' with multiplier proposals, theta are adjusted with sliding window proposals,
#' and the number of shifts is adjusted by splitting and merging, as well as sliding the shifts
#' both within and between branches. Allowed shift locations are specified by the
#' prior function (see make.prior()).
#' @return A bayouMCMCFn object (list).

#model="bd"; tree <- phy; SE=0; ngen=1000; samp=10; chunk=100; control=NULL;tuning=NULL; new.dir=TRUE;plot.freq=100; outname="bayou";ticker.freq=1000; tuning.int=c(0.1,0.2,0.3); startpar=pars; moves=NULL; control.weights=NULL; lik.fn <- bdSplit.lik; plot.fn <- plotSimmap
#tree <- phy; dat <- dat; pred=NULL; SE=0; model="OU"; prior=priorOU; samp=10; chunk=100; control=NULL;tuning=NULL;file.dir=NULL;plot.freq=500;outname="bayou";plot.fn=phenogram; ticker.freq=1000;startpar=NULL; moves=NULL; control.weights=NULL; lik.fn=NULL; perform.checks=TRUE;
bayou.makeMCMC <- function(tree, dat, pred=NULL, SE=0, model="OU", prior, samp=10, chunk=100,
                            control=NULL, tuning=NULL, file.dir=tempdir(), plot.freq=500, outname="bayou",
                            plot.fn=phenogram, ticker.freq=1000, startpar=NULL, moves=NULL,
                            control.weights=NULL, lik.fn=NULL, perform.checks=TRUE, verbose=TRUE){

  #Check if custom or predefined model
  if(is.character(model)){
    model.pars <- switch(model, "OU"=model.OU, "QG"=model.QG, "OUrepar"=model.OUrepar)
  } else {
    model.pars <- model
    model <- "Custom"
  }

  #Perform input checks to see if everything looks like it is specified correctly
  if(perform.checks == TRUE){
    if(model=="Custom") {
      checks <- bayou.checkModel(pars=startpar, tree=tree, dat=dat, pred=pred, SE=SE, prior=prior, model=model.pars, autofix=TRUE, verbose=verbose)
    } else {
      checks <- bayou.checkModel(pars=startpar, tree=tree, dat=dat, pred=pred, SE=SE, prior=prior, model=model, autofix=TRUE, verbose=verbose)
    }
    cache <- checks$autofixed$cache
    startpar <- checks$autofixed$pars
    tree <- checks$autofixed$cache$phy
    dat <- checks$autofixed$cache$dat
    pred <- checks$autofixed$pred
    model.pars <- checks$autofixed$model
    if(all(sapply(checks[-length(checks)], function(x) x)==TRUE)){
      if(verbose){
        cat("seems fine...\n")
      }
    }
  } else {
    # Skip the checks (caution)
    cache <- .prepare.ou.univariate(tree, dat, SE=SE, pred=pred)
    dat <- cache$dat
    if(is.null(startpar)){
      startpar <- priorSim(prior, cache$phy, model, nsim=1, plot=FALSE, shiftpars=model.pars$rjpars, exclude.branches=NULL)$pars[[1]]
    }
  }

  #Check and see if there is anything that is fixed
  fixed <- gsub('^[a-zA-Z]',"",names(attributes(prior)$distributions)[which(attributes(prior)$distributions=="fixed")])
  if("loc" %in% fixed){
    fixed <- c(fixed,"slide")
  }

  #Use predefined moves if specified
  if(is.null(moves)){
    moves <- model.pars$moves
    moves <- moves[which(!(names(moves) %in% fixed))]
  }

  pv <- model.pars$prevalues
  parorder <- model.pars$parorder
  rjpars <- model.pars$rjpars
  outpars <- parorder[which(!(parorder %in% rjpars | parorder %in% model.pars$shiftpars))]

  if(length(fixed)==0 & is.null(control.weights)){
    control.weights <- model.pars$control.weights
    ct <- .buildControl(startpar, prior, control.weights)
  } else {
    if(is.null(control.weights)){
      control.weights <-model.pars$control.weights
      control.weights[fixed[fixed %in% names(control.weights)]] <- 0
    } else {control.weights <- control.weights}
    ct <- .buildControl(startpar, prior, move.weights=control.weights)
  }

  if(is.null(tuning)){
    D <- model.pars$D
  } else {D <- tuning}
  attributes(ct)$splitmergepars <- rjpars

  if(is.null(lik.fn)) lik.fn <- model.pars$lik.fn #.OU.lik

  if(is.null(file.dir)){
    if(verbose){
      cat("\nFile directory specified as NULL. Results will be returned as output.\n")
    }
    oldpar <- startpar
    chunk <- Inf
    files <- NULL
    dir <- NULL
    dir.name <- NULL
    i <- 1
  } else {
    chunk <- 1
    if(file.dir==tempdir()){
      dir.name <- paste(sample(LETTERS,10),collapse="")
      dir <- paste(tempdir(),"/",dir.name,"/",sep="")
      dir.create(dir)
    } else {
      dir.name <- NULL
      dir <- paste(file.dir,"/",sep="")
      if(!dir.exists(dir)){
        dir.create(file.dir)
      }
    }
    filenames <- list(mapsb=paste(dir, outname,".sb",sep=""),
                      mapsloc=paste(dir, outname,".loc",sep=""),
                      mapst2=paste(dir, outname,".t2",sep=""),
                      pars.output=paste(dir, outname,".pars",sep=""),
                      rjpars=paste(dir, outname,".rjpars", sep=""))

    if(any(sapply(filenames, file.exists))){
      warning("Files with given outname already exist in directory '", dir, "'.\n
              Runs will be appended to existing files and starting parameters from last entry used.\n")
      files <- list(mapsb=file(filenames$mapsb,open="a"),
                    mapsloc=file(filenames$mapsloc,open="a"),
                    mapst2=file(filenames$mapst2,open="a"),
                    pars.output=file(filenames$pars.output,open="a"),
                    rjpars=file(filenames$rjpars, open="a"))
      startinf <- .lastpar(files, outpars, startpar, rjpars)
      oldpar <- startinf$pars
      i <- startinf$i
    } else {
      files <- list(mapsb=file(filenames$mapsb,open="w"),
                    mapsloc=file(filenames$mapsloc,open="w"),
                    mapst2=file(filenames$mapst2,open="w"),
                    pars.output=file(filenames$pars.output,open="w"),
                    rjpars=file(filenames$rjpars, open="w"))

      oldpar <- startpar
      oll  <- lik.fn(oldpar, cache, dat, model=model)$loglik
      pr1 <- prior(oldpar,cache)
      i <- 1
      store <- .store.bayou2(1, oldpar, outpars, rjpars, oll, pr1,
                             list("out"=list(), "rjpars"=list(), "sb"=list(), "loc"=list(), "t2"=list()),
                             1, 1, parorder, files)
    }
  }
  oll  <- lik.fn(oldpar, cache, dat, model=model)$loglik
  pr1 <- prior(oldpar,cache)
  store <- list("out"=list(), "rjpars"=list(), "sb"=list(), "loc"=list(), "t2"=list())
  store <- .store.bayou2(i, oldpar, outpars, rjpars, oll, pr1, store, 1, chunk, parorder, files)
  header <- 0
  gbg <- lapply(files, close)
  out <- list('model'=model, 'model.pars'=model.pars, 'dir.name'=dir.name,'dir'=dir, 'outname'=outname, 'tree'=tree, 'dat'=dat, 'pred'=pred, 'SE'=SE, 'tmpdir'=(file.dir==tempdir()), 'startpar'=startpar)

  if(!is.null(plot.freq)){
    environment(plot.fn) <- new.env()
    set.runpars(plot.fn, runpars=list(oldpar=oldpar, cache=cache, i=0))
    tr <- .toSimmap(.pars2map(oldpar, cache),cache)
    tcols <- makeTransparent(rainbow(oldpar$ntheta),alpha=200)
    names(tcols)<- 1:oldpar$ntheta
    plot.fn(tr,dat,colors=tcols,ftype="off")
    plot.dim <- list(par('usr')[1:2],par('usr')[3:4])
  }

  #Define MCMC loop function
  mcmc.loop <- function(ngen, j, oll, pr1, startpar, store, ref1=NULL, type="mcmc", ref=NULL, Bk=NULL, k=NULL, files=NULL, verbose=TRUE){
    if(!is.null(files)){
      closeAllConnections()
      files <- list(mapsb=file(files$mapsb,open="a"),
                    mapsloc=file(files$mapsloc,open="a"),
                    mapst2=file(files$mapst2,open="a"),
                    pars.output=file(files$pars.output,open="a"),
                    rjpars=file(files$rjpars, open="a"))
    }
    if(type=="ss"){
      pB.old <- powerPosteriorFn(k, Bk, oll, pr1, ref1)
    }
    iseq <- (j+1):(j+ngen)
    acceptNames <- unlist(lapply(1:length(moves), function(x) paste(attr(get(moves[[x]]), "types"), names(moves)[x], sep=".")),F,F)
    accept <- acceptN <- data.frame(matrix(rep(0, length(acceptNames)), ncol=length(acceptNames) ))
    names(accept) <- names(acceptN) <- acceptNames
    for (i in iseq){
      ct <- .updateControl(ct, oldpar, fixed)
      u <- runif(1)
      prop <- .proposalFn(u,ct,D,moves,cache,oldpar,prior)
      new.pars <- prop$pars
      propName <- paste(prop$decision,prop$move,sep=".")
      acceptN[propName] <- acceptN[propName]+1
      pr2 <- prior(new.pars,cache)
      hr <- prop$hr
      new <- lik.fn(new.pars, cache, dat, model=model)
      nll <- new$loglik
      nll <- ifelse(is.na(nll), -Inf, nll)
      if(type=="mcmc"){
        aR <- exp(nll-oll+pr2-pr1+hr)
      } else{
        if(type=="ss"){
          ref2 <- ref(new.pars, cache)
          pB.new <- powerPosteriorFn(k, Bk, nll, pr2, ref2)
          bhr <- ifelse(!is.finite(hr), hr, Bk[k]*hr)
          aR <- exp(pB.new-pB.old+bhr)
        } else {
          if(type=="prior.only"){
            aR <- exp(pr2-pr1+hr)
          }
        }

      }

      if (runif(1) < aR){
        oldpar <- new.pars
        pr1 <- pr2
        oll <- nll
        if(type=="ss"){
          ref1 <- ref2
          pB.old <- powerPosteriorFn(k, Bk, oll, pr1, ref1)
        }
        accept[propName] <- accept[propName]+1
      } #else {#accept <- c(accept,0)}

      if(i %% samp == 0){
        store <- .store.bayou2(i, oldpar, outpars, rjpars, oll, pr1, store, 1, chunk, parorder, files, ref=ref1)
      }
      if(!is.null(plot.freq)){
        if(i %% plot.freq==0){
          set.runpars(plot.fn, runpars=list(oldpar=oldpar, cache=cache, i=i))
          tr <- .toSimmap(.pars2map(oldpar, cache),cache)
          tcols <- makeTransparent(rainbow(oldpar$ntheta),alpha=200)
          names(tcols)<- 1:oldpar$ntheta
          plot(plot.dim[[1]],plot.dim[[2]],type="n",xlab="time",ylab="phenotype")
          mtext(paste("gens = ",i," lnL = ",round(oll,2)),3)
          plot.fn(tr,dat,colors=tcols,ftype="off",add=TRUE)
        }
      }
      if(i%%ticker.freq==0 && verbose){
        model.pars$monitor.fn(i, oll, pr1, oldpar, accept, acceptN, header)
        header <- 1
      }
    }
    if(is.null(files)){
      if(type=="ss"){
        chain <- .fchain.local(store, out, ss=TRUE)
      } else {
        chain <- .fchain.local(store, out)
      }
      return(chain)
    } else {
      gbg <- lapply(files,close)
      return(NULL)
    }
  }
  run.mcmc <- function(ngen, prior.only=FALSE, verbose=TRUE){
    if(!is.null(files)){
      closeAllConnections()
      files <- list(mapsb=file(filenames$mapsb,open="a"),
                    mapsloc=file(filenames$mapsloc,open="a"),
                    mapst2=file(filenames$mapst2,open="a"),
                    pars.output=file(filenames$pars.output,open="a"),
                    rjpars=file(filenames$rjpars, open="a"))
      startinf <- .lastpar(files, outpars=outpars, startpar=oldpar, rjpars=rjpars, ss=FALSE)
      oldpar <- startinf$pars
      j <- startinf$i
      oll <- startinf$oll
      pr1 <- startinf$pr1
    } else {
      j <- 1
      files <- NULL
      filenames <- NULL
    }
    if(prior.only==TRUE){
      res <- tryCatch(mcmc.loop(ngen, j=j, oll=oll, pr1=pr1, startpar=oldpar, store=store, ref1 = NULL, type="prior.only", ref=NULL, files=filenames, verbose=verbose))
    } else {
      res <- tryCatch(mcmc.loop(ngen, j=j, oll=oll, pr1=pr1, startpar=oldpar, store=store, ref1 = NULL, type="mcmc", ref=NULL, files=filenames, verbose=verbose))
    }
    #if(!is.null(files)){gbg <- lapply(files, close)}
    return(res)
  }

  run.steppingstone <- function(ngen, chain, Bk, burnin=0.3, plot=TRUE, verbose=TRUE){
    postburn <- floor(max(c(1,burnin*length(chain$gen)))):length(chain$gen)
    ref <- make.refFn(chain, model=model.pars, priorFn=prior, burnin = burnin, plot)
    stores <- ollss <- prss <- refss <- startpars <- list()
    if(!is.null(files)){
      closeAllConnections()
      ssfilenames <- lapply(1:length(Bk), function(y) lapply(filenames, function(x) gsub(paste(outname, ".", sep=""), paste(outname, "_ss", y,".", sep=""), x)))
      ssfiles <- list()
      for(x in 1:length(Bk)){
        ssfiles[[x]] <- list(mapsb=file(ssfilenames[[x]]$mapsb,open="w"),
                             mapsloc=file(ssfilenames[[x]]$mapsloc,open="w"),
                             mapst2=file(ssfilenames[[x]]$mapst2,open="w"),
                             pars.output=file(ssfilenames[[x]]$pars.output,open="w"),
                             rjpars=file(ssfilenames[[x]]$rjpars, open="w"))
        startpars[[x]] <- pull.pars(sample(postburn, 1, replace=FALSE), chain, model=model.pars)
        ollss[[x]] <- lik.fn(startpars[[x]], cache, dat, model=model)$loglik
        prss[[x]] <- prior(startpars[[x]], cache)
        refss[[x]] <- ref(startpars[[x]], cache)
        stores[[x]] <- list("out"=list(), "rjpars"=list(), "sb"=list(), "loc"=list(), "t2"=list())
        stores[[x]] <- .store.bayou2(1, startpars[[x]], outpars, rjpars, ollss[[x]], prss[[x]], stores[[x]], 1, 1, parorder, ssfiles[[x]], ref=refss[[x]])
        gbg <- lapply(ssfiles[[x]], close)
      }
    } else {
      for(x in 1:length(Bk)){
        startpars[[x]] <- pull.pars(sample(postburn, 1, replace=FALSE), chain, model=model.pars)
        ollss[[x]] <- lik.fn(startpars[[x]], cache, dat, model=model)$loglik
        prss[[x]] <- prior(startpars[[x]], cache)
        refss[[x]] <- ref(startpars[[x]], cache)
        stores[[x]] <- list("out"=list(), "rjpars"=list(), "sb"=list(), "loc"=list(), "t2"=list())
        stores[[x]] <- .store.bayou2(1, startpars[[x]], outpars, rjpars, ollss[[x]], prss[[x]], stores[[x]], 1, Inf, parorder, files=NULL, ref=refss[[x]])
      }
      ssfiles <- NULL
      ssfilenames <- NULL
    }
    ssfits <- foreach(j=1:length(Bk), .export = c("mcmc.loop")) %dopar% mcmc.loop(ngen, j=1, ollss[[j]], prss[[j]], startpars[[j]], stores[[j]], ref1=refss[[x]], type="ss", ref=ref, Bk=Bk, k=j, files=ssfilenames[[j]], verbose=verbose) #mcmc.loop <- function(ngen, j, oll, pr1, startpar, store, ref1=NULL, type="mcmc", ref=NULL){
    if(!is.null(ssfiles)){
      outs <- lapply(1:length(Bk), function(x){out$outname <- paste(outname, "_ss", x, sep=""); out})
      chains <- lapply(1:length(Bk), function(x) load.bayou(outs[[x]], saveRDS=FALSE, file=NULL, cleanup=FALSE, ref=TRUE))
    } else {
      chains <- ssfits
      ssfits <- "Results not saved to files"
    }
    postburn <- floor(max(c(1, burnin*length(chains[[1]]$gen)))):( min(sapply(chains, function(x) length(x$gen))) )
    lnr <- .computelnr(chains, Bk, postburn)
    ssres <- list(chains=chains, lnr=lnr$lnr, lnrk=lnr$lnrk, Bk=Bk, fits=ssfits, filenames=ssfilenames, startpars=startpars, refFn=ref)
    class(ssres) <- c("ssMCMC", "list")
    return(ssres)
  }
  out$run <- run.mcmc
  out$steppingstone <- run.steppingstone
  mcmc.load <- function(saveRDS=FALSE, file=NULL, cleanup=FALSE){
    load.bayou(out, saveRDS, file, cleanup)
  }
  out$load <- mcmc.load
  class(out) <- c("bayouMCMCFn", "list")
  return(out)
}

# Function that grabs the last parameters from a set of bayou output files.
.lastpar <- function(files, outpars, startpar, rjpars, ss=FALSE){
  fLs <- sapply(files, function(x) countL(summary(x)$description))
  fL <- min(fLs[fLs > 0])
  if(fL==1){skipL <- 0} else {skipL=fL-1}
  res <- lapply(1:length(files), function(x) try(utils::read.table(summary(files[[x]])$description, skip=skipL), silent=TRUE))
  res <- lapply(1:length(files), function(x) if(inherits(res[[x]], "try-error")){numeric(0)}else{res[[x]]})
  pars <- list()
  parLs <- lapply(startpar, length)[outpars]
  npars <- length(res[[4]])
  j=4 + as.numeric(ss)
  if(length(outpars) > 0){
    for(i in 1:length(outpars)){
      pars[[outpars[i]]] <- unlist(res[[4]][,j:(j+parLs[[i]]-1)],F,F)
      j <- j+1+parLs[[i]]-1
    }
  }
  if(length(rjpars) > 0){
    j <- 1
    for(i in 1:length(rjpars)){
      pars[[rjpars[i]]] <- unlist((res[[5]][j:(j+pars$ntheta-1)]),F,F)
      j <- j+pars$ntheta
    }
  }
  pars$sb <- unlist(res[[1]], F, F)
  pars$loc <- unlist(res[[2]], F, F)
  pars$t2 <- unlist(res[[3]], F, F)
  i <- unlist(res[[4]][1], F, F)
  oll <- unlist(res[[4]][2],F,F)
  pr1 <- unlist(res[[4]][3], F,F)
  if(ss){
    ref1 <- unlist(res[[4]][4], F,F)
    return(list(pars=pars, i=i, oll=oll, pr1=pr1, ref1=ref1))
  } else {
    return(list(pars=pars, i=i, oll=oll, pr1=pr1))
  }

}

# Function to format a locally returned bayou store and format it as a bayouMCMC chain.
.fchain.local <- function(store, bayouFit, ss=FALSE){
  tree <- bayouFit$tree
  dat <- bayouFit$dat
  outname <- bayouFit$outname
  model <- bayouFit$model
  model.pars <- bayouFit$model.pars
  startpar <- bayouFit$startpar
  dir <- bayouFit$dir
  outpars <- model.pars$parorder[!(model.pars$parorder %in% model.pars$rjpars)]
  rjpars <- model.pars$rjpars
  sampled <- sapply(store$out,function(x) !is.null(x))
  chain <- list()
  pars.out <- store$out[sampled]
  chain$gen <- unname(sapply(pars.out,function(x) x[1]))
  chain$lnL <- unname(sapply(pars.out,function(x) x[2]))
  chain$prior <- unname(sapply(pars.out,function(x) x[3]))
  if(ss==TRUE){
    chain$ref <- unname(sapply(pars.out, function(x) x[4]))
  }
  parLs <- lapply(startpar, length)[outpars]
  j=4+as.numeric(ss)
  if(length(outpars) > 0){
    for(i in 1:length(outpars)){
      chain[[outpars[i]]] <- unname(lapply(pars.out, function(x) as.vector(x[j:(j+parLs[[i]]-1)])))
      if(parLs[[i]]==1) chain[[outpars[i]]]=unlist(chain[[outpars[i]]])
      j <- j+1+parLs[[i]]-1
    }
  }
  chain$sb <- store$sb[sampled]
  chain$loc <- store$loc[sampled]
  chain$t2 <- store$t2[sampled]
  rjpars.out <- store$rjpars[sampled]
  if(length(rjpars >0)){
    nrjpars <- length(rjpars)
    for(i in 1:length(rjpars)){
      chain[[rjpars[i]]] <- lapply(rjpars.out, function(x) unlist((x[(1+length(x)/nrjpars*(i-1)):(1+i*length(x)/nrjpars-1)]),F,F))
    }
  }
  attributes(chain)$model <- model
  attributes(chain)$model.pars <- model.pars
  attributes(chain)$tree <- tree
  attributes(chain)$dat <- dat
  class(chain) <- c("bayouMCMC", "list")
  return(chain)
}


countL <- function (file, chunkSize = 5e+07, ...) {
  if (inherits(file, "connection")) {
    con <- file
  }
  else {
    file <- as.character(file)
    con <- gzfile(file, open = "rb")
    on.exit(close(con))
  }
  LF <- as.raw(10)
  CR <- as.raw(13)
  SPC <- as.raw(32L)
  isLastCR <- isLastLF <- FALSE
  isEmpty <- TRUE
  nbrOfLines <- 0L
  while (TRUE) {
    bfr <- readBin(con = con, what = raw(), n = chunkSize)
    if (isLastCR) {
      if (bfr[1L] == LF)
        bfr[1L] <- SPC
    }
    n <- length(bfr)
    if (n == 0L)
      break
    isEmpty <- FALSE
    idxsCR <- which(bfr == CR)
    nCR <- length(idxsCR)
    if (nCR > 0L) {
      idxsCRLF <- idxsCR[(bfr[idxsCR + 1L] == LF)]
      if (length(idxsCRLF) > 0L) {
        bfr <- bfr[-idxsCRLF]
        n <- length(bfr)
        idxsCRLF <- NULL
        nCR <- length(which(bfr == CR))
      }
    }
    nLF <- length(which(bfr == LF))
    nbrOfLines <- nbrOfLines + (nCR + nLF)
    if (n == 0L) {
      isLastCR <- isLastLF <- FALSE
    }
    else {
      bfrN <- bfr[n]
      isLastCR <- (bfrN == CR)
      isLastLF <- (bfrN == LF)
    }
  }
  if (!isEmpty) {
    if (!isLastLF)
      nbrOfLines <- nbrOfLines + 1L
    attr(nbrOfLines, "lastLineHasNewline") <- isLastLF
  }
  nbrOfLines
}
set.runpars <- function(fun, runpars = list()){
  env <- environment(fun)
  for(i in 1:length(runpars)){
    env[[names(runpars)[i]]] <- runpars[[i]]
  }
}
#set.runpars(mymcmc, list(samp=20))
utils::globalVariables(c("j","pv"))
