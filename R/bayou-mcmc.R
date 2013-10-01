#SE=0; model="OU"; ngen=10000; samp=10;chunk=100; control=NULL; tuning=NULL; new.dir=TRUE; plot.freq=500; outname="bayou"; ticker.freq=1000; tuning.int=c(0.1,0.2,0.3); startpar=NULL; moves=NULL; control.weights=NULL

#' Bayesian sampling of multi-optima OU models 
#' 
#' @description Runs a reversible-jump Markov chain Monte Carlo on continuous phenotypic data on a phylogeny, sampling possible shift locations and 
#' shift magnitudes, and shift numbers.
#' 
#' @param tree a phylogenetic tree of class 'phylo'
#' @param dat a named vector of continuous trait values matching the tips in tree
#' @param SE The standard error of the data. Either a single value applied to all the data, or a vector of length(dat).
#' @param model The parameterization of the OU model used. Either "OU" for standard parameterization with \alpha and \sigma^2; "OUrepar" for phylogenetic half-life and \V_y, or 
#' "QG" for the Lande model, with parameters \h^2 (heritability), \P (phenotypic variance), \omega^2 (width of adaptive landscape), and \Ne (effective population size)
#' @param prior A prior function of class 'priorFn' that gives the prior distribution of all parameters
#' @param ngen The number of gneerations to run the Markov Chain
#' @param samp The frequency at which Markov samples are retained
#' @param chunk The number of samples retained in memory before being written to a file
#' @param control A list providing a control object governing how often and which proposals are used
#' @param tuning A named vector that governs how liberal or conservative proposals are that equals the number of proposal mechanisms.
#' @param new.dir If TRUE, then results are stored in a new temporary directory. If FALSE, results are written to the current working directory.
#' @param plot.freq How often plots should be made during the mcmc. If NULL, then plots are not produced
#' @param outname The prefix given to files created by the mcmc
#' @param ticker.freq How often a summary log should be printed to the screen
#' @param tuning.int How often the tuning parameters should be adjusted as a fraction of the total number of generations (currently ignored)
#' @param startpar A list with the starting parameters for the mcmc. If NULL, starting parameters are simulated from the prior distribution
#' @param moves A named list providing the proposal functions to be used in the mcmc. Names correspond to the parameters to be modified in the parameter list. See 'details' for default values.
#' @param control.weights A named vector providing the relative frequency each proposal mechanism is to be used during the mcmc
#' 
#' @details 
#' By default, the alpha, sig2 (and various reparameterizations of these parameters) are adjusted with multiplier proposals, theta are adjusted with sliding window proposals,
#' and the number of shifts is adjusted by splitting and merging, as well as sliding the shifts both within and between branches. Allowed shift locations are specified by the 
#' prior function (see \link(make.prior())). 

bayou.mcmc <- function(tree, dat, SE=0, model="OU", prior, ngen=10000, samp=10, chunk=100, control=NULL, tuning=NULL, new.dir=FALSE, plot.freq=500, outname="bayou", ticker.freq=1000, tuning.int=c(0.1,0.2,0.3), startpar=NULL, moves=NULL, control.weights=NULL){
  if(is.null(moves)){
    moves <- switch(model,"QG"=list(h2=".multiplierProposal",P=".multiplierProposal",w2=".multiplierProposal",Ne=".multiplierProposal",k=".splitmerge",theta=".adjustTheta",slide=".slide"),
                          "OU"=list(alpha=".multiplierProposal",sig2=".multiplierProposal",k=".splitmerge",theta=".adjustTheta",slide=".slide"),
                          "OUrepar"=list(halflife=".multiplierProposal",Vy=".multiplierProposal",k=".splitmerge",theta=".adjustTheta",slide=".slide")) #,"OUcpp"=list(alpha=".multiplierProposal",sig2=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"), "QGcpp"=list(h2=".multiplierProposal",P=".multiplierProposal",w2=".multiplierProposal",Ne=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"),"OUreparcpp"=list(halflife=".multiplierProposal",Vy=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"))
  }
  
  cache <- .prepare.ou.univariate(tree,dat)
  
  if(is.null(startpar)){
    startpar <- priorSim(prior,tree,model,nsim=1,plot=FALSE, exclude.branches=NULL)$pars[[1]]
  }

  if(is.null(control.weights)){
    ct <- .buildControl(startpar, prior, default.weights=model, move.weights=NULL)
  } else {
    ct <- .buildControl(startpar, prior, default.weights=NULL, move.weights=control.weights)
  }
  
  if(is.null(tuning)){
    D <- switch(model, "OU"=list(alpha=1, sig2= 1, k = 4,theta=2,slide=1), "QG"=list(h2=1, P=1, w2=1, Ne=1, k = 4, theta=2, slide=1), "OUrepar"=list(halflife=1, Vy=1, k=4, theta=2, slide=1))#,"OUcpp"=list(alpha=1, sig2= 1,sig2jump=2, k = 4,theta=2,slide=1),"QGcpp"=list(h2=1,P=1,w2=1,Ne=1,sig2jump=2,k=4,theta=2,slide=1),"OUreparcpp"=list(halflife=1,Vy=1,sig2jump=2,k=4,theta=2,slide=1))
  } else {D <- tuning}
  
  if(new.dir){
    dir.name <- paste(sample(LETTERS,10),collapse="")
    dir <- paste(tempdir(),"/",dir.name,"/",sep="")
    dir.create(dir)
    } else {
    dir <- paste(getwd(),"/",sep="")
    dir.name <- sapply(strsplit(dir,'/'), function(x) x[length(x)])
  }
  
  mapsb <<- file(paste(dir, outname,".sb",sep=""),open="w")
  mapsloc <<- file(paste(dir, outname,".loc",sep=""),open="w")
  mapst2 <<- file(paste(dir, outname,".t2",sep=""),open="w")
  pars.output <<- file(paste(dir, outname,".pars",sep=""),open="w")
  
  oldpar <- startpar
  store <- list("out"=NULL, "sb"=NULL, "loc"=NULL, "t2"=NULL)

  
  lik.fn <- .OU.lik
  oll  <- lik.fn(oldpar, cache, dat, SE, model=model)$loglik
  pr1 <- prior(oldpar,cache)
  parorder <- switch(model,"QG"=c("h2","P","w2","Ne","k","ntheta","theta"), "OU"=c("alpha","sig2","k","ntheta","theta"),"OUrepar"=c("halflife","Vy","k","ntheta","theta"),"OUcpp"=c("alpha","sig2","sig2jump","k","ntheta","theta"))#,"QGcpp"=c("h2","P","w2","Ne","sig2jump","k","ntheta","theta"),"OUreparcpp"=c("halflife","Vy","sig2jump","k","ntheta","theta"))
  
  accept.type <- NULL
  accept <- NULL
  if(!is.null(plot.freq)){
    tr <- .toSimmap(.pars2map(oldpar, cache),cache)
    tcols <- makeTransparent(rainbow(oldpar$ntheta),alpha=100)
    names(tcols)<- 1:oldpar$ntheta
    phenogram(tr,dat,colors=tcols,ftype="off")
    plot.dim <- list(par('usr')[1:2],par('usr')[3:4])
  }
  #tuning.int <- round(tuning.int*ngen,0)
  for (i in 1:ngen){
    ct <- .updateControl(ct, oldpar)
    u <- runif(1)
    prop <- .proposalFn(u,ct,D,moves,cache,oldpar)
    new.pars <- prop$pars
    #new.cache <- prop$prop$cache
    accept.type <- c(accept.type,paste(prop$decision,prop$move,sep="."))
    pr2 <- prior(new.pars,cache)
    hr <- prop$hr
    new <- lik.fn(new.pars, cache, dat, SE=SE, model=model)
    nll <- new$loglik
    if (runif(1) < exp(nll-oll+pr2-pr1+hr)){
      oldpar <- new.pars
      pr1 <- pr2
      oll <- nll
      accept <- c(accept,1)
    } else {
      accept <- c(accept,0)
    }
    store <- store.bayOU(i, oldpar, oll, pr1, store, samp, chunk, parorder)
    if(!is.null(plot.freq)){
      if(i %% plot.freq==0){
        tr <- .toSimmap(.pars2map(oldpar, cache),cache)
        tcols <- makeTransparent(rainbow(oldpar$ntheta),alpha=100)
        names(tcols)<- 1:oldpar$ntheta
        plot(plot.dim[[1]],plot.dim[[2]],type="n",xlab="time",ylab="phenotype")
        mtext(paste("gens = ",i," lnL = ",round(oll,2)),3)
        #regime.plot probably doesn't work for simmaps
        #try(regime.plot(oldpar,tr$tree,tcols,type="density",model=model),silent=TRUE)
        phenogram(tr,dat,colors=tcols,ftype="off",add=TRUE)
      }
    }
    #if(i %in% tuning.int){
    #  D <- tune.D(D,accept,accept.type)$D
    #}
    if(i%%ticker.freq==0){
      alpha <- switch(model,"QG"=QG.alpha(oldpar),"OU"=oldpar$alpha,"OUrepar"=OU.repar(oldpar)$alpha,"OUcpp"=oldpar$alpha,"QGcpp"=QG.alpha(oldpar),"OUrepar"=OU.repar(oldpar)$alpha)
      sig2 <- switch(model,"QG"=QG.sig2(oldpar),"OU"=oldpar$sig2,"OUrepar"=OU.repar(oldpar)$sig2,"OUcpp"=oldpar$sig2,"QGcpp"=QG.sig2(oldpar),"OUrepar"=OU.repar(oldpar)$sig2)
      tick <- c(i,oll,pr1,log(2)/alpha,sig2/(2*alpha),oldpar$k,tapply(accept,accept.type,mean))
      tick[-1] <- round(tick[-1],2)
      names(tick)[1:6] <- c('gen','lnL','prior','half.life','Vy','K')
      if(i==ticker.freq){
        cat(c(names(tick),'\n'),sep='\t\t\t')
      }
      cat(c(tick,'\n'),sep='\t\t\t')
    }
  }
  closeAllConnections()
  out <- list('model'=model, 'dir.name'=dir.name,'dir'=dir, 'outname'=outname, 'accept'=accept,'accept.type'=accept.type)
  class(out) <- c("bayouFit", "list")
  return(out)
}

