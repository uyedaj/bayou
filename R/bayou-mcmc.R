#SE=0; ngen=1000; samp=10; chunk=100; control=NULL; tuning=NULL; new.dir=TRUE; plot.freq=NULL; outname="bayou"; ticker.freq=1000; tuning.int=NULL; moves=NULL; control.weights=NULL; lik.fn=NULL; plot.fn=NULL
#model <- model.Impute; plot.fn <- NULL
#startpar=list(alpha=0.1, sig2=3, beta1=1, k=1, ntheta=2, theta=c(4,4), sb=200, loc=0, t2=2)
#' Bayesian sampling of multi-optima OU models 
#' 
#' @description Runs a reversible-jump Markov chain Monte Carlo on continuous phenotypic data on a phylogeny, 
#' sampling possible shift locations and shift magnitudes, and shift numbers.
#' 
#' @param tree a phylogenetic tree of class 'phylo'
#' @param dat a named vector of continuous trait values matching the tips in tree
#' @param SE The standard error of the data. Either a single value applied to all 
#' the data, or a vector of length(dat).
#' @param model The parameterization of the OU model used. Either "OU" for standard parameterization with 
#' alpha and sigma^2; "OUrepar" for phylogenetic half-life and stationary variance (Vy), or "QG" for the 
#' Lande model, with parameters h^2 (heritability), P (phenotypic variance), omega^2 (width of adaptive 
#' landscape), and Ne (effective population size)
#' @param prior A prior function of class 'priorFn' that gives the prior distribution of all parameters
#' @param ngen The number of generations to run the Markov Chain
#' @param samp The frequency at which Markov samples are retained
#' @param chunk The number of samples retained in memory before being written to a file
#' @param control A list providing a control object governing how often and which proposals are used
#' @param tuning A named vector that governs how liberal or conservative proposals are that equals the 
#' number of proposal mechanisms.
#' @param new.dir If TRUE, then results are stored in a new temporary directory. If FALSE, results are 
#' written to the current working directory. If a character string,
#' then results are written to that working directory. 
#' @param plot.freq How often plots should be made during the mcmc. If NULL, then plots are not produced
#' @param outname The prefix given to files created by the mcmc
#' @param plot.fn Function used in plotting, defaults to phytools::phenogram
#' @param ticker.freq How often a summary log should be printed to the screen
#' @param tuning.int How often the tuning parameters should be adjusted as a fraction of the total 
#' number of generations (currently ignored)
#' @param startpar A list with the starting parameters for the mcmc. If NULL, starting parameters are 
#' simulated from the prior distribution
#' @param moves A named list providing the proposal functions to be used in the mcmc. Names correspond to 
#' the parameters to be modified in the parameter list. See 'details' for default values.
#' @param control.weights A named vector providing the relative frequency each proposal mechanism is to 
#' be used during the mcmc
#' @param lik.fn Likelihood function to be evaluated. Defaults to \code{bayou.lik}.
#' 
#' @useDynLib bayou
#' @export
#' @details 
#' By default, the alpha, sig2 (and various reparameterizations of these parameters) are adjusted with 
#' multiplier proposals, theta are adjusted with sliding window proposals,
#' and the number of shifts is adjusted by splitting and merging, as well as sliding the shifts both 
#' within and between branches. Allowed shift locations are specified by the 
#' prior function (see make.prior()). 

#model="bd"; tree <- phy; SE=0; ngen=1000; samp=10; chunk=100; control=NULL;tuning=NULL; new.dir=TRUE;plot.freq=100; outname="bayou";ticker.freq=1000; tuning.int=c(0.1,0.2,0.3); startpar=pars; moves=NULL; control.weights=NULL; lik.fn <- bdSplit.lik; plot.fn <- plotSimmap

bayou.mcmc <- function(tree, dat, SE=0, model="OU", prior, ngen=10000, samp=10, chunk=100, 
                       control=NULL, tuning=NULL, new.dir=FALSE, plot.freq=500, outname="bayou", 
                       plot.fn=phenogram, ticker.freq=1000, tuning.int=c(0.1,0.2,0.3), startpar=NULL, 
                       moves=NULL, control.weights=NULL, lik.fn=NULL){
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
      print(paste("Warning: Fixed parameters '", paste(fixed,collapse=", "), "' not specified, assuming values: ", paste(assumed,collapse=", "),sep="" ))
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
  mcmc.loop <- function(){  
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
        if(i==ticker.freq){
          cat(c(names(tick),'\n'),sep='\t\t\t')
        }
        cat(c(tick,'\n'),sep='\t\t\t')
      }
    }
  }
  mcmc.loop()
  lapply(files, close)
  out <- list('model'=model, 'dir.name'=dir.name,'dir'=dir, 'outname'=outname, 'accept'=accept,'accept.type'=accept.type, 'tree'=tree, 'dat'=dat, 'tmpdir'=ifelse(new.dir==TRUE, TRUE, FALSE), 'startpar'=startpar)
  class(out) <- c("bayouFit", "list")
  return(out)
}

#' Revision of bayou.mcmc that only makes the mcmc loop function, rather than running it itself. 
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
#' @param new.dir If TRUE, then results are stored in a new temporary directory. If FALSE, results 
#' are written to the current working directory. If a character string,
#' then results are written to that working directory. 
#' @param plot.freq How often plots should be made during the mcmc. If NULL, then plots are not 
#' produced
#' @param outname The prefix given to files created by the mcmc
#' @param ticker.freq How often a summary log should be printed to the screen
#' @param plot.fn Function used in plotting, defaults to phytools::phenogram
#' @param tuning.int How often the tuning parameters should be adjusted as a fraction of the total 
#' number of generations (currently ignored)
#' @param startpar A list with the starting parameters for the mcmc. If NULL, starting parameters 
#' are simulated from the prior distribution
#' @param moves A named list providing the proposal functions to be used in the mcmc. Names correspond
#'  to the parameters to be modified in the parameter list. See 'details' for default values.
#' @param control.weights A named vector providing the relative frequency each proposal mechanism is 
#' to be used during the mcmc
#' @param lik.fn Likelihood function to be evaluated. Defaults to \code{bayou.lik}.
#' @param perform.checks A logical indicating whether to use bayou.checkModel to validate model inputs.
#' 
#' @useDynLib bayou
#' @export
#' @details 
#' By default, the alpha, sig2 (and various reparameterizations of these parameters) are adjusted 
#' with multiplier proposals, theta are adjusted with sliding window proposals,
#' and the number of shifts is adjusted by splitting and merging, as well as sliding the shifts 
#' both within and between branches. Allowed shift locations are specified by the 
#' prior function (see make.prior()). 

#model="bd"; tree <- phy; SE=0; ngen=1000; samp=10; chunk=100; control=NULL;tuning=NULL; new.dir=TRUE;plot.freq=100; outname="bayou";ticker.freq=1000; tuning.int=c(0.1,0.2,0.3); startpar=pars; moves=NULL; control.weights=NULL; lik.fn <- bdSplit.lik; plot.fn <- plotSimmap
bayou.makeMCMC <- function(tree, dat, pred=NULL, SE=0, model="OU", prior, samp=10, chunk=100, 
                           control=NULL, tuning=NULL, new.dir=TRUE, plot.freq=500, outname="bayou", 
                           plot.fn=phenogram, ticker.freq=1000, tuning.int=c(0.1,0.2,0.3), 
                           startpar=NULL, moves=NULL, control.weights=NULL, lik.fn=NULL, 
                           perform.checks=TRUE){
  
  if(is.character(model)){
    model.pars <- switch(model, "OU"=model.OU, "QG"=model.QG, "OUrepar"=model.OUrepar)
  } else {
    model.pars <- model
    model <- "Custom"
  }
  
  if(perform.checks == TRUE){
    if(model=="Custom") {
      checks <- bayou.checkModel(pars=startpar, tree=tree, dat=dat, pred=pred, SE=SE, prior=prior, model=model.pars, autofix=TRUE)
    } else {
      checks <- bayou.checkModel(pars=startpar, tree=tree, dat=dat, pred=pred, SE=SE, prior=prior, model=model, autofix=TRUE)
    }
    cache <- checks$autofixed$cache
    startpar <- checks$autofixed$pars
    tree <- checks$autofixed$cache$phy
    dat <- checks$autofixed$cache$dat
    pred <- checks$autofixed$pred
    model.pars <- checks$autofixed$model
    if(all(sapply(checks[-length(checks)], function(x) x)==TRUE)) cat("seems fine...\n")
  } else {
    cache <- .prepare.ou.univariate(tree, dat, SE=SE, pred=pred)
    dat <- cache$dat
    if(is.null(startpar)){
      startpar <- priorSim(prior, cache$phy, model, nsim=1, plot=FALSE, shiftpars=model.pars$rjpars, exclude.branches=NULL)$pars[[1]]
    } 
  }
  
  fixed <- gsub('^[a-zA-Z]',"",names(attributes(prior)$distributions)[which(attributes(prior)$distributions=="fixed")])
  if("loc" %in% fixed){
    fixed <- c(fixed,"slide")
  }
  
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
      dir.name <- NULL
      dir <- paste(new.dir,"/",sep="")
      dir.create(new.dir)
    }
  } 
  
  .lastpar <- function(files){
    fLs <- sapply(files, function(x) countL(summary(x)$description))
    fL <- min(fLs[fLs > 0])
    if(fL==1){skipL <- 0} else {skipL=fL-1}
    res <- lapply(1:length(files), function(x) try(utils::read.table(summary(files[[x]])$description, skip=skipL), silent=TRUE))
    res <- lapply(1:length(files), function(x) if(class(res[[x]])=="try-error"){numeric(0)}else{res[[x]]})
    pars <- list()
    parLs <- lapply(startpar, length)[outpars]
    npars <- length(res[[4]])
    j=4
    if(length(outpars) > 0){
      for(i in 1:length(outpars)){
        pars[[outpars[i]]] <- unlist(res[[4]][,j:(j+parLs[[i]]-1)],F,F)
        j <- j+1+parLs[[i]]-1
      }
    }
    if(length(rjpars) >0){
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
    return(list(pars=pars, i=i, oll=oll, pr1=pr1))
  }
  .lastparSS <- function(files){
    fLs <- sapply(files, function(x) countL(summary(x)$description))
    fL <- min(fLs[fLs > 0])
    if(fL==1){skipL <- 0} else {skipL=fL-1}
    res <- lapply(1:length(files), function(x) try(utils::read.table(summary(files[[x]])$description, skip=skipL), silent=TRUE))
    res <- lapply(1:length(files), function(x) if(class(res[[x]])=="try-error"){numeric(0)}else{res[[x]]})
    pars <- list()
    parLs <- lapply(startpar, length)[outpars]
    npars <- length(res[[4]])
    j=5
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
    ref1 <- unlist(res[[4]][4], F,F)
    return(list(pars=pars, i=i, oll=oll, pr1=pr1, ref1=ref1))
  }
  
  attributes(ct)$splitmergepars <- rjpars
  
  if(is.null(lik.fn)) lik.fn <- model.pars$lik.fn #.OU.lik
  
  filenames <- list(mapsb=paste(dir, outname,".sb",sep=""), 
                    mapsloc=paste(dir, outname,".loc",sep=""),
                    mapst2=paste(dir, outname,".t2",sep=""),
                    pars.output=paste(dir, outname,".pars",sep=""),
                    rjpars=paste(dir, outname,".rjpars", sep=""))
  
  if(any(sapply(filenames, file.exists))){
    warning("Files with given outname already exist in directory '", dir, "'. Runs will be appended to existing files.")
    #ans <- "X"
    #while(!(ans %in% c("Y", "N"))){
    #ans <- readline(cat("Do you want to continue? Continuing will overwrite existing log files. Y/N:","\n"))
    #}
    #if(ans=="N"){stop("MCMC object not made. Either provide a unique 'outname', specify a different directory 'dir' or manually delete files")} else {
    #  cat("Overwriting log files.")
    #}
    files <- list(mapsb=file(filenames$mapsb,open="a"), 
                  mapsloc=file(filenames$mapsloc,open="a"),
                  mapst2=file(filenames$mapst2,open="a"),
                  pars.output=file(filenames$pars.output,open="a"),
                  rjpars=file(filenames$rjpars, open="a"))  
    startinf <- .lastpar(files)
    oldpar <- startinf$pars
    i <- startinf$i
    oll <- startinf$oll
    pr1 <- startinf$pr1
    store <- list("out"=list(), "rjpars"=list(), "sb"=list(), "loc"=list(), "t2"=list())
  } else {
  files <- list(mapsb=file(filenames$mapsb,open="w"), 
                mapsloc=file(filenames$mapsloc,open="w"),
                mapst2=file(filenames$mapst2,open="w"),
                pars.output=file(filenames$pars.output,open="w"),
                rjpars=file(filenames$rjpars, open="w"))  
  
    oldpar <- startpar
    i <- 1
    oll  <- lik.fn(oldpar, cache, dat, model=model)$loglik
    pr1 <- prior(oldpar,cache)
    store <- list("out"=list(), "rjpars"=list(), "sb"=list(), "loc"=list(), "t2"=list())
    store <- .store.bayou2(i, oldpar, outpars, rjpars, oll, pr1, store, 1, 1, parorder, files)
  }

  header <- 0
  gbg <- lapply(files, close)
  
  accept.type <- NULL
  accept <- NULL
  if(!is.null(plot.freq)){
    environment(plot.fn) <- new.env()
    set.runpars(plot.fn, runpars=list(oldpar=oldpar, cache=cache, i=0))
    tr <- .toSimmap(.pars2map(oldpar, cache),cache)
    tcols <- makeTransparent(rainbow(oldpar$ntheta),alpha=200)
    names(tcols)<- 1:oldpar$ntheta
    plot.fn(tr,dat,colors=tcols,ftype="off")
    plot.dim <- list(par('usr')[1:2],par('usr')[3:4])
  }
  #tuning.int <- round(tuning.int*ngen,0)

  mcmc.loop <- function(ngen){
    startinf <- .lastpar(files)
    oldpar <- startinf$pars
    i <- startinf$i
    oll <- startinf$oll
    pr1 <- startinf$pr1
    iseq <- (i+1):(i+ngen)
    for (i in iseq){
      ct <- .updateControl(ct, oldpar, fixed)
      u <- runif(1)
      prop <- .proposalFn(u,ct,D,moves,cache,oldpar,prior)
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
      if(i %% samp == 0){
        store <- .store.bayou2(i, oldpar, outpars, rjpars, oll, pr1, store, 1, 1, parorder, files)
      }
      if(!is.null(plot.freq)){
        if(i %% plot.freq==0){
          set.runpars(plot.fn, runpars=list(oldpar=oldpar, cache=cache, i=i))
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
        model.pars$monitor.fn(i, oll, pr1, oldpar, accept, accept.type, header)
        header <- 1
        #alpha <- switch(model,"QG"=QG.alpha(oldpar),"OU"=oldpar$alpha,"OUrepar"=OU.repar(oldpar)$alpha,"OUcpp"=oldpar$alpha,"QGcpp"=QG.alpha(oldpar),"OUrepar"=OU.repar(oldpar)$alpha, "bd"=oldpar$r, "ffancova"=oldpar$alpha)
        #sig2 <- switch(model,"QG"=QG.sig2(oldpar),"OU"=oldpar$sig2,"OUrepar"=OU.repar(oldpar)$sig2,"OUcpp"=oldpar$sig2,"QGcpp"=QG.sig2(oldpar),"OUrepar"=OU.repar(oldpar)$sig2, "bd"=oldpar$eps, "ffancova"=oldpar$sig2)
        #if(model=="bd"){
        #  tick <- c(i,oll,pr1,median(alpha),median(sig2),oldpar$k,tapply(accept,accept.type,mean))
        #  tick[-1] <- round(tick[-1],2)
        #  names(tick)[1:6] <- c('gen','lnL','prior','r','eps','K')
        #} else {
        #  tick <- c(i,oll,pr1,log(2)/alpha,sig2/(2*alpha),oldpar$k,tapply(accept,accept.type,mean))
        #  tick[-1] <- round(tick[-1],2)
        #  names(tick)[1:6] <- c('gen','lnL','prior','half.life','Vy','K')
        #}
        #if(i==ticker.freq){
        #  cat(c(names(tick),'\n'),sep='\t\t\t')
        #}
        #cat(c(tick,'\n'),sep='\t\t\t')
      }
    }
  }
  mcmc.prior.only <- function(ngen){
    startinf <- .lastpar(files)
    oldpar <- startinf$pars
    i <- startinf$i
    oll <- 0#startinf$oll
    pr1 <- startinf$pr1
    iseq <- (i+1):(i+ngen)
    for (i in iseq){
      ct <- .updateControl(ct, oldpar, fixed)
      u <- runif(1)
      prop <- .proposalFn(u,ct,D,moves,cache,oldpar, prior)
      new.pars <- prop$pars
      new.cache <- prop$prop$cache
      accept.type <- c(accept.type,paste(prop$decision,prop$move,sep="."))
      pr2 <- prior(new.pars,cache)
      hr <- prop$hr
      new <- lik.fn(new.pars, cache, dat, model=model)
      nll <- new$loglik
      nll <- ifelse(is.na(nll), -Inf, nll)
      if (runif(1) < exp(#nll-oll+
                            pr2-pr1+hr)){
        oldpar <- new.pars
        pr1 <- pr2
        oll <- nll
        accept <- c(accept,1)
      } else {
        accept <- c(accept,0)
      }
      if(i %% samp == 0){
        store <- .store.bayou2(i, oldpar, outpars, rjpars, oll, pr1, store, 1, 1, parorder, files)
      }
      if(!is.null(plot.freq)){
        if(i %% plot.freq==0){
          set.runpars(plot.fn, runpars=list(oldpar=oldpar, cache=cache, i=i))
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
        model.pars$monitor.fn(i, oll, pr1, oldpar, accept, accept.type, header)
        header <- 1
        #alpha <- switch(model,"QG"=QG.alpha(oldpar),"OU"=oldpar$alpha,"OUrepar"=OU.repar(oldpar)$alpha,"OUcpp"=oldpar$alpha,"QGcpp"=QG.alpha(oldpar),"OUrepar"=OU.repar(oldpar)$alpha, "bd"=oldpar$r, "ffancova"=oldpar$alpha)
        #sig2 <- switch(model,"QG"=QG.sig2(oldpar),"OU"=oldpar$sig2,"OUrepar"=OU.repar(oldpar)$sig2,"OUcpp"=oldpar$sig2,"QGcpp"=QG.sig2(oldpar),"OUrepar"=OU.repar(oldpar)$sig2, "bd"=oldpar$eps, "ffancova"=oldpar$sig2)
        #if(model=="bd"){
        #  tick <- c(i,oll,pr1,median(alpha),median(sig2),oldpar$k,tapply(accept,accept.type,mean))
        #  tick[-1] <- round(tick[-1],2)
        #  names(tick)[1:6] <- c('gen','lnL','prior','r','eps','K')
        #} else {
        #  tick <- c(i,oll,pr1,log(2)/alpha,sig2/(2*alpha),oldpar$k,tapply(accept,accept.type,mean))
        #  tick[-1] <- round(tick[-1],2)
        #  names(tick)[1:6] <- c('gen','lnL','prior','half.life','Vy','K')
        #}
        #if(i==ticker.freq){
        #  cat(c(names(tick),'\n'),sep='\t\t\t')
        #}
        #cat(c(tick,'\n'),sep='\t\t\t')
      }
    }
  }
  steppingstone.loop <- function(k, Bk, ngen, ssfilenames, ref){
      ssfiles <- list(mapsb=file(ssfilenames[[k]]$mapsb,open="a"), 
                         mapsloc=file(ssfilenames[[k]]$mapsloc,open="a"),
                         mapst2=file(ssfilenames[[k]]$mapst2,open="a"),
                         pars.output=file(ssfilenames[[k]]$pars.output,open="a"),
                         rjpars=file(ssfilenames[[k]]$rjpars, open="a"))
      startinf <- .lastparSS(ssfiles)
      oldpar <- startinf$pars
      i <- startinf$i
      oll <- startinf$oll
      pr1 <- startinf$pr1
      ref1 <- startinf$ref1
      pB.old <- powerPosteriorFn(k, Bk, oll, pr1, ref1)
      iseq <- (i+1):(i+ngen)
      store <- list("out"=list(), "rjpars"=list(), "sb"=list(), "loc"=list(), "t2"=list())
      for (i in iseq){
        ct <- .updateControl(ct, oldpar, fixed)
        u <- runif(1)
        prop <- .proposalFn(u,ct,D,moves,cache,oldpar, prior)
        new.pars <- prop$pars
        #new.cache <- prop$prop$cache
        accept.type <- c(accept.type,paste(prop$decision,prop$move,sep="."))
        pr2 <- prior(new.pars,cache)
        hr <- prop$hr
        new <- lik.fn(new.pars, cache, dat, model=model)
        nll <- new$loglik
        nll <- ifelse(is.na(nll), -Inf, nll)
        ref2 <- ref(new.pars, cache)
        pB.new <- powerPosteriorFn(k, Bk, nll, pr2, ref2)
        bhr <- ifelse(!is.finite(hr), hr, Bk[k]*hr)
        if (runif(1) < exp(pB.new-pB.old+bhr)){
          oldpar <- new.pars
          pr1 <- pr2
          oll <- nll
          ref1 <- ref2
          pB.old <- powerPosteriorFn(k, Bk, oll, pr1, ref1)
          accept <- c(accept,1)
        } else {
          accept <- c(accept,0)
        }
        if(i %% samp == 0){
          store <- .store.bayou2(i, oldpar, outpars, rjpars, oll, pr1, store, 1, 1, parorder, ssfiles, ref=ref1)
        }
        #if(!is.null(plot.freq)){
        #  if(i %% plot.freq==0){
        #    set.runpars(plot.fn, runpars=list(oldpar=oldpar, cache=cache, i=i))
        #    tr <- .toSimmap(.pars2map(oldpar, cache),cache)
        #    tcols <- makeTransparent(rainbow(oldpar$ntheta),alpha=200)
        #   names(tcols)<- 1:oldpar$ntheta
        #    plot(plot.dim[[1]],plot.dim[[2]],type="n",xlab="time",ylab="phenotype")
        #    mtext(paste("gens = ",i," lnL = ",round(oll,2)),3)
        #    #regime.plot probably doesn't work for simmaps
        #    #try(regime.plot(oldpar,tr$tree,tcols,type="density",model=model),silent=TRUE)
        #    plot.fn(tr,dat,colors=tcols,ftype="off",add=TRUE)
        #  }
        #}
        if(i%%ticker.freq==0){
          model.pars$monitor.fn(i, oll, pr1, oldpar, accept, accept.type, header)
          header <- 1
          #alpha <- switch(model,"QG"=QG.alpha(oldpar),"OU"=oldpar$alpha,"OUrepar"=OU.repar(oldpar)$alpha,"OUcpp"=oldpar$alpha,"QGcpp"=QG.alpha(oldpar),"OUrepar"=OU.repar(oldpar)$alpha, "bd"=oldpar$r, "ffancova"=oldpar$alpha)
          #sig2 <- switch(model,"QG"=QG.sig2(oldpar),"OU"=oldpar$sig2,"OUrepar"=OU.repar(oldpar)$sig2,"OUcpp"=oldpar$sig2,"QGcpp"=QG.sig2(oldpar),"OUrepar"=OU.repar(oldpar)$sig2, "bd"=oldpar$eps, "ffancova"=oldpar$sig2)
          #if(model=="bd"){
          #  tick <- c(i,oll,pr1,median(alpha),median(sig2),oldpar$k,tapply(accept,accept.type,mean))
          #  tick[-1] <- round(tick[-1],2)
          #  names(tick)[1:6] <- c('gen','lnL','prior','r','eps','K')
          #} else {
          #  tick <- c(i,oll,pr1,log(2)/alpha,sig2/(2*alpha),oldpar$k,tapply(accept,accept.type,mean))
          #  tick[-1] <- round(tick[-1],2)
          #  names(tick)[1:6] <- c('gen','lnL','prior','half.life','Vy','K')
          #}
          #if(i==ticker.freq){
          #  cat(c(names(tick),'\n'),sep='\t\t\t')
          #}
          #cat(c(tick,'\n'),sep='\t\t\t')
        }
      }
      out <- list('model'=model, 'dir.name'=dir.name,'dir'=dir, 'outname'=paste(outname, "_ss", k,".", sep=""), 'accept'=accept,'accept.type'=accept.type)
      class(out) <- c("ssbayouFit", "list")
      return(out)
    }
  run.mcmc <- function(ngen, prior.only=FALSE){
    closeAllConnections()
    files <- list(mapsb=file(filenames$mapsb,open="a"), 
                  mapsloc=file(filenames$mapsloc,open="a"),
                  mapst2=file(filenames$mapst2,open="a"),
                  pars.output=file(filenames$pars.output,open="a"),
                  rjpars=file(filenames$rjpars, open="a"))  
    if(prior.only==TRUE){
      tryCatch(mcmc.prior.only(ngen))
    } else {
      tryCatch(mcmc.loop(ngen))
    }
    gbg <- lapply(files, close)
  }
  run.steppingstone <- function(ngen, chain, Bk, burnin=0.3, plot=TRUE){
    closeAllConnections()
    ssfilenames <- lapply(1:length(Bk), function(y) lapply(filenames, function(x) gsub(paste(outname, ".", sep=""), paste(outname, "_ss", y,".", sep=""), x)))
    ssfiles <- list(); 
    postburn <- floor(max(c(1,burnin*length(chain$gen)))):length(chain$gen)
    ref <- make.refFn(chain, model=model.pars, priorFn=prior, burnin = burnin, plot)
    for(x in 1:length(Bk)){
      ssfiles[[x]] <- list(mapsb=file(ssfilenames[[x]]$mapsb,open="w"), 
                        mapsloc=file(ssfilenames[[x]]$mapsloc,open="w"),
                        mapst2=file(ssfilenames[[x]]$mapst2,open="w"),
                        pars.output=file(ssfilenames[[x]]$pars.output,open="w"),
                        rjpars=file(ssfilenames[[x]]$rjpars, open="w"))
      #ppFn <-make.powerposteriorFn(Bk, priorFn=prior, refFn = ref, model=model.pars)
      startpars <- pull.pars(sample(postburn, 1, replace=FALSE), chain, model=model.pars)
      olls <- lik.fn(startpars, cache, dat, model=model)$loglik
      prs <- prior(startpars, cache)
      refs <- ref(startpars, cache)
      stores <- list("out"=list(), "rjpars"=list(), "sb"=list(), "loc"=list(), "t2"=list())
      stores <- .store.bayou2(1, startpars, outpars, rjpars, olls, prs, stores, 1, 1, parorder, ssfiles[[x]], ref=refs)
      gbg <- lapply(ssfiles[[x]], close) #lapply(1:length(Bk), function(x) lapply(ssfiles[[x]], close))
    }
    
    ssfits <- foreach(j=1:length(Bk), .export = c("steppingstone.loop")) %dopar% steppingstone.loop(j, Bk, ngen, ssfilenames, ref)
    outs <- lapply(1:length(Bk), function(x){out$outname <- paste(outname, "_ss", x, sep=""); out})
    chains <- lapply(1:length(Bk), function(x) load.bayou(outs[[x]], saveRDS=FALSE, file=NULL, cleanup=FALSE, ref=TRUE))
    postburn <- floor(max(c(1, burnin*length(chains[[1]]$gen)))):( min(sapply(chains, function(x) length(x$gen))) )
    lnr <- .computelnr(chains, Bk, postburn)   
    ssres <- list(chains=chains, lnr=lnr$lnr, lnrk=lnr$lnrk, Bk=Bk, fits=ssfits, filenames=ssfilenames, startpars=startpars, refFn=ref)
    class(ssres) <- c("ssMCMC", "list")
    return(ssres)
  }
  out <- list('run' = run.mcmc, 'steppingstone'=run.steppingstone, 'model'=model, 'model.pars'=model.pars, 'dir.name'=dir.name,'dir'=dir, 'outname'=outname, 'tree'=tree, 'dat'=dat, 'pred'=pred, 'SE'=SE, 'tmpdir'=ifelse(new.dir==TRUE, TRUE, FALSE), 'startpar'=startpar)
  mcmc.load <- function(saveRDS=FALSE, file=NULL, cleanup=FALSE){
    load.bayou(out, saveRDS, file, cleanup)
  }
  out$load <- mcmc.load
  class(out) <- c("bayouMCMCFn", "list")
  return(out)
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
