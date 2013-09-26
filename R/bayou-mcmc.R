bayou.mcmc <- function(tree,dat,SE=0,model="OU",prior,ngen=10000,samp=10,chunk=100,control=NULL,tuning=NULL,new.dir=FALSE,plot.freq=500,outname="bayou",ticker.freq=1000,tuning.int=c(0.1,0.2,0.3),startpar=NULL,moves=NULL,control.weights=NULL){
  if(is.null(moves)){
    moves <- switch(model,"QG"=list(h2=".multiplierProposal",P=".multiplierProposal",w2=".multiplierProposal",Ne=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideBranch",pos=".adjustPos"),
                          "OU"=list(alpha=".multiplierProposal",sig2=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideBranch",pos=".adjustPos"),
                          "OUrepar"=list(halflife=".multiplierProposal",Vy=".multiplierProposal",nb=".splitmergeSimmap",theta=".adjustTheta",slide=".slideBranch",pos=".adjustPos")) #,"OUcpp"=list(alpha=".multiplierProposal",sig2=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"), "QGcpp"=list(h2=".multiplierProposal",P=".multiplierProposal",w2=".multiplierProposal",Ne=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"),"OUreparcpp"=list(halflife=".multiplierProposal",Vy=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"))
  }
  cache <- .prepare.ou.univariate(tree,dat)
  T <- sum(cache$edge.length)
  if(is.null(startpar)){
    startpar <- priorSim(prior,tree,model,nsim=1,plot=FALSE, exclude.branches=NULL)
  }
  cache$maps <- pars2simmap(startpar$pars[[1]],tree,sim.theta=FALSE,theta=startpar$pars[[1]])$tree$maps
  if(is.null(control.weights)){
    ct <- .buildControl(startpar$pars[[1]],default.weights=model,move.weights=NULL,prior)
  } else {ct <- .buildControl(startpar$pars[[1]],default.weights=NULL,move.weights=control.weights,prior)}
  
  if(is.null(tuning)){
    D <- switch(model, "OU"=list(alpha=1, sig2= 1, k = 4,theta=2,slide=1,pos=mean(cache$edge.length)), "QG"=list(h2=1,P=1,w2=1,Ne=1, k = 4,theta=2,slide=1,pos=mean(cache$edge.length)),"OUrepar"=list(halflife=1,Vy=1,k=4,theta=2,slide=1,pos=mean(cache$edge.length)),"OUcpp"=list(alpha=1, sig2= 1,sig2jump=2, k = 4,theta=2,slide=1),"QGcpp"=list(h2=1,P=1,w2=1,Ne=1,sig2jump=2,k=4,theta=2,slide=1),"OUreparcpp"=list(halflife=1,Vy=1,sig2jump=2,k=4,theta=2,slide=1))
  } else { D <- tuning}
  
  if(new.dir){
    dir.name <- paste(sample(LETTERS,10),collapse="")
    dir <- paste(getwd(),"/",dir.name,"/",sep="")
    dir.create(dir)
    setwd(dir)
  } else {
    dir <- getwd()
    dir.name <- sapply(strsplit(dir,'/'),function(x) x[length(x)])
  }
  
  mapsb <- file(paste(outname,".mapsb",sep=""),open="w")
  mapsloc <- file(paste(outname,".maploc",sep=""),open="w")
  mapst2 <- file(paste(outname,".mapst2",sep=""),open="w")
  pars.output <- file(paste(outname,".pars",sep=""),open="w")
  
  oldpar <- startpar$pars[[1]]
  chunk.sb <- list()
  chunk.loc <- list()
  chunk.t2 <- list()
  out <- list()
  
  lik.fn <- switch(model,"OU"=.smOU.lik,"QG"=.smOU.lik,"OUrepar"=.smOU.lik)
  oll  <- lik.fn(oldpar,cache,dat,SE=SE,model=model)$loglik
  pr1 <- prior(oldpar,cache)
  parorder <- switch(model,"QG"=c("h2","P","w2","Ne","k","ntheta","theta"), "OU"=c("alpha","sig2","k","ntheta","theta"),"OUrepar"=c("halflife","Vy","k","ntheta","theta"),"OUcpp"=c("alpha","sig2","sig2jump","k","ntheta","theta"),"QGcpp"=c("h2","P","w2","Ne","sig2jump","k","ntheta","theta"),"OUreparcpp"=c("halflife","Vy","sig2jump","k","ntheta","theta"))
  
  accept.type <- NULL
  accept <- NULL
  if(!is.null(plot.freq)){
    tr <- pars2simmap(oldpar,cache$phy,sim.theta=FALSE,theta=oldpar$theta)
    tcols <- makeTransparent(tr$col,alpha=100)
    phenogram(tr$tree,dat,colors=tr$col,ftype="off")
    plot.dim <- list(par('usr')[1:2],par('usr')[3:4])
  }
  #tuning.int <- round(tuning.int*ngen,0)
  for (i in 1:ngen){
    u <- runif(1)
    prop <- .proposalFn(u,ct,D,moves,cache,oldpar)
    prop$move
    new.pars <- prop$prop$pars
    new.cache <- prop$prop$cache
    accept.type <- c(accept.type,paste(prop$prop$decision,prop$move,sep="."))
    pr2 <- prior(new.pars,cache)
    hr <- prop$prop$hr
    new <- lik.fn(new.cache, new.pars,dat,SE=SE,model=model)
    nll <- new$loglik
    if (runif(1) < exp(nll-oll+pr2-pr1+hr)){
      cache$maps <- new.cache$maps
      oldpar <- new.pars
      pr1 <- pr2
      oll <- nll
      accept <- c(accept,1)
    } else {
      accept <- c(accept,0)
    }
    store.bayOU(i,oldpar,oll,pr1,samp,chunk,parorder)
    if(!is.null(plot.freq)){
      if(i %% plot.freq==0){
        tr <- pars2simmap(oldpar,cache$phy,sim.theta=FALSE,theta=oldpar$theta)
        tcols <- makeTransparent(tr$col,alpha=100)
        plot(plot.dim[[1]],plot.dim[[2]],type="n",xlab="time",ylab="phenotype")
        mtext(paste("gens = ",i," lnL = ",round(oll,2)),3)
        #regime.plot probably doesn't work for simmaps
        #try(regime.plot(oldpar,tr$tree,tcols,type="density",model=model),silent=TRUE)
        phenogram(tr$tree,dat,colors=tr$col,ftype="off",add=TRUE)
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
        cat(c(names(tick),'\n'),sep='\t')
      }
      cat(c(tick,'\n'),sep='\t')
    }
  }
  closeAllConnections()
  return(list('dir.name'=dir.name,'dir'=dir,'accept'=accept,'accept.type'=accept.type))
}
