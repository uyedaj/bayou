make.refFn <- function(chain,.prior,model="OU",burnin=0.3){
  parorder <- switch(model,"QG"=c("h2","P","w2","Ne","nb","ntheta","optima"), "OU"=c("alpha","sig2","nb","ntheta","optima"),"OUrepar"=c("halflife","Vy","nb","ntheta","optima"))
  postburn <- round(0.3*length(chain[[1]]),0):length(chain[[1]])
  fx <- lapply(attributes(.prior)$distributions,function(x) try(get(x),silent=TRUE))
  fx$bprior <- bFUN
  dists <- attributes(.prior)$distributions
  switch.dist <- function(distribution) {try(switch(distribution,"dbeta"="beta", "dcauchy"="cauchy", "dchisq"="chi-squared","dexp"="exponential", "dgamma"="gamma","dlnorm"="lognormal","dnorm"="normal","dpois"="Poisson", "cdpois"="Poisson"),silent=TRUE)}
  dist.names <- lapply(dists,switch.dist)
  names(dist.names) <- sapply(names(dist.names),function(x) try(switch(x,"dalpha"="alpha","dsig2"="sig2","dk"="nb","dnb"="nb","dtheta"="optima","dh2"="h2","dNe"="Ne","dw2"="w2","dP"="P","dVy"="Vy","dhalflife"="halflife"), silent=TRUE))
  tofit <- which(sapply(dist.names,function(x) !is.null(x))==TRUE)
  start <- attributes(.prior)$parameters[tofit]
  start <- lapply(tofit,function(y) if(dist.names[y] %in% c("lognormal","Poisson","normal","exponential","geometric")){NULL} else {start[[y]]})
  fits <- lapply(tofit, function(y) fitdistr(unlist(chain[[names(tofit)[y]]][postburn]),densfun=dist.names[[y]],start=start[[y]])$estimate)
  fits$nb <- c(fits$nb,"kmax"=attributes(.prior)$parameters$nb$kmax)
  setfx <- lapply(tofit,function(y) set.defaults(fx[[y]],defaults=as.list(fits[[y]])))
  setfx <- lapply(tofit,function(y) set.defaults(setfx[[y]],log=TRUE))
  return(function(pars,emap,cache) { dens <- sapply(names(setfx),function(y) sum(setfx[[y]](pars[[y]])));
                                     bdens <- log(1/choose(2*cache$ntips-2,pars$nb));
                                     ldens <- sum(dunif(emap$sh[emap$sh==1],log=TRUE));
                                     return(sum(dens)+bdens+ldens)
  })
}

#make.refFn(chain1,.prior,model="OU",burnin=0.3)->refFn1
#refFn1(pars,emap,cache)
#.prior(pars,emap,cache)+
#.emOU.lik(pars,emap,cache,dat,SE=0,model="OU")$loglik
make.powerPosterior <- function(k,Bk,prior.fn,ref.fn){
  foo <- function(k,Bk,pars,emap,cache,dat,SE,model){
    lik <- .emOU.lik(pars,emap,cache,cache$dat,SE=SE,model=model)$loglik
    prior <- prior.fn(pars,emap,cache)
    ref <- ref.fn(pars,emap,cache)
    coeff <- c(Bk[k],Bk[k],(1-Bk[k]))
    result <- c(lik,prior,ref)
    result[coeff==0] <- 0
    result <- result*coeff
    result <- sum(result)
    return(list(result=result,lik=lik,prior=prior,ref=ref))
  }
  return(foo)
}
#powerPosterior <- make.powerPosterior(1,Bk,.prior,refFn1)
#powerPosterior(101,Bk,pars,emap,cache,dat,SE=0,model="OU")




#control=NULL
#tuning=NULL
#new.dir=FALSE
#plot=TRUE
#plot.freq=500
#outname="bayou.ss"
#ticker.freq=1000
#tuning.int=c(0.1,0.2,0.3)
#startpar=startpar
#moves=NULL
#control.weights=NULL
#model="OU"
#ngen=10^5
#samp=1
#chunk=1
#SE=NA


steppingstoneMCMC <- function(tree,dat,powerPosterior,k,Bk,SE=NA,model="OU",ngen=100000,samp=20,chunk=200,control=NULL,tuning=NULL,new.dir=FALSE,plot=TRUE,plot.freq=500,outname="bayou.ss",ticker.freq=1000,tuning.int=c(0.1,0.2,0.3),startpar=NULL,moves=NULL,control.weights=NULL){
  if(is.null(moves)){
    moves <- switch(model,"QG"=list(h2=".multiplierProposal",P=".multiplierProposal",w2=".multiplierProposal",Ne=".multiplierProposal",nb="splitmerge.emap",optima=".adjustTheta",slide=".slideBranch",pos=".adjustPos"),"OU"=list(alpha=".multiplierProposal",sig2=".multiplierProposal",nb="splitmerge.emap",optima=".adjustTheta",slide=".slideBranch",pos=".adjustPos"),"OUrepar"=list(halflife=".multiplierProposal",Vy=".multiplierProposal",nb="splitmerge.emap",optima=".adjustTheta",slide=".slideBranch",pos=".adjustPos"))
  }
  cache <- .prepare.ou.univariate(tree,dat)
  #FOR NOW PROVIDE STARTPAR
  #  if(is.null(startpar)){
  #    startpar <- .priorsim(.prior,tree,plot=FALSE)
  #  }
  if(is.null(control.weights)){
    ct <- build.control(startpar$pars,startpar$emap,default.weights=model)
  } else {ct <- build.control(startpar$pars,startpar$emap,default.weights=NULL,move.weights=control.weights)}
  
  if(is.null(tuning)){
    D <- switch(model, "OU"=list(alpha=1, sig2= 1, nb = 4,optima=2,slide=1,pos=mean(cache$edge.length)), "QG"=list(h2=1,P=1,w2=1,Ne=1, nb = 4,optima=2,slide=1,pos=mean(cache$edge.length)),"OUrepar"=list(halflife=1,Vy=1,nb=4,optima=2,slide=1,pos=mean(cache$edge.length)))
  } else { D <- tuning}
  
  if(attributes(.prior)$distributions$bprior=="equal"){
    plim <- 0
  } else {plim <- cache$ntips}
  
  if(new.dir){
    dir.name <- paste(.sample(LETTERS,10),collapse="")
    dir <- paste(getwd(),"/",dir.name,"/",sep="")
    dir.create(dir)
    setwd(dir)
  } else {
    dir <- getwd()
    dir.name <- sapply(strsplit(dir,'/'),function(x) x[length(x)])
  }
  
  mapsb <<- file(paste(outname,".",k,".mapsb",sep=""),open="w")
  mapsr2 <<- file(paste(outname,".",k,".mapsr2",sep=""),open="w")
  mapst2 <<- file(paste(outname,".",k,".mapst2",sep=""),open="w")
  pars.output <<- file(paste(outname,".",k,".pars",sep=""),open="w")
  
  oldpar <- startpar$pars
  oldmap <- startpar$emap
  chunk.branch <<- list()
  chunk.r2 <<- list()
  chunk.t2 <<- list()
  out <<- list()
  
  #oll  <- .emOU.lik(oldpar,oldmap,cache,dat,model=model)$loglik
  #pr1 <- .prior(oldpar,oldmap,cache)
  pB.old <- powerPosterior(k,Bk,oldpar,oldmap,cache,dat,SE=0,model=model)
  parorder <- switch(model,"QG"=c("h2","P","w2","Ne","nb","ntheta","optima"), "OU"=c("alpha","sig2","nb","ntheta","optima"),"OUrepar"=c("halflife","Vy","nb","ntheta","optima"))
  
  accept.type <- NULL
  accept <- NULL
  
  tr <- emap2simmap(oldmap,cache$phy)
  tcols <- makeTransparent(tr$col,alpha=100)
  if(plot){
    phenogram(tr,dat,colors=tr$col,ftype="off")
    plot.dim <- list(par('usr')[1:2],par('usr')[3:4])
  }
  #tuning.int <- round(tuning.int*ngen,0)
  Ref <- NULL
  for (i in 1:ngen){
    u <- runif(1)
    prop <- .proposalFn(u,ct,D,moves,cache,oldpar,oldmap)
    new.pars <- prop$prop$pars
    new.emap <- prop$prop$emap
    accept.type <- c(accept.type,paste(prop$prop$decision,prop$move,sep="."))
    #pr2 <- .prior(new.pars,new.emap,cache)
    hr <- prop$prop$hr
    #new <- .emOU.lik(new.pars,new.emap,cache,dat,model=model)
    #nll <- new$loglik
    pB.new <- powerPosterior(k,Bk,new.pars,new.emap,cache,dat,SE=0,model=model)
    if (runif(1) < exp(pB.new$result-pB.old$result+Bk[k]*hr)){
      oldmap <- new.emap
      oldpar <- new.pars
      pr1 <- pB.new$prior
      oll <- pB.new$lik
      pB.old <- pB.new
      accept <- c(accept,1)
      if(i %% samp ==0){
        Ref <- c(Ref,pB.new$ref)
      }
    } else {
      accept <- c(accept,0)
      if(i %% samp==0){
        Ref <- c(Ref,pB.old$ref)
      }
    }
    store.QG(i,oldpar,oldmap,cache,dat,oll,pr1,samp,chunk,parorder)
    if(ct$nb+ct$slide > 0){
      ct <- build.control(oldpar,oldmap,default.weights=model)
    }
    if(plot){
      if(i %% plot.freq==0){
        tr <- emap2simmap(oldmap,cache$phy)
        tcols <- makeTransparent(tr$col,alpha=100)
        plot(plot.dim[[1]],plot.dim[[2]],type="n",xlab="time",ylab="phenotype")
        mtext(paste("gens = ",i," lnL = ",round(oll,2)),3)
        try(regime.plot(oldpar,tr,tcols,type="density",model=model),silent=TRUE)
        phenogram(tr,dat,colors=tr$col,ftype="off",add=TRUE)
      }
    }
    #if(i %in% tuning.int){
    #  D <- tune.D(D,accept,accept.type)$D
    #}
    if(i%%ticker.freq==0){
      alpha <- switch(model,"QG"=QG.alpha(oldpar),"OU"=oldpar$alpha,"OUrepar"=OU.repar(oldpar)$alpha)
      sig2 <- switch(model,"QG"=QG.sig2(oldpar),"OU"=oldpar$sig2,"OUrepar"=OU.repar(oldpar)$sig2)
      tick <- c(i,oll,pr1,log(2)/alpha,sig2/(2*alpha),oldpar$nb,tapply(accept,accept.type,mean))
      tick[-1] <- round(tick[-1],2)
      names(tick)[1:6] <- c('gen','lnL','prior','half.life','Vy','K')
      if(i==ticker.freq){
        cat(c(names(tick),'\n'),sep='\t')
      }
      cat(c(tick,'\n'),sep='\t')
    }
  }
  closeAllConnections()
  return(list('dir.name'=dir.name,'dir'=dir,'accept'=accept,'accept.type'=accept.type,"ref"=Ref))
}

pull.rsample <- function(samp,chain,fit,refFn,model="OU"){
  #pars.list <- lapply(samp,function(y) pull.pars(y,chain,model=model))
  #emap.list <- lapply(samp,function(y) read.emap(chain$branch.shift[[y]],chain$location[[y]],chain$t2[[y]],cache$phy)$emap)
  L <- chain$lnL[samp]+chain$prior[samp]-fit$ref[samp]
  Lmax <- max(L)
  Lfactored <- L-Lmax
  return(list(Lmax=Lmax,Lfactored=Lfactored))
}

computelnr <- function(K.chains,ssfits,Bk,samp,refFn,model="OU"){
  lnr <- list()
  for(i in 1:(length(Bk)-1)){
    Lk <- pull.rsample(samp,K.chains[[i]],ssfits[[i]],refFn,model=model)
    lnr[[i]] <- (Bk[i+1]-Bk[i])*Lk$Lmax+log(1/length(Lk$Lfactored)*sum(exp(Lk$Lfactored)^(Bk[i+1]-Bk[i])))
  }
  return(list("lnr"=sum(unlist(lnr)),"lnrk"=lnr))
} 
