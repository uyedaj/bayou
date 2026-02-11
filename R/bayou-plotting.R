# Utility for getting the starting and ending ages for each regime
.optima.ages <- function(pars,tree){
  nH <- nodeHeights(tree)
  reg <- sapply(tree$maps,function(x) names(x)[length(x)])
  adj <- sapply(tree$maps,function(x) ifelse(length(x)>1,x[1],0))
  abs.age <- nH[,1]+adj
  start <- tapply(abs.age,reg,min)
  end <- rep(max(nH),pars$ntheta)+1
  o <- as.character(1:pars$ntheta)
  names(end) <- names(start)
  return(cbind(start[o],end[o]))
}

#' Make a color transparent (Taken from an answer on StackOverflow by Nick Sabbe)
#'
#' @param someColor A color, either a number, string or hexidecimal code
#' @param alpha The alpha transparency. The maxColorValue is set to 255.
#' @return A character vector of colors in hexadecimal format with the specified transparency applied.

#' @export
makeTransparent <- function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#' Plot a phylogenetic tree with posterior probabilities from a bayouMCMC chain (function adapted from phytools' plotSimmap)
#'
#' @param chain A bayouMCMC chain
#' @param burnin The proportion of runs to be discarded, if NULL, then the value stored in the bayouMCMC chain's attributes is used
#' @param lwd The width of the edges
#' @param edge.type Either "theta" (branches will be colored according to their median value of theta), "regimes" (clades will be assigned to distinct regimes if the posterior probability of a shift
#' on that branch is > pp.cutoff), or "pp" (branches will be colored according to the probability of a shift on that branch). If "none" then edge.color will be assigned to all branches.
#' @param pal A color palette function used to paint the branches (unless edge.type="none")
#' @param pp.cutoff If edge.type=="regimes", the posterior probability above which a shift should be reconstructed on the tree.
#' @param circles a logical value indicating whether or not a circle should be plotted at the base of the node with values that correspond to the posterior probability of having a shift.
#' @param circle.cex.max The cex value of a circle with a posterior probability of 1
#' @param circle.col The color used to fill the circles
#' @param circle.pch the type of symbol used to plot at the node to indicate posterior probability
#' @param circle.lwd the line width of the points plotted at the nodes
#' @param circle.alpha a value between 0 and 255 that indicates the transparency of the circles (255 is completely opaque).
#' @param pp.labels a logical indicating whether the posterior probability for each branch should be printed above the branch
#' @param pp.col The color used for the posterior probability labels
#' @param pp.alpha a logical or numeric value indicating transparency of posterior probability labels. If TRUE, then transparency is ramped from invisible (pp=0), to black (pp=1). If numeric, all labels are given the same transparency. If NULL, then no transparency is given.
#' @param pp.cex the size of the posterior probability labels
#' @param edge.color The color of edges if edge.type="none"
#' @param parameter.sample When edge.type=="theta", the number of samples used to estimate the median "theta" value from each branch. Since this is
#' computationally intensive, this enables you to downsample the chain.
#' @param ... Additional arguments passed to ape's plot.phylo
#'
#' @return **No return value**, called for **side effects**.
#' The function generates a **phylogenetic tree visualization** with branches colored
#' based on posterior probabilities, regimes, or estimated parameters.

#'
#' @export

plotSimmap.mcmc <- function(chain, burnin=NULL, lwd=1, edge.type = c("regimes", "theta", "none", "pp"),
                            pal=rainbow, pp.cutoff=0.3, circles=TRUE, circle.cex.max=3, circle.col="red",
                            circle.pch=21, circle.lwd=0.75, circle.alpha=100, pp.labels=FALSE, pp.col=1,
                            pp.alpha=255, pp.cex=0.75, edge.color = 1, parameter.sample=1000, ...){
  tree <- attributes(chain)$tree
  edge.type <- match.arg(edge.type, c("regimes", "theta", "none", "pp"))
  cache <- .prepare.ou.univariate(tree, attributes(chain)$dat)
  tree <- cache$phy
  if(is.null(burnin)) burnin = attributes(chain)$burnin
  if(is.null(burnin)) burnin = 0
  if(burnin==0) postburn <- 1:length(chain$gen)  else {
    postburn <- round(burnin*length(chain$gen),0):length(chain$gen)
  }
  L <- Lposterior(chain, tree, burnin=burnin)
  if(!is.null(pp.cutoff)){
    pp <- L$pp
    pars <- list()
    pars$sb <- which(pp > pp.cutoff)
    pars$k <- length(pars$sb)
    pars$ntheta <- length(pars$sb)+1
    pars$loc <- L$rel.location[pars$sb]*tree$edge.length[pars$sb]
    pars$t2 <- 2:(length(pars$sb)+1)
    if(length(pars$sb)>0){
      tr <- pars2simmap(pars, tree)$tree
      colors <- NULL
    } else {
      tr <- tree
      colors <- NULL
    }
  } else {
    tr <- tree
    tr$maps <- lapply(tr$edge.length, function(x) setNames(x, 1))
    colors <- setNames(1, 1)
  }
  .colorRamp <- function(trait, .pal, nn){
    strait <- (trait-min(trait))/(max(trait-min(trait)))
    itrait <- floor(strait*(nn-1))+1
    if(!is.null(.pal)){
    return(.pal(nn+1)[itrait])
    } else {
      return(itrait)
    }
  }
  if(edge.type=="none"){
    ape::plot.phylo(tr, edge.color=edge.color, lwd=lwd, ...)
  }
  if(edge.type == "regimes"){
    plotRegimes(tr, col=colors, lwd=lwd, pal=pal, ...)
  }
  if(edge.type == "theta"){
    plotBranchHeatMap(tree, chain, "theta", burnin=burnin, pal=heat.colors, ...)
  }
  if(edge.type == "pp"){
   ape::plot.phylo(tree, edge.color=.colorRamp(L$pp, pal, 100), ...)
  }
  if(circles){
    #theta2 <- L$magnitude.of.theta2
    #root.median <- median(sapply(chain$theta[postburn], function(x) x[1]))
    #theta2[is.na(theta2)] <- root.median
    #theta2 <- theta2 - root.median
    #circle.cols <- sapply(colorRamp(theta2, circle.pal, 100), function(x) makeTransparent(x, circle.alpha))
    circle.cexs <- seq(0, circle.cex.max, length.out=100)[.colorRamp(L$pp, NULL, 100)]
    edgelabels(pch=circle.pch, lwd=circle.lwd, bg=makeTransparent(circle.col, circle.alpha), cex=circle.cexs)
  }
  if(pp.labels){
    edgelabels(round(L$pp,2), col=makeTransparent(pp.col, pp.alpha), cex=pp.cex, frame = "none")
  }
}

#' Adds visualization of regimes to a plot
#'
#' @param pars A bayou formatted parameter list
#' @param tree A tree of class 'phylo'
#' @param cols A vector of colors to give to regimes, in the same order as pars$sb
#' @param type Either "rect", "density" or "lines". "rect" plots a rectangle for the 95\% CI for the stationary
#' distribution of a regime. "density" varies the transparency of the rectangles according to the probability density
#' from the stationary distribution. "lines" plots lines for the mean and 95\% CI's without filling them.
#' @param transparency The alpha transparency value for the maximum density, max value is 255.
regime.plot <- function(pars,tree,cols,type='rect',transparency=100){
  OA <- .optima.ages(pars,tree)
  CIU95 <- pars$theta+2*sqrt(pars$sig2/(2*pars$alpha))
  CIL95 <- pars$theta-2*sqrt(pars$sig2/(2*pars$alpha))
  if(type=="lines"){
    for(i in 1:pars$ntheta){
      lines(c(OA[i,1],OA[i,2]),rep(pars$optima[i],2),col=makeTransparent(cols[i],transparency),lwd=3)
      lines(c(OA[i,1],OA[i,2]),rep(CIU95[i],2),col=makeTransparent(cols[i],transparency),lwd=1.25,lty=2)
      lines(c(OA[i,1],OA[i,2]),rep(CIL95[i],2),col=makeTransparent(cols[i],transparency),lwd=1.25,lty=2)
    }
  }
  if(type=="rect"){
    for(i in 1:pars$ntheta){
      rect(OA[i,1],CIL95[i],OA[i,2],CIU95[i],col=makeTransparent(cols[i],transparency),border=NA)
    }
  }
  if(type=="density"){
    ylim <- par('usr')[3:4]
    for(i in 1:pars$ntheta){
      x <- seq(OA[i,1],OA[i,2],length=10)
      y <- seq(ylim[1],ylim[2],length=100)
      Z <- matrix(nrow=length(x),ncol=length(y))
      for(j in 1:length(x)){
        Z[j,] <- dnorm(y,pars$theta[i],sqrt(pars$sig2/(2*pars$alpha)))
      }
      if(sum(Z)!=0){
        densregion(x,y,Z,colmax=makeTransparent(cols[i],transparency),colmin="transparent")
      }
      lines(c(OA[i,1],OA[i,2]),rep(pars$theta[i],2),col=makeTransparent(cols[i],min(255, 50+(transparency))),lwd=2)
    }
  }
}

#' Plot a pheongram with the posterior density for optima values
#'
#' Plots a phenogram and the posterior density for optima values
#'
#' @param tree A phylogeny of class 'phylo'
#' @param dat A named vector of tip data
#' @param burnin The initial proportion of the MCMC to be discarded
#' @param chain A bayouMCMC object that contains the results of an MCMC chain
#' @param colors An optional named vector of colors to assign to regimes, \code{NULL} results in no regimes being plotted.
#' @param pp.cutoff The posterior probability cutoff value. Branches with posterior probabilities of having a shift above this value
#' will have the average location of the regime shift painted onto the branches.
#' @param K A list with the values of K to be plotted. If \code{NULL} all values of K are combined and a total posterior produced. This
#' allows separate lines to be plotted for different numbers of shifts so that the location of optima can be compared, for example, between
#' all samples that have 1 vs. 2 shifts in the posterior.
#' @param ... Additional parameters passed to \code{phenogram(...)}
#'
#' @return No return value, called for side effects. This function generates a phenogram plot
#' with posterior density overlays for optima values, visualizing the distribution of evolutionary
#' regimes across a phylogenetic tree.

#' @export
phenogram.density <- function(tree, dat, burnin=0, chain ,colors=NULL, pp.cutoff=NULL, K=NULL, ...){
  tree <- reorder(tree,"postorder")
  dat <- dat[tree$tip.label]
  postburn <- round(length(chain$gen)*burnin,0):length(chain$gen)
  chain2 <- lapply(chain,function(x) x[postburn])
  theta <- chain2$theta
  no.theta <- lapply(theta,length)
  min.theta <- min(unlist(theta))
  max.theta <- max(unlist(theta))
  if(is.null(K)){
    K <- list(unique(unlist(no.theta)))
  }
  if(!is.null(pp.cutoff)){
    L <- Lposterior(chain2, tree)
    pp <- L$pp
    pars <- list()
    pars$sb <- which(pp > pp.cutoff)
    pars$k <- length(pars$sb)
    pars$ntheta <- length(pars$sb)+1
    pars$loc <- L$rel.location[pars$sb]*tree$edge.length[pars$sb]
    pars$t2 <- 2:(length(pars$sb)+1)
    if(length(pars$sb)>0){
      tr <- pars2simmap(pars, tree)
      tree <- tr$tree
      colors <- tr$col
      names(colors) <- 1:length(colors)
    } else {
      tr <- tree
      colors <- 1; names(colors) <-1
    }
  }
  if(is.null(colors)){
      ntheta <- length(unique(names(unlist(tree$maps))))
      colors <- rainbow(ntheta)
      names(colors) <- 1:ntheta
    }
  nH <- max(nodeHeights(tree))
  plot(c(0,nH+0.3*nH),c(min(dat)-0.25,max(dat)+0.25),type='n',xlab="Time",ylab="Phenotype")
  phenogram(tree, dat, add=TRUE, colors=colors, spread.labels=FALSE, ...)
  dens.theta <- lapply(1:length(K), function(x) density(unlist(theta[no.theta %in% K[[x]]])))
  tmp <- sapply(1:length(dens.theta),function(Q){lines(nH+dens.theta[[Q]]$y*(0.3*nH)/max(dens.theta[[Q]]$y),dens.theta[[Q]]$x,col=Q+1)})
}

#' S3 method for plotting bayouMCMC objects
#'
#' @param x A mcmc chain of class 'bayouMCMC' produced by the function bayou.mcmc and loaded into the environment using load.bayou
#' @param ... Additional arguments passed to \code{plot.mcmc} from the \code{coda} package
#' @return No return value, called for side effects. This function generates diagnostic
#' trace and density plots for MCMC chains of class `bayouMCMC` to assess convergence
#' and parameter distributions.

#' @export
#' @method plot bayouMCMC
plot.bayouMCMC <- function(x, ...){
  if(is.null(attributes(x)$burnin)){
    start <- 1
  } else {
    start <- round(attributes(x)$burnin*length(x$gen),0)
  }
  postburn <- start:length(x$gen)
  chain2 <- lapply(x,function(x) x[postburn])
  chain.length <- length(chain2$gen)
  univariates <- chain2[sapply(chain2,function(x) length(unlist(x)))==length(chain2$gen)]
  univariates$root <- sapply(chain2$theta, function(x) x[1])
  uni.df <- as.data.frame(univariates)
  uni.df <- uni.df[!duplicated(uni.df$gen),]
  rownames(uni.df) <- uni.df[,1]
  uni.df <- uni.df[,-1]
  plot(mcmc(uni.df), ...)
}


#' Plot parameter list as a simmap tree
#'
#' @param pars A bayou formatted parameter list
#' @param tree A tree of class 'phylo'
#' @param ... Additional arguments passed to plotRegimes
#' @return No return value, called for side effects. This function generates
#' a visualization of a phylogenetic tree with mapped regimes or other parameters.
#' @export
plotBayoupars <- function(pars, tree,...){
  tree <- reorder(tree, 'postorder')
  X <- rep(0, length(tree$tip.label))
  names(X) <- tree$tip.label
  cache <- .prepare.ou.univariate(tree, X)
  tr <- .toSimmap(.pars2map(pars, cache),cache)
  plotRegimes(tr,...)
}

# Experimental function for ancestral state reconstruction for a given OU model
.OU.asr <- function(tree, dat, pars, start=NULL, SE=0){
  phy <- reorder(tree, "postorder")
  dat <- dat[phy$tip.label]
  if(length(SE)>1){
    SE[phy$tip.label]
  }
  if(length(phy$tip.label) > 100) warning("This may take a while for large trees")
  EV <- .vcv.asrOU(phy, dat, pars, SE=SE)
  ntips <- length(phy$tip.label)
  ExpV <- EV$ExpV
  VCV <- EV$VCV
  diag(VCV) <- diag(VCV)+SE^2
  lik.fn <- function(anc){
    -1*dmnorm(as.vector(c(dat, anc)), mean=ExpV[,1], varcov=VCV, log=TRUE)
  }
  if(is.null(start)){
    start = ExpV[(length(dat)+1):(length(phy$edge.length)+1)]
  }
  result <- stats::optim(start, lik.fn, method="L-BFGS-B")
  x <- c(dat, result$par)
  names(x)[(ntips+1):length(x)] <- (ntips+1):length(x)
  return(x)
}


.vcv.asrOU <- function(phy, dat, pars, SE, internal=TRUE){
  cache <- .prepare.ou.univariate(phy, dat, SE=SE)
  phy <- cache$phy
  new.pars <- pars
  ntips <- length(phy$tip.label)
  sig2 <- new.pars$sig2
  alpha <- new.pars$alpha
  D <- dist.nodes(phy)
  Cii <- D[ntips+1,]
  C <- D; C[,] <- 0
  ##Covariance[y_i, y_j]= s2/(2*alpha) * Exp[-alpha*t_ij] *  [1 - Exp(-2*alpha*t_a1)]
  for(i in 1:nrow(D)) for(j in 1:ncol(D))
    C[i,j]<- sig2/(2*alpha)*exp(-alpha*D[i,j])*(1-exp(-2*alpha*(Cii[i]+Cii[j]-D[i,j])))
  ##Calculate expectations
  diag(C) <- diag(C) + 1e-10
  mu <- rep(0, nrow(D))
  mu[1:ntips] <- dat
  W <- .allnodes.W(cache, new.pars)
  ExpV <- W %*% new.pars$theta
  return(list(ExpV=ExpV, VCV=C))
}

.allnodes.W <- function(tree, pars){
  a <- pars$alpha
  s2 <- pars$sig2
  nbranch <- length(tree$edge.length)
  if(inherits(tree,"phylo")){
    X <- rep(NA,length(tree$tip.label))
    names(X) <- tree$tip.label
    cache <- .prepare.ou.univariate(tree,X)
  } else {cache <- tree}
  if(is.null(pars$ntheta)){
    pars$ntheta <- length(pars$theta)
  }
  plook <- function(x){mapply(paste,x[2:length(x)],x[1:(length(x)-1)],sep=",")}
  tB <- cache$desc$anc[1:(cache$n.node+cache$ntips)]
  tB <- mapply(c,1:(cache$n.node+cache$ntips),tB, SIMPLIFY=FALSE)
  lookup <- lapply(tB,plook)
  edge.names <- mapply(paste,cache$edge[,1],cache$edge[,2],sep=",")
  cache$branchtrace <- t(sapply(lookup,function(x) as.numeric(edge.names %in% x)))
  smtree <- pars2simmap(pars, cache$phy)
  maps <- smtree$tree$maps
  allnodes <- cache$n.node+cache$ntips
  W <- matrix(0, ncol=pars$ntheta, allnodes)
  for(i in 1:allnodes){
    m <- maps[as.logical(cache$branchtrace[i,])]
    m <- c(0, rev(unlist(lapply(m, rev))))
    names(m)[1] <- 1
    TH <- sum(m)
    csm <- cumsum(m)
    eT <- exp(-a*TH)*(exp(a*csm[2:length(csm)])-exp(a*csm[1:(length(csm)-1)]))
    w <- tapply(eT, names(csm)[2:length(csm)], sum)
    W[i, as.numeric(names(w))] <- w
    W[i, 1] <- W[i,1] + exp(-a*TH)
  }
  return(W)
}

#' Experimental phenogram plotting function for set of model of model parameters
#'
#' @param pars A bayou formatted parameter list
#' @param tree A tree of class 'phylo'
#' @param dat A named vector of tip data
#' @param regime.col A named vector of colors equal in length to the number of regimes
#' @param SE Standard error of the tip states
#' @param ... Optional arguments passed to \code{phenogram()}
#'
#' @details This is an experimental plotting utility that can plot a phenogram with a given regime painting from
#' a parameter list. Note that it uses optimization of internal node states using matrix inversion, which is very
#' slow for large trees. However, what is returned is the maximum likelihood estimate of the internal node states
#' given the model, data and the parameter values.
#'
#'
#' @examples
#' \donttest{
#' tree <- sim.bdtree(n=50)
#' tree$edge.length <- tree$edge.length/max(branching.times(tree))
#' prior <- make.prior(tree, dists=list(dk="cdpois", dsig2="dnorm",
#'            dtheta="dnorm"), param=list(dk=list(lambda=5, kmax=10),
#'              dsig2=list(mean=1, sd=0.01), dtheta=list(mean=0, sd=3)),
#'                plot.prior=FALSE)
#' pars <- priorSim(prior, tree, plot=FALSE, nsim=1)$pars[[1]]
#' pars$alpha <- 4
#' dat <- dataSim(pars, model="OU", phenogram=FALSE, tree)$dat
#' OUphenogram(pars, tree, dat, ftype="off")
#' }
#'
#' @return No return value. This function generates a phenogram plot as a side effect.
#' @export
OUphenogram <- function(pars, tree, dat, SE=0, regime.col=NULL, ...){
  datanc <- .OU.asr(tree, dat, pars, SE=SE)
  tr <- pars2simmap(pars, reorder(tree,"postorder"))
  OA <- .optima.ages(pars, tr$tree)
  CIU95 <- pars$theta+2*sqrt(pars$sig2/(2*pars$alpha))
  CIL95 <- pars$theta-2*sqrt(pars$sig2/(2*pars$alpha))
  if(is.null(regime.col)){
    regime.cols <- tr$col
  } else {regime.cols <- regime.col}
  ylim <- c(min(c(dat, pars$theta-2*sqrt(pars$sig2/(pars$alpha*2)))), max(c(dat, pars$theta-2*sqrt(pars$sig2/(pars$alpha*2)))))
  phenogram(tr$tree, datanc, colors=regime.cols, ylim=ylim, spread.labels=FALSE, ...)
  for(i in 1:pars$ntheta){
    x <- seq(OA[i,1],OA[i,2],length=10)
    y <- seq(ylim[1],ylim[2],length=100)
    Z <- matrix(nrow=length(x),ncol=length(y))
    for(j in 1:length(x)){
      Z[j,] <- dnorm(y,pars$theta[i],sqrt(pars$sig2/(2*pars$alpha)))
    }
    if(sum(Z)!=0){
      densregion(x,y,Z,colmax=makeTransparent(regime.cols[i]),colmin="transparent")
    }
    lines(c(OA[i,1],OA[i,2]),rep(pars$theta[i],2),col=regime.cols[i],lwd=2)
  }
  phenogram(tr$tree, datanc, , colors=regime.cols, add=TRUE, spread.labels=FALSE,  ...)
}

#' Function to plot the regimes from a simmap tree
#'
#' @param tree A simmap tree of class phylo or simmap with a tree$maps list
#' @param col A named vector of colors to assign to character states, if NULL, then colors are generated from pal
#' @param lwd A numeric value indicating the width of the edges
#' @param pal A color palette function to generate colors if col=NULL
#' @param ... Optional arguments that are passed to plot.phylo
#'
#' @details This function uses plot.phylo to generate coordinates and plot the tree, but plots the
#' 'maps' element of phytools' simmap format. This provides much of the functionality of plot.phylo from
#' the ape package. Currently, only types 'phylogram', 'unrooted', 'radial', and 'cladogram' are allowed. Phylogenies must
#' have branch lengths.
#'
#' @return **No return value**, called for **side effects**.
#' The function **generates a plot** of the phylogenetic tree with color-coded regimes.
#'
#' @export
plotRegimes <- function(tree, col=NULL, lwd=1, pal=rainbow, ...){
  if(is.null(col)){
    regNames <- unique(names(unlist(tree$maps)))
    nreg <- length(regNames)
    col <- setNames(pal(nreg), regNames)
  }
  #nodecols <- col[sapply(tree$maps, function(x) names(x)[1])]
  tmp <- ape::plot.phylo(tree, edge.color="#FFFFFF00", use.edge.length=TRUE, ...)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  #if(lastPP$type != "phylogram") stop("Currently only able to plot phylograms")
  nbranch <- nrow(tree$edge)
  .getBranchCoords <- function(i){
    xx <- lastPP$xx[tree$edge[i,]]
    yy <- lastPP$yy[tree$edge[i,]]
    xdist <- diff(xx)
    ydist <- diff(yy)
    map <- tree$maps[[i]]
    cs <- cumsum(c(0, map))/sum(map)
    colmap <- col[names(map)]
    return(list(xx=xx, yy=yy, xdist=xdist, ydist=ydist, cs=cs, colmap=colmap, nsegs=length(cs)-1, segreg = names(colmap)))
  }
  coords <- lapply(1:nbranch, .getBranchCoords)
  .phylogramLines <- function(x){
    xdist <- x$xdist; ydist <- x$ydist; xx <- x$xx; yy <- x$yy
    cs <- x$cs; nsegs <- x$nsegs; segreg <- x$segreg; colmap <- x$colmap
    if(lastPP$direction %in% c("upwards", "downwards")){
      xcoord <- rbind(xx, matrix(xx[2], nrow=nsegs, ncol=2))
      ycoord <- rbind(rep(yy[1],2), cbind(cs[1:(length(cs)-1)]*ydist+yy[1], cs[2:(length(cs))]*ydist+yy[1]))
      rownames(xcoord) <- rownames(ycoord) <- c(segreg[1], segreg)
      cols <- c(colmap[1], colmap)
      dum <- lapply(1:(nsegs+1), function(i) lines(xcoord[i,], ycoord[i, ], col=cols[i], lwd=lwd))
    }
    if(lastPP$direction %in% c("leftwards", "rightwards")){
      ycoord <- rbind(yy, matrix(yy[2], nrow=nsegs, ncol=2))
      xcoord <- rbind(rep(xx[1],2), cbind(cs[1:(length(cs)-1)]*xdist+xx[1], cs[2:(length(cs))]*xdist+xx[1]))
      rownames(xcoord) <- rownames(ycoord) <- c(segreg[1], segreg)
      cols <- c(colmap[1], colmap)
      dum <- lapply(1:(nsegs+1), function(i) lines(xcoord[i,], ycoord[i, ], col=cols[i], lwd=lwd))
    }
  }
  .cladogramLines <- function(x){
    xdist <- x$xdist; ydist <- x$ydist; xx <- x$xx; yy <- x$yy
    cs <- x$cs; nsegs <- x$nsegs; segreg <- x$segreg; colmap <- x$colmap
    xcoord <- cbind(cs[1:(length(cs)-1)]*xdist+xx[1], cs[2:(length(cs))]*xdist+xx[1])
    ycoord <- cbind(cs[1:(length(cs)-1)]*ydist+yy[1], cs[2:(length(cs))]*ydist+yy[1])
    rownames(xcoord) <- rownames(ycoord) <- segreg
    cols <- colmap
    dum <- lapply(1:nsegs, function(i) lines(xcoord[i,], ycoord[i, ], col=cols[i], lwd=lwd))
  }
  .fanLines <- function(x){
    xdist <- x$xdist; ydist <- x$ydist; xx <- x$xx; yy <- x$yy
    cs <- x$cs; nsegs <- x$nsegs; segreg <- x$segreg; colmap <- x$colmap
    circular.plot(lastPP$edge, lastPP$Ntip, lastPP$Nnode, lastPP$xx, lastPP$yy, )
  }
  if(lastPP$type=="fan") warning("type='fan' not currently supported, plotting a radial cladogram")
  plotfn <- switch(lastPP$type, phylogram=.phylogramLines, cladogram=.cladogramLines, unrooted=.cladogramLines, radial=.cladogramLines, fan=.cladogramLines)
  dum <- lapply(coords, plotfn)
}


#' A function for summarizing the state of a model after a shift
#'
#' @param chain A bayouMCMC chain
#' @param mcmc A bayou mcmc object
#' @param pp.cutoff The threshold posterior probability for shifts to summarize, if 'branches'
#' specified than this is ignored.
#' @param branches The specific branches with shifts to summarize, assuming postordered tree
#'
#' @details shiftSummaries summarizes the immediate parameter values after a shift on a particular
#' branch. Parameters are summarized only for the duration that the particular shift exists. Thus,
#' even global parameters will be different for particular shifts.
#'
#' @return A list with elements:
#' \code{pars} = a bayoupars list giving the location of shifts specified;
#' \code{tree} = The tree;
#' \code{pred} = Predictor variable matrix;
#' \code{dat} = A vector of the data;
#' \code{SE} = A vector of standard errors;
#' \code{PP} = Posterior probabilities of the specified shifts;
#' \code{model} = A list specifying the model used;
#' \code{variables} = The variables summarized;
#' \code{cladesummaries} = A list providing the medians and densities of the distributions of regression
#' variables for each shift;
#' \code{descendents} = A list providing the taxa that belong to each regime
#' \code{regressions} = A matrix providing the regression coefficients for each regime.
#' @export
shiftSummaries <- function(chain, mcmc, pp.cutoff=0.3, branches=NULL){
  cache <- .prepare.ou.univariate(mcmc$tree,mcmc$dat, SE=mcmc$SE, pred=mcmc$pred)
  tree <- cache$phy
  dat <- cache$dat
  pred <- cache$pred
  SE <- cache$SE
  model <- mcmc$model.pars

  if(is.null(attributes(chain)$burnin)){
    L <- Lposterior(chain, tree, burnin=0)
  } else {
    L <- Lposterior(chain, tree, burnin=attributes(chain)$burnin)
  }
  if(is.null(branches)){
    branches <- which(L[,1] > pp.cutoff)
    PP <- L[branches,1]
  } else {
    PP <- L[branches,1]
  }
  if(length(branches) == 0){
    stop("No shifts found with posterior probability above cutoff")
  }

  if(!is.null(model$call)){
    coefs <- paste("beta_", attr(terms(model$call), "term.labels"), sep="")
    coefs <- gsub(":", "x", coefs)
  } else {
    coefs <- NULL
  }
  variables <- c('theta', coefs)
  sumpars <- list(k=length(branches), ntheta=length(branches)+1, sb=branches, t2=2:(length(branches)+1), loc=rep(0, length(branches)))
  .summarizeDerivedState <- function(branch, chain, variables){
    if(branch==0){
      values <- lapply(variables, function(x) sapply(chain[[x]], function(y) y[1]))
    } else {
      SB <- unlist(chain$sb)
      gen <- unlist(lapply(1:length(chain$sb), function(x) rep(x, length(chain$sb[[x]]))))
      ind <- which(SB==branch)
      gen <- gen[ind]
      T2 <- unlist(chain$t2)[ind]
      values <- lapply(variables, function(x) sapply(1:length(T2), function(y)if(length(chain[[x]][[gen[y]]]) > 1) {chain[[x]][[gen[y]]][T2[y]]} else chain[[x]][[gen[y]]]))
    }
    medians <- lapply(values, median)
    densities <- lapply(values, density)
    names(medians) <- names(densities) <- variables
    return(list(medians=medians, densities=densities))
  }
  cladesummaries <- lapply(c(0, branches), function(x) .summarizeDerivedState(x, chain, variables))
  regressions <- do.call(rbind, lapply(cladesummaries, function(x) unlist(x$medians)))
  rownames(regressions) = c("root", sumpars$sb)
  tipregs <- .tipregime(sumpars, tree)
  descendents <- lapply(1:(length(sumpars$sb)+1), function(x) names(tipregs[tipregs==x]))
  out <- list(pars= sumpars, tree=tree, pred=pred, dat=dat, SE=SE, PP=PP, model=model, variables=variables, cladesummaries=cladesummaries, descendents=descendents, regressions=regressions)
  return(out)
}

#' A function to plot a list produced by \code{shiftSummaries}
#'
#' @param summaries A list produced by the function \code{shiftSummaries}
#' @param pal A color palette function
#' @param ask Whether to wait for the user between plotting each shift summary
#' @param single.plot A logical indicating whether to summarize all shifts in a single plot.
#' @param label.pts A logical indicating whether to label the scatter plot.
#' @param ... Additional parameters passed to the function par(...)
#'
#' @details For each shift, this function plots the taxa on the phylogeny that are (usually) in this regime (each taxon
#' is assigned to the specified shifts, thus some descendent taxa may not always be in indicated regime if the shift if
#' they are sometimes in another tipward shift with low posterior probability). The function then plots the distribution
#' of phenotypic states and the predicted regression line, as well as density plots for the intercept and any regression
#' coefficients in the model.
#'
#' @return **No return value**, called for **side effects**.
#' The function **generates visualizations** of shift summaries, including:
#' \itemize{
#'   \item **Phylogenetic tree with shift locations**
#'   \item **Scatter plots of phenotype data**
#'   \item **Density plots for regression coefficients**
#' }
#'
#' @export
plotShiftSummaries <- function(summaries, pal=rainbow, ask=FALSE, single.plot=FALSE, label.pts=TRUE, ...){
  #oldpar <- graphics::par(no.readonly = TRUE)    # code line i
  #on.exit(graphics::par(oldpar))            # code line i + 1
  #px <- par()
  ndens <- length(summaries$cladesummaries[[1]]$densities)
  #par(mfrow=c(2,max(ndens,2)), mar=c(3,3,5,1), bg="black", ask=FALSE, col.axis="white", col.lab="white", col.main="white", ...)
  blank.panels <- prod(par()$mfrow) - (2+ndens)
  #par(ask=ask)
  regressions <- summaries$regressions
  if(ncol(regressions)==1){ regressions <- data.frame(regressions, "slope"=0)}
  dat <- summaries$dat
  tree <- summaries$tree
  sumpars <- summaries$pars
  descendents <- summaries$descendents
  PP <- c("Root",round(summaries$PP,2))
  xlimits <- apply(do.call(rbind,
                           lapply(summaries$cladesummaries, function(x)
                             sapply(x$densities, function(y) range(y$x))
                           )), 2, range)
  xlimits[1,] <- xlimits[1,]-0.1*apply(xlimits, 2, diff)
  xlimits[2,] <- xlimits[2,]+0.1*apply(xlimits, 2, diff)
  if(ndens > 1){
    xint <- setNames(data.frame(summaries$pred)[[1]], names(dat))
    xlimits2 <- range(xint)
  } else {
    xint <- jitter(.tipregime(sumpars, tree))
    xlimits2 <- c(-2, sumpars$ntheta+3)
  }
  if(!single.plot){
    for(i in (1:nrow(regressions))){
      plotBayoupars(sumpars, tree, col=setNames(c(pal(nrow(regressions))[i], rep("gray80", nrow(regressions)-1)), c(i, (1:nrow(regressions))[-i])), cex=0.2)
      plot(xint, dat, pch=21, xlim=xlimits2, bg=makeTransparent("gray80", 100), col =makeTransparent("gray80", 10), main=paste("Posterior prob: ", PP[i], sep=""))
      if(length(descendents[[i]] > 0)){
        if(label.pts) text(xint[descendents[[i]]], dat[descendents[[i]]], labels=names(dat[descendents[[i]]]), col="white", cex=0.4, pos = 2)
        points(xint[descendents[[i]]], dat[descendents[[i]]], pch=21, bg=makeTransparent(pal(nrow(regressions))[i], 100), col =makeTransparent(pal(nrow(regressions))[i], 10))
      } else{
        warnings("No descendents for this shift")
      }
      abline(a=regressions[i,1], b=regressions[i,2], col=pal(nrow(regressions))[i], lwd=2, lty=2)
      dens <- summaries$cladesummaries[[i]]$densities
      gbg <- lapply(1:length(dens), function(y)plot(dens[[y]], col=pal(nrow(regressions))[i], main=names(dens)[y],xlim=c(xlimits[,y])))
      if(blank.panels >0){lapply(1:blank.panels,function(x) plot.new())}
    }
  } else {
    plotBayoupars(sumpars, tree, col=setNames(pal(sumpars$ntheta), 1:sumpars$ntheta), cex=0.2, tip.col="white")
    plot(xint, dat, pch=21, xlim=xlimits2, bg=makeTransparent("gray80", 100), col =makeTransparent("gray80", 10))
    for(i in 1:length(descendents)){
      if(length(descendents[[i]] > 0)){
        if(label.pts) text(xint[descendents[[i]]], dat[descendents[[i]]], labels=names(dat[descendents[[i]]]), col="white", cex=0.4, pos = 2)
        points(xint[descendents[[i]]], dat[descendents[[i]]], pch=21, bg=makeTransparent(pal(nrow(regressions))[i], 100), col =makeTransparent(pal(nrow(regressions))[i], 10))
        } else{
          warnings("No descendents for this shift")
        }
      abline(a=regressions[i,1], b=regressions[i,2], col=pal(nrow(regressions))[i], lwd=2, lty=2)
    }
    varN <- length(summaries$cladesummaries[[1]]$densities)
    for(j in 1:varN){
      dens <- lapply(1:length(summaries$cladesummaries), function(x) summaries$cladesummaries[[x]]$densities[[j]])
      xrange <- quantile((do.call(c, lapply(dens, function(q) q$x))), c(0.01, 0.99))
      ymax <- max(do.call(c, lapply(dens, function(q) q$y)))
      plot(xrange, c(0, ymax*1.1), type="n", main=names(summaries$cladesummaries[[1]]$densities)[j])
      gbg <- lapply(1:length(dens), function(y) lines(dens[[y]], col=pal(nrow(regressions))[y]))
    }

  }
  #px <- px[!(names(px) %in% c("cin", "cra", "cxy", "csi", "din", "page"))]
  #suppressWarnings(par(px))
}

#' A function to plot a heatmap of reconstructed parameter values on the branches of the tree
#'
#' @param tree A phylogenetic tree
#' @param chain A bayou MCMC chain
#' @param variable The parameter to reconstruct across the tree
#' @param burnin The initial proportion of burnin samples to discard
#' @param nn The number of discrete categories to divide the variable into
#' @param pal A color palette function that produces nn colors
#' @param legend_ticks The sequence of values to display a legend for
#' @param legend_settings A list of legend attributes (passed to bayou:::.addColorBar)
#' @param ... Additional options passed to plot.phylo
#'
#' @details legend_settings is an optional list of any of the following:
#'
#' legend - a logical indicating whether a legend should be plotted
#'
#' x - the x location of the legend
#'
#' y - the y location of the legend
#'
#' height - the height of the legend
#'
#' width - the width of the legend
#'
#' n - the number of gradations in color to plot from the palette
#'
#' adjx - an x adjustment for placing text next to the legend bar
#'
#' cex.lab - the size of text labels next to the legend bar
#'
#' text.col - The color of text labels
#'
#' locator - if TRUE, then x and y coordinates are ignored and legend is placed
#' interactively.
#' @return No return value, called for side effects. This function generates
#' a heatmap visualization of reconstructed parameter values on a phylogenetic tree.

#' @export
plotBranchHeatMap <- function(tree, chain, variable, burnin=0, nn=NULL, pal=heat.colors, legend_ticks=NULL, legend_settings=list(plot=TRUE), ...){
  dum <- setNames(rep(1, length(tree$tip.label)), tree$tip.label)
  cache <- .prepare.ou.univariate(tree, dum)
  tree <- cache$phy
  seq1 <- floor(max(seq(burnin*length(chain$gen),1), length(chain$gen), 1))
  if(is.null(legend_ticks)){
    legend_ticks <- seq(min(unlist(chain[[variable]][seq1],F,F)), max(unlist(chain[[variable]][seq1],F,F)), length.out=5)
  }
  if(is.null(nn)) nn <- length(seq1) else { seq1 <- floor(seq(max(burnin*length(chain$gen),1), length(chain$gen), length.out=nn))}
  if(length(nn) > length(chain$gen)) stop("Number of samples greater than chain length, lower nn")
  abranches <- lapply(1:nrow(tree$edge), .ancestorBranches, cache=cache)
  allbranches <- suppressWarnings(sapply(1:nrow(tree$edge), function(x) .branchRegime(x, abranches, chain, variable, seq1, summary=TRUE)))
  ape::plot.phylo(tree, edge.color=.colorRamp(allbranches, pal, 100), ...)
  lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  legend_stuff <- list(x=0.01* lastPP$x.lim[2],
                       y=0,
                       height=0.25*diff(lastPP$y.lim),
                       width=0.01*diff(lastPP$x.lim),
                       n=100,
                       trait=allbranches,
                       ticks=legend_ticks,
                       adjx=0.01*lastPP$x.lim[2],
                       cex.lab=0.5,
                       text.col="black",
                       plot=TRUE,
                       locator=FALSE
  )
  if(length(legend_settings) > 0){
    for(i in 1:length(legend_settings)){
      legend_stuff[[names(legend_settings)[i]]] <- legend_settings[[i]]
    }
  }
  if(legend_stuff$plot) {
    if(legend_stuff$locator){
      lc <- locator(1)
      legend_stuff$x <- lc$x
      legend_stuff$y <- lc$y
      .addColorBar(x=legend_stuff$x, y=legend_stuff$y, height=legend_stuff$height, width=legend_stuff$width, pal=pal, n=legend_stuff$n, trait=allbranches, ticks=legend_ticks, adjx=legend_stuff$adjx, cex.lab=legend_stuff$cex.lab, text.col=legend_stuff$text.col)
    } else .addColorBar(x=legend_stuff$x, y=legend_stuff$y, height=legend_stuff$height, width=legend_stuff$width, pal=pal, n=legend_stuff$n, trait=allbranches, ticks=legend_ticks, adjx=legend_stuff$adjx, cex.lab=legend_stuff$cex.lab, text.col=legend_stuff$text.col)
  }
}


.ancestorBranches <- function(branch, cache){
  ancbranches <- which(sapply(cache$bdesc, function(x) branch %in% x))
  sort(ancbranches, decreasing=FALSE)
}
.branchRegime <- function(branch, abranches, chain, parameter, seqx, summary=FALSE){
  ancs <- c(branch, abranches[[branch]])
  ancshifts <- lapply(1:length(seqx), function(x) chain$t2[[seqx[x]]][which(chain$sb[[seqx[x]]] == ancs[min(which(ancs %in% chain$sb[[seqx[x]]]))])])
  ancshifts <- sapply(ancshifts, function(x) ifelse(length(x)==0, 1, x))
  ests <- sapply(1:length(ancshifts), function(x) chain[[parameter]][[seqx[x]]][ancshifts[x]])
  res <- cbind(ests)
  if(summary){
    return(apply(res, 2, stats::median))
  } else {
    return(res)
  }
}

