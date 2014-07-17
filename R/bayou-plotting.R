#' Utility for getting the starting and ending ages for each regime
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
#' 
#' @export
makeTransparent <- function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#' Plot a phylogenetic tree with posterior probabilities from a bayouMCMC chain (function adapted from phytools' plotSimmap)
#' 
#' @param tree A tree of class phylo
#' @param chain A bayouMCMC chain
#' @param burnin The proportion of runs to be discarded
#' @param colors A named vector of colors corresponding to the regimes in the stored simmap map in the tree
#' @param fsize relative font size for tip labels
#' @param ftype font type- options are "reg", "i" (italics), "b" (bold), or "bi" (bold-italics)
#' @param lwd line width for plotting
#' @param node.numbers a logical value indicating whether or not node numbers should be plotted
#' @param mar vector containing the margins for the plot to be passed to par. If not specified, the default margins are (0.1,0.1,0.1,0.1)
#' @param add a logical value indicating whether or not the tree should be plotted in the current or a new plot
#' @param offset offset for the tip labels
#' @param pp.cutoff the posterior probability above which a shift should be reconstructed on the tree. If this value is provided, it overrides any preexisting simmap associated with the tree.
#' @param circle a logical value indicating whether or not a circle should be plotted at the base of the node with values that correspond to the posterior probability of having a shift.
#' @param circle.pch the type of symbol used to plot at the node to indicate posterior probability
#' @param circle.pal a palette of colors that will be used to color the interior of the circles. This will be varied over the interval proportional to the deviation that occurs at that shift on the phylogeny.
#' @param circle.lwd the line width of the points plotted at the nodes
#' @param circle.alpha a value between 0 and 255 that indicates the transparency of the circles (255 is completely opaque).
#' @param dash a logical value indicating whether or not a dash (or other point) should be plotted on the branches at every location that a shift was present in the posterior chain.
#' @param dash.pal Ignored for now
#' @param dash.cex the relative size of dashes
#' @param dash.pch the plotting symbol for dashes
#' @param dash.lwd the line width for dashes
#' @param dash.alpha the transparency of dashes
#' @param pp.labels a logical indicating whether the posterior probability for each branch should be printed above the branch
#' @param pp.alpha a logical or numeric value indicating transparency of posterior probability labels. If TRUE, then transparency is ramped from invisible (pp=0), to black (pp=1). If numeric, all labels are given the same transparency. If NULL, then no transparency is given. 
#' @param pp.cex the size of the posterior probability labels 
#' @param legend Logical indicating whether or not a legend should be produced
#' @param limits Divergence values assigning the minimum and maximum amounts applied to the extremes of the color palette. Values outside the limits are assigned the most extreme values in the color palette function.
#' 
#' @export
plotSimmap.mcmc <- function (tree, chain, burnin=NULL, colors = NULL, fsize = 1, ftype = "reg", lwd = 0.75, 
                             node.numbers = FALSE, mar = NULL, add = FALSE, offset = NULL, pp.cutoff=NULL,legend=TRUE, limits=NULL, 
                             circle=TRUE, circle.pch=21, circle.pal=cm.colors, circle.lwd=0.75, circle.alpha=200, 
                             dash=FALSE, dash.pal=cm.colors, dash.cex=0.5, dash.pch="|", dash.lwd=0.5, dash.alpha=10, pp.labels=FALSE, pp.alpha=NULL, pp.cex=0.75) {
  dummy <- setNames(rep(1, length(tree$tip.label)), tree$tip.label)
  cache <- .prepare.ou.univariate(tree, dummy)
  tree <- cache$phy
  if(is.null(burnin)) burnin = attributes(chain)$burnin
  if(burnin==0) postburn <- 1:length(chain$gen)  else {
    postburn <- round(burnin*length(chain$gen),0):length(chain$gen)
  }
  L <- Lposterior(chain, tree)
  pull.map <- function(x, chain){
    pars.list <- lapply(x, function(i) list(k=chain$k[[i]], ntheta=chain$ntheta[[i]], theta=chain$theta[[i]], sb=chain$sb[[i]], t2=chain$t2[[i]], loc=chain$loc[[i]]))
    return(pars.list)
  }
  map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal(100)[findInterval(x,seq(limits[1],limits[2],length.out=100+1), all.inside=TRUE)]
  }
  pars.list <- pull.map(postburn, chain)
  Div <- unlist(lapply(pars.list, function(x) .D.from.theta(x, cache)$D), F,F)
  sb <- unlist(lapply(pars.list, function(x) x$sb), F, F)
  sb <- factor(sb, levels=1:nrow(tree$edge))
  ave.Div <- tapply(Div, sb, mean)
  ave.Div[is.na(ave.Div)] <- 0
  
  if(!is.null(pp.cutoff)){
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
  if (is.null(tree$maps)){
    tree$maps <- lapply(tree$edge.length,function(x){names(x) <- 1; x})
  } 
  nb <- length(unique(unlist(sapply(tree$maps,names))))
  ftype <- which(c("off", "reg", "b", "i", "bi") == ftype) - 1
  if (!ftype) 
    fsize = 0
  if (is.null(colors)) {
    colors <- rainbow(nb-1)
    colors <- c("gray50",colors)
    names(colors) <- as.character(1:nb)
  }
  if (class(tree) != "phylo") 
    stop("tree should be object of class 'phylo.'")
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  cw <- reorderSimmap(tree)
  o <- .whichorder(cw$edge[,2],tree$edge[,2])
  ob <- (1:length(tree$edge.length))[o]
  #shifts.o <- shifts[o,]
  pw <- reorderSimmap(tree, "pruningwise")
  n <- length(cw$tip)
  m <- cw$Nnode
  Y <- matrix(NA, m + n, 1)
  Y[cw$edge[cw$edge[, 2] <= length(cw$tip), 2]] <- 1:n
  nodes <- unique(pw$edge[, 1])
  for (i in 1:m) {
    desc <- pw$edge[which(pw$edge[, 1] == nodes[i]), 
                    2]
    Y[nodes[i]] <- (min(Y[desc]) + max(Y[desc]))/2
  }
  root <- length(cw$tip) + 1
  node.height <- matrix(NA, nrow(cw$edge), 2)
  for (i in 1:nrow(cw$edge)) {
    if (cw$edge[i, 1] == root) {
      node.height[i, 1] <- 0
      node.height[i, 2] <- cw$edge.length[i]
    } else {
      node.height[i, 1] <- node.height[match(cw$edge[i, 
                                                     1], cw$edge[, 2]), 2]
      node.height[i, 2] <- node.height[i, 1] + cw$edge.length[i]
    }
  }
  ### Here is stuff I've added
  b.shift <- unlist(chain$sb[postburn])
  r1 <- tree$edge.length[b.shift]-unlist(chain$loc[postburn])
  loc.shift <- r1
  shifts.o <- node.height[match(b.shift,o),1]+loc.shift
  L <- L[o,]
  ####
  if (is.null(mar)) {
    par(mar = c(0.1, 0.1, 0.1, 0.1))
  } else par(mar = mar)
  if (!add) 
    plot.new()
  if (fsize * max(strwidth(cw$tip.label)) < 1) {
    c <- (1 - fsize * max(strwidth(cw$tip.label)))/max(node.height)
    cw$edge.length <- c * cw$edge.length
    cw$maps <- lapply(cw$maps, function(x) x <- c * x)
    node.height <- c * node.height
  } else message("Font size too large to properly rescale tree to window.")
  height <- max(nodeHeights(tree))
  if (!add) 
    plot.window(xlim = c(0, max(node.height) + fsize * 
                           max(strwidth(cw$tip.label))), ylim = c(1, max(Y)))
  for (i in 1:m) lines(node.height[which(cw$edge[, 1] == nodes[i]), 1], Y[cw$edge[which(cw$edge[, 1] == nodes[i]), 2]], col = colors[names(cw$maps[[match(nodes[i], cw$edge[, 1])]])[1]], lwd = lwd)
  for (i in 1:nrow(cw$edge)) {
    x <- node.height[i, 1]
    #  bs <- shifts.o/height
    # points(bs,rep(Y[cw$edge[i,2]],length(bs)),cex=sh.cex,pch="|",lwd=(lwd-1),col=makeTransparent("#000000",alpha=alpha))
    for (j in 1:length(cw$maps[[i]])) {
      lines(c(x, x + cw$maps[[i]][j]), c(Y[cw$edge[i, 
                                                   2]], Y[cw$edge[i, 2]]), col = colors[names(cw$maps[[i]])[j]], 
            lwd = lwd, lend = 2)

      x <- x + cw$maps[[i]][j]
      j <- j + 1
    }
  }
  
  #    if(type=="dashes"){
  #     bs <- shifts.o*c
  #      points(bs,Y[tree$edge[b.shift,2]],cex=sh.cex,pch="|",lwd=(lwd-1),col=makeTransparent("#000000",alpha=alpha))
  #    }
  if(dash){
    bs <- shifts.o*c
    points(bs,Y[tree$edge[b.shift,2]],cex=dash.cex, pch=dash.pch,lwd=dash.lwd, col=makeTransparent("black", dash.alpha))
  }
  if(circle){
    if(all(ave.Div==0)){
     cols <- rep(circle.pal(1), length(ave.Div)) 
    } else {
      cols <- map2color(ave.Div[o], circle.pal, limits=limits)
    }
    points(node.height[,1], Y[cw$edge[, 2]], cex=sqrt(4^2*L[,1]), lwd=circle.lwd,bg=makeTransparent(cols,alpha=circle.alpha),pch=circle.pch)
    if(legend==TRUE){
    legend_image <- as.raster(matrix(rev(circle.pal(100)), ncol=1))
    text(x=1.5, y = round(seq(range(ave.Div)[1],range(ave.Div)[2],l=5),2), labels = seq(range(ave.Div)[1],range(ave.Div)[2],l=5))
    rasterImage(legend_image, 1, 0.25*length(tree$tip.label), 1.01, length(tree$tip.label)-0.25*length(tree$tip.label))
    text(1.02, seq(0.25*length(tree$tip.label), length(tree$tip.label)-0.25*length(tree$tip.label), length.out=5), labels=round(seq(range(ave.Div)[1], range(ave.Div)[2], length.out=5),2), cex=fsize)
    }
    }
  if(pp.labels){
    if(is.null(pp.alpha)){
      pp.col = "black"
    } else {
      if(is.numeric(pp.alpha)){
        pp.col = makeTransparent("black", pp.alpha)
      }
      if(pp.alpha){
        pp.col <- makeTransparent("black",alpha=(0:255)[findInterval(L[,1],seq(0,1,length.out=256))])
      }
    }
    text(node.height[,1]+0.5*apply(node.height,1,diff), Y[cw$edge[, 2]]-(2.5-fsize), labels=round(L[,1],2), cex=pp.cex, pos=3, col=pp.col)
  }
  
  
  if (node.numbers) {
    symbols(0, mean(Y[cw$edge[cw$edge[, 1] == (length(cw$tip) + 
                                                 1), 2]]), rectangles = matrix(c(1.2 * fsize * 
                                                                                   strwidth(as.character(length(cw$tip) + 1)), 1.4 * 
                                                                                   fsize * strheight(as.character(length(cw$tip) + 
                                                                                                                    1))), 1, 2), inches = F, bg = "white", add = T)
    text(0, mean(Y[cw$edge[cw$edge[, 1] == (length(cw$tip) + 
                                              1), 2]]), length(cw$tip) + 1, cex = fsize)
    for (i in 1:nrow(cw$edge)) {
      x <- node.height[i, 2]
      if (cw$edge[i, 2] > length(tree$tip)) {
        symbols(x, Y[cw$edge[i, 2]], rectangles = matrix(c(1.2 * 
                                                             fsize * strwidth(as.character(cw$edge[i, 
                                                                                                   2])), 1.4 * fsize * strheight(as.character(cw$edge[i, 
                                                                                                                                                      2]))), 1, 2), inches = F, bg = "white", add = T)
        text(x, Y[cw$edge[i, 2]], cw$edge[i, 2], cex = fsize)
      }
    }
  }
  if (is.null(offset)) 
    offset <- 0.2 * lwd/3 + 0.2/3
  for (i in 1:n) if (ftype) 
    text(node.height[which(cw$edge[, 2] == i), 2], Y[i], 
         cw$tip.label[i], pos = 4, offset = offset, cex = fsize, 
         font = ftype)
  par(mar = c(5, 4, 4, 2) + 0.1)
}

#' Adds visualization of regimes to a plot
#' 
#' @param pars A bayou formatted parameter list
#' @param tree A tree of class 'phylo'
#' @param cols A vector of colors to give to regimes, in the same order as pars$sb
#' @param type Either "rect", "density" or "lines". "rect" plots a rectangle for the 95% CI for the stationary
#' distribution of a regime. "density" varies the transparency of the rectangles according to the probability density
#' from the stationary distribution. "lines" plots lines for the mean and 95% CI's without filling them. 
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
  phenogram(tree,dat,add=TRUE, colors=colors)#,...)
  dens.theta <- lapply(1:length(K), function(x) density(unlist(theta[no.theta %in% K[[x]]])))
  tmp <- sapply(1:length(dens.theta),function(Q){lines(nH+dens.theta[[Q]]$y*(0.3*nH)/max(dens.theta[[Q]]$y),dens.theta[[Q]]$x,col=Q+1)})
}

#' S3 method for plotting bayouMCMC objects
#' 
#' @param x A mcmc chain of class 'bayouMCMC' produced by the function bayou.mcmc and loaded into the environment using load.bayou
#' @param ... Additional arguments passed to \code{plot.mcmc} from the \code{coda} package
#' 
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
  rownames(uni.df) <- uni.df[,1]
  uni.df <- uni.df[,-1]
  plot(mcmc(uni.df), ...)
}


#' Plot parameter list as a simmap tree
#' 
#' @param pars A bayou formatted parameter list
#' @param tree A tree of class 'phylo'
#' @param ... Additional arguments passed to plotSimmap
#' 
#' @export
plotBayoupars <- function(pars, tree,...){
  mar <- par()$mar
  tree <- reorder(tree, 'postorder')
  X <- rep(0, length(tree$tip.label))
  names(X) <- tree$tip.label
  cache <- .prepare.ou.univariate(tree, X)
  tr <- .toSimmap(.pars2map(pars, cache),cache)
  plotSimmap(tr,...)
  par(mar=mar)
}

#' Experimental function for ancestral state reconstruction for a given OU model
.OU.asr <- function(tree, dat, pars, start=NULL, SE=0){
  phy <- reorder(tree, "postorder")
  dat <- dat[phy$tip.label]
  if(length(SE)>1){
    SE[phy$tip.label]
  }
  if(length(phy$tip.label) > 100) cat("This may take a while for large trees")
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
  result <- optim(start, lik.fn, method="L-BFGS-B")
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
  if(class(tree)=="phylo"){
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
#' @examples
#' \dontrun{
#' tree <- sim.bdtree(n=50)
#' tree$edge.length <- tree$edge.length/max(branching.times(tree))
#' prior <- make.prior(tree, dists=list(dk="cdpois", dsig2="dnorm", dtheta="dnorm"), 
#'              param=list(dk=list(lambda=5, kmax=10), dsig2=list(mean=1, sd=0.01), 
#'                    dtheta=list(mean=0, sd=3)), plot.prior=FALSE)
#' pars <- priorSim(prior, tree, plot=FALSE, nsim=1)$pars[[1]]
#' pars$alpha <- 4
#' dat <- dataSim(pars, model="OU", phenogram=FALSE, tree)$dat
#' OUphenogram(pars, tree, dat, ftype="off")
#' }
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
  phenogram(tr$tree, datanc, colors=regime.cols, ylim=ylim, ...)
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
  phenogram(tr$tree, datanc, , colors=regime.cols, add=TRUE, ...)
}