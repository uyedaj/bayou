optima.ages <- function(pars,tree){
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

regime.plot <- function(pars,tree,cols,type='rect',model="OU"){
  if(model=="QG"){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model=="OUrepar"){
    repar <- OU.repar(pars)
    pars$alpha <- repar$alpha
    pars$sig2 <- repar$sig2
  }
  OA <- optima.ages(pars,tree)
  CIU95 <- pars$optima+2*sqrt(pars$sig2/(2*pars$alpha))
  CIL95 <- pars$optima-2*sqrt(pars$sig2/(2*pars$alpha))
  if(type=="lines"){
    for(i in 1:pars$ntheta){
      lines(c(OA[i,1],OA[i,2]),rep(pars$optima[i],2),col=cols[i],lwd=3)
      lines(c(OA[i,1],OA[i,2]),rep(CIU95[i],2),col=cols[i],lwd=1.25,lty=2)
      lines(c(OA[i,1],OA[i,2]),rep(CIL95[i],2),col=cols[i],lwd=1.25,lty=2)
    }
  }
  if(type=="rect"){
    for(i in 1:pars$ntheta){
      rect(OA[i,1],CIL95[i],OA[i,2],CIU95[i],col=cols[i],border=NA)
    }
  }
  if(type=="density"){
    ylim <- par('usr')[3:4]
    for(i in 1:pars$ntheta){
      x <- seq(OA[i,1],OA[i,2],length=10)
      y <- seq(ylim[1],ylim[2],length=100)
      Z <- matrix(nrow=length(x),ncol=length(y))
      for(j in 1:length(x)){
        Z[j,] <- dnorm(y,pars$optima[i],sqrt(pars$sig2/(2*pars$alpha)))
      }
      if(sum(Z)!=0){
        densregion(x,y,Z,colmax=cols[i],colmin="transparent")
      }
      lines(c(OA[i,1],OA[i,2]),rep(pars$optima[i],2),col=cols[i],lwd=2)
    }
  }
}

makeTransparent <- function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#' @export
plotSimmap.mcmc <- function (tree,chain,burnin=NULL,colors = NULL, fsize = 1, ftype = "reg", lwd = 2, 
                             pts = TRUE, node.numbers = FALSE, mar = NULL, add = FALSE, 
                             offset = NULL,alpha=10,sh.cex=1,type="dashes",circle.col=NULL,pch=21) {
  tree <- reorder.phylo(tree,order="postorder")
  if(type=="circles"){
    L <- Lposterior(chain,tree)
  }
  if (is.null(tree$maps)){
    tree$maps <- lapply(tree$edge.length,function(x){names(x) <- 1; x})
  } 
  nb <- length(unique(unlist(sapply(tree$maps,names))))
  if (class(tree) == "multiPhylo") {
    par(ask = TRUE)
    for (i in 1:length(tree)) plotSimmap(tree[[i]], colors = colors, 
                                         fsize = fsize, ftype = ftype, lwd = lwd, pts = pts, 
                                         node.numbers = node.numbers)
  } else {
    
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
    o <- phytools:::whichorder(cw$edge[,2],tree$edge[,2])
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
      }
      else {
        node.height[i, 1] <- node.height[match(cw$edge[i, 
                                                       1], cw$edge[, 2]), 2]
        node.height[i, 2] <- node.height[i, 1] + cw$edge.length[i]
      }
    }
    #shifts.o <- node.height[,1]+shifts.o
    postburn <- round(burnin*length(chain$branch.shift),0):length(chain$branch.shift)
    b.shift <- unlist(chain$branch.shift[postburn])
    r1 <- tree$edge.length[b.shift]-unlist(chain$location[postburn])
    loc.shift <- r1
    shifts.o <- node.height[match(b.shift,o),1]+loc.shift
    L <- L[o,]
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
    for (i in 1:m) lines(node.height[which(cw$edge[, 1] == 
                                             nodes[i]), 1], Y[cw$edge[which(cw$edge[, 1] == nodes[i]), 
                                                                      2]], col = colors[names(cw$maps[[match(nodes[i], 
                                                                                                             cw$edge[, 1])]])[1]], lwd = lwd)
    for (i in 1:nrow(cw$edge)) {
      x <- node.height[i, 1]
      #  bs <- shifts.o/height
      # points(bs,rep(Y[cw$edge[i,2]],length(bs)),cex=sh.cex,pch="|",lwd=(lwd-1),col=makeTransparent("#000000",alpha=alpha))
      for (j in 1:length(cw$maps[[i]])) {
        lines(c(x, x + cw$maps[[i]][j]), c(Y[cw$edge[i, 
                                                     2]], Y[cw$edge[i, 2]]), col = colors[names(cw$maps[[i]])[j]], 
              lwd = lwd, lend = 2)
        if (pts)
          points(c(x, x + cw$maps[[i]][j]), c(Y[cw$edge[i, 
                                                        2]], Y[cw$edge[i, 2]]), pch = 20, lwd = (lwd - 
                                                                                                   1))
        x <- x + cw$maps[[i]][j]
        j <- j + 1
      }
    }
    
    if(type=="dashes"){
      bs <- shifts.o*c
      points(bs,Y[tree$edge[b.shift,2]],cex=sh.cex,pch="|",lwd=(lwd-1),col=makeTransparent("#000000",alpha=alpha))
    }
    if(type=="circles"){
      if(is.null(circle.col)){circle.col="#000000"} else {circle.col=circle.col[o]}
      points(node.height[,1], Y[cw$edge[, 2]],cex=L[,1]*5,lwd=1.2,bg=makeTransparent(circle.col,alpha=alpha),pch=pch)
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
  }
  par(mar = c(5, 4, 4, 2) + 0.1)
}

regime.plot2 <- function(pars,tree,cols,type='rect',alpha=255){
  OA <- optima.ages(pars,tree)
  CIU95 <- pars$optima+2*sqrt(pars$sig2/(2*pars$alpha))
  CIL95 <- pars$optima-2*sqrt(pars$sig2/(2*pars$alpha))
  if(type=="lines"){
    for(i in 1:pars$ntheta){
      lines(c(OA[i,1],OA[i,2]),rep(pars$optima[i],2),col=cols[i],lwd=3)
      lines(c(OA[i,1],OA[i,2]),rep(CIU95[i],2),col=cols[i],lwd=1.25,lty=2)
      lines(c(OA[i,1],OA[i,2]),rep(CIL95[i],2),col=cols[i],lwd=1.25,lty=2)
    }
  }
  if(type=="rect"){
    for(i in 1:pars$ntheta){
      rect(OA[i,1],CIL95[i],OA[i,2],CIU95[i],col=cols[i],border=NA)
    }
  }
  if(type=="density"){
    ylim <- par('usr')[3:4]
    for(i in 1:pars$ntheta){
      x <- seq(OA[i,1],OA[i,2],length=10)
      y <- seq(ylim[1],ylim[2],length=100)
      Z <- matrix(nrow=length(x),ncol=length(y))
      for(j in 1:length(x)){
        Z[j,] <- dnorm(y,pars$optima[i],sqrt(pars$sig2/(2*pars$alpha)))
      }
      if(sum(Z)!=0){
        densregion(x,y,Z,colmax=makeTransparent(cols[i],alpha),colmin="transparent")
      }
      lines(c(OA[i,1],OA[i,2]),rep(pars$optima[i],2),col=cols[i],lwd=2)
    }
  }
}

#' Plot a pheongram with the posterior density for optima values
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
#' @export
#' @method plot bayouMCMC
plot.bayouMCMC <- function(chain, ...){
  if(is.null(attributes(chain)$burnin)){
    start <- 1
  } else {
    start <- round(attributes(chain)$burnin*length(chain$gen),0)
  }
  postburn <- start:length(chain$gen)
  chain2 <- lapply(chain,function(x) x[postburn])
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