optima.ages <- function(pars,tree){
  reg <- sapply(tree$maps,function(x) names(x)[length(x)])
  adj <- sapply(tree$maps,function(x) ifelse(length(x)>1,x[1],0))
  abs.age <- nodeHeights(tree)[,1]+adj
  start <- tapply(abs.age,reg,min)
  end <- rep(max(abs.age),pars$ntheta)+1
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

plotSimmap.posterior <- function(chain,i,tree,dat,alpha=255){
  emap <- read.emap(chain$branch.shift[[i]],chain$location[[i]],chain$t2[[i]],tree)$emap
  pars <- list("alpha"=chain$alpha[[i]],"sig2"=chain$sig2[[i]],nb=chain$nb[[i]],optima=chain$optima[[i]],ntheta=chain$ntheta[[i]])
  tr <- emap2simmap(emap,tree)
  phenogram(tr,dat,colors=tr$col,ftype="off")
  regime.plot2(pars,tr,cols=tr$col,type="density",alpha=100)
  phenogram(tr,dat,colors=tr$col,add=TRUE,ftype="off") 
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
    pars$loc <- L$rel.location[pars$sb]
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

##Not working yet. Attempt to do ancestral state reconstruction under a multi-optimum OU model. 
#acemOU <- function(pars, tree, dat){
#  ntips <- length(tree$tip.label)
#  names(X) <- 1:(2*ntips-1)
#  cache <- .prepare.ou.univariate(tree, dat)
#  W <- .mOU.W(cache, pars)
#  E.X <- (W%*%pars$theta)[,1]
#  ouMatrix <- function(vcvMatrix, alpha)
##  {  vcvDiag<-diag(vcvMatrix)
#     diagi<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag))
#     diagj<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag), byrow=T)
#     Tij = diagi + diagj - (2 * vcvMatrix)
##     vcvRescaled = (1 / (2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcvMatrix))
##     return(vcvRescaled)
#  }
#  C <- vcvPhylo(cache$phy,anc.nodes=TRUE)
#  Sigma <- ouMatrix(C, pars$alpha)
#  likancFx <- function(px){
##    x <- c(dat, px)
#    res <- -dmnorm(x, mean=E.X, varcov=Sigma, log=TRUE)
#    return(res)
#  }
#  ip <- E.X[(ntips+1):(2*ntips-2)]
#  np <- length(ip)
#  x <- ip
#  x[419-226-1]<- ip[419-226-1]+0.1
#  likancFx(x)
#  res.opt <- nlminb(ip, likancFx, lower=rep(0,np), upper=rep(max(dat)*2,np))
#  res.opt <- nlm(likancFx, ip)
  
#}

##Return a weight matrix for all tips and internal nodes.
.mOU.W <- function(cache, pars){
    ntips <- cache$ntips
    plook <- function(x){mapply(paste,x[2:length(x)],x[1:(length(x)-1)],sep=",")}
    tB <- cache$desc$anc
    tB <- mapply(c,1:(2*ntips-1),tB)
    lookup <- lapply(tB,plook)
    edge.names <- mapply(paste,cache$edge[,1],cache$edge[,2],sep=",")
    branchtrace <- t(sapply(lookup,function(x) as.numeric(edge.names %in% x)))
    nH <- nodeHeights(cache$phy)
    nbranch <- length(cache$edge.length)
    #create a vector called shifts that indicates the number of shifts for each branch
    nshifts <- table(pars$sb)
    shifts <- rep(0,nbranch)
    shifts[as.numeric(attributes(nshifts)$dimnames[[1]])]<- nshifts
    #Create an index equal to the number of segments that identifies the branch on which each segment is found
    irow <- rep(1:nbranch,shifts+1)
    #For now, starting height is just the height of the node
    csbase <- cache$nH[irow]
    #Calculate the ending height by sorting the edge.length and the location of shifts by their branch identity and location
    csadd <- c(cache$edge.length, pars$loc)
    tmp.o <- c(1:nbranch, pars$sb)
    names(csadd) <- tmp.o
    add.o <- order(tmp.o,csadd)
    csadd <- csadd[add.o]
    #Ending height of the segment
    csmaps <- csadd + csbase
    #We need to know what the ending theta is for each segment, so we sort pars$t2 as we did for pars$loc, but +1 because t2 is the ending regime
    t2index <- add.o[which(add.o > nbranch)]
    t2b <- c(rep(1,length(csmaps)))
    t2b[match(t2index,add.o)+1] <- pars$t2[t2index-nbranch]
    #Now we need to cascade these regime down the tree. We won't need to cascade sandwiches, as they are trapped on the branch they occur. So we find them below:
    loc.o <- order(pars$loc,decreasing=TRUE)
    sandwiches <- tmp.o[duplicated(pars$sb[loc.o])]
    # And remove them:
    if(length(sandwiches)>0){
      sb.down <- pars$sb[-sandwiches]
      t2.down <- pars$t2[-sandwiches]
    } else {sb.down <- pars$sb; t2.down <- pars$t2}
    #Now we order the sb's and t2's to prepare for a postorder tree traversal
    sb.o <- order(sb.down)
    sb.down <- sb.down[sb.o]
    t2.down <- t2.down[sb.o]
    sb.desc <- cache$bdesc[sb.down]
    #Loop traveling down the tree, saving all descendents that are from that shift into the vector censored. These branches cannot be modified by shifts further down the tree.
    censored <- NULL
    name.o <- names(csmaps)
    names(t2b) <- name.o
    for(i in 1:length(sb.desc)){
      sb.desc[[i]] <- sb.desc[[i]][!(sb.desc[[i]] %in% censored)]
      censored <- c(censored, sb.desc[[i]])
      t2b[name.o[name.o %in% sb.desc[[i]]]] <- t2.down[i]
    }
    names(csmaps) <- t2b
    multips <- which(irow[2:length(irow)]==irow[1:(length(irow)-1)])
    #Set segments with more than one shift per branch to start at end of last shift
    csbase[multips+1] <- csmaps[multips]
    mOU.W.nH <- function(nh, cache, pars, csbase, csmaps, irow, nbranch){
      #Exponential term 1
      oW <- pars$alpha*(csbase-nh)
      #Exponential term 2
      nW <- (csmaps-csbase)*pars$alpha
    #If value of expnential term is too large (resulting in overflow), then use approximation
      if(any(nW>500)){
        tmp <- ifelse(nW>500, exp(nW+oW), exp(oW)*(exp(nW)-1))
      } else {
        tmp <- exp(oW)*(exp(nW)-1)
      }
    #Set up branch weight matrix
      bW <- matrix(0,nrow=nbranch,ncol=pars$ntheta)
    #Set up index over matrix, so that values go to right row index and column, based on the name of the segment in maps
      index <- irow + (as.integer(names(tmp))-1)*nbranch
      if(any(duplicated(index))){
        tmp <- tapply(tmp,index,sum)
        bW[as.numeric(names(tmp))] <- tmp
      } else {bW[index] <- tmp}
      bW[449:450,1] <- bW[449:450,1]+exp(-nh*pars$alpha)
      return(bW)
    }
    W.nn <- NULL
    for(i in (ntips+2):(2*ntips-1)){
      W.nn <- rbind(W.nn,branchtrace[i,]%*%mOU.W.nH(nH[which(cache$edge[,2]==i),2], cache, pars, csbase, csmaps, irow, nbranch))
    }
    W <- .parmap.W(cache, pars)
    W.all <- rbind(W, W.nn)
    return(W.all)
}

#' S3 method for plotting bayouMCMC objects
#' 
#' @export
#' @method plot bayouMCMC
plot.bayouMCMC <- function(chain, burnin=0, ...){
  postburn <- round(length(chain$gen)*burnin,0):length(chain$gen)
  chain2 <- lapply(chain,function(x) x[postburn])
  chain.length <- length(chain2$gen)
  univariates <- chain2[sapply(chain2,function(x) length(unlist(x)))==length(chain2$gen)]
  univariates$root <- sapply(chain2$theta, function(x) x[1])
  uni.df <- as.data.frame(univariates)
  rownames(uni.df) <- uni.df[,1]
  uni.df <- uni.df[,-1]
  plot(mcmc(uni.df), ...)
}


#' S3 method for plotting bayoupars objects
#' 
#' @export
#' @method plot bayoupars
plot.bayoupars <- function(pars, tree,...){
  tree <- reorder(tree, 'postorder')
  X <- rep(0, length(tree$tip.label))
  names(X) <- tree$tip.label
  cache <- .prepare.ou.univariate(tree, X)
  tr <- .toSimmap(.pars2map(pars, cache),cache)
  plotSimmap(tr,...)
}