.sim.OU <- function(ns, totaltime, alpha, sig2, theta, x0){
  time.seq <- seq(0, totaltime, length.out=ns)
  time.step <- (totaltime - 0)/(ns-1)
  if(length(alpha)==1){
    alpha <- rep(alpha, ns)
  } 
  if(length(theta)==1){
    theta <- rep(theta, ns)
  }
  if(length(sig2)==1){
    sig2 <- rep(sig2, ns)
  }
  MM <- array(dim=ns)
  Xi <- x0
  for(i in 1:ns){
    Xi <- .ou.next(Xi, alpha[i], sig2[i], theta[i], time.step)
    MM[i] <- Xi
  }
  return(MM)
}

.simOU.onmap <- function(mapi, alpha, sig2, theta, ptsperunit, x0){
  regimes <- as.numeric(names(mapi))
  nopts <- round(mapi*ptsperunit,0)
  totalpts <- sum(nopts)
  theta <- unlist(lapply(1:length(mapi), function(x) rep(theta[regimes[x]], nopts[x])))
  .sim.OU(totalpts, sum(mapi), alpha, sig2, theta, x0)
}

.ou.next <- function(x, alpha, sig2, theta, tt){
  Exp <- x*exp(-alpha*tt) + theta*(1-exp(-alpha*tt))
  stats::rnorm(1, Exp, sqrt(sig2*tt))
}

#' A function to visualize a multi-optimum OU process evolving on a phylogeny
#' 
#' @param pars A bayou parameter list to simulate the OU process from
#' @param tree A phylogenetic tree
#' @param ptsperunit A number giving the number of points to simulate per unit time
#' @param pal A color palette function
#' @param aph The alpha value for transparency of the lines
#' @param lwd The width of the lines
#' 
#' @export
plotOUtreesim <- function(pars, tree, ptsperunit=100, pal=rainbow, aph=255, lwd=1){
  simmapTree <- pars2simmap(pars, tree)
  maps <- simmapTree$tree$maps
  alpha <- pars$alpha
  sig2 <- pars$sig2
  theta <- pars$theta
  x0 <- pars$theta[1]
  cols <- setNames(pal(pars$ntheta), 1:pars$ntheta)
  brcols <- cols[sapply(1:length(maps), function(x) names(utils::tail(maps[[x]], 1)))]
  nH <- nodeHeights(tree)
  history <- list()
  Xi <- x0
  for(i in length(maps):1){
    if(i < length(maps) -1) Xi <- history[[which(tree$edge[,2]==tree$edge[i,1])]][length(history[[which(tree$edge[,2]==tree$edge[i,1])]])]
    history[[i]] <- .simOU.onmap(maps[[i]], alpha, sig2, theta, ptsperunit, Xi)
  }
  plot(c(min(nH), max(nH)), c(min(unlist(history)), max(unlist(history))), type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
  invisible(lapply(1:length(history), function(x) lines(seq(nH[x,1], nH[x,2], length.out=length(history[[x]])), history[[x]], col=makeTransparent(brcols[[x]], alpha=aph),lwd=lwd)))
}