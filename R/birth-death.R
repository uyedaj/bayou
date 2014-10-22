#pars <- new.pars
#splitNode <- pars$sb
#r <- pars$r
#eps <- pars$eps
#index <- c(0, pars$t2-1)
#phy <- cache$phy

getSplitLikelihood<-function(phy, splitNode, r, eps, index) {

	richness<-NULL
	richness <- geiger:::.check.richness(phy = phy, richness = richness);
	phyData <- geiger:::.treedata.medusa(phy = phy, richness = richness, warnings = FALSE); ## modified prune.tree.merge.data for multiple trees (jme)
	
	obj <- geiger:::.make.cache.medusa(phy = phy, richness = richness, fx = geiger:::.get.parallel(1), shiftCut = shiftCut);

	zSplit<-obj

	for(nextNode in splitNode) {
			zSplit<-.split.z.at.node.medusa(node = nextNode, z = zSplit$z, desc = list(stem = obj$desc.stem, node = obj$desc.node), shiftCut = shiftCut)
	}
	
	maxPartition<-max(zSplit$z[,"partition"])
	
	res<-0
	
	for(i in 1: maxPartition) {
		
		whichParam<-which(index==i-1)
		
		part<-zSplit$z[,"partition"]==i
		newLik<-.lik.partition.medusa(zSplit$z[part,], model="bd")
		res<-res+newLik(c(r[whichParam], eps[whichParam]))
	}
		
	res
}		

.split.z.at.node.medusa <- function (node, z, desc, shiftCut, extract = FALSE) {
  descendants <- NULL
  if (shiftCut == "stem") {
    descendants <- desc$stem
  } else {
    descendants <- desc$node
  }
  part <- z[, "partition"]
  base <- min(part[z[, 1] == node | z[, 2] == node])
  tag <- max(part) + 1
  i <- descendants[[node]]
  idx <- i[part[i] == base]
  z[idx, "partition"] <- tag
  if (extract) {
    z <- z[idx, , drop = FALSE]
  }
  return(list(z = z, affected = c(unique(part[idx]), tag)))
}
.lik.partition.medusa <- function (partition, model) {
  is.int <- is.na(partition[, "n.t"])
  is.pend <- !is.int
  n.int <- sum(is.int)
  n.pend <- sum(is.pend)
  if (n.int + n.pend != length(partition[, 1])) {
    stop("You messed up, yo.")
  }
  int <- partition[is.int, , drop = FALSE]
  pend <- partition[is.pend, , drop = FALSE]
  sum.int.t.len <- sum(int[, "t.len"])
  int.t.0 <- int[, "t.0"]
  pend.n.0 <- pend[, "n.0"]
  pend.n.t <- pend[, "n.t"]
  pend.t.len <- pend[, "t.len"]
  f <- function(pars) {
    if (model == "bd") {
      r <- pars[1]
      epsilon <- pars[2]
      if (r < 0 || epsilon < 0 || epsilon >= 1) {
        return(-Inf)
      }
    }
    else {
      r <- pars[1]
      epsilon <- 0
      if (r < 0) {
        return(-Inf)
      }
    }
    l.int <- numeric()
    l.pend <- numeric()
    if (n.int == 0) {
      l.int <- 0
    }
    else {
      l.int <- n.int * log(r) - r * sum.int.t.len - sum(log(1 - 
                                                              (epsilon * exp(-r * int.t.0))))
    }
    if (n.pend == 0) {
      l.pend <- 0
    }
    else {
      ert <- exp(r * pend.t.len)
      B <- (ert - 1)/(ert - epsilon)
      l.pend <- sum(log(1 - B) + (pend.n.t - 1) * log(B))
    }
    return(l.int + l.pend)
  }
}

medusa2bayou <- function(phy, splitNode, r, eps, index,...){
  phy <- reorder(phy, "postorder")  
  bid <- sapply(splitNode, function(x) which(phy$edge[,2]==x))
  pars <- list(r=r, eps=eps, sb = bid, t2=index[which(index!=0)]+1, loc=phy$edge.length[bid]-.Machine$double.eps*10)
  class(pars) <- "bayoupars"
  plotBayoupars(pars, phy, ...)
  return(pars)
}

