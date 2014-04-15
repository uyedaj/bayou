emap2simmap <- function(emap,tree){
  foo <- function(x){
    tmp <- unlist(emap[x,c('r1','r2')])
    names(tmp) <- c(emap$t1[x],emap$t2[x])
    tmp[tmp>0]
  }
  nb <- sum(emap$sh)
  if(nb>0){
    col <- c("#000000",rainbow(nb))
  } else {col <- 1}
  names(col) <- 1:(nb+1)
  tree$maps <- lapply(1:dim(emap)[1],foo)
  tree$col <- col
  tree
}


#' Convert a bayou parameter list into a simmap formatted phylogeny
#' 
#' @export
pars2simmap <- function(pars,tree,theta=NULL,root.theta=0){
  sb <- pars$sb
  loc <- pars$loc
  t2 <- pars$t2
  Th <- theta
  nbranch <- length(tree$edge.length)
  maps <- lapply(tree$edge.length,function(x){y <- x; names(y) <- 1; y})
  dup <- which(duplicated(sb))
  if(length(dup)>0){
    maps[sb[-dup]] <- lapply(1:length(sb[-dup]),addshift2map,maps=maps,sb=sb[-dup],loc=loc[-dup],t2=t2[-dup])
  } else {
    maps[sb] <- lapply(1:length(sb),addshift2map,maps=maps,sb=sb,loc=loc,t2=t2)
  }
  for(i in dup){
    maps[[sb[i]]] <-addshift2map(i,maps=maps,sb=sb,loc=loc,t2=t2)
  }
  nopt <- rep(1,nbranch)
  for(i in nbranch:1){
    if(i %in% sb){
      opt <- as.integer(names(maps[[i]])[length(maps[[i]])])
      nopt[tree$edge[i,2]] <- opt
      names(maps[[i]])[1] <- nopt[tree$edge[i,1]]
    } else {
      names(maps[[i]])[1] <- nopt[tree$edge[i,1]] 
      nopt[tree$edge[i,2]] <- nopt[tree$edge[i,1]]
    }
  }
  shiftdown <- nopt[tree$edge[,1]]
  new.maps <- lapply(1:nbranch,function(x){names(maps[[x]])[1] <- shiftdown[x]; maps[[x]]})
  new.maps <- maps
  for(j in 1:nbranch){
    names(new.maps[[j]])[1] <-shiftdown[j]
  }
  anc.theta <- unlist(lapply(new.maps[sb],function(x) as.integer(names(x)[length(x)-1])),F,F)
  o <- rev(order(sb,loc*-1))
  shifted.maps <- new.maps[sb[o]]
  t1 <- rep(NA,length(t2))
  for(i in 1:length(t2)){
    nm <- as.integer(names(maps[[sb[o][i]]]))
    t1[nm[2:length(nm)]-1] <- nm[1:(length(nm)-1)]
    Th[t2[o[i]]] <- Th[t1[o[i]]]
  }
  new.tree <- tree
  new.tree$maps <- new.maps
  new.pars <- pars
  col <- c(1,rainbow(pars$k))
  names(col) <- 1:(pars$k+1)
  return(list(tree=new.tree,pars=new.pars,col=col))
}

.pars2map <- function(pars, cache){
  nbranch <- length(cache$edge.length)
  nshifts <- table(pars$sb)
  shifts <- rep(0,nbranch)
  shifts[as.numeric(attributes(nshifts)$dimnames[[1]])]<- nshifts
  irow <- rep(1:nbranch,shifts+1)
  segs <- c(cache$edge.length, pars$loc)
  tmp.o <- c(1:nbranch, pars$sb)
  names(segs) <- tmp.o
  add.o <- order(tmp.o,segs)
  segs <- segs[add.o]
  ind <- names(segs)
  t2index <- add.o[which(add.o > nbranch)]
  t2b <- c(rep(1,length(segs)))
  t2b[match(t2index,add.o)+1] <- pars$t2[t2index-nbranch]
  loc.o <- order(pars$loc,decreasing=TRUE)
  sandwiches <- loc.o[duplicated(pars$sb[loc.o])]
  if(length(sandwiches)>0){
    sb.down <- pars$sb[-sandwiches]
    t2.down <- pars$t2[-sandwiches]
  } else {sb.down <- pars$sb; t2.down <- pars$t2}
  sb.o <- order(sb.down)
  sb.down <- sb.down[sb.o]
  t2.down <- t2.down[sb.o]
  sb.desc <- cache$bdesc[sb.down]
  desc.length <- unlist(lapply(sb.desc, length),F,F)
  sb.desc <- sb.desc[desc.length>0]
  names(t2b) <- names(segs)
  sb.desc2 <- unlist(sb.desc,F,F)
  sb.dup <- duplicated(sb.desc2)
  sb.desc3 <- sb.desc2[!sb.dup]
  t2.names <- rep(t2.down[desc.length>0], unlist(lapply(sb.desc,length),F,F))
  t2.names <- t2.names[!sb.dup]
  t2b[as.character(unlist(sb.desc3,F,F))] <- t2.names
  base <- duplicated(names(segs))*c(0,segs[1:(length(segs)-1)])
  segs <- segs-base
  #maps <- lapply(1:nbranch, function(x) segs[ind==x])
  #maps <- lapply(maps, function(x) if(length(x) >1) {c(x[1],diff(x[1:length(x)]))} else x)
  return(list(segs=segs,theta=t2b))
}

QG.alpha <- function(pars){
  pars$h2*pars$P/(pars$P+pars$w2*pars$P)
}
QG.sig2 <- function(pars){
  (pars$h2*pars$P)/pars$Ne
}

OU.repar <- function(pars){
  alpha <- log(2)/pars$halflife
  sig2 <- (2*log(2)/(pars$halflife))*pars$Vy
  return(list(alpha=alpha,sig2=sig2))
}


#' Return an edge map for a given generation from a bayOU mcmc list
pull.emap <- function(i,chain,cache){
  sb <- chain$branch.shift[[i]]
  t2 <- chain$t2[[i]]
  sl <- chain$location[[i]]
  phy <- cache$phy
  shifts <- rep(0,length(cache$edge.length))
  shifts[sb] <- 1
  nopt <- rep(1,length(shifts)+1)
  opt <- 1
  j <- length(t2)
  phy$maps <- lapply(phy$edge.length,function(x){names(x) <- 1; x})
  if(length(sb)>0){
    phy$maps[sb] <- lapply(1:length(sb),function(x){y <- c(sl[x],phy$maps[[sb[x]]]-sl[x]);names(y)[2] <- t2[x] ;y})
    for(i in length(shifts):1){
      if(shifts[i]==1){
        opt <- t2[j]
        nopt[phy$edge[i,2]] <- opt
        j <- j-1
      } else {
        nopt[phy$edge[i,2]] <- nopt[phy$edge[i,1]]
      }
    }
  }
  phy$maps <- lapply(1:length(shifts),function(x){ names(phy$maps[[x]])[1] <- nopt[phy$edge[x,1]];phy$maps[[x]] })
  phy$maps <- lapply(1:length(shifts),function(x){ names(phy$maps[[x]])[length(phy$maps[[x]])] <- nopt[phy$edge[x,2]];phy$maps[[x]] })
  S <- phy$edge.length
  S[sb] <- phy$edge.length[sb]-sl
  S2 <- rep(0,length(shifts))
  S2[sb] <- sl
  edge.map <- data.frame(phy$edge,nopt[phy$edge[,1]],nopt[phy$edge[,2]],shifts,phy$tip.label[phy$edge[,2]],S,S2,phy$edge.length)
  colnames(edge.map)= c("e1","e2","t1","t2","sh","tip","r1","r2","r")
  return(list(phy=phy,emap=edge.map))
}

#' Generate an edge map from a vector of branches, shift locations and optima assignments
read.emap <- function(sb,sl,t2,phy){
  shifts <- rep(0,length(phy$edge.length))
  shifts[sb] <- 1
  nopt <- rep(1,length(shifts))
  opt <- 1
  j <- 1
  phy$maps <- lapply(phy$edge.length,function(x){names(x) <- 1; x})
  if(length(sb)>0){
    phy$maps[sb] <- lapply(1:length(sb),function(x){y <- c(sl[x],phy$maps[[sb[x]]]-sl[x]);names(y)[2] <- x+1 ;y})
    for(i in length(shifts):1){
      if(shifts[i]==1){
        opt <- t2[j]
        nopt[phy$edge[i,2]] <- opt
        j <- j+1
      } else {
        nopt[phy$edge[i,2]] <- nopt[phy$edge[i,1]]
      }
    }
  }
  phy$maps <- lapply(1:length(shifts),function(x){ names(phy$maps[[x]])[1] <- nopt[phy$edge[x,1]];phy$maps[[x]] })
  phy$maps <- lapply(1:length(shifts),function(x){ names(phy$maps[[x]])[length(phy$maps[[x]])] <- nopt[phy$edge[x,2]];phy$maps[[x]] })
  S <- phy$edge.length
  S[sb] <- sl
  S2 <- rep(0,length(shifts))
  S2[sb] <- phy$edge.length[sb]-sl
  edge.map <- data.frame(phy$edge,nopt[phy$edge[,1]],nopt[phy$edge[,2]],shifts,phy$tip.label[phy$edge[,2]],S,S2,phy$edge.length)
  colnames(edge.map)= c("e1","e2","t1","t2","sh","tip","r1","r2","r")
  return(list(phy=phy,emap=edge.map))
}

.toSimmap <- function(map, cache){
  maps <- lapply(1:length(cache$edge.length), function(x){ y <- map$segs[names(map$segs)==x]; names(y) <- map$theta[names(map$theta)==x]; y })  
  tree <- cache$phy
  tree$maps <- maps
  return(tree)
}

#' Converts OUwie data into bayou format
#' 
#' \code{OUwie2bayou} calculates the probability density of a value 
#' 
#' @param tree A phylogenetic tree with states at internal nodes as node labels
#' @param trait A data frame in OUwie format
#' @export
OUwie2bayou <- function(tree, trait){
  tree <- reorder(tree, 'postorder')
  tip.states <- trait[,2]
  names(tip.states) <- trait[,1]
  states <- c(tip.states[tree$tip.label], tree$node.label)
  states <- unname(states)
  e1 <- states[tree$edge[,1]]
  e2 <- states[tree$edge[,2]]
  sb <- which(e1 != e2)
  loc <- 0.5*tree$edge.length[sb]
  t2 <- as.numeric(factor(e2[sb]))+1
  k <- length(sb)
  ntheta <- length(unique(t2))+1
  pars <- list(k=k, ntheta=ntheta, sb=sb, loc=loc, t2=t2)
  class(pars) <- c("bayoupars","list")
  return(pars)
}

#' Converts bayou data into OUwie format
#' 
#' \code{OUwie2bayou} calculates the probability density of a value 
#' 
#' @param pars A list with parameter values specifying \code{sb} = the branches with shifts,
#' \code{loc} = the location on branches where a shift occurs and \code{t2} = the optima to which
#' descendants of that shift inherit
#' @param tree A phylogenetic tree
#' @param dat A vector of tip states
#' @export
bayou2OUwie <- function(pars, tree, dat){
  if(is.null(names(dat))){
    warning("No labels on trait data, assuming the same order as the tip labels")
  } else {dat <- dat[tree$tip.label]}
  ntips <- length(tree$tip.label)
  cache <- .prepare.ou.univariate(tree, dat)
  tr <- .toSimmap(.pars2map(pars, cache),cache)
  tips <- which(tr$edge[,2] <= ntips)
  node.states <- sapply(tr$maps, function(x) names(x)[1])
  names(node.states) <- tr$edge[,1]
  node.states <- rev(node.states[unique(names(node.states))])
  tr$node.label <- as.numeric(node.states)
  tip.states <- sapply(tr$maps[tips], function(x) names(x)[length(x)])
  names(tip.states) <- tr$tip.label[tr$edge[tips,2]]
  tip.states <- as.numeric(tip.states[tr$tip.label])
  OUwie.dat <- data.frame("Genus_species"=tr$tip.label, "Reg"= tip.states, "X"= dat)
  rownames(OUwie.dat) <- NULL
  return(list(tree=tr, dat=OUwie.dat))
}
   
  
