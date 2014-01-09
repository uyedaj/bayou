#' bayOU internal function. 
#' 
#' \code{.repars} is an internal function and not generally called by the user
#' 
#' This is an internal function borrowed from geiger.
.repars <- function (pars, expected) 
{
  if (!length(pars) == length(expected)) 
    stop(paste("The following 'pars' are expected:\n\t", 
               paste(expected, collapse = "\n\t", sep = ""), sep = ""))
  if (all(!is.null(nm <- names(pars)))) {
    if (!all(nm %in% expected)) 
      stop(paste("The following 'pars' are unexpected:\n\t", 
                 paste(nm[!nm %in% expected], collapse = "\n\t", 
                       sep = ""), sep = ""))
    if (length(unique(nm)) != length(expected)) 
      stop(paste("The following 'pars' are expected:\n\t", 
                 paste(expected, collapse = "\n\t", sep = ""), 
                 sep = ""))
    mm = match(expected, nm)
    return(pars[mm])
  }
  else {
    return(pars)
  }
}
#' bayOU internal function. 
#' 
#' \code{.set.defaults} is an internal function and not generally called by the user
#' 
#' This is an internal function borrowed from diversitree.
.set.defaults <- function (f, ..., defaults = NULL) {
  dots <- match.call(expand.dots = FALSE)[["..."]]  
  if (missing(defaults)) 
    defaults <- dots  
  else if (is.list(defaults)) 
    defaults <- c(dots, defaults)
  else stop("'defaults' must be a list")
  if (is.null(defaults)) 
    return(f)
  if (!all(names(defaults) %in% names(formals(f)))) 
    stop("Unknown defaults")
  att <- attributes(f)
  formals(f)[names(defaults)] <- defaults
  attributes(f) <- att
  f
}

#' bayOU internal function. 
#' 
#' \code{.prepare.ou.univariate} is an internal function and not generally called by the user
#' 
#' This is an internal function modified from geiger's function .prepare.bm.univariate for use with OU models.
.prepare.ou.univariate <- function(tree,X, SE=0, ...){
  ntips <- length(tree$tip.label)
  rownames(tree$edge) <- 1:(length(tree$edge[,1]))
  cache <- .prepare.bm.univariate(tree, X, SE=SE)#, ...)
  ind <- as.numeric(rownames(cache$edge))
  cache$nH <- nodeHeights(tree)[ind,1]
  cache$maps <- tree$maps[ind]
  cache$mapped.edge <- tree$mapped.edge[ind,]
  cache$height <- max(nodeHeights(tree))
  cache$ntips <- length(X)
  cache$ind <- ind
  cache$ordering <- "postorder"
  cache$ht <- geiger:::.heights.cache(cache)
  cache$edge <- unname(cache$edge)
  plook <- function(x){mapply(paste,x[2:length(x)],x[1:(length(x)-1)],sep=",")}
  tB <- cache$desc$anc[1:ntips]
  tB <- mapply(c,1:ntips,tB, SIMPLIFY=FALSE)
  lookup <- lapply(tB,plook)
  edge.names <- mapply(paste,cache$edge[,1],cache$edge[,2],sep=",")
  cache$branchtrace <- t(sapply(lookup,function(x) as.numeric(edge.names %in% x)))
  cache$bdesc <- lapply(edge.names,function(branch) which(edge.names %in% unique(unlist(sapply(lookup,function(look) if(branch %in% look) look[1:which(branch==look)])))))
  cache$bdesc <- lapply(cache$bdesc,function(x) x[-length(x)])
  cache$lookup <- lookup
  rownames(cache$edge)=NULL
  return(cache)
}


#' bayOU internal function. 
#' 
#' \code{.prepare.bm.univariate} is an internal function and not generally called by the user
#' 
#' This is an internal function modified from geiger's function .prepare.bm.univariate for use with OU models.
.prepare.bm.univariate <- geiger:::.prepare.bm.univariate

#' bayOU internal function. 
#' 
#' \code{.sample} is an internal function and not generally called by the user
#' 
#' This is an internal function modified from the base function \code{sample()} \\
#' that provides consistent results with variable sample size.
.sample <- function (x, size, replace = FALSE, prob = NULL, ...) {
  x[sample(length(x), size, replace, prob, ...)]
  #if (missing(size)) 
  #  size <- length(x)
  #x[.Internal(sample(length(x), size, replace, prob))]
}

#' bayOU internal function. 
#' 
#' \code{.heights.cache} is an internal function and not generally called by the user
#' 
#' This is an internal function taken from geiger.
.heights.cache <- function (cache) {
  if (is.null(cache$ordering) || cache$ordering != "postorder") {
    stop("'cache' should be postordered")
  }
  n <- cache$n.tip
  n.node <- cache$n.node
  xx <- numeric(n + n.node)
  for (i in nrow(cache$edge):1) xx[cache$edge[i, 2]] <- xx[cache$edge[i, 
                                                     1]] + cache$edge.length[i]
  root = ifelse(is.null(cache$root.edge), 0, cache$root.edge)
  depth = max(xx)
  tt = depth - xx
  idx = 1:length(tt)
  dd = cache$edge.length[idx]
  mm = match(1:length(tt), c(cache$edge[, 2], n + 1))
  dd = c(cache$edge.length, root)[mm]
  ss = tt + dd
  res = cbind(ss, tt)
  rownames(res) = idx
  colnames(res) = c("start", "end")
  res = data.frame(res)
  res
}

print.priorFn <- function(x, ...){
  cat("prior function for bayou\n")
  cat(paste("expecting ", attributes(x)$model, " model\n", sep=""))
  cat("'pars' should be a list with named parameter values: list(", paste(gsub('^[a-zA-Z]',"",names(attributes(x)$param)),collapse=", "),")\n",sep="")
  cat("prior distribution functions used:\n")
  print(unlist(attributes(x)$dist))
  cat("\n")
  cat("definition:\n")
  attributes(x) <- NULL
  print(x)
}

print.refFn <- function(x, ...){
  cat("reference function for bayou\n")
  cat(paste("expecting ", attributes(x)$model, " model\n", sep=""))
  cat("'pars' should be a list with named parameter values: list(", paste(names(attributes(x)$dist),collapse=", "),")\n",sep="")
  cat("prior distribution functions used:\n")
  print(unlist(attributes(x)$dist))
  cat("\n")
  cat("definition:\n")
  attributes(x) <- NULL
  print(x)
}


.filldown.emap <- function(emap, tree){
  shifts <- emap$sh
  K <- sum(shifts)
  nopt <- rep(1,length(shifts)+1)
  opt <- 1
  for(i in length(shifts):1){
    if(shifts[i]==1){
      opt <- opt+1
      nopt[emap$e2[i]] <- opt
    } else {
      nopt[emap$e2[i]] <- nopt[emap$e1[i]]
    }
  }
  edge.map <- data.frame(tree$edge,nopt[tree$edge[,1]],nopt[tree$edge[,2]],shifts,tree$tip.label[tree$edge[,2]],emap$r1,emap$r2,tree$edge.length)
  names(edge.map) <- c("e1","e2","t1","t2","sh","tip","r1","r2","r")
  return(edge.map)
}

ouMatrix <- function(vcvMatrix, alpha)
{  vcvDiag<-diag(vcvMatrix)
   diagi<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag))
   diagj<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag), byrow=T)
   Tij = diagi + diagj - (2 * vcvMatrix)
   vcvRescaled = (1 / (2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcvMatrix))
   return(vcvRescaled)
}

identify.branches <- function(tree, n, fixed.loc=TRUE, plot.simmap=TRUE){
  mar.old <- par('mar')
  par(mfrow=c(1,1), mar=c(0.1,0.1,0.1,0.1))
  tree <- reorder(tree,"postorder")
  plot(tree, cex=0.5)
  L <- get("last_plot.phylo",envir=.PlotPhyloEnv)
  nH <- nodeHeights(tree)
  xx <- L$xx
  yy <- L$yy
  yH <- yy[tree$edge[,2]]
  sb <- numeric()
  loc.x <- numeric()
  for(i in 1:n){
    lx <- locator(1)    
    xmatch <- which(nH[,1] < lx$x & nH[,2] > lx$x)
    ymatch <- which(abs(yH[xmatch]-lx$y)==min(abs(yH[xmatch]-lx$y)))
    sb[i] <- xmatch[ymatch]
    loc.x[i] <- lx$x
    points(lx$x,yH[xmatch][ymatch],pch=21,bg="gray50")
    text(lx$x,yH[xmatch][ymatch],labels=sb[i],pos=4)
  }
  loc <- loc.x - nH[sb,1]
  if(plot.simmap){
    dum <- numeric(length(tree$tip.label))
    names(dum) <- tree$tip.label
    cache <- .prepare.ou.univariate(tree,dum)
    pars <- list(sb=sb, loc=loc, t2=1:n+1)
    tr <- .toSimmap(.pars2map(pars, cache),cache)
    cols <- rainbow(length(sb)+1)
    names(cols) <- 1:(length(sb)+1)
    plotSimmap(tr, pts=FALSE, fsize=0.5,colors=cols)
  }
  par(mar=mar.old)
  out <- list(sb=sb)
  if(fixed.loc) out$loc <- loc
  return(out)
}
