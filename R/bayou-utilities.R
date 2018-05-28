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
##Merging .prepare.ou.phylolm and .prepare.ou.univariate
##Clean this up to make sure I actually need all this!!
.prepare.ou.univariate <- function(tree,X, SE=0, pred=NULL, ...){
  ntips <- length(tree$tip.label)
  rownames(tree$edge) <- 1:(length(tree$edge[,1]))
  cache <- .prepare.bm.univariate(tree, X, SE=SE)#, ...)
  ind <- as.numeric(rownames(cache$edge))
  cache$n <- ntips
  cache$N <- nrow(cache$phy$edge)
  cache$nH <- phytools::nodeHeights(tree)[ind,1]
  if(is.null(pred)){
    pred <- cbind(rep(0, ntips))
    rownames(pred) <-  cache$tip.label
  }
  cache$maps <- tree$maps[ind]
  cache$mapped.edge <- tree$mapped.edge[ind,]
  cache$height <- max(phytools::nodeHeights(tree))
  cache$anc <- cache$phy$edge[,1]
  cache$des <- cache$phy$edge[,2]
  cache$ntips <- length(X)
  cache$ind <- ind
  cache$ordering <- "postorder"
  cache$ht <- .heights.cache(cache)
  cache$edge <- unname(cache$edge)
  o <- match(tree$tip.label, cache$tip.label)
  cache$pred <- cbind(pred[o,])
  colnames(cache$pred) <- colnames(pred)
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
  cache$distFromRoot <- .pruningwise.distFromRoot(cache$phy, n=cache$n, N=cache$N)
  cache$tipFromRoot <- cache$distFromRoot[1:cache$n]
  cache$ultrametric <- as.numeric(is.ultrametric(cache$phy))
  D <- max(cache$distFromRoot[1:cache$n]) - cache$distFromRoot[1:cache$n]
  cache$D <- D - mean(D)
  phy2 <- cache$phy
  for (i in 1:cache$n) {
    tmp = phy2$edge.length[which(cache$phy$edge[,2]==i)]
    phy2$edge.length[which(cache$phy$edge[,2]==i)] = tmp + cache$D[i]      
  }
  cache$times <- .pruningwise.branching.times(phy2, n=cache$n, des=phy2$edge[,2], anc=phy2$edge[,1])
  names(cache$times) <- (cache$n+1):(cache$n+cache$phy$Nnode) 
  cache$Tmax <- max(cache$times)
  cache$branches.anc <- lapply(1:cache$N, function(x) which(names(cache$times)==cache$phy$edge[,1][x]))
  cache$branches.des <- lapply(1:cache$N, function(x) which(names(cache$times)==cache$phy$edge[,2][x]))
  cache$branches.des2 <- cache$branches.des
  cache$branches.des2[which(sapply(cache$branches.des,length)==0)] <- -1
  cache$branches.des2 <- unlist(cache$branches.des2)
  cache$branches.anc2 <- unlist(cache$branches.anc)
  cache$externalEdge <- cache$phy$edge[,2]<=cache$n
  return(cache)
}

#' bayOU internal function. 
#' 
#' \code{.prepare.bm.univariate} is an internal function and not generally called by the user
#' 
#' This is an internal function modified from geiger's function .prepare.bm.univariate for use with OU models.
.prepare.bm.univariate <- function (phy, dat, nodes = NULL, SE = NA, control = list(binary = TRUE, ultrametric = FALSE)) {
  ct = list(binary = TRUE, ultrametric = FALSE)
  ct[names(control)] = control
  td = treedata(phy, dat, sort = TRUE, warnings = FALSE)
  phy = reorder(td$phy, "postorder")
  if (ct$binary) 
    if (!is.binary.tree(phy)) 
      stop("'phy' should be a binary tree")
  if (ct$ultrametric) 
    if (!is.ultrametric(phy)) 
      stop("'phy' should be an ultrametric tree")
  if (is.null(phy$edge.length)) 
    stop("'phy' must include branch lengths in units of time")
  if (ncol(td$data) > 1) 
    stop("'dat' should be univariate")
  dat = td$data[, 1]
  seTMP = structure(rep(NA, length(dat)), names = names(dat))
  if (is.null(SE)) 
    SE = NA
  if (length(SE) > 1) {
    if (is.null(names(SE))) 
      stop("'SE' should be a named vector")
    if (!all(names(dat) %in% names(SE))) 
      stop("names in 'SE' must all occur in names of 'dat'")
    seTMP[names(SE[names(dat)])] = SE[names(dat)]
    SE = seTMP
  }
  else {
    if (is.numeric(SE)) {
      seTMP[] = SE
      SE = seTMP
    }
    else {
      SE = seTMP
    }
  }
  if (!all(is.na(SE) | SE >= 0)) 
    stop("'SE' values should be positive (including 0) or NA")
  cache = .cache.tree(phy)
  N = cache$n.tip
  n = cache$n.node
  m <- s <- g <- numeric(N + n)
  g[1:N] = 1
  m[] = NA
  m[1:N] = dat
  s[1:N] = SE
  if (!is.null(nodes)) {
    nn = (N + 1):(N + n)
    vec = .cache.y.nodes(m, s, g, nn, phy, nodes = nodes)
  }
  else {
    vec = rbind(m = m, s = s)
    attr(vec, "given") = g
    attr(vec, "adjse") = as.numeric(is.na(s))[1:N]
  }
  cache$SE = SE
  cache$dat = dat[match(phy$tip.label, names(dat))]
  cache$phy = phy
  cache$y = vec
  return(cache)
}

#' Internal function from geiger
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

#' Internal function from geiger
.cache.tree <- function (phy) {
  ordxx = function(children, is.tip, root) {
    todo <- list(root)
    i <- root
    repeat {
      kids <- children[i, ]
      i <- kids[!is.tip[kids]]
      if (length(i) > 0) 
        todo <- c(todo, list(i))
      else break
    }
    as.vector(unlist(rev(todo)))
  }
  edge <- phy$edge
  edge.length <- phy$edge.length
  idx <- seq_len(max(edge))
  n.tip <- Ntip(phy)
  tips <- seq_len(n.tip)
  root <- n.tip + 1
  is.tip <- idx <= n.tip
  desc = .cache.descendants(phy)
  children <- desc$fdesc
  if (!max(sapply(children, length) == 2)) {
    children = NULL
    order = NULL
    binary = FALSE
  }
  else {
    children <- rbind(matrix(NA, n.tip, 2), t(matrix(unlist(children), 
                                                     nrow = 2)))
    order <- ordxx(children, is.tip, root)
    binary = TRUE
  }
  len <- edge.length[mm <- match(idx, edge[, 2])]
  ans <- list(tip.label = phy$tip.label, node.label = phy$node.label, 
              len = len, children = children, order = order, root = root, 
              n.tip = n.tip, n.node = phy$Nnode, tips = tips, edge = edge, 
              edge.length = edge.length, nodes = phy$edge[, 2], binary = binary, 
              desc = desc)
  ans
}

#' Internal function from geiger
.cache.descendants <- function (phy) {
  N = as.integer(Ntip(phy))
  n = as.integer(Nnode(phy))
  phy = reorder(phy, "postorder")
  zz = list(N = N, MAXNODE = N + n, ANC = as.integer(phy$edge[, 
                                                              1]), DES = as.integer(phy$edge[, 2]))
  res = .Call("cache_descendants", phy = zz, package = "geiger")
  return(res)
}

#' Internal function from geiger
.cache.y.nodes <- function (m, s, g, nn, phy, nodes) {
  if (is.numeric(nodes) & is.vector(nodes)) {
    if (!all(names(nodes) %in% nn)) 
      stop("'nodes' must have (integer) names corresponding to the internal nodes of 'phy'")
    nodes = data.frame(cbind(node = as.integer(names(nodes)), 
                             mean = nodes, SE = 0), stringsAsFactors = FALSE)
  }
  else {
    if (!all(c("taxon1", "taxon2", "mean", "SE") %in% colnames(nodes))) {
      flag = FALSE
      if (!all(c("mean", "SE") %in% colnames(nodes)) | 
            is.null(rownames(nodes))) {
        flag = TRUE
      }
      else if (!all(rr <- as.integer(rownames(nodes)) %in% 
                      nn)) {
        flag = TRUE
      }
      if (flag) 
        stop("'nodes' must minimally have column names: 'taxon1', 'taxon2', 'mean', and 'SE'")
      nodes = as.data.frame(nodes)
      nodes$node = as.integer(rownames(nodes))
    }
    else {
      nodes = as.data.frame(nodes)
      if (!is.numeric(nodes$mean) | !is.numeric(nodes$SE)) {
        stop("'nodes' must have numeric vectors for 'mean' and 'SE'")
      }
      if (!all(zz <- unique(c(as.character(nodes$taxon1), 
                              as.character(nodes$taxon2))) %in% phy$tip.label)) {
        stop(paste("Some taxa appear missing from 'phy':\n\t", 
                   paste(zz[!zz %in% phy$tip.label], collapse = "\n\t", 
                         sep = ""), sep = ""))
      }
      nodes$node = apply(nodes[, c("taxon1", "taxon2")], 
                         1, .mrca, phy = phy)
    }
    if (!length(unique(nodes$node)) == nrow(nodes)) {
      stop("Some nodes multiply constrained:\n\t", paste(nodes$node[duplicated(nodes$node)], 
                                                         collapse = "\n\t", sep = ""), sep = "")
    }
  }
  nidx = nodes$node
  if (any(zz <- g[nidx] == 1)) 
    stop("Some nodes already constrained:\n\t", paste(nidx[which(zz)], 
                                                      collapse = "\n\t", sep = ""), sep = "")
  m[nidx] = as.numeric(nodes$mean)
  s[nidx] = as.numeric(nodes$SE)
  g[nidx] = 1
  vec = rbind(m = m, s = s)
  attr(vec, "given") = g
  attr(vec, "adjse") = as.numeric(is.na(s))
  vec
}

#' Internal function from geiger
.mrca <- function (labels, phy)  {
mm = labels
if (all(is.character(labels))) {
  ll = c(phy$tip.label, phy$node.label)
  mm = match(labels, ll)
  if (any(is.na(mm))) 
    stop("Some 'labels' not encountered in 'phy'")
}
if (!all(is.numeric(mm))) 
  stop("Supply 'labels' as a character or integer vector")
if (length(u <- unique(mm)) == 1) 
  return(u)
aa = unlist(lapply(mm, function(x) .get.ancestors.of.node(x, 
                                                          phy)))
tt = table(aa)
max(as.integer(names(tt[tt == length(labels)])))
}

.get.ancestor.of.node <- function (node, phy) {
  return(phy$edge[which(phy$edge[, 2] == node), 1])
}

.get.ancestors.of.node <- function (node, phy) {
  a = c()
  if (node == (root <- Ntip(phy) + 1)) 
    return(NULL)
  f = .get.ancestor.of.node(node, phy)
  a = c(a, f)
  if (f > root) 
    a = c(a, .get.ancestors.of.node(f, phy))
  return(a)
}


#geiger:::.prepare.bm.univariate

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

#' S3 method for printing priorFn objects
#' 
#' @param x A function of class 'priorFn' produced by \code{make.prior}
#' @param ... Additional arguments passed to \code{print}
#' 
#' @export
#' @method print priorFn
print.priorFn <- function(x, ...){
  cat("prior function for bayou\n")
  cat(paste("expecting ", attributes(x)$model, " model\n", sep=""))
  cat("'pars' should be a list with named parameter values: list(", paste(gsub('^[a-zA-Z]',"",names(attributes(x)$param)),collapse=", "),")\n",sep="")
  cat("prior distribution functions used:\n")
  print(unlist(attributes(x)$dist))
  cat("\n")
  cat("definition:\n")
  attributes(x) <- NULL
  print(x, ...)
}

#' S3 method for printing refFn objects
#' 
#' @param x A function of class 'refFn' produced by make.refFn
#' @param ... Additional arguments passed to \code{print}
#' 
#' @export
#' @method print refFn
print.refFn <- function(x, ...){
  cat("reference function for bayou\n")
  cat(paste("expecting ", attributes(x)$model, " model\n", sep=""))
  cat("'pars' should be a list with named parameter values: list(", paste(names(attributes(x)$dist),collapse=", "),")\n",sep="")
  cat("prior distribution functions used:\n")
  print(unlist(attributes(x)$dist))
  cat("\n")
  cat("definition:\n")
  attributes(x) <- NULL
  print(x, ...)
}

#' bayOU internal function. 
#' 
#' \code{.heights.cache} is an internal function and not generally called by the user
#' 
#' This function calculates the change in optima at each shift point on the tree
.D.from.theta <- function(pars, cache, sort=TRUE){
  map <- .pars2map(pars, cache)
  reg.shifts <- cbind(map$theta[which(duplicated(names(map$theta)))-1], map$theta[which(duplicated(names(map$theta)))])
  D <- pars$theta[reg.shifts[,2]] - pars$theta[reg.shifts[,1]]
  D <- D[order(reg.shifts[,2])]
  return(list(D=D, map=reg.shifts))
}

# Transforms a BM variance-covariance matrix into an OU VCV with a alpha parameter
.ouMatrix <- function(vcvMatrix, alpha)
{  vcvDiag<-diag(vcvMatrix)
   diagi<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag))
   diagj<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag), byrow=T)
   Tij = diagi + diagj - (2 * vcvMatrix)
   vcvRescaled = (1 / (2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcvMatrix))
   return(vcvRescaled)
}

#' Identify shifts on branches of a phylogenetic tree
#' 
#' \code{identifyBranches} opens an interactive phylogeny plot that allows the user to specify the location
#' of shifts in a phylogenetic tree.
#' @param tree An object of class 'phylo'
#' @param n The number of shifts to map interactively onto the phylogeny
#' @param fixed.loc A logical indicating whether the exact location on the branch should be returned, or the shift will be free to move along the branch
#' @param plot.simmap A logical indicating whether the resulting painting of regimes should be plotted following the selection shift location.
#' 
#' @description This is a convenience function for mapping regimes interactively on the phylogeny. The method locates the nearest branch to where the
#' cursor is clicked on the plot and records the branch number and the location selected on the branch.
#' 
#' @return Returns a list with elements "sb" which contains the branch numbers of all selected branches with length "n". If "fixed.loc=TRUE", then the list also
#' contains a vector "loc" which contains the location of the selected shifts along the branch. 
#' 
#' @export
identifyBranches <- function(tree, n, fixed.loc=TRUE, plot.simmap=TRUE){
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

#' Internal function taken from phytools
.whichorder <- function (x, y){ 
  sapply(x, function(x, y) which(x == y), y = y)
}


#' Internal function to determine tip regimes:
.tipregime <- function(pars, tree){
  ntips <- length(tree$tip.label)
  tree <- reorder(tree, "postorder")
  tr <- pars2simmap(pars, tree)
  tip.reg <- as.numeric(sapply(tr$tree$map, function(x) names(x)[length(x)]))
  tip.reg <- tip.reg[which(tr$tree$edge[,2] <= ntips)]
  o <- order(tr$tree$edge[,2][which(tr$tree$edge[,2] <= ntips)])
  tip.reg <- tip.reg[o]
  tip.reg <- setNames(tip.reg, tr$tree$tip.label)
  return(tip.reg)
}


#' Function for checking parameter lists, prior and models are consistent and error-free
#' 
#' @param pars A list of parameters that will be specified as starting parameter
#' @param tree An object of class ``phylo''
#' @param  dat A named data vector that matches the tip lables in the provided tree
#' @param pred A matrix or data frame with named columns with predictor data represented 
#' in the specified formula
#' @param SE The standard error of the data. Either a single value applied to all the data, 
#' or a vector of length(dat).
#' @param prior A prior function made using make.prior
#' @param model Either one of c("OU", "QG" or "OUrepar") or a list specifying the model
#' to be used.
#' @param autofix A logical that indicates whether certain errors should be automatically fixed.
#' 
#' @details A series of checks are performed, run internally within bayou.makeMCMC, but can also
#' be run on provided inputs prior to this. Errors are reported. 
#' 
#' If autofix == TRUE, then the following errors will be automatically corrected:
#' 
#' Branch lengths == 0; any branches of length 0 will be given length .Machine$double.eps
#' is.binary(tree) == FALSE; runs multi2di
#' pars do not match prior$fixed; parameters are resimulated from prior
#' 
#' @return A list of results of the checks and if 'autofix==TRUE', then ..$autofixed returns a
#' list of all the input elements, with corrections.
#' 
#' @export  
bayou.checkModel <- function(pars=NULL, tree, dat, pred=NULL, SE=0, 
                             prior, model="OU", autofix=TRUE){
  checks <- list()
  
  ## Step 1. Check to see if the model already exists in bayou or is custom.
  if(class(model)=="character"){
    if(model %in% c("OU", "QG", "OUrepar")){
      model.name <- model
      model <- switch(model, "OU"=model.OU, "QG"=model.QG, "OUrepar"=model.OUrepar)
    } else {stop("Specified model not recognized by bayou. Choose 'OU', 'QG', 'OUrepar'
                 or provide a bayou model object (list) that fully specifies the 
                 model to be run.")}
  } else {
    ## If custom, perform checks of each of the expected elements in the model specification list
    model.name <- "custom"
    expected_model_elements <- c("moves", "control.weights", "D", "parorder", "fixedpars", "rjpars",
                        "shiftpars", "call" ,"expFn","prevalues" , "monitor.fn", "lik.fn")
    checks$model.expected_elements <- assertthat::validate_that(all(names(model) %in% expected_model_elements))
    checks$model.moves_type <- assertthat::validate_that(class(model$moves)=="list")
    moves <- sapply(model$moves, function(x) ifelse(x=="fixed", "function", class(get(x))))
    checks$model.moves_available <- assertthat::validate_that(all(moves == "function"))
    checks$model.weights_list <- assertthat::validate_that(class(model$control.weights)=="list")
    control_weights <- sapply(model$control.weights, function(x) class(x))
    checks$model.weights_type <- assertthat::validate_that(all(control_weights=="numeric"))
    checks$model.weights_moves <- assertthat::validate_that(length(model$control.weights[model$control.weights>0])==length(model$moves))
    checks$model.weights_D <- assertthat::validate_that(length(model$control.weights)==length(model$D))
    checks$model.D_list <- assertthat::validate_that(class(model$D)=="list")
    D.type <- lapply(model$D, function(x) class(x))
    checks$model.D_type <- assertthat::validate_that(all(D.type == "numeric"))
    checks$model.parorder_type <- assertthat::validate_that(class(model$parorder)=="character")
    checks$model.parorder_length <- assertthat::validate_that(length(model$parorder) >= length(model$moves))
    if(!is.null(model$rjpars)) checks$model.rjpars_type <- assertthat::validate_that(class(model$rjpars)=="character")
    checks$model.shiftpars_type <- assertthat::validate_that(identical(model$shiftpars,c('sb', 'loc', 't2')))
    if(!is.null(model$monitor.fn)){
      checks$model.monitor_type <- assertthat::validate_that(class(model$monitor.fn)=="function")
      monitor_argnames <- names(formals(model$monitor.fn))
      checks$model.monitor_args <- assertthat::validate_that(identical(monitor_argnames,
                                                  c("i", "lik", "pr", "pars", "accept", "accept.type", "j")))
    }
    checks$model.likfn_type <- assertthat::validate_that(class(model$lik.fn)=="function")
    likfn_argnames <- names(formals(model$lik.fn))
    checks$model.likfn_args <- assertthat::validate_that(identical(likfn_argnames,
                                                       c("pars", "cache", "X", "model")))
    
  }
  ## Step 2. Check pars if provided.
  expected_pars <- c(model$parorder, model$shiftpars)
  if(is.null(pars)){
    #checks$pars <- "No starting parameters provided, skipping parameter checks"
  } else {
    
    rjpars <- model$rjpars
    checks$pars.expected <- assertthat::validate_that(all(names(pars) %in% expected_pars))
    checks$pars.k_sb <- assertthat::validate_that(pars$k == length(pars$sb))
    checks$pars.k_loc <- assertthat::validate_that(pars$k == length(pars$loc))
    checks$pars.k_t2 <- assertthat::validate_that(pars$k == length(pars$t2))
    checks$pars.k_ntheta <- assertthat::validate_that(pars$k >= pars$ntheta-1)
    checks$pars.ntheta_theta <- assertthat::validate_that(pars$ntheta == length(pars$theta))
    if(length(rjpars)>0) checks$pars.rjnumber <- assertthat::validate_that(all(sapply(pars[rjpars], length) == pars$ntheta))
  }
  
  ## Step 3. Check tree
    checks$tree_binary <- assertthat::validate_that(is.binary(tree)==TRUE)
    if(checks$tree_binary != TRUE & autofix==TRUE){
      tree <- multi2di(tree)
      checks$tree_binary <- paste(checks$tree_binary, "autofixed", sep="; ")
    } 
    checks$tree_branchLengths <- assertthat::validate_that(all(tree$edge.length > 0))
  if(checks$tree_branchLengths != TRUE & autofix==TRUE){
    tree$edge.length[tree$edge.length==0] <- .Machine$double.eps
    checks$tree_branchLengths <- paste(checks$tree_branchLengths, "autofixed", sep="; ")
  } 

  ## Step 4. Check data
      if(is.null(dat)){
        checks$dat <- "No data provided, creating dummy dataset for further checks"
        dat <- setNames(rep(0, length(tree$tip.label)), tree$tip.label)
      }
      checks$dat_type <- assertthat::validate_that(class(dat)=="numeric")
      checks$dat_names <- assertthat::validate_that(all(names(dat) %in% tree$tip.label))
  
  ## Step 5. Check to see if cache object can be made
      cache <- try(.prepare.ou.univariate(tree=tree, X=dat, SE=SE, pred=pred))
      checks$cache_creation <- assertthat::validate_that(class(cache)=="list")
      
  ## Step 6. Check prior function
      priorAttributes <- attributes(prior)
      checks$prior.can_simulate <- assertthat::validate_that(class(try(priorSim(prior, tree, plot=FALSE)))=="list")
      if(is.null(pars)){
        pars <- priorSim(prior,tree,plot=FALSE)$pars[[1]]
      } else {
        if(length(priorAttributes$fixed)>0){
          checks$prior.pars_fixed <- assertthat::validate_that(identical(priorAttributes$fixed, pars[names(priorAttributes$fixed)]))
          if(checks$prior.pars_fixed != TRUE & autofix==TRUE){
            checks$prior.pars_fixed <- paste0(checks$prior.pars_fixed, 
                  "; starting parameters resimulated from the prior")
            pars <- priorSim(prior, tree, plot=FALSE)$pars[[1]]
          }
        }
      }
      checks$prior.distFx_available <- assertthat::validate_that(all(sapply(priorAttributes$functions, class)=="function"))
      if("dk" %in% names(priorAttributes$distributions)){
        if(priorAttributes$distributions$dk=="fixed"){
          checks$prior.k_fixed <- assertthat::validate_that(all(c("t2", "sb", "k") %in% names(priorAttributes$fixed)))
        } else {
          checks$prior.k_t2 <- assertthat::validate_that(sum(duplicated(pars$t2))==0, msg = "Cannot have convergent regimes with reversible-jump model")
        }
      } 
      checks$prior.prior_finite <- assertthat::validate_that(is.finite(try(prior(pars))))
      priorNames <- c(priorAttributes$parnames, names(priorAttributes$fixed))
      checks$prior.model_prior <- assertthat::validate_that(all(priorNames %in% expected_pars))
      
  ## Step 7. Check likelihood function
      checks$likfn_finite <- assertthat::validate_that(is.finite(try(model$lik.fn(pars, cache, dat, model=model.name)$loglik)))
  
  cat("Checking inputs for errors:")
  sapply(checks, function(x) if(x==TRUE) cat(".") else cat("!"))
  cat("\n")
  warns <- which(sapply(checks, function(x) x!=TRUE))
  if(length(warns)>0) {
    cat("The following errors were found, proceed with caution\n")
    dum <- lapply(1:length(warns), function(x) cat(paste(x,". ", checks[[warns[x]]],"\n", sep="")))
  }
  checks$autofixed <- list(pars=pars, tree=tree, dat=dat, pred=pred, SE=SE, prior=prior, model=model)
  checks$autofixed$cache <- cache
  return(checks)
  
  
}







