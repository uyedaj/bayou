#' bayOU internal function. 
#' 
#' \code{.fix.root.bm} is an internal function and not generally called by the user
#' 
#' This is an internal function from geiger.
.fix.root.bm <- geiger:::.fix.root.bm

#' bayOU internal function. 
#' 
#' \code{.ou.cache} is an internal function and not generally called by the user
#' 
#' This is an internal function that modifies the internal function \code{.ou.cache} in geiger for efficiency.
.ou.cache.fast <- function (cache) 
{
  ht = cache$ht
  N = cache$n.tip
  Tmax = ht$start[N + 1]
  mm = match(1:nrow(ht), cache$edge[, 2])
  ht$t1 = Tmax - ht$end[cache$edge[mm, 1]]
  ht$t2 = ht$start - ht$end + ht$t1
  z = function(alpha) {
    if (alpha < 0) 
      stop("'alpha' must be positive valued")
    if (alpha == 0){
      bl = ht$t2-ht$t1
    } else {
      bl = (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - ht$t2)) * 
        -(expm1(-2 * alpha * ht$t2)) - (1/(2 * alpha)) * 
        exp(-2 * alpha * (Tmax - ht$t1)) * -(expm1(-2 * 
                                                      alpha * ht$t1))
    }
    cache$len = bl
    cache
  }
  attr(z, "argn") = "alpha"
  return(z)
}
#' bayOU internal function. 
#' 
#' \code{fastbm.lik} is an internal function and not generally called by the user
#' 
#' This is an internal function that modifies the internal function \code{bm.lik} in geiger for efficiency.
.fastbm.lik <- function (cache, dat,SE = NA, model = "OU", ...) {
  cache$dat <- dat
  cache$y[1,][1:cache$ntips] <- dat
  #cache = .prepare.bm.univariate(phy, dat, SE = SE, ...)
  cache$ordering = attributes(cache$phy)$order
  cache$N = cache$n.tip
  cache$n = cache$n.node
  cache$nn = (cache$root + 1):(cache$N + cache$n)
  cache$intorder = as.integer(cache$order[-length(cache$order)])
  cache$tiporder = as.integer(1:cache$N)
  cache$z = length(cache$len)
  FUN = switch(model, OU = .ou.cache.fast(cache))
  ll.bm.direct = function(cache, sigsq, q = NULL, drift = NULL, 
                          se = NULL) {
    n.cache = cache
    given = attr(n.cache$y, "given")
    if (is.null(q)) {
      llf = FUN()
    }
    else {
      llf = FUN(q)
    }
    ll = llf$len
    dd = 0
    if (!is.null(drift)) 
      dd = drift
    adjvar = as.integer(attr(n.cache$y, "adjse"))
    adjSE = any(adjvar == 1)
    .xxSE = function(cache) {
      vv = cache$y["s", ]^2
      ff = function(x) {
        if (any(ww <- adjvar == 1)) {
          vv[which(ww)] = x^2
          return(vv)
        }
        else {
          return(vv)
        }
      }
      return(ff)
    }
    modSE = .xxSE(n.cache)
    vv = as.numeric(modSE(se))
    datC = list(len = as.numeric(ll), intorder = as.integer(n.cache$intorder), 
                tiporder = as.integer(n.cache$tiporder), root = as.integer(n.cache$root), 
                y = as.numeric(n.cache$y["m", ]), var = as.numeric(vv), 
                n = as.integer(n.cache$z), given = as.integer(given), 
                descRight = as.integer(n.cache$children[, 1]), descLeft = as.integer(n.cache$children[, 
                                                                                                      2]), drift = as.numeric(dd))
    parsC = as.numeric(rep(sigsq, n.cache$z))
    out = .Call("bm_direct", dat = datC, pars = parsC, PACKAGE = "geiger")
    loglik <- sum(out$lq)
    if (is.na(loglik)) 
      loglik = -Inf
    attr(loglik, "ROOT.MAX") = out$initM[datC$root]
    class(loglik) = c("glnL", class(loglik))
    return(loglik)
  }
  class(ll.bm.direct) <- c("bm", "dtlik", "function")
  fx_exporter = function() {
    attb = c()
    if (!is.null(qq <- argn(FUN))) {
      adjQ = TRUE
      attb = c(attb, qq)
    }
    else {
      adjQ = FALSE
    }
    attb = c(attb, "sigsq")
    if (any(attr(cache$y, "adjse") == 1)) {
      attb = c(attb, "SE")
    }
    if (model == "drift") {
      attb = c(attb, "drift")
    }
    cache$attb = attb
    lik <- function(pars, ...) {
      recache = function(nodes = NULL, root = ROOT.MAX, 
                         cache) {
        r.cache = cache
        if (root == ROOT.MAX) {
          rtmx = TRUE
        }
        else if (root %in% c(ROOT.OBS, ROOT.GIVEN)) {
          rtmx = FALSE
          r.cache$attb = c(cache$attb, "z0")
        }
        else {
          stop("unusable 'root' type specified")
        }
        r.cache$ROOT.MAX = rtmx
        if (!is.null(nodes)) {
          m = r.cache$y["m", ]
          s = r.cache$y["s", ]
          g = attr(r.cache$y, "given")
          nn = r.cache$nn
          r.cache$y = geiger:::.cache.y.nodes(m, s, g, nn, r.cache$phy, 
                                     nodes = nodes)
        }
        r.cache
      }
      rcache = recache(..., cache = cache)
      attb = rcache$attb
      if (missing(pars)) 
        stop(paste("The following 'pars' are expected:\n\t", 
                   paste(attb, collapse = "\n\t", sep = ""), sep = ""))
      pars = .repars(pars, attb)
      names(pars) = attb
      if (adjQ) 
        q = pars[[qq]]
      else q = NULL
      sigsq = pars[["sigsq"]]
      if ("SE" %in% attb) 
        se = pars[["SE"]]
      else se = NULL
      if ("drift" %in% attb) 
        drift = -pars[["drift"]]
      else drift = 0
      if ("z0" %in% attb) 
        rcache = .fix.root.bm(pars[["z0"]], rcache)
      ll = ll.bm.direct(cache = rcache, sigsq = sigsq, 
                        q = q, drift = drift, se = se)
      return(ll)
    }
    attr(lik, "argn") = attb
    attr(lik, "cache") <- cache
    class(lik) = c("bm", "function")
    lik
  }
  likfx = fx_exporter()
  return(likfx)
}
#' bayOU internal function. 
#' 
#' \code{.emOU.lik} is an internal function and not generally called by the user
#' 
#' This is an internal function that calculates the likelihood of a multioptima OU model.
.emOU.lik <- function(pars,emap,cache,X,SE=0,model="OU",fast=TRUE){
  if(model=="QG"){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model=="OUrepar"){
    repar <- OU.repar(pars)
    pars$alpha <- repar$alpha
    pars$sig2 <- repar$sig2
  }
  TotExp=exp(-cache$height*pars$alpha)
  W <- .edgemap.W(cache,pars,emap,TotExp)
  if(pars$ntheta>1){
    E.th=W%*%pars$theta
  } else {E.th=W*pars$theta}
  X.c<-X-as.vector(E.th)
  if(fast){
    lnL.fx<-.fastbm.lik(cache,X.c,SE=TRUE,model="OU")
    #lnL.fx<-bm.lik(cache$phy,X.c,SE=NA,model="OU")
    loglik <- lnL.fx(pars=c(pars$alpha,pars$sig2,0),root=ROOT.GIVEN)
  } else {
    lnL.fx <- bm.lik(cache,X.c,SE=TRUE,model="OU")
    loglik <- lnL.fx(pars=c(pars$alpha,pars$sig2,0),root=ROOT.GIVEN)
  }
  list(loglik=loglik,W=W,theta=pars$theta,resid=X.c,Exp=E.th)
}

emOU.lik <- function(pars,emap,tree,X,SE=0, method="pruning",model="OU"){
  if(model=="QG"){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model=="OUrepar"){
    repar <- OU.repar(pars)
    pars$alpha <- repar$alpha
    pars$sig2 <- repar$sig2
  }
  if(class(tree)=="phylo"){
    cache = .prepare.ou.univariate(tree, X, SE = SE)
    #cache <- .prepare.ou.univariate(tree, X)
  } else {cache <- tree}
  TotExp=exp(-cache$height*pars$alpha)
  W <- .edgemap.W(cache,pars,emap,TotExp,X)
  E.th=as.matrix(W)%*%pars$theta
  if(method=="pruning"){system.time({
    ##Convert the data into residuals from the expected values
    X.c<-X-as.vector(E.th)
    ##Make the pruning algorithm function for a single-optimum OU model
    lnL<-.fastbm.lik(cache,X.c,SE=TRUE,model=c("OU"))
    ##Call the pruning algorithm function we just made
    loglik=lnL(pars=c(pars$alpha,pars$sig2, 0),root=ROOT.GIVEN)
  })->time}
  if(method=="invert"){system.time({
    ##Standard calculation of the likelihood by inverting the VCV matrix
    ouMatrix <- function(vcvMatrix, alpha)
    {  vcvDiag<-diag(vcvMatrix)
       diagi<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag))
       diagj<-matrix(vcvDiag, nrow=length(vcvDiag), ncol=length(vcvDiag), byrow=T)
       Tij = diagi + diagj - (2 * vcvMatrix)
       vcvRescaled = (1 / (2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcvMatrix))
       return(vcvRescaled)
    }
    Sigma <- ouMatrix(vcv.phylo(cache$phy),pars$alpha)*pars$sig2
    diag(Sigma) <- diag(Sigma)+SE^2
    X.c <- X-as.vector(E.th)
    loglik <- dmnorm(as.vector(X.c),rep(0,length(X.c)),Sigma,log=TRUE)
  })->time}
  list(loglik=loglik,W=W,theta=pars$theta,resid=X.c,Exp=E.th,time=time)
}

#' Calculate the likelihood of a multi-optima OU model from simmap format regimes
#' 
#' Calculates the likelihood of an OU model with regimes stored in tree$maps. 
#' 
#' @rdname sdOU.lik
#' 
smOU.lik <- function(pars,tree,X,SE=0,model="OU"){
  if(class(tree)=="phylo"){
    cache <- .prepare.ou.univariate(tree, X, SE=SE)
  } else {cache <- tree}
  if(model=="QG"){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model=="OUrepar"){
    repar <- OU.repar(pars)
    pars$alpha <- repar$alpha
    pars$sig2 <- repar$sig2
  }
  W <- .simmap.W(cache,pars)
  if(pars$ntheta>1){
    E.th <- W%*%pars$theta
  } else {E.th <- W*pars$theta}
  X.c<-X-as.vector(E.th)
  lnL.fx<-.fastbm.lik(cache,X.c,SE=TRUE,model="OU")
  #lnL.fx<-bm.lik(cache$phy,X.c,SE=NA,model="OU")
  loglik <- lnL.fx(pars=c(pars$alpha,pars$sig2,0),root=ROOT.GIVEN)
  list(loglik=loglik,W=W,theta=pars$theta,resid=X.c,Exp=E.th)
}

#' 
#' 
.smOU.lik <- function(pars,cache,X,SE=0,model="OU"){
  if(model=="QG"){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model=="OUrepar"){
    repar <- OU.repar(pars)
    pars$alpha <- repar$alpha
    pars$sig2 <- repar$sig2
  }
  W <- .simmap.W(cache,pars)
  if(pars$ntheta>1){
    E.th=W%*%pars$theta
  } else {E.th=W*pars$theta}
  X.c<-X-as.vector(E.th)
  lnL.fx<-.fastbm.lik(cache,X.c,SE=TRUE,model="OU")
  loglik <- lnL.fx(pars=c(pars$alpha,pars$sig2,0),root=ROOT.GIVEN)
  list(loglik=loglik,W=W,theta=pars$theta,resid=X.c,Exp=E.th)
}

#' Calculate the likelihood of a multi-optima OU model from parameter list that includes shifts branches (sb),
#' shift locations (loc) and shift optima (t2)
#' 
#' Calculates the likelihood of an OU model with regimes specified by a parameter list
#' 
#' @rdname sdOU.lik
#' @export
OU.lik <- function(pars,tree,X,SE=0,model="OU"){
  if(class(tree)=="phylo"){
    cache <- .prepare.ou.univariate(tree, X, SE=SE)
  } else {cache <- tree}
  if(model=="QG"){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model=="OUrepar"){
    repar <- OU.repar(pars)
    pars$alpha <- repar$alpha
    pars$sig2 <- repar$sig2
  }
  W <- .parmap.W(cache,pars)
  if(pars$ntheta>1){
    E.th <- W%*%pars$theta
  } else {E.th <- W*pars$theta}
  X.c<-X-as.vector(E.th)
  lnL.fx<-.fastbm.lik(cache,X.c,SE=TRUE,model="OU")
  #lnL.fx<-bm.lik(cache$phy,X.c,SE=NA,model="OU")
  loglik <- lnL.fx(pars=c(pars$alpha,pars$sig2,0),root=ROOT.GIVEN)
  list(loglik=loglik,W=W,theta=pars$theta,resid=X.c,Exp=E.th)
}

#' 
#' 
.OU.lik <- function(pars,cache,X,SE=0,model="OU"){
  if(model=="QG"){
    pars$alpha <- QG.alpha(pars)
    pars$sig2 <- QG.sig2(pars)
  }
  if(model=="OUrepar"){
    repar <- OU.repar(pars)
    pars$alpha <- repar$alpha
    pars$sig2 <- repar$sig2
  }
  W <- .parmap.W(cache,pars)
  if(pars$ntheta>1){
    E.th=W%*%pars$theta
  } else {E.th=W*pars$theta}
  X.c<-X-as.vector(E.th)
  lnL.fx<-.fastbm.lik(cache,X.c,SE=TRUE,model="OU")
  loglik <- lnL.fx(pars=c(pars$alpha,pars$sig2,0),root=ROOT.GIVEN)
  list(loglik=loglik,W=W,theta=pars$theta,resid=X.c,Exp=E.th)
}