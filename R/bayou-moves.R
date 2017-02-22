## General proposal function caller, all should have these inputs
.proposalFn <- function(u, ct, D, moves, cache, pars, prior=NULL){
  #ct <- .updateControl(ct,oldpar)
  ctM <-ct[sapply(ct,length)==1]
  move <- names(ct)[u < cumsum(unlist(ctM))][1]
  .moveFn <- get(moves[[move]])
  prop <- .moveFn(move=move,ct=ct,pars=pars,cache=cache,d=D[[move]], prior=prior)
  prop$move <- move
  prop$u <- u
  return(prop)
}


.addshift2map <- function(x,maps=maps,sb=sb,loc=loc,t2=t2){
  m <- maps[[sb[x]]]
  cs.m <- cumsum(m)
  o <- min(which(cs.m>loc[x]))
  if(o==1){
    m[o] <- loc[x]
  } else {
    m[o] <- loc[x]-cs.m[o-1]
  }
  new.m <- cs.m[o]-loc[x]
  names(new.m) <- t2[x]
  M <- c(m[1:o],new.m,m[1:length(m) > o])
  return(M)
}

###Generic multiplier proposal
.multiplierProposal <- function(move,cache,pars,d,ct=NULL, prior=NULL){
  m <- exp(d*(stats::runif(1)-0.5))
  prop <- pars[[move]]*m
  lnHastingsRatio <- log(m)
  pars.new <- pars
  pars.new[[move]] <- prop
  #prior1 <- .prior(pars,emap,cache)
  #prior2 <- .prior(pars.new,emap,cache)
  return(list("pars"=pars.new, "hr"=lnHastingsRatio))
}

.slidingWindowProposal <- function(cache, pars, d, move, ct=NULL, prior=NULL){
  prop <- d*(stats::runif(1)-0.5)+pars[[move]]
  lnHastingsRatio <- 0
  pars.new <- pars
  pars.new[[move]] <- prop
  #pr <- .prior(pars.new,emap,cache)-.prior(pars,emap,cache)
  return(list("pars" = pars.new, "hr"=lnHastingsRatio))
}

##Adjust a randomly selected theta parameter
.adjustTheta <- function(cache, pars, d, type="slidingwindow",move=NULL,ct=NULL, prior=NULL){
  j <- sample(1:pars$ntheta,1)
  if(type=="slidingwindow"){
    ##Generate sliding-window proposal
    prop <- d*(stats::runif(1)-0.5)+pars$theta[j]
    lnHastingsRatio <- 0
    pars.new <- pars
    pars.new$theta[j] <- prop
    #pr <- .prior(pars.new,emap,cache)-.prior(pars,emap,cache)
    return(list("pars" = pars.new, "hr"=lnHastingsRatio, "theta" = j))
  }
  if(type=="multiplier"){
    ##Generate multiplier proposal
    m <- exp(d*(stats::runif(1)-0.5))
    prop <- pars$theta[j]*m
    lnHastingsRatio <- log(m)
    pars.new <- pars
    pars.new$theta[j] <- prop
    #pr <- .prior(pars.new,emap,cache)-.prior(pars,emap,cache)
    return(list("pars" = pars.new, "hr"=lnHastingsRatio, "theta" = j))
  }
}

##Adjust a randomly selected theta parameter
.vectorMultiplier <- function(cache, pars, d, move,ct=NULL, prior=NULL){
  j <- sample(1:length(pars[[move]]),1)
  ##Generate multiplier proposal
  m <- exp(d*(stats::runif(1)-0.5))
  prop <- pars[[move]][j]*m
  lnHastingsRatio <- log(m)
  pars.new <- pars
  pars.new[[move]][j] <- prop
  #pr <- .prior(pars.new,emap,cache)-.prior(pars,emap,cache)
  return(list("pars" = pars.new, "hr"=lnHastingsRatio, "theta" = j))
}

.vectorSlidingWindow <- function(cache, pars, d, move,ct=NULL, prior=NULL){
  j <- sample(1:length(pars[[move]]),1)
  prop <- d*(stats::runif(1)-0.5)+pars[[move]][j]
  lnHastingsRatio <- 0
  pars.new <- pars
  pars.new[[move]][j] <- prop
  #pr <- .prior(pars.new,emap,cache)-.prior(pars,emap,cache)
  return(list("pars" = pars.new, "hr"=lnHastingsRatio, "j" = j))
}

#' MCMC move for sliding a shift up or down to neighboring branches, or within a branch
.slidespace <- function(j, pars, cache, ct, map){
  sb <- rep(NA, 5)
  t2 <- pars$t2[j]
  #map.i <- which(names(map$segs)==(pars$sb[j]))
  map.i <- which(map$branch==pars$sb[j])
  m.sb <- map$segs[map.i]
  m.t2 <- map$theta[map.i]
  m.i <- which(m.t2==t2)
  sb[1:2] <- map$branch[map.i][1]
  U0 <- m.sb[m.i]
  D0 <- m.sb[m.i-1]
  nodes <- cache$edge[pars$sb[j],]
  R1 = FALSE
  if(m.i-1 == 1){
    if(nodes[1]==cache$ntips+1){
      root.branch <- setdiff(c(2*cache$ntips-3, 2*cache$ntips-2), pars$sb[j])
      #D1 <- map$segs[as.character(root.branch)]
      ind <- match(root.branch, map$branch)
      D1 <- map$segs[ind]
      sb[5] <- map$branch[ind]
      R1 = TRUE
    } else {
      anc <- which(cache$edge[,2] == nodes[1])
      ind <- which(map$branch==anc)
      anc.m <- map$segs[ind]
      D1 <- anc.m[length(anc.m)]
      sb[5] <- map$branch[ind][1]
    }
  } else {
    D1 <- 0
  }
  if(m.i == length(map.i)){
    desc <- which(cache$edge[,1] == nodes[2])
    if(length(desc)>0){
      #desc.m <- map$segs[as.character(desc)]
      ind <- match(desc, map$branch)
      desc.m <- map$segs[match(desc, map$branch)]
      sb[3:4] <- map$branch[ind]
      U1 <- desc.m[1]
      U2 <- desc.m[2]
    } else {
      U1 <- U2 <- 0
      } 
    } else {
      U1 <- U2 <- 0
    }
  pp = c(U0, D0, U1, U2, D1)
 # sb = as.numeric(names(pp))
  if(any(ct$sb$bmax[sb]==0, na.rm=TRUE)){
    pp[sb %in% which(ct$sb$bmax==0)] <- 0
  }
  names(pp) <- c("U0", "D0", "U1", "U2", "D1")
  if(R1) names(pp)[5] <- "R1"
  if(any(ct$sb$bmax!=Inf)){
    tb.sb <- table(pars$sb[!(pars$sb==pars$sb[j])])
    nm.sb <- as.numeric(names(tb.sb))
    full <- nm.sb[ct$sb$bmax[nm.sb] <= tb.sb]
    if(length(full)>0) pp[sb %in% full] <- 0
  }
  return(list(sb=sb,pp=pp))
}

.slide <- function(pars, cache, d, ct, move=NULL, prior=NULL){
  map <- .pars2map(pars,cache)
  j <- sample(1:pars$k,1)
  t2 <- pars$t2[j]
  space <- .slidespace(j, pars, cache, ct, map)
  pars.new <- pars
  mv <- .sample(1:5,1,prob=space$pp/sum(space$pp))
  type <- names(space$pp)[mv]
  l <- stats::runif(1,0,space$pp[type])
  if(type=="U0"){
    pars.new$loc[j] <- l + pars$loc[j]
  } else {
    if(type=="D0"){
      pars.new$loc[j] <- pars$loc[j] - l
    } else {
      pars.new$sb[j] <- space$sb[mv]
      if(type %in% c("U1","U2")){
        pars.new$loc[j] <- l
      } else {
        if(type == "D1"){
          pars.new$loc[j] <- cache$edge.length[space$sb[mv]] - l
        } else {
          if (type=="R1"){
            pars.new$loc[j] <- l
            pars.new$theta[c(1,t2)] <- pars.new$theta[c(t2,1)]
          }
        }
      }
    }
  }
  map.new <- .pars2map(pars.new, cache)
  space.new <- .slidespace(j,pars.new, cache, ct, map.new)
  hr <- log(1/sum(space.new$pp))-log(1/(sum(space$pp)))
  return(list(pars=pars.new, hr=hr, decision = type))
}


#' MCMC move for splitting or collapsing a shift on phylogeny
.splitmerge <- function(pars, cache, d, ct, move=NULL, prior=NULL){
  nbranch <- length(cache$edge.length)
  TH <- sum(cache$edge.length)
  v <- stats::runif(1)
  sb.max <- ct$sb$bmax
  sb.taken <- rep(0,2*cache$ntips-2)
  sb.table <- table(pars$sb)
  sb.taken[as.numeric(names(sb.table))] <- sb.table
  sb.prob <- ct$sb$prob
  sb.prob[sb.max <= sb.taken] <- 0
  if(v < ct$bk[pars$ntheta]/(ct$bk[pars$ntheta]+ct$dk[pars$ntheta])){
    decision <- "birth"
    sb.j <- sample(1:(2*cache$ntips-2),1,prob=sb.prob)
    loc.j <- stats::runif(1,min=0,max=cache$edge.length[sb.j])
    t2.j <- pars$ntheta+1
    pars.new <- pars
    pars.new$sb <- c(pars$sb, sb.j)
    pars.new$loc <- c(pars$loc, loc.j)
    pars.new$t2 <- c(pars$t2, t2.j)
    map.new <- .pars2map(pars.new, cache)
    t1 <- map.new$theta[max(which(map.new$theta==t2.j))-1]
    t2W <- sum(map.new$segs[map.new$theta==t2.j])
    t1W <- sum(map.new$segs[map.new$theta==t1])
    r <- t2W/(t1W+t2W)
    u <- stats::runif(1,-0.5,0.5)*d
    pars.new$theta[t1] <- pars$theta[t1]-u*r
    pars.new$theta[t2.j] <- pars$theta[t1]+u*(1-r)   
    pars.new$k <- pars$k+1
    pars.new$ntheta <- pars$ntheta + 1
    pars.new$sb <- c(pars$sb,sb.j)
    pars.new$loc <- c(pars$loc,loc.j)
    pars.new$t2 <- c(pars$t2,t2.j)
    hr <- log(ct$dk[pars.new$ntheta]*1/pars.new$k*d)-log(ct$bk[pars$ntheta]*sb.prob[sb.j]/sum(sb.prob))
  } else {
    decision <- "death"
    j <- sample(1:pars$k,1)
    pars.new <- pars
    pars.new$k <- pars$k-1
    pars.new$ntheta <- pars$ntheta-1
    pars.new$sb <- pars$sb[-j]
    pars.new$loc <- pars$loc[-j]
    pars.new$t2 <- pars$t2[-pars$k]
    pars.new$theta <- pars$theta[-(j+1)]
    map <- .pars2map(pars, cache)
    t2.j <- pars$t2[j]
    sb.j <- pars$sb[j]
    t1 <- map$theta[max(which(map$theta==t2.j))-1]
    t2W <- sum(map$segs[map$theta==t2.j])
    t1W <- sum(map$segs[map$theta==t1])
    r <- t2W/(t1W+t2W)
    pars.new$theta[t1-(t1>t2.j)] <- pars$theta[t1]*(1-r)+pars$theta[t2.j]*r
    sb.prob[sb.j] <- ct$sb$prob[sb.j]
    hr <- log(ct$bk[pars.new$ntheta]*sb.prob[sb.j]/sum(sb.prob))-log(ct$dk[pars$ntheta]*1/pars$k*d)
  }
  return(list(pars=pars.new, hr=hr, decision=decision, sb.prob=sb.prob))
}
#.add2map <- function(map, cache, pars, sb.j, loc.j, t2.j){
#    j <- which(names(map$segs)==sb.j)
#    m <- map$segs[j]
#    new.m <- c(m,loc.j)
#    o <- order(new.m)
#    new.m <- new.m[o]
#    names(new.m) <- rep(sb.j,length(o))
#    t2 <- map$theta[j]
#    o.t2 <- which(o==length(o))
#    new.t2 <- c(t2[1:o.t2],t2.j)
#    if(o.t2 != length(t2)){
#      new.t2 <- c(new.t2, t2[(o.t2+1):(length(t2))])
#    }
#    names(new.t2) <- rep(sb.j,length(o))
#    new.segs <- c(map$segs[1:(min(j)-1)], new.m)
#    new.theta <- c(map$theta[1:(min(j)-1)], new.t2)
#    if(length(new.segs)!=(length(map$segs)+1)){
#      new.segs <- c(new.segs, map$segs[(max(j)+1):length(map$segs)])
#      new.theta <- c(new.theta, map$theta[(max(j)+1):length(map$segs)])
#    }
#    if(o.t2 == length(t2)){
#      desc <- cache$bdesc[sb.j][[1]]
#      desc.desc.sb <- cache$bdesc[pars$sb[pars$sb %in% desc]]
##      relabel <- setdiff(desc,unlist(desc.desc.sb,F,F))
#        new.theta[as.character(relabel)] <- t2.j
#      }
#    }
 # return(list(segs=new.segs, theta=new.theta))
#}

#' MCMC move for splitting or collapsing a shift on phylogeny
.splitmergebd <- function(pars, cache, d, ct, move=NULL, prior=NULL){
  splitmergepars <- attributes(ct)$splitmergepars
  nbranch <- length(cache$edge.length)
  TH <- sum(cache$edge.length)
  v <- stats::runif(1)
  sb.max <- ct$sb$bmax
  sb.taken <- rep(0,2*cache$ntips-2)
  sb.table <- table(pars$sb)
  sb.taken[as.numeric(names(sb.table))] <- sb.table
  sb.prob <- ct$sb$prob
  sb.prob[sb.max <= sb.taken] <- 0
  if(v < ct$bk[pars$ntheta]/(ct$bk[pars$ntheta]+ct$dk[pars$ntheta])){
    decision <- "birth"
    sb.j <- sample(1:(2*cache$ntips-2),1,prob=sb.prob)
    loc.j <- stats::runif(1,min=0,max=cache$edge.length[sb.j])
    t2.j <- pars$ntheta+1
    pars.new <- pars
    pars.new$sb <- c(pars$sb, sb.j)
    pars.new$loc <- c(pars$loc, loc.j)
    pars.new$t2 <- c(pars$t2, t2.j)
    map.new <- .pars2map(pars.new, cache)
    t1 <- map.new$theta[max(which(map.new$theta==t2.j))-1]
    t2W <- sum(map.new$segs[map.new$theta==t2.j])
    t1W <- sum(map.new$segs[map.new$theta==t1])
    r <- t2W/(t1W+t2W)
    for(i in 1:length(splitmergepars)){
      u <- stats::runif(1,-0.5,0.5)*d[i]
      pars.new[[splitmergepars[i]]][t1] <- pars[[splitmergepars[i]]][t1]-u*r
      pars.new[[splitmergepars[i]]][t2.j] <- pars[[splitmergepars[i]]][t1]+u*(1-r)
    }
    pars.new$k <- pars$k+1
    pars.new$ntheta <- pars$ntheta + 1
    pars.new$sb <- c(pars$sb,sb.j)
    pars.new$loc <- c(pars$loc,loc.j)
    pars.new$t2 <- c(pars$t2,t2.j)
    hr <- log(ct$dk[pars.new$ntheta]*1/pars.new$k*prod(d))-log(ct$bk[pars$ntheta]*sb.prob[sb.j]/sum(sb.prob))
  } else {
    decision <- "death"
    j <- sample(1:pars$k,1)
    pars.new <- pars
    pars.new$k <- pars$k-1
    pars.new$ntheta <- pars$ntheta-1
    pars.new$sb <- pars$sb[-j]
    pars.new$loc <- pars$loc[-j]
    pars.new$t2 <- pars$t2[-pars$k]
    for(i in 1:length(splitmergepars)){
      pars.new[[splitmergepars[i]]] <- pars[[splitmergepars[i]]][-(j+1)]
    }
    map <- .pars2map(pars, cache)
    t2.j <- pars$t2[j]
    sb.j <- pars$sb[j]
    t1 <- map$theta[max(which(map$theta==t2.j))-1]
    t2W <- sum(map$segs[map$theta==t2.j])
    t1W <- sum(map$segs[map$theta==t1])
    r <- t2W/(t1W+t2W)
    for(i in 1:length(splitmergepars)){
      pars.new[[splitmergepars[i]]][t1-(t1>t2.j)] <- pars[[splitmergepars[i]]][t1]*(1-r)+pars[[splitmergepars[i]]][t2.j]*r
    }
    sb.prob[sb.j] <- ct$sb$prob[sb.j]
    hr <- log(ct$bk[pars.new$ntheta]*sb.prob[sb.j]/sum(sb.prob))-log(ct$dk[pars$ntheta]*1/pars$k*prod(d))
  }
  return(list(pars=pars.new, hr=hr, decision=decision, sb.prob=sb.prob))
}

.slide2 <- function(pars, cache, d, ct, move=NULL, prior=NULL){
  splitmergepars <- attributes(ct)$splitmergepars
  map <- .pars2map(pars,cache)
  j <- sample(1:pars$k,1)
  t2 <- pars$t2[j]
  space <- .slidespace(j, pars, cache, ct, map)
  pars.new <- pars
  mv <- .sample(1:5,1,prob=space$pp/sum(space$pp))
  type <- names(space$pp)[mv]
  l <- stats::runif(1,0,space$pp[type])
  if(type=="U0"){
    pars.new$loc[j] <- l + pars$loc[j]
  } else {
    if(type=="D0"){
      pars.new$loc[j] <- pars$loc[j] - l
    } else {
      pars.new$sb[j] <- space$sb[mv]
      if(type %in% c("U1","U2")){
        pars.new$loc[j] <- l
      } else {
        if(type == "D1"){
          pars.new$loc[j] <- cache$edge.length[space$sb[mv]] - l
        } else {
          if (type=="R1"){
            pars.new$loc[j] <- l
            for(i in 1:length(splitmergepars)){
              pars.new[[splitmergepars[i]]][c(1,t2)] <- pars.new[[splitmergepars[i]]][c(t2,1)]
            }
            
          }
        }
      }
    }
  }
  map.new <- .pars2map(pars.new, cache)
  space.new <- .slidespace(j,pars.new, cache, ct, map.new)
  hr <- log(1/sum(space.new$pp))-log(1/(sum(space$pp)))
  return(list(pars=pars.new, hr=hr, decision = type))
}

.jointHalflifeVyProposal <- function(move, cache, pars, d, ct = NULL, prior=NULL){
  dSig <- matrix(c(d[1], d[2]*d[1], d[2]*d[1], d[1]), ncol=2)
  m <- mnormt::rmnorm(1, c(1,1), dSig)
  if(move == "halflife"){
    move2 <- "Vy"
  }
  if(move == "Vy"){
    move2 <- "halflife"
  }
  prop1 <- pars[[move]]*m[1]
  prop2 <- pars[[move2]]*m[2]
  pars.new <- pars
  pars.new[[move]] <- prop1
  pars.new[[move2]] <- prop2
  lnHastingsRatio <- mnormt::dmnorm(m, c(1,1), dSig, log = TRUE) - mnormt::dmnorm(1/m, c(1,1), dSig, log=TRUE) 
  return(list(pars = pars.new, hr=lnHastingsRatio))
}

.vectorSlidingWindowSplit <- function(move, cache, pars, d, ct=NULL, prior=NULL){
  j <- sample(1:length(pars[[move]]),1)
  if(j==1){
    prop <- d[1]*(stats::runif(1)-0.5)+pars[[move]][j]
    lnHastingsRatio <- 0
    pars.new <- pars
    pars.new[[move]][j] <- prop
    #pr <- .prior(pars.new,emap,cache)-.prior(pars,emap,cache)
    return(list("pars" = pars.new, "hr"=lnHastingsRatio, "j" = j))
  } else {
    hc <- .sample(1:length(d), 1, replace=FALSE)
    prop <- d[1]*(stats::runif(1)-0.5)+pars[[move]][j]
    lnHastingsRatio <- 0
    pars.new <- pars
    pars.new[[move]][j] <- prop
    #pr <- .prior(pars.new,emap,cache)-.prior(pars,emap,cache)
    return(list("pars" = pars.new, "hr"=lnHastingsRatio, "j" = j))
  }
}

## Split merge proposal from Green 1995, only works for positive parameters
.splitmergeGreen <- function(pars, cache, d, ct, move=NULL, prior=NULL){
  splitmergepars <- attributes(ct)$splitmergepars
  nbranch <- length(cache$edge.length)
  TH <- sum(cache$edge.length)
  v <- stats::runif(1)
  sb.max <- ct$sb$bmax
  sb.taken <- rep(0,2*cache$ntips-2)
  sb.table <- table(pars$sb)
  sb.taken[as.numeric(names(sb.table))] <- sb.table
  sb.prob <- ct$sb$prob
  sb.prob[sb.max <= sb.taken] <- 0
  if(v < ct$bk[pars$ntheta]/(ct$bk[pars$ntheta]+ct$dk[pars$ntheta])){
    decision <- "birth"
    sb.j <- sample(1:(2*cache$ntips-2),1,prob=sb.prob)
    loc.j <- 0 #runif(1,min=0,max=cache$edge.length[sb.j])
    t2.j <- pars$ntheta+1
    pars.new <- pars
    pars.new$sb <- c(pars$sb, sb.j)
    pars.new$loc <- c(pars$loc, loc.j)
    pars.new$t2 <- c(pars$t2, t2.j)
    map.new <- .pars2map(pars.new, cache)
    t1 <- map.new$theta[max(which(map.new$theta==t2.j))-1]
    t2W <- sum(map.new$segs[map.new$theta==t2.j])
    t1W <- sum(map.new$segs[map.new$theta==t1])
    r <- t2W/(t1W+t2W)
    .hr <- NULL
    for(i in 1:length(splitmergepars)){
      u <- stats::runif(1,0,1)*d[i]
      pars.new[[splitmergepars[i]]][t1] <- pars[[splitmergepars[i]]][t1]*(u/(1-u))^(t1W/(t1W+t2W))
      pars.new[[splitmergepars[i]]][t2.j] <- pars[[splitmergepars[i]]][t1]*((1-u)/(u))^(t2W/(t1W+t2W))
      .hr <- c(.hr, log((pars.new[[splitmergepars[i]]][t1] + pars.new[[splitmergepars[i]]][t2.j])^2/pars[[splitmergepars[i]]][t1]))
    }
    pars.new$k <- pars$k+1
    pars.new$ntheta <- pars$ntheta + 1
    pars.new$sb <- c(pars$sb,sb.j)
    pars.new$loc <- c(pars$loc,loc.j)
    pars.new$t2 <- c(pars$t2,t2.j)
    hr <- sum(.hr)+log(d)
  } else {
    decision <- "death"
    j <- sample(1:pars$k,1)
    pars.new <- pars
    pars.new$k <- pars$k-1
    pars.new$ntheta <- pars$ntheta-1
    pars.new$sb <- pars$sb[-j]
    pars.new$loc <- pars$loc[-j]
    pars.new$t2 <- pars$t2[-pars$k]
    .hr <- NULL
    map <- .pars2map(pars, cache)
    t2.j <- pars$t2[j]
    sb.j <- pars$sb[j]
    t1 <- map$theta[max(which(map$theta==t2.j))-1]
    t2W <- sum(map$segs[map$theta==t2.j])
    t1W <- sum(map$segs[map$theta==t1])
    r <- t2W/(t1W+t2W)
    for(i in 1:length(splitmergepars)){
      pars.new[[splitmergepars[i]]] <- pars[[splitmergepars[i]]][-(j+1)]
    }
    
    for(i in 1:length(splitmergepars)){
      pars.new[[splitmergepars[i]]][t1-(t1>t2.j)] <- pars[[splitmergepars[i]]][t1]*(1-r)+pars[[splitmergepars[i]]][t2.j]*r
      .hr <- c(.hr, log(pars.new[[splitmergepars[i]]][t1-(t1>t2.j)]/(pars[[splitmergepars[i]]][t1] + pars[[splitmergepars[i]]][t2.j])^2))
    }
    sb.prob[sb.j] <- ct$sb$prob[sb.j]
    hr <- sum(.hr)-1*log(d)
  }
  return(list(pars=pars.new, hr=hr, decision=decision, sb.prob=sb.prob))
}

## Split merge proposal where new value is drawn from the prior
.splitmergePrior <- function(pars, cache, d, ct, move=NULL, prior){
  splitmergepars <- attributes(ct)$splitmergepars
  nbranch <- length(cache$edge.length)
  TH <- sum(cache$edge.length)
  v <- stats::runif(1)
  sb.max <- ct$sb$bmax
  sb.taken <- rep(0,2*cache$ntips-2)
  sb.table <- table(pars$sb)
  sb.taken[as.numeric(names(sb.table))] <- sb.table
  sb.prob <- ct$sb$prob
  sb.prob[sb.max <= sb.taken] <- 0
  if(v < ct$bk[pars$ntheta]/(ct$bk[pars$ntheta]+ct$dk[pars$ntheta])){
    decision <- "birth"
    sb.j <- sample(1:(2*cache$ntips-2),1,prob=sb.prob)
    loc.j <- 0 #runif(1,min=0,max=cache$edge.length[sb.j])
    t2.j <- pars$ntheta+1
    pars.new <- pars
    pars.new$sb <- c(pars$sb, sb.j)
    pars.new$loc <- c(pars$loc, loc.j)
    pars.new$t2 <- c(pars$t2, t2.j)
    map.new <- .pars2map(pars.new, cache)
    t1 <- map.new$theta[max(which(map.new$theta==t2.j))-1]
    t2W <- sum(map.new$segs[map.new$theta==t2.j])
    t1W <- sum(map.new$segs[map.new$theta==t1])
    r <- t2W/(t1W+t2W)
    .hr <- NULL
    for(i in 1:length(splitmergepars)){
      #u <- runif(1,0,1)*d[i]
      rfx <- attributes(prior)$rfunctions[[paste("d", splitmergepars[i], sep="")]]
      dfx <- attributes(prior)$functions[[paste("d", splitmergepars[i], sep="")]]
      u <- rfx(1)
      #pars.new[[splitmergepars[i]]][t1] <- pars[[splitmergepars[i]]][t1]
      #pars.new[[splitmergepars[i]]][t2.j] <- pars[[splitmergepars[i]]][t1]*((1-u)/(u))^(t2W/(t1W+t2W))
      pars.new[[splitmergepars[i]]][t2.j] <- u
      .hr <- c(.hr, -1*dfx(u))
    }
    pars.new$k <- pars$k+1
    pars.new$ntheta <- pars$ntheta + 1
    pars.new$sb <- c(pars$sb,sb.j)
    pars.new$loc <- c(pars$loc,loc.j)
    pars.new$t2 <- c(pars$t2,t2.j)
    hr <- sum(.hr)
  } else {
    decision <- "death"
    j <- sample(1:pars$k,1)
    pars.new <- pars
    pars.new$k <- pars$k-1
    pars.new$ntheta <- pars$ntheta-1
    pars.new$sb <- pars$sb[-j]
    pars.new$loc <- pars$loc[-j]
    pars.new$t2 <- pars$t2[-pars$k]
    .hr <- NULL
    map <- .pars2map(pars, cache)
    t2.j <- pars$t2[j]
    sb.j <- pars$sb[j]
    t1 <- map$theta[max(which(map$theta==t2.j))-1]
    t2W <- sum(map$segs[map$theta==t2.j])
    t1W <- sum(map$segs[map$theta==t1])
    r <- t2W/(t1W+t2W)
    for(i in 1:length(splitmergepars)){
      pars.new[[splitmergepars[i]]] <- pars[[splitmergepars[i]]][-(j+1)]
    }
    
    for(i in 1:length(splitmergepars)){
      rfx <- attributes(prior)$rfunctions[[paste("d", splitmergepars[i], sep="")]]
      u <- rfx(1)
      dfx <- attributes(prior)$functions[[paste("d", splitmergepars[i], sep="")]]
      pars.new[[splitmergepars[i]]][t1-(t1>t2.j)] <- pars[[splitmergepars[i]]][t1]*(1-r)+pars[[splitmergepars[i]]][t2.j]*r
      .hr <- c(.hr, dfx(pars[[splitmergepars[i]]][t2.j]))
    }
    sb.prob[sb.j] <- ct$sb$prob[sb.j]
    hr <- sum(.hr)
  }
  return(list(pars=pars.new, hr=hr, decision=decision, sb.prob=sb.prob))
}

