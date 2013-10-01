.proposalFn <- function(u,ct,D,moves,cache,oldpar){
  #ct <- .updateControl(ct,oldpar)
  if(oldpar$k==0){
    ct <- .updateControl(ct, pars)
  }
  ctM <-ct[sapply(ct,length)==1]
  move <- names(ct)[u < cumsum(unlist(ctM))][1]
  .moveFn <- get(moves[[move]])
  prop <- .moveFn(move=move,ct=ct,pars=oldpar,cache=cache,d=D[[move]])
  prop$move <- move
  prop$u <- u
  return(prop)
}


.splitmerge.emap <- function(emap,pars,cache,plim=plim,d,move=NULL){
  j <- round(runif(1,0.5,2*cache$ntips-2+0.5),0)
  emap.new <- emap
  pars.new <- pars
  b <- dim(emap)[1]
  if(emap$sh[j]==1){
    decision <- "death"
    emap.new$sh[j] <- 0
    emap.new$r1[j] <- emap.new$r[j]
    emap.new$r2[j] <- 0
    old.t <- emap.new$t2[j]
    new.t <- emap.new$t1[j]
    t1W <- sum(emap.new$r1[emap.new$t1==new.t]+emap.new$r2[emap.new$t1==new.t])
    t2W <- sum(emap.new$r1[emap.new$t2==old.t]+emap.new$r2[emap.new$t2==old.t])
    r <- t2W/(t1W+t2W)
    emap.new$t1[emap.new$t1==old.t] <- new.t
    emap.new$t2[emap.new$t2==old.t] <- new.t
    if(sum(emap.new$t1>old.t)>0){
      emap.new$t1[emap.new$t1>old.t] <-  emap.new$t1[emap.new$t1>old.t]-1
    }
    if(sum(emap.new$t2>old.t)>0){
      emap.new$t2[emap.new$t2>old.t] <-  emap.new$t2[emap.new$t2>old.t]-1
    }
    pars.new$optima[new.t] <- (1-r)*pars.new$optima[new.t]+(r)*pars.new$optima[old.t]
    pars.new$optima <- pars.new$optima[-old.t]
    pars.new$ntheta <- pars$ntheta-1
    pars.new$nb <- pars$nb-1
    hr <- ifelse(pars.new$nb==pars$nb-1,log(1/d),log(0))
  }
  if(emap$sh[j]==0){
    decision <- "birth"
    emap.new$sh[j] <- 1
    emap.new$t2[j] <- pars$ntheta+1
    u <- runif(1)
    emap.new$r1[j] <- u*emap.new$r[j]
    emap.new$r2[j] <- (1-u)*emap.new$r[j]
    pars.new$ntheta <- pars$ntheta+1
    pars.new$nb <- pars$nb+1
    t1W <- sum(emap.new$r1[emap.new$t1==emap.new$t1[j]]+emap.new$r2[emap.new$t1==emap.new$t1[j]])
    t2W <- sum(emap.new$r1[emap.new$t2==emap.new$t2[j]]+emap.new$r2[emap.new$t2==emap.new$t2[j]])
    r <- t2W/(t1W+t2W)
    u <- runif(1,-0.5,0.5)*d
    pars.new$optima[emap.new$t1[j]] <- pars$optima[emap.new$t1[j]]-u*r
    pars.new$optima[emap.new$t2[j]] <- pars$optima[emap.new$t1[j]]+u*(1-r)   
    if(emap.new$e2[j]>cache$ntips){
      nopt <- rep(1,b)
      nopt[emap.new$e2] <- emap.new$t2
      for(i in (j-1):1){
        if(emap.new$sh[i]==0){
          nopt[emap.new$e2[i]] <- nopt[emap.new$e1[i]]
        }
      }
      emap.new$t1 <- nopt[emap.new$e1]
      emap.new$t2 <- nopt[emap.new$e2]
    }
    hr <-  ifelse(pars.new$nb==pars$nb+1,log(d),log(0))
  }
  list(pars=pars.new,emap=emap.new,decision=decision,branch=j,hr=hr)
}

addshift2map <- function(x,maps=maps,sb=sb,loc=loc,t2=t2){
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

#maps <- lapply(tree$edge.length,function(x){y <- x; names(y) <- 1; y})
#maps[sb] <- lapply(1:pars$k,function(x){ y <- c(maps[[sb[x]]],loc[x]); names(y)[length(y)] <- t2[x]; y})
.splitmergeSimmap <- function(pars,cache,ct,d,move=NULL,maps=NULL){
  if(is.null(maps)){
    maps <- pars2simmap(pars,cache$phy,sim.theta=FALSE,theta=pars$theta)$tree$maps
  }
  T <- sum(cache$edge.length)
  v <- runif(1)
  if(v < ct$bk[pars$ntheta]/(ct$bk[pars$ntheta]+ct$dk[pars$ntheta])){
    decision <- "birth"
    sb.j <- sample(1:(2*ntips-2),1,prob=cache$edge.length/T)
    loc.j <- runif(1,min=0,max=cache$edge.length[sb.j])
    t2.j <- max(pars$t2)+1
    maps.new <- maps
    maps.new[[sb.j]] <- addshift2map(1,maps.new,sb.j,loc.j,t2.j)
    t1 <- as.integer(names(maps.new[[sb.j]])[which(names(maps.new[[sb.j]])==t2.j)-1])
    nopt <- rep(1,nbranch)
    for(i in nbranch:1){
      if(i %in% c(pars$sb,sb.j)){
        opt <- as.integer(names(maps.new[[i]])[length(maps.new[[i]])])
        nopt[tree$edge[i,2]] <- opt
        names(maps.new[[i]])[1] <- nopt[tree$edge[i,1]]
      } else {
        names(maps.new[[i]])[1] <- nopt[tree$edge[i,1]] 
        nopt[tree$edge[i,2]] <- nopt[tree$edge[i,1]]
      }
    }
    shiftdown <- nopt[tree$edge[,1]]
    for(j in 1:nbranch){
      names(maps.new[[j]])[1] <-shiftdown[j]
    }
    segs <- unlist(maps.new,FALSE,TRUE)
    t2W <- sum(segs[names(segs)==t2.j])
    t1W <- sum(segs[names(segs)==t1])
    r <- t2W/(t1W+t2W)
    u <- runif(1,-0.5,0.5)*d
    pars.new <- pars
    pars.new$theta[t1] <- pars$theta[t1]-u*r
    pars.new$theta[t2.j] <- pars$theta[t1]+u*(1-r)   
    pars.new$k <- pars$k+1
    pars.new$ntheta <- pars$ntheta + 1
    pars.new$sb <- c(pars$sb,sb.j)
    pars.new$loc <- c(pars$loc,loc.j)
    pars.new$t2 <- c(pars$t2,t2.j)
    hr <- log(ct$dk[pars.new$ntheta]*T*d)-log(ct$bk[pars$ntheta]*pars.new$k)
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
    t2 <- pars$t2[j]
    t1 <- as.integer(names(maps[[pars$sb[j]]])[which(names(maps[[pars$sb[j]]])==t2)-1])
    segs <- unlist(maps,FALSE,TRUE)
    t2W <- sum(segs[names(segs)==t2])
    t1W <- sum(segs[names(segs)==t1])
    r <- t2W/(t1W+t2W)
    pars.new$theta[t1-(t1>t2)] <- pars$theta[t1]*(1-r)+pars$theta[t2]*r
    hr <- log(ct$bk[pars.new$ntheta]*pars.new$k)-log(ct$dk[pars$ntheta]*T*d)
    maps.new <- pars2simmap(pars.new,cache$phy,sim.theta=FALSE,theta=pars.new$theta)$tree$maps
  }
  cache$maps <- maps.new
  return(list(pars=pars.new,cache=cache,decision=decision,hr=hr))
}


###Generic multiplier proposal
.multiplierProposal <- function(move,cache,pars,d,ct=NULL){
  m <- exp(d*(runif(1)-0.5))
  prop <- pars[[move]]*m
  lnHastingsRatio <- log(m)
  pars.new <- pars
  pars.new[[move]] <- prop
  #prior1 <- .prior(pars,emap,cache)
  #prior2 <- .prior(pars.new,emap,cache)
  return(list("pars"=pars.new, "hr"=lnHastingsRatio))
}

##Adjust a randomly selected theta parameter
.adjustTheta <- function(cache, pars, d, type="slidingwindow",move=NULL,ct=NULL){
  j <- sample(1:pars$ntheta,1)
  if(type=="slidingwindow"){
    ##Generate sliding-window proposal
    prop <- d*(runif(1)-0.5)+pars$theta[j]
    lnHastingsRatio <- 0
    pars.new <- pars
    pars.new$theta[j] <- prop
    #pr <- .prior(pars.new,emap,cache)-.prior(pars,emap,cache)
    return(list("pars" = pars.new, "hr"=lnHastingsRatio, "theta" = j))
  }
  if(type=="multiplier"){
    ##Generate multiplier proposal
    m <- exp(d*(runif(1)-0.5))
    prop <- pars$theta[j]*m
    lnHastingsRatio <- log(m)
    pars.new <- pars
    pars.new$theta[j] <- prop
    #pr <- .prior(pars.new,emap,cache)-.prior(pars,emap,cache)
    return(list("pars" = pars.new, "hr"=lnHastingsRatio, "theta" = j))
  }
}


#' MCMC move for sliding a shift up or down to neighboring branches, or within a branch
.slidespace <- function(j, pars, cache, ct, map){
  t2 <- pars$t2[j]
  map.i <- which(names(map$segs)==(pars$sb[j]))
  m.sb <- map$segs[map.i]
  m.t2 <- map$theta[map.i]
  m.i <- which(m.t2==t2)
  U0 <- m.sb[m.i]
  D0 <- m.sb[m.i-1]
  nodes <- cache$edge[pars$sb[j],]
  R1 = FALSE
  if(m.i-1 == 1){
    if(nodes[1]==cache$ntips+1){
      root.branch <- setdiff(c(2*cache$ntips-3, 2*cache$ntips-2), pars$sb[j])
      D1 <- map$segs[as.character(root.branch)]
      R1 = TRUE
    } else {
      anc <- which(cache$edge[,2] == nodes[1])
      anc.m <- map$segs[names(map$segs)==anc]
      D1 <- anc.m[length(anc.m)]
    }
  } else {
    D1 <- 0
  }
  if(m.i == length(map.i)){
    desc <- which(cache$edge[,1] == nodes[2])
    if(length(desc)>0){
      desc.m <- map$segs[as.character(desc)]
      U1 <- desc.m[1]
      U2 <- desc.m[2]
    } else {
      U1 <- U2 <- 0
      } 
    } else {
      U1 <- U2 <- 0
    }
  pp = c(U0, D0, U1, U2, D1)
  sb = as.numeric(names(pp))
  names(pp) <- c("U0", "D0", "U1", "U2", "D1")
  if(R1) names(pp)[5] <- "R1"
  if(any(ct$sb$bmax!=Inf)){
    tb.sb <- table(pars$sb[-j])
    nm.sb <- as.numeric(names(tb.sb))
    full <- nm.sb[ct$sb$bmax[nm.sb] <= tb.sb]
    if(length(full)>0) pp[sb %in% full] <- 0
  }
  return(list(sb=sb,pp=pp))
}

.slide <- function(pars, cache, d, ct, move=NULL){
  map <- .pars2map(pars,cache)
  j <- sample(1:pars$k,1)
  t2 <- pars$t2[j]
  space <- .slidespace(j, pars, cache, ct, map)
  pars.new <- pars
  mv <- .sample(1:5,1,prob=space$pp/sum(space$pp))
  type <- names(space$pp)[mv]
  l <- runif(1,0,space$pp[type])
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
.splitmerge <- function(pars, cache, d, ct, move=NULL){
  nbranch <- length(cache$edge.length)
  TH <- sum(cache$edge.length)
  v <- runif(1)
  sb.max <- ct$sb$bmax
  sb.taken <- rep(0,2*cache$ntips-2)
  sb.table <- table(pars$sb)
  sb.taken[as.numeric(names(sb.table))] <- sb.table
  sb.prob <- ct$sb$prob
  sb.prob[sb.max <= sb.taken] <- 0
  if(v < ct$bk[pars$ntheta]/(ct$bk[pars$ntheta]+ct$dk[pars$ntheta])){
    decision <- "birth"
    sb.j <- sample(1:(2*cache$ntips-2),1,prob=sb.prob)
    loc.j <- runif(1,min=0,max=cache$edge.length[sb.j])
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
    u <- runif(1,-0.5,0.5)*d
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
