
###Internal function for calculating the allowable space surrounding a shift for a slide move. Moves can only occur within a branch, or to immediately adjacent neighboring branches. Furthermore, a shift cannot be moved past another preexisting shift, and is also constrained by the tips and the root of the tree. If a shift moves within a branch, the proposal ratio is 1. If it shifts from one branch to a neighboring branch, however, it will change because the allowable space for that move will change.
.slidespace.dep <- function(j,cache,pars){
  m <- cache$maps[[pars$sb[j]]]
  t2 <- pars$t2[j]
  nodes <- cache$edge[pars$sb[j],]
  U0 <- unname(m[names(m)==t2])
  D0 <- unname(m[which(names(m)==t2)-1])
  if(which(names(m)==t2)==length(m) & nodes[2] > cache$ntips){
    sbU <- which(cache$edge[,1]==nodes[2])
    mU <- cache$maps[sbU]
    U1 <- unname(mU[[1]][1])
    U2 <- unname(mU[[2]][1])
  } else {
    sbU <- c(0,0)
    U1 <- 0
    U2 <- 0
  }
  if(which(names(m)==t2)<3 & nodes[1]!=cache$ntips+1){
    sbD <- unname(which(tree$edge[,2]==nodes[1]))
    mD <- cache$maps[[sbD]]
    D1 <- unname(mD[length(mD)])
  } else {D1 <- 0; sbD <- 0}
  sbj <- c(U0=pars$sb[j],D0=pars$sb[j],U1=sbU[1],U2=sbU[2],D1=sbD)
  return(list(sb=sbj,pp=c(U0=U0,D0=D0,U1=U1,U2=U2,D1=D1)))
}

.slideCPP <- function(cache,pars,d,move=NULL,ct=NULL){
  j <- sample(1:pars$k,1)
  t2 <- pars$t2[j]
  space <- .slidespace(j,cache,pars)
  mv <- sample(1:5,1,prob=space$pp/sum(space$pp))
  type <- names(space$sb)[mv]
  l <- runif(1,0,space$pp[type])
  maps.new <- cache$maps
  pars.new <- pars
  if(type=="U0"){
    m <- cache$maps[[space$sb[mv]]]
    m[which(names(m)==t2)-1] <- m[which(names(m)==t2)-1] + (m[names(m)==t2]-l)
    m[names(m)==t2] <- l
    maps.new[[space$sb[mv]]] <- m
    cs.m <- cumsum(m)
    pars.new$loc[j] <- unname(cs.m[which(names(cs.m)==t2)-1])
  } else {
    if(type=="D0"){
      m <- cache$maps[[space$sb[mv]]]
      m[names(m)==t2] <- (m[which(names(m)==t2)-1]-l)+m[names(m)==t2]
      m[which(names(m)==t2)-1] <- l
      maps.new[[space$sb[mv]]] <- m
      cs.m <- cumsum(m)
      pars.new$loc[j] <- unname(cs.m[which(names(cs.m)==t2)-1])
    } else {
      pars.new$sb[j] <- space$sb[mv]
      if(type %in% c("U1","U2")){
        m0 <- cache$maps[[pars$sb[j]]]
        t1 <- names(m0)[which(names(m0)==t2)-1]
        m0[length(m0)-1] <- m0[length(m0)-1]+m0[length(m0)]
        m0 <- m0[-length(m0)]
        m <- cache$maps[[space$sb[mv]]]
        m <- c(l,m[1]-l,m[-1])
        names(m)[1] <- t1
        maps.new[[space$sb[mv]]] <- m
        maps.new[[pars$sb[j]]] <- m0
        cs.m <- cumsum(m)
        pars.new$loc[j] <- unname(cs.m[which(names(cs.m)==t2)-1])
      } else {
        if(type=="D1"){
          m0 <- cache$maps[[pars$sb[j]]]
          t1 <- names(m0)[which(names(m0)==t2)-1]
          m0[2] <- m0[2]+m0[1]
          m0 <- m0[-1]
          m <- cache$maps[[space$sb[mv]]]
          m <- c(m[-length(m)],m[length(m)]-l,l)
          names(m)[length(m)] <- t2
          maps.new[[space$sb[mv]]] <- m
          maps.new[[pars$sb[j]]] <- m0
          cs.m <- cumsum(m)
          pars.new$loc[j] <- unname(cs.m[which(names(cs.m)==t2)-1])
        }
      }
      maps.new <- pars2simmap(pars.new,cache$phy,sim.theta=FALSE,theta=pars.new$theta)$tree$maps
    }
  }
  cache$maps <- maps.new
  new.space <- .slidespace(j,cache,pars.new)
  hr <- log(1/sum(new.space$pp))-log(1/(sum(space$pp)))
  return(list(pars=pars.new,cache=cache,hr=hr))
}




