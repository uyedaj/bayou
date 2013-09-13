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