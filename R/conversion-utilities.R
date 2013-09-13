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
