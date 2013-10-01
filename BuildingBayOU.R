require(devtools)
require(roxygen2)
require(testthat)
load_all()
document()
test_dir("inst/tests/")
install("~/bayOU_1.0")
phytools.funs <- c("nodeHeights","plotSimmap","phenogram")
cat(paste0("importFrom(phytools, ", paste(phytools.funs, collapse=", "), ")"),
    paste0("export(", paste(phytools.funs, collapse=", "), ")"),
    file = "NAMESPACE",
    sep = "\n")
lnL <- bm.lik(tree,dat,model="BM")

data(chelonia.simmap)
tree <- chelonia.simmap$tree
dat <- chelonia.simmap$dat
emap <- chelonia.simmap$emap
#tree <- sim.bdtree(n=1000,stop="taxa")
#tree <- reorder.phylo(tree,"postorder")
prior <- make.prior(tree,model="OU",dists=list(dk="dnbinom"),param=list(dk=list(size=1,prob=0.002),dsb=list(bmax=c(0,0,0,0,0,rep(1,444),0),prob=1)),type="pars",plot.prior=TRUE)
startpar <- priorSim(prior,tree,plot=TRUE,nsim=1,exclude.branches=NULL, ftype="off",pts=F)$pars[[1]]
mapped.tree <- pars2simmap(startpar,tree,theta=startpar$theta,root.theta=startpar$theta[1])
plotSimmap(mapped.tree$tree,col=mapped.tree$col,ftype="off",pts=F)
simdat<-dataSim(startpar,"OU",tree,map.type="pars",phenogram=TRUE,ftype="off")
cache <- .prepare.ou.univariate(tree, simdat$dat)
cache$maps <- mapped.tree$tree$maps
.parmap.W(cache, startpar)
#Testing moves:
ct <- .buildControl(startpar, prior, default.weights="OU")
map<-.pars2map(startpar,cache)
maps <- lapply(1:length(cache$edge.length), function(x){ y <- map$segs[names(map$segs)==x]; names(y) <- map$theta[names(map$theta)==x]; y })
mapped.tree$tree$maps <- maps
plotSimmap(mapped.tree$tree,col=mapped.tree$col,ftype="off",pts=F)
#system.time(for(i in 1:500){
  map<-.pars2map(startpar,cache)
#})
system.time(
  for(i in 1:50){
 prop <- .splitmerge(startpar, cache, 1, ct)$pars
 map <- .pars2map(prop, cache)
 maps <- lapply(1:length(cache$edge.length), function(x){ y <- map$segs[names(map$segs)==x]; names(y) <- map$theta[names(map$theta)==x]; y })
 mapped.tree$tree$maps <- maps
 plotSimmap(mapped.tree$tree,col=mapped.tree$col,ftype="off",pts=F)
}
  )
system.time(
#  startpar -> pars
  for(i in 1:100){
slide.prop <- .slide(pars, cache, 1, ct)
#pars <- slide.prop$pars
#map <- .pars2map(pars, cache)
#maps <- lapply(1:length(cache$edge.length), function(x){ y <- map$segs[names(map$segs)==x]; names(y) <- map$theta[names(map$theta)==x]; y })
#mapped.tree$tree$maps <- maps
#plotSimmap(mapped.tree$tree,col=mapped.tree$col,ftype="off",pts=F)
}
  )
rm(list=ls(all=TRUE))
load_all()
data(chelonia.simmap)
tree <- chelonia.simmap$tree
dat <- chelonia.simmap$dat
SE <- 0.1
cache <- .prepare.ou.univariate(tree,dat)
prior <- make.prior(tree,model="OU",param=list(dtheta=list(mean=3.5, sd=2),dk=list(lambda=15, kmax=2*226-2), dsb=list(bmax=1,prob=1)),type="pars",plot.prior=TRUE)
fit1 <- bayou.mcmc(tree,dat,SE=0,model="OU",prior,ngen=5000,samp=10,chunk=100,control=NULL,tuning=NULL,new.dir=TRUE,plot.freq=500,outname="bayou",ticker.freq=1000,tuning.int=c(0.1,0.2,0.3),startpar=NULL,moves=NULL,control.weights=NULL)
prior <- make.prior(tree,model="OU",param=list(dtheta=list(mean=3.5, sd=2),dk=list(lambda=15, kmax=2*226-2), dsb=list(bmax=Inf,prob=tree$edge.length)),type="pars",plot.prior=TRUE)
fit2 <- bayou.mcmc(tree,dat,SE=0,model="OU",prior,ngen=10000,samp=10,chunk=100,control=NULL,tuning=NULL,new.dir=TRUE,plot.freq=5000,outname="bayou",ticker.freq=1000,tuning.int=c(0.1,0.2,0.3),startpar=NULL,moves=NULL,control.weights=NULL)
chain <- load.bayou(fit2, save.Rdata=TRUE, cleanup=FALSE)
summary(chain)
