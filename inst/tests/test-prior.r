context("testing prior functions")

test_that("testing prior functions", {
  data(chelonia)
  tree <- chelonia$phy
  dat <- chelonia$dat
  cache <- .prepare.ou.univariate(tree, dat)
  pars <- list(alpha=0.1, sig2=1, k=16, theta=c(3,4,5,6), sb=c(411,400,47), loc=c(23, 21, 33))
  QGpars <- list(h2=0.1,P=1,w2=0.9,Ne=1,k=16,theta=c(3,4,5,6), sb=c(411,400,47), loc=c(23, 21, 33))
  prior <- make.prior(tree,dists=list(dalpha="dunif",dsig2="dunif"),param=list(dalpha=list(min=0,max=1),dsig2=list(min=0,max=1),dsb=list(bmax=1,prob=1)),plot.prior=FALSE)
  QGprior <- make.prior(tree,dists=list(dh2="dunif",dP="dunif",dw2="dunif",dNe="dunif"),param=list(dh2=list(min=0,max=1),dP=list(min=0,max=1),dw2=list(min=0,max=1),dNe=list(min=0,max=1)),model="QG", plot.prior=FALSE)
  expect_that(class(try(make.prior(tree, plot.prior=FALSE),silent=TRUE))[1],equals("priorFn"))
  expect_that(class(try(make.prior(tree,model="QG", plot.prior=FALSE),silent=TRUE))[1],equals("priorFn"))
  expect_that(class(try(make.prior(tree,model="OUrepar",plot.prior=FALSE),silent=TRUE))[1],equals("priorFn"))
  expect_that(prior(pars),equals(-94.87692,tolerance=0.0001))
  expect_that(QGprior(QGpars),equals(prior(pars)))
  f1tree <- tree
  f1tree$edge.length[pars$sb] <- tree$edge.length[pars$sb]*100
  priorf1 <- make.prior(tree,dists=list(dalpha="dunif",dsig2="dunif"),param=list(dalpha=list(min=0,max=1),dsig2=list(min=0,max=1),dsb=list(bmax=Inf,prob=tree$edge.length)),plot.prior=FALSE)
  priorf2 <- make.prior(f1tree,dists=list(dalpha="dunif",dsig2="dunif"),param=list(dalpha=list(min=0,max=1),dsig2=list(min=0,max=1),dsb=list(bmax=Inf,prob=f1tree$edge.length)),plot.prior=FALSE)
  expect_that(priorf1(pars)< priorf2(pars), is_true())  
  expect_that(QGprior(pars),throws_error("Missing parameters:  h2 P w2 Ne"))
})
