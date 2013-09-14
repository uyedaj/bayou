context("weight matrix can be calculated")
test_that("weight matrix can be calculated", {
  data(chelonia.simmap)
  tree <- chelonia.simmap$tree
  dat <- chelonia.simmap$dat
  emap <- chelonia.simmap$emap
  cache <- .prepare.ou.univariate(tree, dat)
  pars <- list(alpha=0.1, sig2=1, k=16, optima=c(3,4,5,6), ntheta=4)
  expect_that(apply(simmap.W(tree, pars),1,sum),equals(rep(1,226)))
  expect_that(apply(simmap.W(cache, pars),1,sum),equals(rep(1,226)))
  TotExp <- exp(-cache$height*pars$alpha)
  expect_that(apply(.simmap.W(cache,pars,TotExp),1,sum),equals(rep(1,226)))
  expect_that(apply(edgemap.W(tree,pars,emap),1,sum),equals(rep(1,226)))
  expect_that(apply(.edgemap.W(cache,pars,emap,TotExp),1,sum),equals(rep(1,226)))
  expect_that(simmap.W(tree,pars),equals(edgemap.W(tree,pars,emap)))
})
