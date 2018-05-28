context("weight matrix can be calculated")
test_that("weight matrix can be calculated", {
  data(chelonia)
  tree <- chelonia$phy
  dat <- chelonia$dat
  cache <- .prepare.ou.univariate(tree, dat)
  pars <- list(alpha=0.01, sig2=1, k=3, theta=c(3,4,5), ntheta=3, sb=c(411, 400, 47), loc=c(25, 17, 33), t2=c(2,3,3))
  TotExp <- exp(-cache$height*pars$alpha)
  stree <- pars2simmap(pars, tree)
  sW <- simmap.W(stree$tree, pars)
  expect_that(apply(sW, 1, sum), equals(rep(1,226)))
  bW <- bayou:::C_weightmatrix(cache, pars)$W
  expect_that(bW, equals(sW))
  expect_that(apply(bW,2,sum), equals(c(79.58365, 117.46516, 28.95119), tolerance=0.0001))
})
