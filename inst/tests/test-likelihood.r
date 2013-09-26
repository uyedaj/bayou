context("can calculate likelihoods")
test_that("can calculate likelihoods", {
  data(chelonia.simmap)
  tree <- chelonia.simmap$tree
  dat <- chelonia.simmap$dat
  emap <- chelonia.simmap$emap
  cache <- .prepare.ou.univariate(tree, dat)
  pars <- list(alpha=0.1, sig2=1, k=16, theta=c(3,4,5,6), ntheta=4)
  QGpars <- list(h2=1,P=1,w2=9,Ne=1,k=16,theta=c(3,4,5,6), ntheta=4)
  TotExp <- exp(-cache$height*pars$alpha)
  expect_that(is.finite(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU")$loglik[1]),is_true())
  expect_that(emOU.lik(pars,emap,tree,dat,SE=0.1,model="OU")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  expect_that(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  expect_that(.emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  expect_that(emOU.lik(QGpars,emap,cache,dat,SE=0.1,model="QG")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  expect_that(.emOU.lik(QGpars,emap,cache,dat,SE=0.1,model="QG")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  expect_that(smOU.lik(pars,tree,dat,SE=0.1,model="OU")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  expect_that(smOU.lik(QGpars,tree,dat,SE=0.1,model="QG")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  expect_that(.smOU.lik(pars,cache,dat,SE=0.1,model="OU")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  expect_that(.smOU.lik(QGpars,cache,dat,SE=0.1,model="QG")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  pars$sb <- which(emap$sh==1)
  pars$loc <- emap$r1[emap$sh==1]
  pars$t2 <- emap$t2[emap$sh==1]
  QGpars$sb <- which(emap$sh==1)
  QGpars$loc <- emap$r1[emap$sh==1]
  QGpars$t2 <- emap$t2[emap$sh==1]
  expect_that(.OU.lik(pars,cache,dat,SE=0.1,model="OU")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  expect_that(OU.lik(pars,cache,dat,SE=0.1,model="OU")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  expect_that(.OU.lik(QGpars,cache,dat,SE=0.1,model="QG")$loglik[1],
              equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  })
