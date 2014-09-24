context("Medusa likelihoods can be calculated")
test_that("Medusa likelihoods can be calculated and match Medusa", {
  require(devtools)
  require(geiger)
  require(testthat)
  #Make sure you're in the bayou package working directory
  #setwd("~/repos/bayou/bayou_1.0/")
  load_all()
  dat <- get(data(chelonia))
  phy <- dat$phy
  phy <- reorder(phy, "cladewise")
  richness <- NULL
  richness <- geiger:::.check.richness(phy = phy, richness = richness);
  
  ##Try lowering the threshold to 1, it breaks:
  resMedusa <- medusa(phy, richness, cut="node", threshold = 3, ncores = 1)
  #newo <- c(2,1,3)
  K <- length(resMedusa$summary$Shift.Node)-1
  ##Shuffling doesn't work...
  #newo <- sample(1:K, K, replace=FALSE)
  newo <- 1:K 
  splitNode <- as.numeric(resMedusa$summary$Shift.Node)[-1][newo] #c(234, 386, 452, 91, 282)
  r <- resMedusa$summary$r
  eps <- resMedusa$summary$epsilon#c(0.4, 0.5, 0.6)
  eps <- ifelse(is.na(eps), 0, eps)
  index <- c(0, newo)
  shiftCut <- "node"
  getSplitLikelihood(phy, splitNode, r, eps, index)
  resMedusa$model$lnLik
  sum(as.numeric(resMedusa$summary$Ln.Lik.part))
  expect_equal(getSplitLikelihood(phy, splitNode, r, eps, index), resMedusa$model$lnLik, tolerance=1e-7)
  #Plot shifts and get parameter in bayou format:
  #bparsMedusa <- medusa2bayou(phy, splitNode, r, eps, index)
})