#context("data can be loaded")
testthat::test_that("data can be loaded", {
  data(chelonia, package="geiger")
  testthat::expect_equal(length(chelonia$phy$tip.label),226)
})
