#context("data can be loaded")
testthat::test_that("data can be loaded", {
  data(chelonia)
  expect_that(length(chelonia$phy$tip.label),equals(226))
})
