#' bayou-package
#' 
#' @name bayou-package
#' @aliases bayou-package bayou
#' @title Bayesian Fitting of Ornstein-Uhlenbeck Models to Phylogenies
#' @description A package for inferring adaptive evolution to phylogenetic 
#' comparative data using Bayesian reversible-jump estimation of 
#' multi-optima Ornstein-Uhlenbeck models.
#' @author Josef C Uyeda
#' @useDynLib bayou
#' @import ape geiger phytools coda Rcpp MASS mnormt fitdistrplus denstrip grDevices graphics stats assertthat foreach
#' @importFrom utils globalVariables read.table tail
NULL
# @importFrom denstrip densregion 
# @importFrom MASS mvrnorm
# @importFrom mnormt dmnorm
# @importFrom fitdistrplus fitdist
# @importFrom grDevices as.raster col2rgb heat.colors rainbow rgb
# @importFrom graphics abline curve lines locator mtext par plot plot.new points rasterImage rect text
# @importFrom stats density dmultinom dnorm dpois median model.frame model.matrix na.pass optim quantile reorder rnorm rpois runif sd setNames terms

# @importFrom assertthat validate_that
