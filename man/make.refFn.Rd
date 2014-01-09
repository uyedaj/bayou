\name{make.refFn}
\alias{make.refFn}
\title{Make a reference function in bayou}
\usage{
  make.refFn(chain, prior, burnin = 0.3, plot = TRUE)
}
\arguments{
  \item{chain}{An mcmc chain produced by
  \code{bayou.mcmc()} and loaded with \code{load.bayou()}}

  \item{prior}{The prior function used to generate the mcmc
  chain}

  \item{burnin}{The proportion of the mcmc chain to be
  discarded when generating the reference function}

  \item{plot}{Logical indicating whether or not a plot
  should be created}
}
\value{
  Returns a reference function of class "refFn" that takes
  a parameter list and returns the log density given the
  reference distribution. If \code{plot=TRUE}, a plot is
  produced showing the density of variable parameters and
  the fitted distribution from the reference function (in
  red).
}
\description{
  This function generates a reference function from a mcmc
  chain for use in marginal likelihood estimation.
}
\details{
  Distributions are fit to each mcmc chain and the
  best-fitting distribution is chosen as the reference
  distribution for that parameter using the method of Fan
  et al. (2011). For positive continuous parameters
  \code{alpha, sigma^2, halflife, Vy, w2, Ne}, Log-normal,
  exponential, gamma and weibull distributions are fit. For
  continuous distributions \code{theta}, Normal, Cauchy and
  Logistic distributions are fit. For discrete
  distributions, \code{k}, negative binomial, poisson and
  geometric distributions are fit. Best-fitting
  distributions are determined by AIC.
}

