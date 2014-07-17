\name{dhalfcauchy}
\alias{dhalfcauchy}
\alias{phalfcauchy}
\alias{qhalfcauchy}
\alias{rhalfcauchy}
\title{Half cauchy distribution taken from the R package LaplacesDemon (Hall, 2012).}
\usage{
  dhalfcauchy(x, scale = 25, log = FALSE)

  phalfcauchy(q, scale = 25)

  qhalfcauchy(p, scale = 25)

  rhalfcauchy(n, scale = 25)
}
\arguments{
  \item{x}{A parameter value for which the density should
  be calculated}

  \item{scale}{The scale parameter of the half-Cauchy
  distributoin}

  \item{log}{A logical indicating whether the log density
  should be returned}

  \item{q}{A vector of quantiles}

  \item{p}{A vector of probabilities}

  \item{n}{The number of observations}
}
\description{
  \code{dhalfcauchy} returns the probability density for a
  half-Cauchy distribution
}

