\name{OU.lik}
\alias{OU.lik}
\title{Function for calculating likelihood of an OU model in bayou using pruning algorithm or matrix inversion}
\usage{
  OU.lik(pars, tree, X, SE = 0, model = "OU",
    invert = FALSE)
}
\arguments{
  \item{pars}{A list of parameters to calculate the
  likelihood}

  \item{tree}{A phylogenetic tree of class 'phylo'}

  \item{X}{A named vector giving the tip data}

  \item{SE}{A named vector or single number giving the
  standard errors of the data}

  \item{model}{Parameterization of the OU model. Either
  "OU", "QG" or "OUrepar".}

  \item{invert}{A logical indicating whether the likelihood
  should be solved by matrix inversion, rather than the
  pruning algorithm. This is primarily present to test that
  calculation of the likelihood is correct.}
}
\value{
  A list returning the log likelihood ("loglik"), the
  weight matrix ("W"), the optima ("theta"), the residuals
  ("resid") and the expected values ("Exp").
}
\description{
  Function for calculating likelihood of an OU model in
  bayou using pruning algorithm or matrix inversion
}
\details{
  This function can be used for calculating single
  likelihoods using previously implemented methods. It is
  likely to become deprecated and replaced by
  \code{bayou.lik} in the future, which is based on
  \code{phylolm}'s threepoint algorithm, which works on
  non-ultrametric trees and is substantially faster.
}

