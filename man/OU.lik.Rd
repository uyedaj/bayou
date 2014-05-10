\name{OU.lik}
\alias{OU.lik}
\title{Function for calculating likelihood of an OU model in bayou using the threepoint algorithm}
\usage{
  OU.lik(pars, tree, X, SE = 0, model = "OU")
}
\arguments{
  \item{pars}{A list of parameters to calculate the
  likelihood}

  \item{cache}{A bayou cache object generated using
  .prepare.ou.univariate}

  \item{X}{A named vector giving the tip data}

  \item{SE}{A named vector or single number giving the
  standard errors of the data}

  \item{model}{Parameterization of the OU model. Either
  "OU", "QG" or "OUrepar".}
}
\value{
  A list returning the log likelihood ("loglik"), the
  weight matrix ("W"), the optima ("theta"), the residuals
  ("resid") and the expected values ("Exp").
}
\description{
  Function for calculating likelihood of an OU model in
  bayou using the threepoint algorithm
}
\details{
  This function will likely be deprecated in the near
  future, and replaced by bayou.lik. Currently, it is
  included because it is capable of handling standard
  errors.
}

