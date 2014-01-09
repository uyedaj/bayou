\name{make.powerposteriorFn}
\alias{make.powerposteriorFn}
\title{Makes a power posterior function in bayou}
\usage{
  make.powerposteriorFn(k, Bk, priorFn, refFn)
}
\arguments{
  \item{k}{The step in the sequence being estimated}

  \item{Bk}{The sequence of steps to be taken from the
  reference function to the posterior}

  \item{priorFn}{The prior function to be used in marginal
  likelihood estimation}

  \item{refFn}{The reference function generated using
  \code{make.refFn()} from a preexisting mcmc chain}
}
\value{
  A function of class "powerposteriorFn" that returns a
  list of four values: \code{result} (the log density of
  the power posterior), \code{lik} (the log likelihood),
  \code{prior} (the log prior), \code{ref} the log
  reference density.
}
\description{
  This function generates a power posterior function for
  estimation of marginal likelihood using the stepping
  stone method
}
\details{
  For use in stepping stone estimation of the marginal
  likelihood using the method of Fan et al. (2011).
}

