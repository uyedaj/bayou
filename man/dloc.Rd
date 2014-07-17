\name{dloc}
\alias{dloc}
\alias{rloc}
\title{Probability density function for the location of the shift along the branch}
\usage{
  dloc(loc, min = 0, max = 1, log = TRUE)

  rloc(k, min = 0, max = 1)
}
\arguments{
  \item{loc}{The location of the shift along the branch}

  \item{min}{The minimum position on the branch the shift
  can take}

  \item{max}{The maximum position on the branch the shift
  can take}

  \item{log}{A logical indicating whether the log density
  should be returned}

  \item{k}{The number of shifts to return along a branch}
}
\description{
  Since unequal probabilities are incorporated in
  calculating the density via \code{dsb}, all branches are
  assumed to be of unit length. Thus, the \code{dloc}
  function simply returns 0 if \code{log=TRUE} and 1 if
  \code{log=FALSE}.
}
\details{
  \code{dloc} calculates the probability of a shift
  occuring at a given location along the branch assuming a
  uniform distribution of unit length \code{rloc} randomly
  generates the location of a shift along the branch
}

