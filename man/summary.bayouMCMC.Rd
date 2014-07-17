\name{summary.bayouMCMC}
\alias{summary.bayouMCMC}
\title{S3 method for summarizing bayouMCMC objects}
\usage{
  \method{summary}{bayouMCMC} (object, ...)
}
\arguments{
  \item{object}{A bayouMCMC object}

  \item{...}{Additional arguments passed to \code{print}}
}
\value{
  An invisible list with two elements: \code{statistics}
  which provides summary statistics for a bayouMCMC chain,
  and \code{branch.posteriors} which summarizes branch
  specific data from a bayouMCMC chain.
}
\description{
  S3 method for summarizing bayouMCMC objects
}

