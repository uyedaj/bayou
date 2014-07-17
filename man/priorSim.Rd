\name{priorSim}
\alias{priorSim}
\title{Simulates parameters from bayou models}
\usage{
  priorSim(prior, tree, plot = TRUE, nsim = 1, ...)
}
\arguments{
  \item{prior}{A prior function created by
  \code{bayou::make.prior}}

  \item{tree}{A tree of class 'phylo'}

  \item{plot}{A logical indicating whether the simulated
  parameters should be plotted}

  \item{nsim}{The number of parameter sets to be simulated}

  \item{...}{Parameters passed on to
  \code{plotSimmap(...)}}
}
\value{
  A list of bayou parameter lists
}
\description{
  \code{priorSim} Simulates parameters from the prior
  distribution specified by \code{make.prior}
}

