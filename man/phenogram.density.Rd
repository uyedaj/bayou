\name{phenogram.density}
\alias{phenogram.density}
\title{Plot a pheongram with the posterior density for optima values}
\usage{
  phenogram.density(tree, dat, burnin = 0, chain,
    colors = NULL, pp.cutoff = NULL, K = NULL, ...)
}
\arguments{
  \item{tree}{A phylogeny of class 'phylo'}

  \item{dat}{A named vector of tip data}

  \item{burnin}{The initial proportion of the MCMC to be
  discarded}

  \item{chain}{A bayouMCMC object that contains the results
  of an MCMC chain}

  \item{colors}{An optional named vector of colors to
  assign to regimes, \code{NULL} results in no regimes
  being plotted.}

  \item{pp.cutoff}{The posterior probability cutoff value.
  Branches with posterior probabilities of having a shift
  above this value will have the average location of the
  regime shift painted onto the branches.}

  \item{K}{A list with the values of K to be plotted. If
  \code{NULL} all values of K are combined and a total
  posterior produced. This allows separate lines to be
  plotted for different numbers of shifts so that the
  location of optima can be compared, for example, between
  all samples that have 1 vs. 2 shifts in the posterior.}

  \item{...}{Additional parameters passed to
  \code{phenogram(...)}}
}
\description{
  Plots a phenogram and the posterior density for optima
  values
}

