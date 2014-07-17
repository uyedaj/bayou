\name{regime.plot}
\alias{regime.plot}
\title{Adds visualization of regimes to a plot}
\usage{
  regime.plot(pars, tree, cols, type = "rect",
    transparency = 100)
}
\arguments{
  \item{pars}{A bayou formatted parameter list}

  \item{tree}{A tree of class 'phylo'}

  \item{cols}{A vector of colors to give to regimes, in the
  same order as pars$sb}

  \item{type}{Either "rect", "density" or "lines". "rect"
  plots a rectangle for the 95% CI for the stationary
  distribution of a regime. "density" varies the
  transparency of the rectangles according to the
  probability density from the stationary distribution.
  "lines" plots lines for the mean and 95% CI's without
  filling them.}

  \item{transparency}{The alpha transparency value for the
  maximum density, max value is 255.}
}
\description{
  Adds visualization of regimes to a plot
}

