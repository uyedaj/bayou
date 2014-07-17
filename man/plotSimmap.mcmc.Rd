\name{plotSimmap.mcmc}
\alias{plotSimmap.mcmc}
\title{Plot a phylogenetic tree with posterior probabilities from a bayouMCMC chain (function adapted from phytools' plotSimmap)}
\usage{
  plotSimmap.mcmc(tree, chain, burnin = NULL,
    colors = NULL, fsize = 1, ftype = "reg", lwd = 0.75,
    node.numbers = FALSE, mar = NULL, add = FALSE,
    offset = NULL, pp.cutoff = NULL, legend = TRUE,
    limits = NULL, circle = TRUE, circle.pch = 21,
    circle.pal = cm.colors, circle.lwd = 0.75,
    circle.alpha = 200, dash = FALSE, dash.pal = cm.colors,
    dash.cex = 0.5, dash.pch = "|", dash.lwd = 0.5,
    dash.alpha = 10, pp.labels = FALSE, pp.alpha = NULL,
    pp.cex = 0.75)
}
\arguments{
  \item{tree}{A tree of class phylo}

  \item{chain}{A bayouMCMC chain}

  \item{burnin}{The proportion of runs to be discarded}

  \item{colors}{A named vector of colors corresponding to
  the regimes in the stored simmap map in the tree}

  \item{fsize}{relative font size for tip labels}

  \item{ftype}{font type- options are "reg", "i" (italics),
  "b" (bold), or "bi" (bold-italics)}

  \item{lwd}{line width for plotting}

  \item{node.numbers}{a logical value indicating whether or
  not node numbers should be plotted}

  \item{mar}{vector containing the margins for the plot to
  be passed to par. If not specified, the default margins
  are (0.1,0.1,0.1,0.1)}

  \item{add}{a logical value indicating whether or not the
  tree should be plotted in the current or a new plot}

  \item{offset}{offset for the tip labels}

  \item{pp.cutoff}{the posterior probability above which a
  shift should be reconstructed on the tree. If this value
  is provided, it overrides any preexisting simmap
  associated with the tree.}

  \item{circle}{a logical value indicating whether or not a
  circle should be plotted at the base of the node with
  values that correspond to the posterior probability of
  having a shift.}

  \item{circle.pch}{the type of symbol used to plot at the
  node to indicate posterior probability}

  \item{circle.pal}{a palette of colors that will be used
  to color the interior of the circles. This will be varied
  over the interval proportional to the deviation that
  occurs at that shift on the phylogeny.}

  \item{circle.lwd}{the line width of the points plotted at
  the nodes}

  \item{circle.alpha}{a value between 0 and 255 that
  indicates the transparency of the circles (255 is
  completely opaque).}

  \item{dash}{a logical value indicating whether or not a
  dash (or other point) should be plotted on the branches
  at every location that a shift was present in the
  posterior chain.}

  \item{dash.pal}{Ignored for now}

  \item{dash.cex}{the relative size of dashes}

  \item{dash.pch}{the plotting symbol for dashes}

  \item{dash.lwd}{the line width for dashes}

  \item{dash.alpha}{the transparency of dashes}

  \item{pp.labels}{a logical indicating whether the
  posterior probability for each branch should be printed
  above the branch}

  \item{pp.alpha}{a logical or numeric value indicating
  transparency of posterior probability labels. If TRUE,
  then transparency is ramped from invisible (pp=0), to
  black (pp=1). If numeric, all labels are given the same
  transparency. If NULL, then no transparency is given.}

  \item{pp.cex}{the size of the posterior probability
  labels}

  \item{legend}{Logical indicating whether or not a legend
  should be produced}

  \item{limits}{Divergence values assigning the minimum and
  maximum amounts applied to the extremes of the color
  palette. Values outside the limits are assigned the most
  extreme values in the color palette function.}
}
\description{
  Plot a phylogenetic tree with posterior probabilities
  from a bayouMCMC chain (function adapted from phytools'
  plotSimmap)
}

