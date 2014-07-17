\name{pars2simmap}
\alias{pars2simmap}
\title{Convert a bayou parameter list into a simmap formatted phylogeny}
\usage{
  pars2simmap(pars, tree)
}
\arguments{
  \item{pars}{A list that contains \code{sb} (a vector of
  branches with shifts), \code{loc} (a vector of shift
  locations), \code{t2} (a vector of theta indices
  indicating which theta is present after the shift).}

  \item{tree}{A tree of class 'phylo'}
}
\value{
  A list with elements: \code{tree} A simmap formatted
  tree, \code{pars} bayou formatted parameter list, and
  \code{cols} A named vector of colors that can be passed
  to \code{plotSimmap}.
}
\description{
  This function converts a bayou formatted parameter list
  specifying regime locations into a simmap formatted tree
  that can be plotted using \code{plotSimmap} from
  phytools.
}
\details{
  \code{pars2simmap} takes a list of parameters and
  converts it to simmap format
}
\examples{
tree <- reorder(sim.bdtree(n=100), "postorder")

pars <- list(k=5, sb=c(195, 196, 184, 138, 153), loc=rep(0, 5), t2=2:6)
tr <- pars2simmap(pars, tree)
plotSimmap(tr$tree, col=tr$col)
}

