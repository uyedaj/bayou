\name{pull.pars}
\alias{pull.pars}
\title{Utility function for retrieving parameters from an MCMC chain}
\usage{
  pull.pars(i, chain, model = "OU")
}
\arguments{
  \item{i}{An integer giving the sample to retrieve}

  \item{chain}{A bayouMCMC chain}

  \item{model}{The parameterization used, either "OU", "QG"
  or "OUrepar"}
}
\value{
  A bayou formatted parameter list
}
\description{
  Utility function for retrieving parameters from an MCMC
  chain
}
\examples{
\dontrun{
tree <- sim.bdtree(n=30)
tree$edge.length <- tree$edge.length/max(branching.times(tree))
prior <- make.prior(tree, dists=list(dk="cdpois", dsig2="dnorm", dtheta="dnorm"),
             param=list(dk=list(lambda=15, kmax=32), dsig2=list(mean=1, sd=0.01),
                 dtheta=list(mean=0, sd=3)), plot.prior=FALSE)
pars <- priorSim(prior, tree, plot=FALSE, nsim=1)$pars[[1]]
dat <- dataSim(pars, model="OU", phenogram=FALSE, tree)$dat
fit <- bayou.mcmc(tree, dat, model="OU", prior=prior, new.dir=TRUE, ngen=5000, plot.freq=NULL)
chain <- load.bayou(fit, save.Rdata=TRUE, cleanup=TRUE)
plotBayoupars(pull.pars(300, chain), tree)
}
}

