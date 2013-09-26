\name{make.prior}
\alias{make.prior}
\title{Make a prior function for bayOU}
\usage{
  make.prior(tree, dists = list(), param = list(),
    plot.prior = TRUE, model = "OU", type = "pars")
}
\arguments{
  \item{tree}{A tree object of class "phylo"}

  \item{dists}{A list providing the function names of the
  distribution functions describing the prior distributions
  of parameters (see details). If no distributions are
  provided for a parameter, default values are given. Note
  that the names are provided as text strings, not the
  functions themselves.}

  \item{param}{A list providing the parameter values of the
  prior distributions (see details).}

  \item{plot.prior}{A logical indicating whether the prior
  distributions should be plotted.}

  \item{model}{One of three specifications of the OU
  parameterization used. Takes values \code{"OU"} (alpha &
  sig2), \code{"QG"} (h2, P, w2, Ne), or \code{"OUrepar"}
  (halflife,Vy)}

  \item{type}{Specifies the type of input for regime
  mapping used. Takes values \code{"pars"} (regime shits
  specified in the parameter list by \code{pars$sb} and
  \code{pars$loc}), \code{"emap"} (regime shifts specified
  by a provided edge map), or \code{"simmap"} (regime
  shifts obtained from \code{tree$maps})}
}
\value{
  returns a prior function that calculates the log prior
  density for a set of parameter values provided in a list
  with correctly named values.
}
\description{
  This function generates a prior function to be used for
  bayOU according to user specifications.
}
\details{
  Default distributions and parameter values are given as
  follows: OU:
  \code{list(dists=list("dalpha"="dlnorm","dsig2"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="equal","dloc"="dunif"),
  param=list("dalpha"=list(),"dsig2"=list(),"dtheta"=list(),"dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list()))}
  QG:
  \code{list(dists=list("dh2"="dbeta","dP"="dlnorm","dw2"="dlnorm","dNe"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="equal","dloc"="dunif"),
  param=list("dh2"=list(shape1=1,shape2=1),"dP"=list(),"dw2"=list(),"dNe"=list(),"dtheta"=list(),"dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list()))}
  OUrepar:
  \code{list(dists=list("dhalflife"="dlnorm","dVy"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="equal","dloc"="dunif"),
  param=list("dhalflife"=list("meanlog"=0.25,"sdlog"=1.5),"dVy"=list("meanlog"=1,"sdlog"=2),"dk"=list(lambda=1,kmax=2*ntips-2),"dtheta"=list(),"dloc"=list(min=0,max=1)),"dsb"=list())}

  \code{dalpha, dsig2, dh2, dP, dw2, dNe, dhalflife}, and
  \code{dVy} must be positive continuous distributions and
  provide the parameters used to calculate alpha and
  sigma^2 of the OU model. \code{dtheta} must be continuous
  and describes the prior distribution of the optima. dk is
  the prior distribution for the number of shifts. For
  Poisson and conditional Poisson (cdpois) are provided the
  parameter \code{lambda}, which provides the total number
  of shifts expected on the tree (not the rate per unit
  branch length). Otherwise, \code{dk} can take any
  positive, discrete distribution. dsb indicates the prior
  probability of a given set of branches having shifts.
  Three options are possible:
  \code{"equal"},\code{"equalnotips"} or \code{"free"}.
  \code{"equal"} assigns equal probability to each branch
  sampled without replacement. \code{"equalnotips"} but
  restricts shifts to only internal branches. \code{"free"}
  allows multiple shifts per branch, and determines the
  prior probability of a set of shift locations using a
  multinomial distribution with probabilities proportional
  to branch lengths taken from the tree. \code{"dloc"}
  indicates the prior probability of a shift on a given
  branch. Currently, all locations are given equal density
  of 1. All distributions are set to return log-transformed
  probability densities.
}
\examples{
## Load data
data(chelonia.simmap)
tree <- chelonia.simmap$tree
dat <- chelonia.simmap$dat

##Make prior function
prior <- make.prior(tree,dists=list(dalpha="dunif",dsig2="dunif"),param=list(dalpha=list(min=0,max=10),dsig2=list(min=0,max=10)),type="pars")

##Define a set of parameter values
emap <- chelonia.simmap$emap
pars <- list(alpha=0.1, sig2=1, k=16, theta=c(3,4,5,6), sb=which(emap$sh==1), loc=emap$r1[which(emap$sh==1)])

##Calculate prior probability
prior(pars)

##Return a list of functions used to calculate prior
attributes(prior)$functions

##Return parameter values used in prior distribution
attributes(prior)$parameters

##Alternative parameterization using phylogenetic half-life and stationary variance
pars <- list(halflife=20, Vy=1, k=16, theta=c(3,4,5,6), sb=which(emap$sh==1), loc=emap$r1[which(emap$sh==1)])
prior2 <- make.prior(tree,dists=list(dhalflife="dlnorm",dVy="dlnorm"),param=list(dhalflife=list(meanlog=3,sdlog=1),dVy=list(meanlog=0,sdlog=2)),model="OUrepar",type="pars")
prior2(pars)

##Specify different types of inputs
prior.pars <- make.prior(tree, model="OUrepar", type="pars")
prior.emap <- make.prior(tree, model="OUrepar", type="emap")
prior.simmap <- make.prior(tree, model="OUrepar", type="simmap")

prior.pars(pars)
prior.emap(pars,emap=emap)
prior.simmap(pars,tree)
}

