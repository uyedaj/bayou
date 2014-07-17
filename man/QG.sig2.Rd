\name{QG.sig2}
\alias{QG.sig2}
\title{Calculates the sigma^2 parameter from a QG model}
\usage{
  QG.sig2(pars)
}
\arguments{
  \item{pars}{A bayou formatted parameter list with
  parameters h2 (heritability), P (phenotypic variance) and
  Ne (Effective population size)}
}
\value{
  An sig2 value according to the equation \code{alpha =
  h2*P/(Ne)}.
}
\description{
  Calculates the sigma^2 parameter from a QG model
}

