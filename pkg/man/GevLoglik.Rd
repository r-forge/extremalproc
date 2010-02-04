\name{GevLoglik}
\alias{GevLoglik}
\title{Generalised Extreme Value log-Likelihood}
\description{
  the function returns the value of the log-likelihood for the
  generalised extreme value distribution.
}
\usage{
GevLoglik(data, param)
} 
\arguments{
  \item{data}{A vector with the observed extreme values.}
  \item{param}{A vector with model parameters. The order have to be the
    following: location, scale, shape.}
}

\value{
  The returned value is the log-likelihood computed for the data and
  parameter values passed. 
}


\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}}


\keyword{spatial}
