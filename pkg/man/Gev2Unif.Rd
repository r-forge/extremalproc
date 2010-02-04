\name{Gev2Unif}
\alias{Gev2Unif}
\title{Generalised Extreme Value To Uniform Distribution}
\description{
  the function transform the data with generalized extreme value (gev)
  to data with uniform distribution.
}
\usage{
Gev2Unif(data)
} 
\arguments{
  \item{data}{A vector or a (\eqn{n\times d}{n x d})-matrix of extreme
    values. Columns represent the \eqn{d} different variables, rows represent
    the \eqn{n} observations per each variable.}
}

\details{
The \code{data} (eventually each column) is assumed to be gev
distributed. The routine estimates the parameters of the gev
distribution by Maximum-likelihood. Then the \code{data} are
transform using the gev cumulative distribution function
evaluated at the estimated parameters.
}

\value{
  The returned object is a vector or a \eqn{n\times d}{n x d})-matrix of
  uniformly distributed data.
}

\seealso{\code{\link{GevLoglik}}, \code{\link{fGev}},
  \code{\link{optim}}}


\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}}


\keyword{spatial}
