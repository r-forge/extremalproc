\name{fGev}
\alias{fGev}
\title{Maximum-likelihood Fitting of the Generalized Extreme Value
  Distribution}
\description{
  the function returns the parameters' estimates and the variances of
  the estimates for the Generalized Extreme Value Distribution.
}
\usage{
fGev(data, method='Nelder-Mead', start, std.err=TRUE)
} 
\arguments{
  \item{data}{A vector of extreme values.}
  \item{method}{The optimization method (see \code{\link{optim}} for
    details). 'Nelder-Mead' is the default.}
  \item{start}{A named list with the initial values for the
    parameters over which the likelihood is to be maximized.}
  \item{std.err}{Logical; if \code{TRUE} (the default), the standard
    errors are returned.}
}

\details{
If \code{start} is omitted the routine is computing the
    starting values using moment estimators.
}

\value{
  The returned object is a list with:
  \itemize{
    \item \code{param}, the vector of  parameters' estimates;
    \item \code{varcov} the matrix of the variance-covariance of
    the estimates.
    }
}

\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.
}

\seealso{\code{\link{GevLoglik}}, \code{\link{optim}}}


\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}}


\keyword{spatial}
