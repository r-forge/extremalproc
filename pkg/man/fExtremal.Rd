\name{fExtremal}
\alias{fExtremal}
\title{Maximum Composite Likelihood Fitting
  of Extremal Processes}
\description{
  The function estimates parameters of
  extremal processes (Gaussian and \emph{t}) with composite (Pseudo) likelihood method.
 Model parameters are those of the undelrying correlation function. The function allows  any parameters to be held fixed if desired.
  }
\usage{
fExtremal(coord, corrmodel, data, fixed=list(nugget=0),
          method='L-BFGS-B', model, parscale=TRUE, start=NULL)
}

\arguments{
  \item{coord}{A (\eqn{n\times2}{n x 2})-matrix of coordinates (where n is the numbe of per each coordinate), or vector (\eqn{n\times1}{n x 1}) of coordinates.}
  \item{coormodel}{String; the name of a correlation model, see for a complete description Section \bold{Details} and \command{\link{CovarianceFct}} from the package \code{\link{RandomFields}}.}
  \item{data}{Matrix of values measured at coordinates \code{coord};
    columns represent the data collected at the different locations, rows represent the data replications observed at the sites for different time steps (see \bold{Details}).}
  \item{fixed}{A named list giving the values for the parameters to be considered fixed. For the selected correlation function the listed parameters will be not estimated. If \code{list(nugget=0)} (the default), the nugget effect is ignored.}
  \item{method}{String; the name of the optimization method. If \code{L-BFGS-B} (the default) uses the \emph{box constraints} version of a quasi-Newton method (see \command{\link{optim}} for details).} 
  \item{model}{String; the name of a Extremal process. The options are: \code{Extremal_g} for the Extremal Gaussian process and \code{Extremal_t} for the Extremal \emph{t} process, see \bold{Details}.}
  \item{parscale}{A named list giving values for the scale parameters over which the composite likelihood is to be maximized. If \code{TRUE} (the default) the routine attempts to set opportune scaling values (see \command{\link{optim}} for details).}
  \item{start}{A named list giving the initial values for the parameters over which the composite likelihood is to be maximized. If \code{NULL} (the default) the routine attempts to find good starting values.}
} 
  

\details{
  The implemented Extremal processes are the Gaussian and the (Student)
  \emph{t}, see \command{\link{rExtremal}} for details. With the argument  
  \code{model} should be given one of the options: \code{Extremal_g} for the former process and
  \code{Extremal_t} for the later.

  With an Extremal processes a correlation function associated with the undelrying
  stochastic process (Gaussian or \emph{t}) need to be defined, see 
  \command{\link{CovarianceFct}}
  for details. The implemented correlation models are: 
  \code{exponential}, \code{gauss} (Gaussian), 
  \code{gencauchy} (generalised Cauchy) and \code{stable} for the Extremal Gaussian 
  process. Instead for the Extremal \emph{t} are: \code{cauchy}, \code{exponential}, 
  \code{gauss} (Gaussian), \code{gencauchy} (generalised Cauchy), \code{stable} and 
  \code{whittlematern} (Whittle-Matern).  

  The argument \code{data} expects extreme observations formatted as
  a \eqn{n\times d}{n x d} matrix. The \emph{d} columns are interpreted
  as extreme values observed at specific locations of some region.
  The \emph{n} rows are interpreted as \emph{n} independent and indentically distribuited
  extreme replications over time. For instance, the raw data (i.e. daily
  rainfall, temperatures, ect.) should be aggregated temporarly in blocks of
  successive observations (monthly, quarterly, yearly, etc) for then compute
  the maximum in each block. In this way one obtain a vector of block maxima
  per each site, for example if \emph{n} years of data are available one
  can compute \emph{d} vector of \emph{n} annual maxima. Note, that for
  theoretical reasons each marginal distribution (column) is assumed to
  belong to the Generalised Extreme Value class. 
  
  The optimisations are performed using
  \command{\link[stats]{optim}}. One can choose different optimization
  algorithms setting the argument \code{method} and
  provide by \code{parscale} a vector of scaling values associated with
  the vector of parameters (see \command{\link{optim}}). The latter if
  is set to \code{TRUE} then the routine attempts to set opportune
  scaling values, if \code{FALSE} then the routine set the scaling
  values all equal to one. Otherwise an arbitrary vector of scaling
  values can be given as a list.   

}
\value{
  Returns an object of class \code{"Extremal"}.

  An object of class \code{"Extremal"} is a list
  containing at most the following components (in alphabetical order):
  
  \item{CLIC}{The composite likelihood information criterio at the maximum composite likelihood estimates.}
  \item{coord}{The coordinates where the extremes are observed.}
  \item{convergence}{Response returned by the \code{\link{optim}}.}
  \item{corrmodel}{The correlation function of the underlying process.}
  \item{data}{The extreme observations passed as argument \code{data}.}
  \item{fixed}{A vector containing the fixed parameter.}
  \item{iterations}{Number of iterations used by the maximizer algorithm (see \code{\link{optim}}).}
  \item{logLik}{The composite log-likelihood value at the maximum.}
  \item{message}{Error message returned by \code{\link{optim}}.}
  \item{model}{The type of Extremal process.}
  \item{param}{A vector containing all optimized parameters.}
  \item{std.err}{A vector containing the standard errors.}
  \item{varcov}{The variance covariance matrix.}
}
\references{
Davison, A. C. and Gholamrezaee, M. (2009), Geostatistics of extremes. \emph{Journal of the Royal Statistical
Society, Series B}. To appear.

Davison, A. C., Padoan, S. A. and Ribatet, M. (2010). Modelling of Spatial Extremes. \emph{Submitted}.

Padoan, S. A. (2008) Computational Methods for Complex Problems in Extreme Value Theory. PhD Thesis. University of Padova.
\url{http://paduaresearch.cab.unipd.it/1047/}

Padoan, S. A. Ribatet, M and Sisson, S. A. (2010). Likelihood-Based Inference for
Max-Stable Processes. \emph{Journal of the American Statistical Association, Theory & Methods}.
In press.
}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}}

\seealso{
  \command{\link{CovarianceFct}},
  \command{\link{rExtremal}}.
}
\examples{
## Maximum-Composite likelihood fitting for Extremal processes.
## The fitting procedure is implemented for underlying 
## Gaussian and t processes.
## In both cases the parameters are those of the underlying
## correlation function and plus the degree of freedom for the
## t process.

library(RandomFields)

## Define first a (irregular) grid of points where to simulate 
## synthetic replication from the process.

# Setting simulation design:
numcoord <- 30  # number of points
numsim <- 60    # number of observations per point
lower <- 0      # lower point simulation region
upper <- 40     # upper point simulation region

# Points randomly located:
x <- runif(numcoord, lower, upper)
y <- runif(numcoord, lower, upper)
coord <- cbind(x, y)


##################################################################
## Example n. 1, Extremal t process.
## Stable correlation function.
##
##################################################################

# Setting process parameters

numblock <- 1000       # number of observations per block

model <- 'Extremal_t'  # Type of Extremal process
corrmodel <- 'stable'  # Type of correlation function
mean <- 0              # Location of the underlying t model
nugget <- 0            # Nugget effect parameter
scale <- 20            # Scale of the stable covariance function
power <- 1.5           # Power of the stable correlation function
df <- 3                # Degree of freedom of the underlying t model

param <- list(df = df, mean=mean, nugget=nugget,
              scale=scale, power=power)

# Simulation of the Extremal t process in the specified points:
sim <- rExtremal(coord, corrmodel, model=model, numsim=numsim,
                 numblock=numblock, param=param)

# Fitting Extremal t process by composite likelihood:
fit <- fExtremal(coord, corrmodel, sim, model=model)

# Results:
summary(fit)


##################################################################
## Example n. 2, Extremal Gaussian process.
## Stable correlation function
##
##################################################################

# Setting process parameters

numblock <- 1000       # number of observations per block

model <- 'Extremal_g'  # Type of Extremal process
corrmodel <- 'stable'  # Type of correlation function
mean <- 0              # Mean of the underlying Gaussian model
nugget <- 0            # Nugget effect parameter
scale <- 40            # Scale of the stable correlation function
power <- 1.8           # Power of the stable correlation function

param <- list(mean=mean, nugget=nugget, scale=scale, power=power)

# Simulation of the Extremal Gaussian process in the specified points:
sim <- rExtremal(coord, corrmodel, model=model, numsim=numsim,
                 numblock=numblock, param=param)

# Fitting Extremal Gaussian process by composite likelihood:
fit <- fExtremal(coord, corrmodel, sim, model=model)

# Results:
summary(fit)


}
\keyword{spatial}
