\name{rExtremal}
\alias{rExtremal}
\title{Simulation of extremal processes}
\description{
  the function returns simulations of an Extremal (Gaussian or \emph{t}) random
  field. The Extremal Gaussian process is with Gumbel marginals and the
  Extremal \emph{t} process is with Frechet marginals.
}
\usage{
rExtremal(coord, corrmodel, fitted=NULL, grid=FALSE, model,
          numblock=300, numsim=1, param)
} 
\arguments{
  \item{coord}{A (\eqn{d \times 2}{d x 2})-matrix of coordinates of points
    where the extremes are simulated. The column size represents the two
    plane dimensions and the row size the number of points}
  \item{coormodel}{correlation model of the underlying random field
    (see \bold{Details} and \command{\link{CovarianceFct}} for a
    complete description).}
  \item{fitted}{An object of class \code{"Extremal"}; all the
    informations are taken from the list of components of the fitted
    model (see \bold{Details}).}
  \item{grid}{Logical; if \code{FALSE} (the default) the vector \code{coord},
    is interpreted as a two dimensional vector of arbitrary coordinates
    on the plane (see \bold{Details}).}
  \item{model}{String; the name of the Extremal process (see \bold{Details}).}
  \item{numblock}{The block maxima size; the numeric value is
    interpreted as the number of succesive replications of the underlying
    random field from which the maximum is computed (see
    \bold{Details}).}
  \item{numsim}{Number of independent realisations to generate.}
  \item{param}{A named list giving the values for the parameters;
    the parameters are given depending by the parametrised class of
    correlation functions (see \bold{Details} and \command{\link{CovarianceFct}}).
    An example may be \code{param=list(mean=0, nugget=0, scale=10,...)}.}
}
\details{
  Let \eqn{Y(x)}, \eqn{x \in R^d}, be a Gaussian process with finite
  mean and unit variance and isotropic correlation function \eqn{\rho(h)} with
  \eqn{h \in R}.
  Consider \eqn{Y(x)_i} with \eqn{i=1,..,n},
  \eqn{n} independent copies of it. Then the Extremal Gaussian process (\code{Extremal_g})
  is define by
  \deqn{Z(x) = max_i=1,..n (Y(x/s_n)_i - b_n ) / a_n,}
  where \eqn{a_n > 0}, \eqn{b_n \in R} are sequences of
  constants and \eqn{s_n} are scaling factor that scale the domain of
  \eqn{Y(x)}. For opportune sequences \eqn{a_n}, \eqn{b_n}, \eqn{s_n}
  and for large \eqn{n} the process \eqn{Z(x)} has univariate Gumbel marginals and
  all the \eqn{d}-dimensional marginal distribution functions (\eqn{d>1}) are known as the
  Husler-Reiss models.

  Let \eqn{Q} be a Chi-square random variable with \eqn{\nu} degree of
  freedom. Consider \eqn{Q(x)} as independent copies of \eqn{Q}
  depending by \eqn{x}. Define also \eqn{W(x)=\sqrt{Q(x) / \nu}}.
  Then the Extremal \emph{t} process (\code{Extremal_t}) is defined by
  \deqn{Z(x) = max_i=1,..n (Y(x)_i / W(x) - d_n ) / c_n,}
  where \eqn{c_n > 0} and \eqn{d_n \in R} are sequences of
  constants. For opportune sequences \eqn{c_n}, \eqn{d_n}
  and for large \eqn{n} the process \eqn{Z(x)} has univariate Frechet marginals
  with parameter \eqn{\nu} and all
  the \eqn{d}-dimensional marginal distribution functions (\eqn{d>1}) have known closed forms.
  We call the multivariate distributions, Extremal \emph{t} models. More details soon will
  be available with the user guide (see also \bold{References}).

  The number of independent replications of \eqn{Y(x)} is
  given by the parameter \code{numblock}. Thus, one realization of \eqn{Z(x)} is computed
  componetwise over \code{numblock} independent replications of
  \eqn{Y(x)}. Note, that all the \eqn{d}-dimensional marginal
  distribution functions are asymptotically justified, therefore
  the blockwise of successive replications should be as much large as
  its possible. In the other hand increasing the value of
  \code{numblock} decreases the simulation speed. A right choice is a
  tradeoff between accuracy of the asymptotic approximation and the
  simulation speed. We have found that a \code{numblock} of 500 is a
  good compromise for most of the cases.   
  
  The \eqn{d \times 2}{d x 2}-matrix of coordinates (\code{coord})
  are interpreted differently depending on the argument \code{grid}.
  In particular:
  \itemize{
    \item if \code{grid=TRUE} then \code{coord} are the coordinates of
    a `regular' grid of points on the plane. Regular means that each
    coordinate is an increasing equi-spaced sequence.
    The resulting grid is a \eqn{d \times d}{d x d} matrix.  
    \item if \code{grid=FALSE} then \code{coord}
    is interpreted as a two dimensional vector of `arbitrary' coordinates
    of the plane. Arbitrary means that each coordinate is not an
    increasing equi-spaced sequence but instead is random sequence.}
  
  Classes of correlation functions (\code{corrmodel}) for the underlying
  random fields are defined in \command{\link{CovarianceFct}}. For
  the Extremal Gaussian process the implemented models are:
  \itemize{
    \item \code{exponential};
    \item \code{gauss} (gaussian);
    \item \code{gencauchy} (generalised Cauchy);
    \item \code{stable}.
    }
    Instead for the Extremal \emph{t} process are:
    \itemize{
      \item \code{cauchy};
      \item \code{exponential};
      \item \code{gauss} (gaussian);
      \item \code{gencauchy} (generalised Cauchy);
      \item \code{stable};
      \item \code{whittlematern} (Whittle-Matern).
    }
  
    The parameters (\code{param}) are those of the selected correlation
  model (see \command{\link{CovarianceFct}}). Extremal processes work
  with correlation instead of covariance functions. Hence,
  this corresponds in the covariance parametrization to set equal to one
  the \code{variance} parameter. Parameter values need to be given as a
  named list with the following syntax: \code{param=list(mean=.,
    nugget=., scale=.,...)}, where the \code{mean} represents the mean
    of the underlying Gaussian process. Depending on the selected
    correlation model, extra parameters need eventually to be passed (\code{...}). In
    case of an Extremal \emph{t} process also the degree of freedom (\code{df}) need
    to be passed (see \bold{Examples}).  
 
}
\value{
  The returned object is a matrix of simulated values from a selected Extremal
  process. In particular depending on \code{numsim} and \code{grid}:\cr
    if \code{numsim=1} and \cr
    \code{grid=FALSE}, then an one-dimensional vector of length \eqn{d} is returned,\cr
    \code{grid=TRUE},  then a \eqn{d \times d}{d x d} matrix is returned.\cr
    
    if \code{numsim>1} and \cr
    \code{grid=FALSE}, then a \code{numsim} \eqn{\times d}{x d} matrix is returned; the rows
    contain the realisations,\cr
    \code{grid=TRUE}, then a \eqn{d \times d \times}{d x d x} \code{numsim}
    is returned; the last dimension contains the realisations.
}
\references{
  Demarta, S. and McNeil, A. J. (2005). The t copula and related copulas
 \emph{International Statistical Review}, \bold{73}, 111-129.
  
  Husler, J. and Reiss, R. (1989). Maxima of normal random vectors:
  between independence and complete dependence
  \emph{Statistics and Probability Letters}, \bold{7}, 283--286.

   Kabluchko, Z., Schlather, M., de Haan, L. (2009). Stationary
   max-stable fields associated 
   to negative definite functions. \emph{Annals of Probability}, \bold{37}, 2042-2065.

  Nikoloulopoulos, A. K., H. Joe, and H. Li (2009). Extreme value properties of multivariate t
 copulas. \emph{Extremes}, \bold{12}, 129--148.

 for Gaussian process simulations see \command{\link{RFMethods}}.
}
\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}}
%\note{}


\seealso{
  \command{\link{fExtremal}},
  \command{\link{CovarianceFct}},
  \command{\link{GaussRF}},
  \code{\link{RandomFields}}.
}

\examples{
## Simulation of Extremal processes on a `regular' grid
## of points (i.e. rectangles).
## However, simulations can be performed when points
## have arbitrary locations.

library(RandomFields)

# Definition of a grid of points:

step <- .2     # steps
numsim <- 1    # number of realizations per point
lower <- 0     # lower point of the region
upper <- 20    # upper point of the region

x <- y <- seq(lower, upper, by=step)
coord <- cbind(x, y)

colour <- rainbow(100)

##################################################################
## Example n. 1, Extremal t process.
## Stable correlation function
##
##################################################################

# Setting model parameters:

numblock <- 500        # number of replications per block

model <- "Extremal_t"  # Type of process
corrmodel <- "stable"  # Type of correlation function
mean <- 0              # Mean of the underlying Gaussian model
nugget <- 0            # Nugget effect
scale <- 5             # Scale of the stable covariance function
power <- 1.5           # Smoothing parameter
df <- 3                # Degree of freedom of the underlying t model

param <- list(df=df, mean=mean, nugget=nugget, scale=scale, power=power)

# Simulation:
simET <- rExtremal(coord, corrmodel, model=model, numsim=numsim,
                   grid = TRUE, numblock=numblock, param=param)

# Visualization:
image(x, y, simET, col=colour)


##################################################################
## Example n. 2, Extremal Gaussian process.
## Stable correlation function
##
##################################################################

# Setting model parameters:

numblock <- 500        # number of replications per block

model <- "Extremal_g"  # Type of process
corrmodel <- "stable"  # Type of correlation function
mean <- 0              # Mean of the underlying Gaussian model
nugget <- 0            # Nugget effect
scale <- 5             # Scale of the stable covariance function
power <- 1.5           # Smoothing parameter

param <- list(mean=mean, nugget=nugget, scale=scale, power=power)

# Simulation:
simEG <- rExtremal(coord, corrmodel, model=model, numsim=numsim,
                   grid = TRUE, numblock=numblock, param=param)

# Visualization:
image(x, y, simEG, col=colour)


# Visualization of both processes:
par(mfrow=c(1,2))
image(x, y, simEG, col=colour)
image(x, y, simET, col=colour)

}
\keyword{spatial}
