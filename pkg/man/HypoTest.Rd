\name{HypoTest}
\alias{HypoTest}
\title{Statistical Hypothesis Tests for Nested Models}
\description{
  The function performs statistical hypothesis tests for nested models
  based on composite likelihood versions of: Wald-type, score-type and
  likelihood ratio (Wilks-type) statistics.
}

\usage{
HypoTest(object1, object2, \dots, statistic)
}

\arguments{
  \item{object1}{An object of class \code{'extremal'}.}
  \item{object2}{An object of class \code{'extremal'} that is a nested
    model within \code{object1}.}
  \item{\dots}{Further successively nested objects.}
  \item{statistic}{String; the name of the statistic used within the
  hypothesis test (see \bold{Details}).}
} 
  

\details{
 Consider an \eqn{d}-dimensional random vector \eqn{\mathbf{Y}}{Y} with
 probability density function \eqn{f(\mathbf{y};\mathbf{\theta})}{f(y;theta)},
where \eqn{\mathbf{\theta} \in \Theta}{theta in Theta} is a
\eqn{p}-dimensional vector of parameters. Suppose that
\eqn{\mathbf{\theta}=(\mathbf{\psi},\mathbf{\tau})}{theta=(psi, tau)}
 can be partitioned in a \eqn{q}-dimensional subvector \eqn{\psi}{psi}
 and \eqn{p-q}-dimensional subvector \eqn{\tau}{tau}. Assume also be
 interested in testing the specific values of the vector
 \eqn{\psi}{psi}. Then one can use some statistical
 hypothesis tests for testing the null hypothesis \eqn{H_0:
 \psi=\psi_0}{H_0: psi=psi_0} against the alternative \eqn{H_1:
 \psi \neq \psi_0}{H_1: psi <> psi_0}. Composite likelihood versions
of 'Wald' and 'score' statistics have the usual asymptotic
chi-square distribution with \eqn{q} degree of freedom. The Wald-type statistic is
\deqn{W=(\hat{\psi}-\psi_0)^T G_{\psi \psi}(\hat{\theta})(\hat{\psi}-\psi_0),}{%
W=(hat{psi}-psi_0)^T G_{psi psi} (hat{theta})(hat{psi}-psi_0),}
where \eqn{G_{\psi \psi}}{G_{psi psi}} is the \eqn{q \times q}{q x q}
 submatrix of the Godambe information pertaining to \eqn{\psi}{psi} and
 \eqn{\hat{\theta}}{hat{theta}} is the maximum likelihood estimator from
 the full model. The score-type statistic (Rao test) is
 \deqn{W=s_{\mathbf{\psi}}\{\mathbf{\psi}_0,
 \hat{\mathbf{\tau}}(\mathbf{\psi}_0)\}^T H^{\psi
   \psi}(\hat{\theta}_\psi) \{G^{\psi \psi}(\hat{\theta}_\psi)\}^{-1}
 H^{\psi \psi}(\hat{\theta}_\psi) 
 s_{\mathbf{\psi}}\{\mathbf{\psi}_0, \hat{\mathbf{\tau}}(\mathbf{\psi}_0)\},}{%
   W=s_{psi}{psi_0, hat{tau}(psi_0)}^T H^{psi psi}(hat{theta}_psi)
   \{G^{psi psi}(hat{theta}_psi)\}^{-1} H^{psi psi}(hat{theta}_psi)
   s_{psi}{psi_0, hat{tau}(psi_0)},}
 where \eqn{H^{\psi \psi}}{H^{psi psi}} is the \eqn{q \times q}{q x q}
 submatrix of the inverse of \eqn{H(\theta)}{H(theta)} pertaining to
 \eqn{\psi}{psi} (the same for \eqn{G}) and
 \eqn{\hat{\theta}_\psi}{hat{theta}_psi} is the constrained maximum
 likelihood estimate of \eqn{\theta}{theta} for fixed \eqn{\psi}{psi}.
 These two statistics can be called from the
 routine \code{HypoTest} assigning at the argument \code{statistic}
 respectively the values: \code{Wald} and \code{Rao}.
 

 Alternatively to the 'Wald' and 'score' statistics one can use the
 composite likelihood ratio statistic (Wilks) that is given by
 \deqn{W=2[C \ell(\hat{\mathbf{\theta}};\mathbf{y}) - C \ell\{\mathbf{\psi}_0,
   \hat{\mathbf{\tau}}(\mathbf{\psi}_0);\mathbf{y}\}].}{%
 W=2[Cl(hat{theta};y) - Cl{psi_0, hat{tau}(psi_0);y}].}
 The asymptotic distribution of the composite likelihood ratio statistic
 is given by
 \deqn{W \dot{\sim} \sum_{i} \lambda_i \chi^2,}{%
   W ~ sum_i lambda_i Chi^2_i,}
 for \eqn{i=1,\ldots,q}{i=1,...,q}, where \eqn{\chi^2_i}{Chi^2_i} are
 \eqn{q} iid copies of a chi-square one random variable and
 \eqn{\lambda_1,\ldots,\lambda_q}{lambda_1,...,lambda_q}
 are the eigenvalues of the matrix \eqn{(H^{\psi \psi})^{-1} G^{\psi
 \psi}}{(H^{psi psi})^-1 G^{psi psi}}. There exist several adjustments
 to the composite likelihood ratio statistic in order to get an
 approximate \eqn{\chi^2_q}{Chi^2_q}. For example, Rotnitzky and Jewell
 (1990) proposed the adjustment \eqn{W'= W / \bar{\lambda}}{W'= W /
   bar{lambda}} where \eqn{\bar{\lambda}}{bar{lambda}} is the average
 of the eigenvalues \eqn{\lambda_i}{lambda_i}. This statistic can be
 called within the routine by the value: \code{Wilks-RJ}. A better
 solution is proposed by Satterhwaite (1946) defining \eqn{W''=\nu W /
 (q \bar{\lambda})}{W''= nu W / (q bar{lambda})}, where \eqn{\nu=\sum_i
 \lambda / \sum_i \lambda^2_i}{nu = sum_i lambda / sum_i lambda^2_i} for
\eqn{i=1\ldots,q}{i=1...,q}, is the effective number of the degree of
freedom. Note that in this case the distribution of the likelihood ratio
 statistic is a chi-square random variable with \eqn{\nu}{nu} degree of
 freedom. This statistic can be called from the routine assigning the
 value: \code{Wilks-S}.

}

\value{
 An object of class \code{c("data.frame")}. The object contain a table
 with the results of the tested models. The rows represent the
 responses for each model and the columns the following results:
 \item{Num.Par}{The number of the model's parameters.}
 \item{Diff.Par}{The number of parameters of the model in the previous row
   minus the number of parameters in the actual row.}
 \item{Df}{The effective number of degree of freedom of the chi-square distribution.}
 \item{Chisq}{The observed value of the statistic.}
 \item{Pr(>chisq)}{The p-value calculated by comparing the quantile
 \code{Chisq} with a chi-squared distribution on \code{Df} degrees
    of freedom.}
}

\references{
Chandler, R. E., and Bate, S. (2007). Inference for Clustered Data Using the
Independence log-likelihood, \emph{Biometrika}, \bold{94}, 167--183.
  
Varin, C., Reid, N. and Firth, D. (2010), An Overview of Composite
Likelihood Methods. Submitted.

Padoan, S. A. Ribatet, M and Sisson, S. A. (2010). Likelihood-Based Inference for
Max-Stable Processes. \emph{Journal of the American Statistical Association, Theory & Methods}.
In press.

Satterthwaite, F. E. (1946). An Approximate Distribution of Estimates of
Variance Components. \emph{Biometrics Bulletin}, \bold{2}, 110--114.
}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}}

\seealso{
  \command{\link{fExtremal}}.
}
\examples{
library(RandomFields)
set.seed(2010)

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

# Simulation of the first Extremal t process:
sim <- rExtremal(coord, corrmodel, model=model, numsim=numsim,
                 numblock=numblock, param=param)

# Fitting Extremal t process estimating all the parameters:
fit1 <- fExtremal(coord, corrmodel, sim, model=model)

# Fitting Extremal t process fixing the power parameter:
fit2 <- fExtremal(coord, corrmodel, sim, model=model,
                  fixed=list(nugget=0, power=1))

# Fitting Extremal t process fixing power and scale parameters:
fit3 <- fExtremal(coord, corrmodel, sim, model=model,
                  fixed=list(nugget=0, power=1, scale=30))

# Results:
# Composite likelihood ratio (Wilks) statistic - RJ-type:
HypoTest(fit1, fit2, fit3, statistic='Wilks-RJ')

# Composite likelihood ratio (Wilks) statistic - S-type:
HypoTest(fit1, fit2, fit3, statistic='Wilks-S')

# Composite Wald-type statistic:
HypoTest(fit1, fit2, fit3, statistic='Wald')

# Composite score-type (Rao) statistic:
HypoTest(fit1, fit2, fit3, statistic='Rao')

}
\keyword{spatial}
