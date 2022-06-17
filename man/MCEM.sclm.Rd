\name{MCEM.sclm}
\alias{MCEM.sclm}
\title{
ML estimation of spatial censored linear models via MCEM algorithm
}
\description{
This function fits a spatial censored linear model using the Monte Carlo EM (MCEM) algorithm.
It supports left, right, interval, or missing values in the dependent variable.
It provides estimates and standard errors (SEs) of parameters. The SEs are computed through the observed information
matrix using the method developed by \insertCite{louis1982finding;textual}{RcppCensSpatial}.
}
\usage{
MCEM.sclm(y, x, ci, lcl=NULL, ucl=NULL, coords, init.phi, init.nugget,
  type="exponential", kappa=NULL, lower=c(0.01,0.01), upper=c(30,30),
  MaxIter=500, nMin=20, nMax=5000, error=1e-05, show.SE=TRUE)
}
\arguments{
\item{y}{vector of responses.}

\item{x}{design matrix.}

\item{ci}{vector of censoring indicators. For each observation: \code{1} if censored/missing
and \code{0} otherwise.}

\item{lcl, ucl}{Vectors of length \eqn{n} that represent the lower and upper bounds of the interval, which contains the true value of the censored observation. Default =\code{NULL}, indicating no-censored data. For left censoring, \code{lcl=-Inf} and \code{ucl=c}. For right censoring, \code{lcl=c} and \code{ucl=Inf}. For interval censoring, \code{lcl} and \code{ucl} must be finite, while missing data could be defined by setting \code{lcl=-Inf} and \code{ucl=Inf}.}

\item{coords}{2D spatial coordinates.}

\item{init.phi}{initial value for the spatial scaling parameter.}

\item{init.nugget}{initial value for the nugget effect parameter.}

\item{type}{type of spatial correlation function: \code{'exponential'}, \code{'gaussian'},
\code{'matern'}, and \code{'pow.exp'} for exponential, gaussian, matern, and power exponential, respectively.}

\item{kappa}{parameter for all spatial correlation functions. See \code{\link{CovMat}}.}

\item{lower, upper}{vectors of lower and upper bounds for the optimization method. If unspecified, the default is
\code{c(0.01,0.01)} for lower and \code{c(30,30)} for upper.}

\item{MaxIter}{maximum number of iterations for the MCEM algorithm. By default \code{=500}.}

\item{nMin}{initial sample size for Monte Carlo integration. By default \code{=20}.}

\item{nMax}{maximum sample size for Monte Carlo integration. By default \code{=5000}.}

\item{error}{maximum convergence error. By default \code{=1e-5}.}

\item{show.SE}{\code{TRUE} or \code{FALSE}. It indicates if the standard errors should be estimated. By default \code{=TRUE}.}
}
\value{
Return an object of class \code{sclm}. Generic functions such as \code{print} and \code{summary}
have methods to show the results of the fit. The function \code{plot} can be used to extract convergence
graphs for the parameters estimates.

Specifically, the following components are returned:

\item{Theta}{estimated parameters in all iterations, \eqn{\theta = (\beta, \sigma^2, \phi, \tau^2)}.}
\item{theta}{final estimation of \eqn{\theta = (\beta, \sigma^2, \phi, \tau^2)}.}
\item{beta}{estimated \eqn{\beta}.}
\item{sigma2}{estimated \eqn{\sigma^2}.}
\item{phi}{estimated \eqn{\phi}.}
\item{tau2}{estimated \eqn{\tau^2}.}
\item{EY}{MC approximation of the first moment for the truncated normal distribution.}
\item{EYY}{MC approximation of the second moment for the truncated normal distribution.}
\item{SE}{vector of standard errors of \eqn{\theta = (\beta, \sigma^2, \phi, \tau^2)}.}
\item{InfMat}{observed information matrix.}
\item{loglik}{log-likelihood for the MCEM method.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{Iterations}{number of iterations needed to converge.}
\item{time}{processing time.}
\item{range}{the effective range.}
}
\details{
The spatial Gaussian model is given by

\eqn{Y = X\beta + \xi},

where \eqn{Y} is the \eqn{n}x\eqn{1} response vector, \eqn{X} is the \eqn{n}x\eqn{q} design matrix,
\eqn{\beta} is the \eqn{q}x\eqn{1} vector of regression coefficients to be estimated, and \eqn{\xi}
is the error term. Which is normally distributed with zero-mean and covariance matrix
\eqn{\Sigma=\sigma^2 R(\phi) + \tau^2 I_n}. We assume that \eqn{\Sigma} is non-singular and
\eqn{X} has a full rank \insertCite{diggle2007springer}{RcppCensSpatial}.

The estimation process is performed via the MCEM algorithm, initially proposed by
\insertCite{wei1990monte;textual}{RcppCensSpatial}. The Monte Carlo (MC) approximation starts with a
sample of size \code{nMin}; at each iteration, the sample size increases (\code{nMax-nMin})/\code{MaxIter},
and at the last iteration, the sample size is \code{nMax}. The random observations are sampled
through the slice sampling algorithm available in package \code{relliptical}.
}
\note{
The MCEM final estimates correspond to the mean of the estimates obtained at each iteration after deleting
the half and applying a thinning of 3.

To fit a regression model for non-censored data, just set \code{ci} as a vector of zeros.
}
\examples{
# Example 1: left censoring data
set.seed(1000)
n = 50   # Test with another values for n
coords = round(matrix(runif(2*n,0,15),n,2), 5)
x = cbind(rnorm(n), rnorm(n))
data = rCensSp(c(2,-1), 2, 3, 0.70, x, coords, "left", 0.08, 0, "matern", 1)

fit = MCEM.sclm(y=data$y, x=x, ci=data$ci, lcl=data$lcl, ucl=data$ucl, coords,
                init.phi=2.50, init.nugget=0.75, type="matern", kappa=1,
                MaxIter=20, nMax=1000, error=1e-4)
print(fit)

\donttest{
# Example 2: missing data
yMiss = data$y
yMiss[20] = NA
ci = data$ci
ci[20] = 1
ucl = data$ucl
ucl[20] = Inf

fit1 = MCEM.sclm(y=yMiss, x=x, ci=ci, lcl=data$lcl, ucl=ucl, coords,
                 init.phi=2.50, init.nugget=0.75, type="matern", kappa=1,
                 MaxIter=20, nMax=1000, error=1e-4)
summary(fit1)
rbind(fit$theta, fit1$theta)}
}
\references{
\insertAllCited
}
\seealso{
\code{\link{EM.sclm}}, \code{\link{SAEM.sclm}}, \code{\link{predict.sclm}}
}
\author{
Katherine L. Valeriano, Alejandro Ordoñez, Christian E. Galarza, and Larissa A. Matos.
}
