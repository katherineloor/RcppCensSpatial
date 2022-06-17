\name{SAEM.sclm}
\alias{SAEM.sclm}
\title{
ML estimation of spatial censored linear models via SAEM algorithm
}
\description{
This function fits a spatial censored linear model using the Stochastic Approximation EM (SAEM) algorithm.
It supports left, right, interval, or missing values in the dependent variable.
It provides estimates and standard errors (SEs) of parameters. The SEs are computed through the observed information
matrix using the method developed by \insertCite{louis1982finding;textual}{RcppCensSpatial}.
}
\usage{
SAEM.sclm(y, x, ci, lcl=NULL, ucl=NULL, coords, init.phi, init.nugget,
  type="exponential", kappa=NULL, lower=c(0.01, 0.01), upper=c(30,30),
  MaxIter=300, M=20, pc=0.25, error=1e-05, show.SE=TRUE)
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

\item{type}{type of spatial correlation function: '\code{exponential}', '\code{gaussian}',
'\code{matern}', and '\code{pow.exp}' for exponential, gaussian, matern, and power exponential, respectively.}

\item{kappa}{parameter for all spatial correlation functions. See \code{\link{CovMat}}.}

\item{lower, upper}{vectors of lower and upper bounds for the optimization method. If unspecified, the default is
\code{c(0.01,0.01)} for lower and \code{c(30,30)} for upper.}

\item{MaxIter}{maximum number of iterations of the SAEM algorithm. By default \code{=300}.}

\item{M}{number of Monte Carlo samples for stochastic approximation. By default \code{=20}.}

\item{pc}{percentage of iterations of the SAEM algorithm with no-memory. By default \code{=0.25}.}

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
\item{EY}{stochastic approximation of the first moment for the truncated normal distribution.}
\item{EYY}{stochastic approximation of the second moment for the truncated normal distribution.}
\item{SE}{vector of standard errors of \eqn{\theta = (\beta, \sigma^2, \phi, \tau^2)}.}
\item{InfMat}{observed information matrix.}
\item{loglik}{log-likelihood for the SAEM method.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{Iterations}{number of iterations needed to converge.}
\item{time}{processing time.}
\item{range}{the effective range.}
}
\details{
The spatial Gaussian model is given by

\eqn{Y = X\beta + \xi},

where \eqn{Y} is the \eqn{n x 1} response vector, \eqn{X} is the \eqn{n x q} design matrix,
\eqn{\beta} is the \eqn{q x 1} vector of regression coefficients to be estimated, and \eqn{\xi}
is the error term which is normally distributed with zero-mean and covariance matrix
\eqn{\Sigma=\sigma^2 R(\phi) + \tau^2 I_n}. We assume that \eqn{\Sigma} is non-singular and
\eqn{X} has full rank \insertCite{diggle2007springer}{RcppCensSpatial}.

The estimation process is performed via the SAEM \insertCite{delyon1999convergence}{RcppCensSpatial} algorithm.
The spatial SAEM algorithm was previously proposed by \insertCite{lachos2017influence;textual}{RcppCensSpatial}
and \insertCite{ordonez2018geostatistical;textual}{RcppCensSpatial} and is available in package \code{CensSpatial}.
The difference between this package to \code{CensSpatial} is that the random observations are sampled
through the slice sampling algorithm available in package \code{relliptical} and the optimization procedure
by the \code{roptim} package.

This model is also a particular case of the Spatio-temporal model defined by \insertCite{valeriano2021likelihood;textual}{RcppCensSpatial},
when the number of temporal observations is equal to one. The computing codes of the Spatio-temporal
SAEM algorithm are available in the package \code{StempCens}.
}
\note{
The SAEM final estimates correspond to the estimates obtained at the last iteration of the
algorithm.

To fit a regression model for non-censored data, just set \code{ci} as a vector of zeros.
}
\examples{
# Example 1: 10\% of right-censored observations
set.seed(1000)
n = 50   # Test with another values for n
coords = round(matrix(runif(2*n,0,15),n,2), 5)
x = cbind(rbinom(n,1,0.50), rnorm(n), rnorm(n))
data = rCensSp(c(1,4,-2), 2, 3, 0.50, x, coords, "right", 0.10, 0, "matern", 2)

fit = SAEM.sclm(y=data$y, x=x, ci=data$ci, lcl=data$lcl, ucl=data$ucl, coords,
                init.phi=2, init.nugget=1, type="matern", kappa=2,
                MaxIter=20, error=1e-4)
fit
\donttest{
# Example 2: censored and missing observations
set.seed(123)
n = 200
coords = round(matrix(runif(2*n,0,20),n,2), 5)
x = cbind(1, rnorm(n), rexp(n))
data = rCensSp(c(1,4,-1), 2, 4, 0.50, x, coords, "left", 0.10)
data$y[c(10,20)] = NA
data$ci[c(10,20)] = 1
data$lcl[c(10,20)] = -Inf
data$ucl[c(10,20)] = Inf

fit2 = SAEM.sclm(y=data$y, x=x, ci=data$ci, lcl=data$lcl, ucl=data$ucl, coords,
                 init.phi=2, init.nugget=1, MaxIter=300, error=1e-4)
fit2$theta  # estimates
fit2$SE     # standard error
plot(fit2)}
}
\references{
\insertAllCited
}
\seealso{
\code{\link{EM.sclm}}, \code{\link{MCEM.sclm}}, \code{\link{predict.sclm}}
}
\author{
Katherine L. Valeriano, Alejandro Ordoñez, Christian E. Galarza, and Larissa A. Matos.
}
