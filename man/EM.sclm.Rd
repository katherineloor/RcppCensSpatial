\name{EM.sclm}
\alias{EM.sclm}
\title{
ML estimation of spatial censored linear models via EM algorithm
}
\description{
This function fits a spatial censored linear model using the Expectation-Maximization (EM) algorithm.
It supports left, right, interval, or missing values in the dependent variable.
It provides estimates and standard errors (SEs) of parameters. The SEs are computed through the observed information
matrix using the method developed by \insertCite{louis1982finding;textual}{RcppCensSpatial}.
}
\usage{
EM.sclm(y, x, ci, lcl=NULL, ucl=NULL, coords, init.phi, init.nugget,
  type="exponential", kappa=NULL, lower=c(0.01,0.01), upper=c(30,30),
  MaxIter=300, error=1e-05, show.SE=TRUE)
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

\item{MaxIter}{maximum number of iterations for the EM algorithm. By default \code{=300}.}

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
\item{EY}{first moment for the truncated normal distribution computed in the last iteration.}
\item{EYY}{second moment for the truncated normal distribution computed in the last iteration.}
\item{SE}{vector of standard errors of \eqn{\theta = (\beta, \sigma^2, \phi, \tau^2)}.}
\item{InfMat}{observed information matrix.}
\item{loglik}{log-likelihood for the EM method.}
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

The estimation process is performed via the EM algorithm, initially proposed by
\insertCite{dempster1977maximum;textual}{RcppCensSpatial}. The conditional expectations are
computed using the function \code{meanvarTMD} available in the package \code{MomTrunc}.
}
\note{
The EM final estimates correspond to the estimates obtained at the last iteration of the
EM algorithm.

To fit a regression model for non-censored data, just set \code{ci} as a vector of zeros.
}
\examples{
# Simulated example: 10\% of left-censored observations
set.seed(1000)
n = 50   # Test with another values for n
coords = round(matrix(runif(2*n,0,15),n,2), 5)
x = cbind(rnorm(n), runif(n))
data = rCensSp(c(-1,3), 2, 4, 0.5, x, coords, "left", 0.10, 0, "gaussian")

fit = EM.sclm(y=data$y, x=x, ci=data$ci, lcl=data$lcl, ucl=data$ucl,
              coords=coords, init.phi=3, init.nugget=1, type="gaussian",
              error=1e-4)
fit
}
\references{
\insertAllCited
}
\seealso{
\code{\link{MCEM.sclm}}, \code{\link{SAEM.sclm}}, \code{\link{predict.sclm}}
}
\author{
Katherine L. Valeriano, Alejandro Ordoñez, Christian E. Galarza, and Larissa A. Matos.
}
