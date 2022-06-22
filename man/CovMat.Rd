\name{CovMat}
\alias{CovMat}
\title{
Covariance matrix for spatial models
}
\description{
It computes the spatial variance-covariance matrix considering exponential,
gaussian, matérn, or power exponential correlation function.
}
\usage{
CovMat(phi, tau2, sigma2, dist, type="exponential", kappa=NULL)
}
\arguments{
\item{phi}{spatial scaling parameter.}

\item{tau2}{nugget effect parameter.}

\item{sigma2}{partial sill parameter.}

\item{dist}{\eqn{n \times n} distance matrix.}

\item{type}{type of spatial correlation function: \code{'exponential'},
\code{'gaussian'}, \code{'matern'}, and \code{'pow.exp'} for exponential,
gaussian, matérn, and power exponential, respectively.}

\item{kappa}{parameter for some spatial correlation functions. For exponential
and gaussian \code{kappa=NULL}, for power exponential \code{0 < kappa <= 2}, and for matérn correlation function \code{kappa > 0}.}
}
\value{
An \eqn{n \times n} spatial covariance matrix.
}
\details{
The spatial covariance matrix is given by

\eqn{\Sigma = [Cov(s_i, s_j )] = \sigma^2 R(\phi) + \tau^2 I_n},

where \eqn{\sigma^2 > 0} is the partial sill, \eqn{\phi > 0} is the spatial
scaling parameter, \eqn{\tau^2 > 0} is known as the nugget effect in the
geostatistical framework, \eqn{R(\phi)} is the \eqn{n \times n} correlation
matrix computed from a correlation function, and \eqn{I_n} is the \eqn{n \times n} identity matrix.

The spatial correlation functions available are:
\describe{
\item{\strong{Exponential}:}{\eqn{Corr(d) = exp(-d/\phi)},}

\item{\strong{Gaussian}:}{\eqn{Corr(d) = exp(-(d/\phi)^2)},}

\item{\strong{Matern}:}{\eqn{Corr(d) = \frac{1}{2^{(\kappa-1)}\Gamma(\kappa)}\left(\frac{d}{\phi}\right)^\kappa K_\kappa \left( \frac{d}{\phi} \right)},}

\item{\strong{Power exponential}:}{\eqn{Corr(d) = exp(-(d/\phi)^\kappa)},}

where \eqn{d \geq 0} is the Euclidean distance between two observations,
\eqn{\Gamma(.)} is the gamma function, \eqn{\kappa} is the smoothness parameter, and \eqn{K_\kappa(.)} is the modified Bessel function of the second kind of
order \eqn{\kappa}.
}}
\examples{
set.seed(1000)
n = 20
x = round(runif(n,0,10), 5)     # X coordinate
y = round(runif(n,0,10), 5)     # Y coordinate
Ms = dist2Dmatrix(cbind(x, y))

Cov = CovMat(phi=5, tau2=0.8, sigma2=2, dist=Ms, type="exponential")
}
\seealso{
\code{\link{dist2Dmatrix}}, \code{\link{EM.sclm}}, \code{\link{MCEM.sclm}}, \code{\link{SAEM.sclm}}
}
\author{
Katherine L. Valeriano, Alejandro Ordoñez, Christian E. Galarza, and Larissa A. Matos.
}
