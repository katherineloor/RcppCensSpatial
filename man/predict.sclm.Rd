\name{predict.sclm}
\alias{predict.sclm}
\title{Prediction in spatial models with censored/missing responses}
\usage{
\method{predict}{sclm}(object, locPre, xPre, ...)
}
\arguments{
\item{object}{object of class \code{'sclm'} given as output of \code{\link{EM.sclm}}, \code{\link{MCEM.sclm}},
or \code{\link{SAEM.sclm}} function.}

\item{locPre}{matrix of coordinates for which prediction is performed.}

\item{xPre}{matrix of covariates for which prediction is performed.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
The function returns a data frame with:
\item{xcoord}{x coordinates.}
\item{ycoord}{y coordinates.}
\item{predValues}{predicted values.}
\item{sdPred}{predicted standard deviations.}
}
\description{
This function performs spatial prediction in a set of new \code{S} spatial locations.
}
\details{
This function predicts using the Mean Squared Error (MSE) criterion, which takes the conditional
expectation E(Y|X) as the best linear predictor.
}
\examples{\donttest{
set.seed(1000)
n = 120
coords = round(matrix(runif(2*n,0,15),n,2), 5)
x = cbind(rbinom(n,1,0.50), rnorm(n), rnorm(n))
data = rCensSp(c(1,4,-1), 2, 3, 0.50, x, coords, "left", 0.10,
               20, "exponential")
data1 = data$Data

# Estimation: EM algorithm
fit1 = EM.sclm(y=data1$y, x=data1$x, ci=data1$ci, lcl=data1$lcl, ucl=data1$ucl,
               coords=data1$coords, init.phi=2.50, init.nugget=1,
               show.SE=TRUE, error=1e-4)

# Estimation: SAEM algorithm
fit2 = SAEM.sclm(y=data1$y, x=data1$x, ci=data1$ci, lcl=data1$lcl, ucl=data1$ucl,
                 coords=data1$coords, init.phi=2.50, init.nugget=1,
                 show.SE=TRUE, error=1e-4)

# Estimation: MCEM algorithm
fit3 = MCEM.sclm(y=data1$y, x=data1$x, ci=data1$ci, lcl=data1$lcl, ucl=data1$ucl,
                 coords=data1$coords, init.phi=2.50, init.nugget=1,
                 MaxIter=300, show.SE=TRUE, error=1e-4)
rbind(fit1$theta, fit2$theta, fit3$theta)

# Prediction
data2 = data$TestData
pred1 = predict(fit1, data2$coords, data2$x)
pred2 = predict(fit2, data2$coords, data2$x)
pred3 = predict(fit3, data2$coords, data2$x)

# Cross-validation
mean((data2$y - pred1$predValues)^2)
mean((data2$y - pred2$predValues)^2)
mean((data2$y - pred3$predValues)^2)}
}
\seealso{
\code{\link{EM.sclm}}, \code{\link{MCEM.sclm}}, \code{\link{SAEM.sclm}}
}
\author{
Katherine L. Valeriano, Alejandro Ordoñez, Christian E. Galarza, and Larissa A. Matos.
}
