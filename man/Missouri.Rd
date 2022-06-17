\docType{data}
\name{Missouri}
\alias{Missouri}
\title{
TCDD concentration data
}
\description{
The level of dioxin (2,3,7,8-tetrachlorodibenzo-p-dioxin or TCDD) data was
collected in November 1983 by the U.S. Environmental Protection Agency (EPA) in several areas of a highway
from Missouri, USA. Only the locations used in the geostatistical analysis by
\insertCite{zirschky1986geostatistical}{RcppCensSpatial} are showed.
The TCDD measurement was subject to a limit of detection (\code{cens}), thereby the
\code{TCDD} data are left-censored.
}
\format{
A data frame with 127 observations and 5 variables:
\describe{
     \item{xcoord}{x coordinate of start of each transect (ft).}
     \item{ycoord}{y coordinate of start of each transect (ft).}
     \item{TCDD}{TCDD concentrations (mg/kg).}
     \item{transect}{transect length (ft).}
     \item{cens}{indicator of censoring (left-censored observations).}
}}
\source{
\insertRef{zirschky1986geostatistical}{RcppCensSpatial}
}
\usage{
data("Missouri")
}
\examples{\donttest{
set.seed(83789)
data("Missouri")
y = log(Missouri$TCDD)
cc = Missouri$cens
coord = cbind(Missouri$xcoord/100, Missouri$ycoord)
x = matrix(1, length(y), 1)
lcl = rep(-Inf, length(y))
ucl = y

fit2 = MCEM.sclm(y=y, x=x, ci=cc, lcl=lcl, ucl=ucl, coords=coord, init.phi=5,
              init.nugget=1, type="exponential", lower=c(1e-5,1e-5),
              upper=c(50,50), MaxIter=200, nMax=1000, error=1e-5)
summary(fit2)
plot(fit2)}
}
\keyword{datasets}
