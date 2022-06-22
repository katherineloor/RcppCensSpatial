\name{dist2Dmatrix}
\alias{dist2Dmatrix}
\title{Distance matrix computation}
\usage{
dist2Dmatrix(coords)
}
\arguments{
\item{coords}{2D spatial coordinates.}
}
\value{
An \eqn{n \times n} distance matrix.
}
\description{
It computes the Euclidean distance matrix for a set of coordinates.
}
\examples{
n = 100
set.seed(1000)
x = round(runif(n,0,10), 5)     # X coordinate
y = round(runif(n,0,10), 5)     # Y coordinate
Mdist = dist2Dmatrix(cbind(x, y))
}
\author{
Katherine L. Valeriano, Alejandro Ordoñez, Christian E. Galarza, and Larissa A. Matos.
}
