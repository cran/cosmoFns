\name{D.M}
\alias{D.M}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comoving distance
}
\description{
Function computes comoving distance in a flat cosmology.
}
\usage{
D.M(z, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71)
}
\arguments{
  \item{z}{Redshift}
  \item{omega.m}{Omega matter parameter}
  \item{omega.lambda}{Omega lambda parameter}
  \item{H.0}{Hubble constant in km/s/Mpc}
}

\value{
  Comoving distance in Mpc
}
\references{
Hogg (2000), arXiv:astro-ph/9905116, equations (16) and (15)
}
\author{
A. Harris
}
\note{
For flat universe, \code{omega.k = 0}, so transverse and line-of-sight
comoving distances are equal.
}


\examples{
D.M(2.3)
}

\keyword{ misc }
