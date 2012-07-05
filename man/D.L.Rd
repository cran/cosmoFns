\name{D.L}
\alias{D.L}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Luminosity distance
}
\description{
Function computes luminosity distance in a flat cosmology.
}
\usage{
D.L(z, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71)
}
\arguments{
  \item{z}{Redshift}
  \item{omega.m}{Omega matter parameter}
  \item{omega.lambda}{Omega lambda parameter}
  \item{H.0}{Hubble constant in km/s/Mpc}
}

\value{
  Luminosity distance in Mpc
}
\references{
Hogg (2000), arXiv:astro-ph/9905116, equation (21)
}
\author{
A. Harris
}

\examples{
D.L(2.3)
}
\keyword{ misc }
