\name{dComovVol}
\alias{dComovVol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential comoving volume
}
\description{
Function computes differential comoving volume in a flat cosmology.
}
\usage{
dComovVol(z, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71)
}

\arguments{
  \item{z}{Redshift}
  \item{omega.m}{Omega matter parameter}
  \item{omega.lambda}{Omega lambda parameter}
  \item{H.0}{Hubble constant in km/s/Mpc}
}

\value{
Differential comoving volume in Mpc^3
}
\references{
Hogg (2000), arXiv:astro-ph/9905116, equation (28)
}
\author{
A. Harris
}


\examples{
dComovVol(2.3)
}

\keyword{misc}
