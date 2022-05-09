\name{dimmingFactor}
\alias{dimmingFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Flux dimming factor
}
\description{
Function computes flux dimming factor in a flat cosmology.
}
\usage{
dimmingFactor(z, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71)
}

\arguments{
  \item{z}{Redshift}
  \item{omega.m}{Omega matter parameter}
  \item{omega.lambda}{Omega lambda parameter}
  \item{H.0}{Hubble constant in km/s/Mpc}
}

\value{
  Flux dimming factor, unnormalied.  Mathematically, it is \code{(1+z)/D.L^2}.
  This is the factor that scales luminosity densitiy in the observed
  frame to flux density in the observed frame.
}
\references{
Hogg (2000), arXiv:astro-ph/9905116: section 7, part of equation (22)
}
\author{
A. Harris
}

\seealso{
\code{\link{D.L}}
}

\examples{
z <- seq(0.1, 5, 0.1)
df <- dimmingFactor(z)
plot(z, df/max(df), t='l', xlab='z', ylab='Normalized dimming factor')
}

\keyword{ misc }

