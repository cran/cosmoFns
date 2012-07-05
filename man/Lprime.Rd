\name{Lprime}
\alias{Lprime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Line luminosity, L'
}
\description{
Compute L' line luminosity
}
\usage{
Lprime(intInt, z, f.rest = 115.27, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71)
}

\arguments{
  \item{intInt}{Integrated intensity in Jy km/s}
  \item{z}{Redshift}
  \item{f.rest}{Line rest frequency in GHz}
  \item{omega.m}{Omega matter parameter}
  \item{omega.lambda}{Omega lambda parameter}
  \item{H.0}{Hubble constant in km/s/Mpc}
}

\value{
Rest-frame line luminosity in \code{K km/s pc^-2}.
}
\references{
  Solomon, Downes & Radford (1992), ApJ 398, L29, equation (3)
}
\author{
A. Harris
}
\note{
For flat universe, \code{omega.k = 0}.  Useful for empirical mass
estimates. L' is proportional to the brightness temperature of the
transition.
}

\seealso{
\code{\link{lineLum}}, \code{\link{mass.CO}}
}
\examples{
snu <- 1.e-3 # 1 mJy peak
wid <- 400   # 400 km/s wide
intInt <- 1.06*snu*wid  # Gaussian line
z <- 2.3
Lprime(intInt, z)
}

\keyword{ misc }

