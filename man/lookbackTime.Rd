\name{lookbackTime}
\Rdversion{1.1}
\alias{lookbackTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cosmic lookback time
}
\description{
Compute cosmic lookback time given z and cosmological parameters
}
\usage{
lookbackTime(z, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{Redshift}
  \item{omega.m}{Omega matter parameter}
  \item{omega.lambda}{Omega lambda parameter}
  \item{H.0}{Hubble constant in km/s/Mpc}
}
\details{
Defaults for \code{omega.m}, \code{omega.lambda}, and \code{omega.m},
are from WMAP cosmology; \code{omega.k} (curvature term) is computed
from relationship between \code{omegas} in flat cosmology
(\code{omega.k = 0}).
}
\value{
  Lookback time in Gyr.
}
\references{
"Principles of Physical Cosmology," P.J. Peebles, Princeton c. 1993,
(5.63); "Distance Measures in Cosmology," Hogg (2000),
arXiv:astro-ph/9905116, equation (30);
"First-year WMAP observations...", Spergel et al., ApJS 148:175 (2003)
}
\author{
  A. Harris
}
\examples{
# lookback time for z = 2
lookbackTime(2)
# Inverse problem, age of Earth (4.6 Gyr) example:
uniroot(function(x) lookbackTime(x) - 4.6, c(0,2))$root
}
\keyword{ misc }

