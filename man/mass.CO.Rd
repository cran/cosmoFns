\name{mass.CO}
\alias{mass.CO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Molecular mass
}
\description{
  Compute molecular mass (default CO J = 1-0) from L' and
  empirical conversion factor.
}
\usage{
mass.CO(intInt, z, alpha = 0.8, f.rest = 115.27, omega.m = 0.27,
omega.lambda = 0.73, H.0 = 71)
}

\arguments{
  \item{intInt}{Integrated intensity in Jy km/s}
  \item{z}{Redshift}
  \item{alpha}{Empirical mass conversion factor, see details}
  \item{f.rest}{Line rest frequency in GHz}
  \item{omega.m}{Omega matter parameter}
  \item{omega.lambda}{Omega lambda parameter}
  \item{H.0}{Hubble constant in km/s/Mpc}
}
\details{
\code{alpha} is an empirical mass conversion factor.  The exact value is
a topic of considerable debate.  For CO, see
Solomon and Vanden Bout (2005), also Tacconi et al. (2008) for reviews.
}
\value{
Gas mass in solar masses.
}
\references{
Solomon, Downes & Radford (1992), ApJ 398, L29, equations (3) and (4);
Solomon & Vanden Bout (2005) ARA&A 43, 677; Tacconi et al. (2008) ApJ
680, 246.
}
\author{
A. Harris
}

\seealso{
\code{\link{Lprime}}
}
\examples{
snu <- 1.e-3 # 1 mJy peak
wid <- 400   # 400 km/s wide
intInt <- 1.06*snu*wid  # Gaussian line
z <- 2.3
mass.CO(intInt, z)
}

\keyword{ misc }

