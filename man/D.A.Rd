\name{D.A}
\alias{D.A}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Angular diameter distance
}
\description{
  Function computes angular diameter distance
}
\usage{
D.A(z, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71)
}

\arguments{
  \item{z}{Redshift}
  \item{omega.m}{Omega matter parameter}
  \item{omega.lambda}{Omega lambda parameter}
  \item{H.0}{Hubble constant in km/s/Mpc}
}

\value{
Angular distance in Mpc
}
\references{
Hogg (2000), arXiv:astro-ph/9905116, equation (18)
}
\author{
A. Harris
}
\note{
For flat universe, \code{omega.k = 0}.
}

\examples{
D.A(2.3)

z <- seq(0.1, 5, 0.1)
d <- D.A(z)
plot(z, d/max(d), t='l', xlab='z', ylab='Normalized D.A')

}

\keyword{ misc }

