\name{GaussWS}
\alias{GaussWS}

\title{Convenience function to look-up values from a object returned by GenerateGauss.}

\description{
Dependent on the resolution at which \code{\link{GenerateGauss}} was set to use, the resolution coordinates (array indices) cannot immediately be interpreted as being measured in meters. This function returns wind speeds from the Gauss object accoring to x, y, z simply provided in meters.
}

\usage{
GaussWS(Gauss, x, y, z)
}

\arguments{
	\item{Gauss}{
		must be an object returned by \code{\link{GenerateGauss}}.
	}
	\item{x}{
		must be a single value. Provide desired distance in meters.
	}
	\item{y}{
		must be a single value. Provide desired distance in meters. May be negative.
	}
	\item{z}{
		must be a single value. Provide desired distance in meters.
	}
}

\details{
The Gaussian wake model is loosely based on the initial contribution by Bastankhah & Porte-Agel (2014).
}

\value{
\code{GaussWS} returns a single number which can be considered a wind speed in the wake of a turbine at location x, y, and z.
}

\references{
Bastankhah, M., & Porte-Agel, F. (2014). A new analytical model for wind-turbine wakes. Renewable Energy, 70, 116-123.
}

\author{
Carsten Croonenbroeck
}

\seealso{
Use \code{\link{GenerateGauss}} to compute the three-dimensional tensor array object containing the wind speed data. See \code{\link{QuickGauss3D}} for the same algorithm, immediately returning the wind speed at one single point only.
}

\examples{
\dontrun{
GaussWS(Gauss, 100, 1, 100)
GaussWS(Gauss, 80, -40, 90)
GaussWS(Gauss, 200, 40, 150)
}
}

\keyword{GenerateGauss}
\keyword{GaussWS}
\keyword{QuickGauss3D}
