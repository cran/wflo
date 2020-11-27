\name{QuickGauss3D}
\alias{QuickGauss3D}

\title{For an incoming wind speed at reference height, this function computes a 3D Gaussian model based wind speed.}

\description{
This function computes Gaussian wind speeds.
}

\usage{
QuickGauss3D(x, y, z, u = 8, refHeight = 10)
}

\arguments{
	\item{x}{
		must be a single value. Provide distance in meters.
	}
	\item{y}{
		must be a single value. Provide distance in meters.
	}
	\item{z}{
		must be a single value. Provide distance in meters.
	}
	\item{u}{
		must be a single value. Provide wind speed in meters per second.
	}
	\item{refHeight}{
		must be a single value. Provide reference height (the height at which wind speed u was measured) in meters.
	}
}

\details{
The Gaussian wake model is loosely based on the initial contribution by Bastankhah & Porte-Agel (2014).
}

\value{
\code{QuickGauss3D} returns a single number which can be considered a wind speed in the wake of a turbine at location x, y, and z.
}

\references{
Bastankhah, M., & Porte-Agel, F. (2014). A new analytical model for wind-turbine wakes. Renewable Energy, 70, 116-123.
}

\author{
Carsten Croonenbroeck
}

\note{
Note that the model assumes that along the x axis, x = 0 is the turbine location. x expands along the wind direction downwind. y denoted whether a point is 'left' or 'right' the x axis. Thus, the x-z plane is the plane along the x axis and perpendicular to the ground. The z axis is hight, starting at 0 = ground level.
}

\seealso{
Use \code{\link{GenerateGauss}} to compute the three-dimensional tensor array object containing the wind speed data. Use \code{\link{GaussWS}} for a convenience function to look-up the values from the returned array.
}

\examples{
QuickGauss3D(100, 1, 100)
QuickGauss3D(200, -40, 120)
QuickGauss3D(50, 40, 70)
}

\keyword{GenerateGauss}
\keyword{GaussWS}
\keyword{QuickGauss3D}
