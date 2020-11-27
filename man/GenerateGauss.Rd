\name{GenerateGauss}
\alias{GenerateGauss}

\title{For an incoming wind speed at reference height, this function computes a 3D tensor object containing Gaussian model based wind speeds.}

\description{
This function pre-computes Gaussian wind speeds and stores them in a 3D array, similar to voxels. Using that 'table', wind speeds can be looked up very quickly, which makes Gaussian wake feasible during WFLO runs.
}

\usage{
GenerateGauss(u = 8, refHeight = 10, maxX = 500, resY = 100,
resZ = 100, Verbose = TRUE)
}

\arguments{
	\item{u}{
		measured wind speed at reference height. Will mostly be measured in meters per second.
	}
	\item{refHeight}{
		reference height in meters. This is the height at which the incoming wind speed u is measured.
	}
	\item{maxX}{
		the number of steps down the x axis for which to compute the model.
	}
	\item{resY}{
		the number of steps along the y axis for which to compute the model. Note that as y may take negative values, the resolution space should be chosen not too small, here. If, e.g., resY = 100, this means that y may take values from -50 to 50, which may be too low a resolution in some cases.
	}
	\item{resZ}{
		the number of steps up the z axis for which to compute the model.
	}
	\item{Verbose}{
		selectes whether the function displays status reports during computation, as it may take some time, dependent on the resolution setting.
	}
}

\details{
Users may choose to compute a rather fine resolution run over night and then save the returned object so it can be loaded in future sessions. The Gaussian wake model is loosely based on the initial contribution by Bastankhah & Porte-Agel (2014).
}

\value{
\code{GenerateGauss} returns the three-dimensional array containing wind speeds.
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
Use \code{\link{GaussWS}} for a convenience function to look-up the values from the returned array. See \code{\link{QuickGauss3D}} for the same algorithm, immediately returning the wind speed at one single point only.
}

\examples{
\dontrun{
Gauss <- GenerateGauss(maxX = 500, resY = 1000, resZ = 1000)
}
}

\keyword{GenerateGauss}
\keyword{GaussWS}
\keyword{QuickGauss3D}
