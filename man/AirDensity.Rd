\name{AirDensity}
\alias{AirDensity}

\title{Provides convenient computations for air density.}

\description{
Implements a set of computations and constants to compute air density as a function of altitude (a.s.l.), temperature and relative humidity.
}

\usage{
AirDensity(Altitude = 0, Temperature = 20, phi = 0.76)
}

\arguments{
	\item{Altitude}{
		altitude in meters above sea level.
	}
	\item{Temperature}{
		temperature in degrees Celsius.
	}
	\item{phi}{
		relative humidity, must be within [0, 1].
	}
}

\details{
The function first computes air pressure at user provided target altitude based on the barometric formula. Then, saturation vapor pressure is computed using the Magnus formula. With both, the density of air is computed afterward.
}

\value{
\code{AirDensity} returns a single value representing air density in kg per cubic meter.
}

\author{
Carsten Croonenbroeck
}

\note{
This function returns valid values for negative altitudes, but only for temperature values between -45 and +60 degrees Celsius.
}

\seealso{
\code{\link{WindspeedLog}} and \code{\link{WindspeedHellmann}} to compute wind speed at target heights using two slightly different established models.
}

\examples{
AirDensity()
AirDensity(300, 25, 0.6)
AirDensity(0, 0, 0) # Standard conditions according to DIN 1343:1990,
                    # returns the expected value of 1.292.
}

\keyword{AirDensity}
