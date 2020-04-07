\name{GetAngle}
\alias{GetAngle}

\title{Returns the angle between two turbines.}

\description{
As seen from a point in the wind farm, computes the angle to another point in that farm.
}

\usage{
GetAngle(x1, y1, x2, y2)
}

\arguments{
	\item{x1}{
		must be a single value. Provide the \code{x} location of the first turbine.
	}
	\item{y1}{
		must be a single value. Provide the \code{y} location of the first turbine.
	}
	\item{x2}{
		must be a single value. Provide the \code{x} location of the second turbine.
	}
	\item{y2}{
		must be a single value. Provide the \code{y} location of the second turbine.
	}
}

\details{
From point 2's point of view, computes the angle to point 1.
}

\value{
\code{GetAngle} returns a single number between 0 and 360. If both points are identical, the return value is 0.
}

\author{
Carsten Croonenbroeck
}

\note{
Note that this function returns an angle in arithmetic degrees notation. To convert to azimuth notation, use \code{\link{Geo2Ari}}.
}

\seealso{
Use \code{\link{JensenAngle}} to compute the wake cone and with it, use \code{\link{JensenTrapezoid}} to see if another turbine B is in turbine A's wake.
}

\examples{
GetAngle(0.2, 0.2, 0.1, 0.1)
## Looking from point (0.1, 0.1) at point (0.2, 0.2), the angle is 45° (arithmetic).
}

\keyword{Jensen}
\keyword{Profit}
