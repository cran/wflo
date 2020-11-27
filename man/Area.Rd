\name{Area}
\alias{Area}

\title{Computes the overlap area for the partial Jensen wake.}

\description{
The partial Jensen wake model requires the overlap area of the rotor disc at the downwind location with the wake disc. This function computes it.
}

\usage{
Area(a, b, d)
}

\arguments{
	\item{a}{
		must be a single value. Provide the radius of the first circle in meters.
	}
	\item{b}{
		must be a single value. Provide the radius of the second circle in meters.
	}
	\item{d}{
		must be a single value. Provide the distance between the two circles' centers.
	}	
}

\details{
If a turbine is downwind another turbine, the wake cone of that upwind turbine may only partially cover the rotor disc of the second turbine. For the partial wake model it is necessary to compute the covered area.
}

\value{
\code{Area} returns covered area in square meters.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{JensenTrapezoid}} to check whether there are wake effects present. \code{\link{FarmVars}} for the data object.
}

\examples{
Area(60, 40, 50)
# Returns 2930.279.
}

\keyword{Area}
\keyword{Profit}
