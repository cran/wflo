\name{GK2LonLat}
\alias{GK2LonLat}

\title{Converts Gauss-Kruger coordinates to longitude/latitude coordinates.}

\description{
This function accepts a pair of Gauss-Kruger (GK) coordinates and converts them to longitude (decimal system, east) and latitude (decimal system, north) coordinates by performing re-projection.
}

\usage{
GK2LonLat(X, Y)
}

\arguments{
	\item{X}{
		a numeric value for the 'right' value in the GK coordinate system.
	}
	\item{Y}{
		a numeric value for the 'top' value in the GK coordinate system.
	}
}

\details{
Assumes that the GK zone is EPSG 31467.
}

\value{
Returns a vector of two values where longitude is first, latitude second.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{LonLat2GK}} for the inverse.
}

\examples{
LonLat <- GK2LonLat(3702793, 5998319) # Will return c(12.09750, 54.07548).
}

\keyword{Geolocation}
\keyword{GK2LonLat}
