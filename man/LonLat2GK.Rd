\name{LonLat2GK}
\alias{LonLat2GK}

\title{Converts longitude/latitude coordinates to Gauss-Kruger coordinates.}

\description{
This function accepts a pair longitude (decimal system, east) and latitude (decimal system, north) coordinates and converts them to Gauss-Kruger (GK) coordinates by performing re-projection.
}

\usage{
LonLat2GK(Lon, Lat)
}

\arguments{
	\item{Lon}{
		longitude value (decimal system, east).
	}
	\item{Lat}{
		latitude value (decimal system, north).
	}
}

\details{
The resulting values are in GK zone EPSG 31467.
}

\value{
Returns a vector of two values where X is first, Y second.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{GK2LonLat}} for the inverse.
}

\examples{
LonLat2GK(12.09750, 54.07548) # Will return c(3702793, 5998320).
}

\keyword{Geolocation}
\keyword{LonLat2GK}
