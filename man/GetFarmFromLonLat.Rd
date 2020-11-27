\name{GetFarmFromLonLat}
\alias{GetFarmFromLonLat}

\title{Accesses full FarmData set and returns a compatible list object for the requested location.}

\description{
This function accepts a pair longitude (decimal system, east) and latitude (decimal system, north) coordinates as well a a desired edge length in m and if valid, returns a list object similar to \code{\link{FarmData}}, but at the specified location and with the specified size.
}

\usage{
GetFarmFromLonLat(Top, Left, EdgeLen, DoPlot = FALSE)
}

\arguments{
	\item{Top}{
		longitude value (decimal system, east).
	}
	\item{Left}{
		latitude value (decimal system, north).
	}
\item{EdgeLen}{
		edge length in meters. Note that the returned farm data object is always a square area and thus, contains square matrices.
	}
\item{DoPlot}{
		optionally plots the annual energy production (AEP) 'landscape' in the returned object using a \code{\link{terrain.colors}} coloring scheme, together with a color key.
	}
}

\details{
Requires that the full \code{\link{FarmData}} dataset is loaded. See \code{\link{AcquireData}} on how to obtain it.
}

\value{
Returns a list object following the structure of \code{\link{FarmData}}, but only containing the farm area starting from the top-left point specified and as many 'tiles' required to meet the desired edge length at 200 m tiles resolution.
}

\author{
Carsten Croonenbroeck
}

\seealso{
See \code{\link{Index2GK}} to convert index coordinates to Gauss-Kruger coordinates, and \code{\link{GK2LonLat}} to convert Gauss-Kruger coordinates to longitude/latitude coordinates.
}

\examples{
# This will return a farm at the specified location, edge length 5,000 m (5 km).
# Requires full data set to be loaded.
\dontrun{
MyFarm <- GetFarmFromLonLat(51.49594, 11.58818, 5000, DoPlot = FALSE)
}
}

\keyword{Geolocation}
\keyword{GetFarmFromLonLat}
