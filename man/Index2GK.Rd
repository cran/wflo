\name{Index2GK}
\alias{Index2GK}

\title{Converts FarmData indices to Gauss-Kruger coordinates.}

\description{
This function accepts indices that are within the FarmData space and converts them to a pair of Gauss-Kruger (GK) coordinates.
}

\usage{
Index2GK(X, Y)
}

\arguments{
	\item{X}{
		a numeric value for the 'X' value in the FarmData space.
	}
	\item{Y}{
		a numeric value for the 'Y' value in the FarmData space.
	}
}

\details{
Assumes that the GK zone is EPSG 31467. References to indices are valid if the full data set is loaded, see \code{\link{FarmData}} and \code{\link{AcquireData}}.
}

\value{
For index values between 1 and 3,250 for X and between 1 and 4,400 for Y, returns Gauss-Kruger coordinates valid in the EPSG 31467 coordinate system.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{GK2Index}} for the inverse.
}

\examples{
Index2GK(1, 4400) # Will return c(3280000, 5230000), the lower left point.
Index2GK(3250, 1) # Will return c(3929800, 6109800), the top right point.
}

\keyword{Geolocation}
\keyword{Index2GK}
