\name{GK2Index}
\alias{GK2Index}

\title{Converts Gauss-Kruger coordinates to FarmData indices.}

\description{
This function accepts Gauss-Kruger (GK) coordinates that are within the FarmData space and converts them to a pair of indices that point to locations within the FarmData matrices.
}

\usage{
GK2Index(X, Y)
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
Assumes that the GK zone is EPSG 31467. References to indices valid if the full data set is loaded, see \code{\link{FarmData}} and \code{\link{AcquireData}}.
}

\value{
Returns index values between 1 and 3,250 for X and between 1 and 4,400 for Y.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{Index2GK}} for the inverse.
}

\examples{
GK2Index(3280000, 5230000) # Will return c(1, 4400), the lower left point.
GK2Index(3929800, 6109800) # Will return c(3250, 1), the top right point.
}

\keyword{Geolocation}
\keyword{GK2Index}
