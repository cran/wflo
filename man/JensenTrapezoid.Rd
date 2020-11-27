\name{JensenTrapezoid}
\alias{JensenTrapezoid}

\title{Computes the four corner points of a Jensen trapezoid (or cone).}

\description{
Provided a wind direction, a point and a downwind length, computes the corner points of a Jensen trapezoid (in the literature oftentimes called a 'cone').
}

\usage{
JensenTrapezoid(WindDir, Point, x, Margin = TRUE)
}

\arguments{
	\item{WindDir}{
		unique value, a wind direction in degrees.
	}
	\item{Point}{
		a vector of two containing x and y coordinates of a point. Usually, both will be between 0 and 1.
	}
	\item{x}{
		downwind distance of the cone in meters.
	}
	\item{Margin}{
		specified whether to add a small margin to the distance parameter (i.e. internally, compute x * 1.1). Defaults to TRUE.
	}
}

\value{
\code{JensenTrapezoid} returns a matrix of four columns (the four corner points) and two rows (x and y coordinates).
}

\references{
Jensen, N. O. (1983). A note on wind generator interaction. Roskilde: Riso National Laboratory. Riso-M, No. 2411
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{JensenAngle}} to compute 'cone' angle information, \code{\link{PointInPolygon}} for a test whether a point is inside a polygon.
}

\examples{
MyTrapezoid <- JensenTrapezoid(45, c(0.5, 0.5), 500)
}

\keyword{Jensen}
\keyword{Trapezoid}
