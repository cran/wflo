\name{GetDirInfo}
\alias{GetDirInfo}

\title{Returns average wind direction and direction standard deviation for a turbine's location.}

\description{
For a turbine's location represented by \code{x} and \code{y}, looks up the wind direction from the matrix \code{Dirs} and the corresponding standard deviation from matrix \code{SDs}. Internally transforms coordinates of \code{x} and \code{y} from problem space (usually unit square) to the matrix space of \code{Dirs} and \code{SDs}, respectively.
}

\usage{
GetDirInfo(x, y, Dirs, SDs)
}

\arguments{
	\item{x}{
		must be a single value containing the 'x' location of a turbine in problem space.
	}
	\item{y}{
		must be a single value containing the 'y' location of a turbine in problem space.
	}
	\item{Dirs}{
		a matrix containing average yearly wind directions. Usually, the third element of the list object \code{\link{FarmData}} will be used as this matrix.
	}
	\item{SDs}{
		a matrix containing average yearly wind direction standard deviations. Usually, the fourth element of the list object \code{\link{FarmData}} will be used as this matrix.
	}
}

\value{
\code{GetDirInfo} returns a vector of two elements, the first being the average wind direction in degrees and the second being the corresponding standard deviation. Note that degrees are meant in the arithmetic degrees system (0° being east, ascending counterclockwise). To transform into an azimuth system (0° being north, ascending clockwise), use function \code{\link{Geo2Ari}}. Also note that wind directions are meant to denote 'where the wind is going to' rather than 'where the wind is coming from'.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{Profit}} to see where to use \code{GetDirInfo}, \code{\link{Yield}} for a similar function for adjusted yield. See \code{\link{FarmVars}} for the data object.
}

\examples{
Dirs <- FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
e$FarmVars$StartPoint:e$FarmVars$EndPoint]
SDs <- FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
e$FarmVars$StartPoint:e$FarmVars$EndPoint]
GetDirInfo(0.5, 0.7, Dirs, SDs)
## Returns wind direction and standard deviation for the given location
## and provided the matrices given.
}

\keyword{Wind direction}
\keyword{Standard deviation}
