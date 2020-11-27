\name{Height}
\alias{Height}

\title{Returns the elevation of a turbine's location.}

\description{
For a turbine's location represented by \code{x} and \code{y}, looks up the elevation from the matrix \code{Elev}. Internally transforms coordinates of \code{x} and \code{y} from problem space (usually unit square) to the matrix space of \code{Elev}.
}

\usage{
Height(x, y, Elev)
}

\arguments{
	\item{x}{
		must be a single value containing the 'x' location of a turbine in problem space.
	}
	\item{y}{
		must be a single value containing the 'y' location of a turbine in problem space.
	}
	\item{Elev}{
		a matrix containing heights. Usually, the fifth element of the list object \code{\link{FarmData}} will be used as this matrix.
	}
}

\details{
Height is a convenience function for looking up heights as required e.g. for a partial Jensen wake model, independent from the actual size of the area under investigation.
}

\value{
\code{Height} returns a single value, the elevation in meters.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{Profit}} to see where to use \code{Height}. See \code{\link{FarmData}} for the data set.
}

\examples{
## Returns adjusted yield for the given location.
P <- c(0.5868695, 0.9722714)
Height(P[1], P[2], FarmData[[5]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
    e$FarmVars$StartPoint:e$FarmVars$EndPoint])
}

\keyword{Height}
\keyword{Profit}
