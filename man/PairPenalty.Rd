\name{PairPenalty}
\alias{PairPenalty}

\title{Returns the Jensen wake penalty factor for a pair of turbines.}

\description{
As seen from a turbine in the wind farm, computes the wake penalty factor for another turbine in that farm.
}

\usage{
PairPenalty(x1, y1, x2, y2, Dirs, SDs)
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
	\item{Dirs}{
		a matrix containing average yearly wind directions. Usually, the third element of the list object \code{\link{FarmData}} will be used as this matrix.
	}
	\item{SDs}{
		a matrix containing average yearly wind direction standard deviations. Usually, the fourth element of the list object \code{\link{FarmData}} will be used as this matrix.
	}
}

\details{
First, this function uses \code{\link{GetAngle}} to compute the angle between the two points provided, as seen from point 2's point of view. It then obtains the wind direction at point 2 using \code{\link{GetDirInfo}}. After that, the distance between the two points is computed. With it, the wake cone is computed using \code{\link{JensenAngle}} to check whether point 2 is in point 1's wake using \code{\link{IsInWake}}. If that is the case, \code{\link{JensenFactor}} is used to compute the penalty factor.
}

\value{
\code{PairPenalty} returns a single number between 0 and 1. If point 2 is not in the wake of point 1, the function returns 1.
}

\author{
Carsten Croonenbroeck
}

\seealso{
Use \code{\link{JensenFactor}} to see how this function operates. See \code{\link{FarmVars}} for the data object.
}

\examples{
Dirs <- FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
e$FarmVars$StartPoint:e$FarmVars$EndPoint]
SDs <- FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
e$FarmVars$StartPoint:e$FarmVars$EndPoint]
PairPenalty(0.9, 0.8, 0.6, 0.9, Dirs, SDs)
## Weak wake penalty
PairPenalty(0.1, 0.1, 0.6, 0.9, Dirs, SDs)
## No wake penalty
}

\keyword{Jensen}
\keyword{Profit}
