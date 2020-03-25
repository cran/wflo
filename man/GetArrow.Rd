\name{GetArrow}
\alias{GetArrow}

\title{Simple helper function for \code{\link{PlotResult}}.}

\description{
Given a point, an angle information (in arithmetic degrees), and a length information, computes start and end points of an arrow usable via the \code{\link{arrows}} function.
}

\usage{
GetArrow(BaseX, BaseY, Degrees, Frac = 25)
}

\arguments{
	\item{BaseX}{
		must be a single value containing the x value of a point which later will be the center of the arrow.
	}
	\item{BaseY}{
		must be a single value containing the y value of a point which later will be the center of the arrow.
	}
	\item{Degrees}{
		must be a single value containing the desired rotation degree of the arrow.
	}
	\item{Frac}{
		must be a single value containing the length of the arrow. Default is 25 and for convenience, this parameter should in most cases be identical to FarmVars$Width.
	}
}

\details{
This function will be used internally by \code{\link{PlotResult}}.
}

\value{
\code{GetArrow} returns a vector of four values representing x and y for the start point and x and y for the end point of an arrow (in that order).
}

\author{
Carsten Croonenbroeck
}

\seealso{
Use \code{\link{PlotResult}} to visualize the optimization result. See \code{\link{FarmVars}} for the data object.
}

\examples{
GetArrow(0.5, 0.5, 45)
#At c(0.5, 0.5), generates an arrow pointing in north-eastern direction.
}

\keyword{FarmVars}
\keyword{PlotResult}
