\name{ImposeVectorField}
\alias{ImposeVectorField}

\title{Simple helper function for \code{\link{PlotResult}}.}

\description{
Draws a set of arrows over an existing plot. Will be used internally by \code{\link{PlotResult}} to visualize the vector field of wind directions over a plot of the wind farm.
}

\usage{
ImposeVectorField(xNum = 5, yNum = 5, ColMain = "black",
  ColBand = "white", Frac = 25, DoSDs = TRUE)
}

\arguments{
	\item{xNum}{
		must be a single value containing the desired number of arrows horizontally.
	}
	\item{yNum}{
		must be a single value containing the desired number of arrows vertically.
	}
	\item{ColMain}{
		must be a single value containing the desired color for the main wind direction arrow.
	}
	\item{ColBand}{
		must be a single value containing the desired color for the secondary arrows representing the standard deviation arrows. Only used if \code{DoSDs = TRUE}.
	}
	\item{Frac}{
		must be a single value containing the desired length information for the arrows. Will be passed to \code{\link{GetArrow}} internally.
	}
	\item{DoSDs}{
		must be \code{TRUE} or \code{FALSE}. If \code{TRUE}, will not only draw one arrow for the main wind direction per point, but also two additional arrows (using color \code{ColBand}) representing the main wind direction plus/minus half the standard deviation.
	}	
}

\details{
This function will be used internally by \code{\link{PlotResult}}.\cr
\code{ImposeVectorField} requires \code{\link{FarmData}} to be present and an existing plot to impose the vector field over.
}

\value{
\code{ImposeVectorField} returns nothing.
}

\author{
Carsten Croonenbroeck
}

\seealso{
Use \code{\link{PlotResult}} to visualize the optimization result.
}

\examples{
plot(c(0, 1), c(0, 1))
ImposeVectorField(ColBand = "red")
}

\keyword{FarmVars}
\keyword{PlotResult}
