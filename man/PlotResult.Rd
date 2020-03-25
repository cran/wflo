\name{PlotResult}
\alias{PlotResult}

\title{Visualizes the wind farm layout optimization result.}

\description{
This function draws the adjusted yields of the wind farm under investigation using \code{\link{image}}, superimposes a contour plot using \code{\link{contour}} and arrows for the wind directions using \code{\link{ImposeVectorField}} and then draws the points for the turbines' locations (aligned to their respective raster grid centers).
}

\usage{
PlotResult(Result)
}

\arguments{
	\item{Result}{
		must be an optimization result as returned by an optimizer such as \code{\link{optim}}. Usually, this will be a \code{list} at least containing '$par'. Most optimizers comply with this R standard. If your optimizer does not, wrap its result into a list containing an object '$par' containing the optimization result (i.e., a vector of points).
	}
}

\details{
For maximum convenience and compatibility with numeric optimizers of most kinds, this function expects nothing but the usual optimization result \code{list}. This, however, requires that the \code{\link{FarmData}} dataset as well as the \code{\link{FarmVars}} object is present defining additional settings.
}

\value{
\code{PlotResult} returns nothing.
}

\author{
Carsten Croonenbroeck
}

\seealso{
Use \code{\link{PlotResult}} to visualize the optimization result.
}

\examples{
#Will not provide a very good result
NumTurbines <- 4
set.seed(1357)
Result <- optim(par = runif(NumTurbines * 2), fn = Profit,
  method = "L-BFGS-B", lower = rep(0, NumTurbines * 2),
  upper = rep(1, NumTurbines * 2))
Result
PlotResult(Result)
}

\keyword{Jensen}
\keyword{Profit}
