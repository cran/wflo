\name{Yield}
\alias{Yield}

\title{Returns yearly yield for a turbine's location.}

\description{
For a turbine's location represented by \code{x} and \code{y}, looks up the (adjusted) yield from the matrix \code{Adj}. Internally transforms coordinates of \code{x} and \code{y} from problem space (usually unit square) to the matrix space of \code{Adj}.
}

\usage{
Yield(x, y, Adj)
}

\arguments{
	\item{x}{
		must be a single value containing the 'x' location of a turbine in problem space.
	}
	\item{y}{
		must be a single value containing the 'y' location of a turbine in problem space.
	}
	\item{Adj}{
		a matrix containing adjusted yields. Usually, the first element of the list object \code{\link{FarmData}} will be used as this matrix.
	}
}

\details{
Adjusted yields are the projected yearly average yields dependent on wind speed, hub height and other settings at each point in the raster data. Annual Energy Production (AEP) at a specific location, weighted by a location quality correction factor, produces adjusted yields. This adjustment returns a better guess on the marketable yield at a specific point. For details on the data, see the data set description to this package.\cr
Note that \code{\link{Profit}} internally multiplies the outcome of \code{Yield} by \code{e$FarmVars$Price} to obtain revenue. Users who replace the function by \code{e$Yield} need to provide that manually, if revenue is desired.
}

\value{
\code{Yield} returns a single value.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{Profit}} to see where to use \code{Yield}, \code{\link{Cost}} for a similar function for yearly cost. \code{\link{FarmData}} for the data set.
}

\examples{
## Returns adjusted yield for the given location.
Adj <- FarmData[[1]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
e$FarmVars$StartPoint:e$FarmVars$EndPoint]
Yield(0.5, 0.7, Adj)

## Replace the function by another function
## also called 'Yield', embedded in environment e.
## Also, see the vignette.
\donttest{
e$Yield <- function(x, y, AEP) #x, y \in R
{
	return(x + y)
}
set.seed(1357)
Result <- pso::psoptim(par = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
Result
rm(Yield, envir = e)
}
}

\keyword{Yield}
\keyword{Profit}
