\name{Cost}
\alias{Cost}

\title{Stub for a turbine's cost function.}

\description{
A function that returns the yearly installation costs for a given set of turbines (provide x and y for the turbines' locations). In its present form it only returns the 'UnitCost' value from the \code{\link{FarmVars}} settings object per turbine.
}

\usage{
Cost(x, y)
}

\arguments{
	\item{x}{
		can be a single value or a numeric vector of values, contains the 'x' location(s) of turbines.
	}
	\item{y}{
		can be a single value or a numeric vector of values, contains the 'y' location(s) of turbines.
	}
}

\details{
Note that \code{x} and \code{y} should both be of length \code{n}, i.e. the numbers of values they contain should match the number of turbines in the current wind farm layout problem.\cr
This function is a stub and can and should be replaced by something reasonable in an actual wind farm layout problem.
}

\value{
\code{Cost} returns a vector of values. The number of elements matches the length of \code{x} and \code{y}.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{Profit}} to see where to use \code{Cost}, \code{\link{Yield}} for a similar function for yearly yield.
}

\examples{
## Returns a vector of two, c(100000, 100000).
Cost(c(0.5, 0.7), c(0.2, 0.3))

## Replace the function by another function
## also called 'Cost', embedded in environment e.
## Also, see the vignette.
\donttest{
if(requireNamespace("pso")){
library(pso)
e$Cost <- function(x, y) #x, y \in R^n
{
	retVal <- rep(e$FarmVars$UnitCost, min(length(x), length(y)))
	retVal[x > 0.5] <- retVal[x > 0.5] * 2
	return(retVal)
}
set.seed(1357)
Result <- psoptim(par = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
Result
rm(Cost, envir = e)
}
}
}

\keyword{Cost}
\keyword{Profit}
