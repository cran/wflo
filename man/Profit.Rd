\name{Profit}
\alias{Profit}

\title{Computes the economic profit for a given wind farm layout configuration}

\description{
This is the bread-and-butter function to this package. It takes a set of points (turbine locations) and based on the adjusted yields in the field specified via the \code{\link{FarmVars}} settings object, checks whether the layout is valid, computes the Jensen penalties, generates the total marketable power production of this layout, takes cost and sale price into consideration and finally computes the farm's economic profit. Note, however, that since most numeric optimizers by default operate as a minimizer, \code{Profit} returns the negative profit, i.e. a negative number for positive profit values and a positive number if cost is greater than revenue.
}

\usage{
Profit(X)
}

\arguments{
	\item{X}{
		must be a numeric vector containing an even number of values, at least two. If the length of the vector is N, then n = N / 2 is the number of points. The first values 1,...,n are interpreted as x coordinates and the subsequent n + 1,...,N values are the y coordinates.
	}
}

\details{
For maximum convenience and compatibility with numeric optimizers of most kinds, this function expects nothing but the problem vector (the set of points) as a parameter. This, however, requires that the \code{\link{FarmData}} dataset as well as the \code{\link{FarmVars}} object is present defining additional settings.\cr
As \code{Profit} returns the negative profit for compatibility with minimizers, reverse the sign for actual profit values.\cr
\code{Profit} requires the optimizer to accept box constraints, as \code{Profit} can not compute values outside the wind farm boundary (the function's domain). If your optimizer does not accept box constraints, embed \code{Profit} into a wrapper function that returns the sum of costs if it least one point is outside the boundary.
}

\value{
\code{Profit} returns a single number. The result is the negative profit for any valid setting of points and consequently, the sum of all costs for invalid settings.
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
##########################################
#Will provide a somewhat better result
#Necessary to install pso
\donttest{
if(requireNamespace("pso")){
library(pso)
NumTurbines <- 4
Result <- psoptim(par = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
Result
PlotResult(Result)}}
##########################################
#Simple wrapper function for optimizers not accepting box constraints
NumTurbines <- 4
lower <- rep(0, NumTurbines * 2)
upper <- rep(1, NumTurbines * 2)
Wrapper <- function(X)
{
	xSel <- seq(from = 1, to = length(X) - 1, by = 2)
	x <- X[xSel]
	y <- X[xSel + 1]

	if (any(x < lower) | any(x > upper) | any(y < lower) | any(y > upper))
	{
		return(sum(rep(FarmVars$UnitCost, length(x))))
	}

	return(Profit(X))
}
\dontrun{
set.seed(1357)
Result <- optim(par = runif(NumTurbines * 2), fn = Wrapper, method = "SANN")
Result
PlotResult(Result)}
}

\keyword{Jensen}
\keyword{Profit}
