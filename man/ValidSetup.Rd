\name{ValidSetup}
\alias{ValidSetup}

\title{Checks whether all turbine locations provided satisfy the minimum distance criterion.}

\description{
For a set of turbine locations represented by \code{x} and \code{y}, checks whether all possible pairs satisfy the minimum distance criterion.
}

\usage{
ValidSetup(x, y)
}

\arguments{
	\item{x}{
		must be a numeric vector of at least two values, contains the 'x' location(s) of turbines.
	}
	\item{y}{
		must be a numeric vector of at least two values, contains the 'y' location(s) of turbines.
	}
}

\value{
\code{ValidSetup} returns \code{TRUE} if all pairs of turbines are at least as far away from each other as 'MinDist' from the \code{\link{FarmVars}} settings object requests, or \code{FALSE}, else.
}

\author{
Carsten Croonenbroeck
}

\examples{
ValidSetup(c(0.5, 0.7), c(0.7, 0.9))
## Returns TRUE if FarmVars$MinDist is at its standard value (0.1).
}

\keyword{Profit}
