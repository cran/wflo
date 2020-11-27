\name{MosettiTurbineCost}
\alias{MosettiTurbineCost}

\title{Returns Mosetti's cost model for a number of turbines.}

\description{
This function returns the dimensionless cost based on Mosetti et al. (1994).
}

\usage{
MosettiTurbineCost(N)
}

\arguments{
	\item{N}{
		number of turbines for which cost is to be computed.
	}
}

\details{
Wind farm cost does not increase linearly with the number of turbines. Economies of scale lead to disproportionate cost increases. This function returns the dimensionless cost based on Mosetti et al. (1994).
}

\value{
\code{MosettiTurbineCost} returns a dimensionless cost factor.
}

\references{
G. Mosetti, C. Poloni, B. Diviacco (1994). Optimization of wind turbine positioning in large windfarms by means of a genetic algorithm, Journal of Wind Engineering and Industrial Aerodynamics, 51(1), 105-116.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{Cost}} for a plug-in cost function for profit.
}

\examples{
MosettiTurbineCost(1)
MosettiTurbineCost(20)
}

\keyword{MosettiTurbineCost}
