\name{ProfitContributors}
\alias{ProfitContributors}

\title{Computes profit contributions for all points in a setup solution.}

\description{
This function takes an optimizer result (or simple vector object) and computes the profit contribution for each point.
}

\usage{
ProfitContributors(Result)
}

\arguments{
	\item{Result}{
		must be an optimization result as returned by an optimizer such as \code{\link{optim}} or a vector. Usually, this will be a \code{list} at least containing '$par'.
	}
}

\details{
Neglecting a possibly invalid setup (caused by at least one point), computes the profit contributions. If the setup provided is valid, the sum of contributions is identical to the (absulute value of the) returned value of \code{\link{Profit}}.
}

\value{
\code{ProfitContributors} returns a matrix of two columns and n rows, whith n being the number of turbines. The first column is a sequence of 1:n, representing the turbine IDs, while the second column contains the actual profit contributions.
}

\author{
Carsten Croonenbroeck
}

\seealso{
See \code{\link{ValidSetup}} to see how a setup is categorized as valid or not. Use \code{\link{PlotResult}} to visualize the optimization result.
}

\examples{
#Prints a result and uses the profit contributions as labels.
NumTurbines <- 4
set.seed(1235)
Result <- optim(par = runif(NumTurbines * 2), fn = Profit,
  method = "L-BFGS-B", lower = rep(0, NumTurbines * 2),
  upper = rep(1, NumTurbines * 2))
MyLabels <- ProfitContributors(Result)
MyLabels
#PlotResult(Result, DoLabels = TRUE, Labels = MyLabels[, 2])

# Given a valid setup, this should be TRUE.
#identical(abs(Profit(Result$par)), sum(MyLabels[, 2]))
}

\keyword{Jensen}
\keyword{Profit}
