\name{ShowWakePenalizers}
\alias{ShowWakePenalizers}

\title{Visualizes the points causing/'suffering' from wake effects.}

\description{
After optimization, this function draws a reduced 'field' and shows which points cause wake penalties on which points. 'Causers' are displayed in grey, 'sufferers', i.e. points that are in (possibly multiple) wakes of (several) causers are displayed in gold.
}

\usage{
ShowWakePenalizers(Result, Cones = TRUE, VectorField = TRUE)
}

\arguments{
	\item{Result}{
		must be an optimization result as returned by an optimizer such as \code{\link{optim}}. Usually, this will be a \code{list} at least containing '$par'. Most optimizers comply with this R standard. If your optimizer does not, wrap its result into a list containing an object '$par' containing the optimization result (i.e., a vector of points).
	}
\item{Cones}{
		triggers whether the Jensen cones are to be displayed as well. Defaults to \code{TRUE}.
	}
\item{VectorField}{
		controls whether or not to impose the wind directions vector field (arrows). Defaults to \code{TRUE}.
	}
}

\details{
Enables the researcher to inspect the optimization result in an a-posteriori analysis of Jensen's wake model.
}

\value{
\code{ShowWakePenalizers} invisibly returns a square matrix (dimension: 'number of turbines') of wake penalty pairs: Columns are sufferers, rows are causers. For each row/column pair unequal to one, the turbine in row i sheds wake on the turbine in column j. The actual numbers are those returned by \code{\link{JensenFactor}}.
}

\author{
Carsten Croonenbroeck
}

\seealso{
Use \code{\link{PlotResult}} to visualize the optimization result.
}

\examples{
\dontrun{
#Will show that turbine 1 sheds wake on turbine 5.
NumTurbines <- 8
set.seed(1357)
Result <- optim(par = runif(NumTurbines * 2), fn = Profit,
  method = "L-BFGS-B", lower = rep(0, NumTurbines * 2),
  upper = rep(1, NumTurbines * 2))
ShowWakePenalizers(Result)
}}

\keyword{Jensen}
\keyword{Wake}
