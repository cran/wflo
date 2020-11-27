\name{ShowWakePenalizers}
\alias{ShowWakePenalizers}

\title{Visualizes the points causing/'suffering' from wake effects.}

\description{
After optimization, this function draws a reduced 'field' and shows which points cause wake penalties on which points. 'Causers' are displayed in grey, 'sufferers', i.e. points that are in (possibly multiple) wakes of (several) causers are displayed in gold.
}

\usage{
ShowWakePenalizers(Result, Cones = TRUE, All = FALSE, VectorField = TRUE,
ImposeMST = FALSE, NS = FALSE, NSOnly = FALSE, Exaggerate = TRUE,
DoubleLine = TRUE, Frames = 100, Alpha = 0.5, MaxContrast = TRUE, Soften = TRUE)
}

\arguments{
	\item{Result}{
		must be an optimization result as returned by an optimizer such as \code{\link{optim}}. Usually, this will be a \code{list} at least containing '$par'. Most optimizers comply with this R standard. If your optimizer does not, wrap its result into a list containing an object '$par' containing the optimization result (i.e., a vector of points).
	}
\item{Cones}{
		triggers whether the Jensen cones are to be displayed as well. Defaults to \code{TRUE}.
	}
\item{All}{
		triggers whether all Jensen cones are to be displayed, instead of just causers and sufferers. Defaults to \code{FALSE}, as \code{TRUE} results in a mess most of the times.
	}
\item{VectorField}{
		controls whether or not to impose the wind directions vector field (arrows). Defaults to \code{TRUE}.
	}
\item{ImposeMST}{
		if \code{TRUE}, computes and imposes a minimum spanning tree for the vertices, useful for cable installation purposes. Defaults to \code{FALSE}.
	}
\item{NS}{
		if \code{TRUE}, imposes a computationally intensive Navier-Stokes based wake simulation image over the plot. Defaults to \code{FALSE}.
	}
\item{NSOnly}{
		if \code{TRUE}, imposes just the Navier-Stokes wake simulation (no Jensen wakes). Defaults to \code{FALSE}.
	}
\item{Exaggerate}{
		if \code{TRUE}, displays turbines as squares in the Navier-Stokes visualization. The turbine is then still to scale, but in reality the rotor is a disc instead of a cube, of course. This overly emphasizes wakes, which may be desired in most cases. Thus, defaults to \code{TRUE}. Set to \code{FALSE} for a more realistic look of the rotor discs, but this will make the wakes appear less severe.
	}
\item{DoubleLine}{
		applies only if \code{Exaggerate == TRUE}. The rotor disc is then a disc instead of a cube (or line instead of square), but the line will be two instead of one pixel wide. This option is a trade-off between \code{Exaggerate == TRUE} and \code{Exaggerate == FALSE}. Set to \code{TRUE} (the default) if \code{Exaggerate == FALSE} to apply. Note that only if \code{Exaggerate == FALSE} and \code{DoubleLine == FALSE}, the rotor disc is a one pixel wide line, the most realistic setting which however shows only weak wakes.
	}
\item{Frames}{
		number of frames to pre-compute for the Navier-Stokes simulation. 100 is a good value. Less makes the computation faster, but more inaccurate. More than 100 in only necessary for very low wind speeds. Thus, defaults to \code{100}.
	}
\item{Alpha}{
		The Navier-Stokes wake simulation is imposed semi-transparent over the usual image. This controls the transparency. Alpha = 0 means transparent, Alpha = 1 means opaque. Defaults to \code{0.5}.
	}
\item{MaxContrast}{
		increases the contrast to maximum (this is equivalent to histogram stretching). Defaults to \code{TRUE}.
	}
\item{Soften}{
		applies a Gaussian soften filter to the Navier-Stokes image to have a more homogeneous visualization. Defaults to \code{TRUE}.
	}
}

\details{
Enables the researcher to inspect the optimization result in an a-posteriori analysis of Jensen's wake model. Parts of the Navier-Stokes code are inspired by Stam (2003).
}

\value{
\code{ShowWakePenalizers} invisibly returns a square matrix (dimension: 'number of turbines') of wake penalty pairs: Columns are sufferers, rows are causers. For each row/column pair unequal to one, the turbine in row i sheds wake on the turbine in column j. The actual numbers are those returned by \code{\link{JensenFactor}}.
}

\references{
Stam, S. (2003). Real-Time Fluid Dynamics for Games. \url{https://www.researchgate.net/publication/2560062_Real-Time_Fluid_Dynamics_for_Games}, 1-17
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

###############

#This may take some time. A progress bar is displayed.
Result = list(par = e$FarmVars$BenchmarkSolution)
ShowWakePenalizers(Result, All = TRUE, NS = TRUE, VectorField = FALSE, Alpha = 0.7)

###############

#Generates a few frames of an animation.
Result <- list(par = e$FarmVars$BenchmarkSolution)
MakeFrames <- function(Frames)
{
	for (i in 100:(100 + Frames - 1))
	{
		bmp(file = paste(i, ".bmp", sep = ""))
		ShowWakePenalizers(Result, All = TRUE, NS = TRUE, Frames = i)
		dev.off()
	}
}
system.time(MakeFrames(30)) #May take an hour.

}
}

\keyword{Jensen}
\keyword{Wake}
