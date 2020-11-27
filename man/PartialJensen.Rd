\name{PartialJensen}
\alias{PartialJensen}

\title{Computes the partial Jensen wake effect.}

\description{
If a rotor disc of a first turbine is only partially covered by a wake cone of a second turbine, the penalty must be adjusted accordingly. This function returns partial penalties, if partial wake applies, or full penalty, if not.
}

\usage{
PartialJensen(x2, y2, x1, y1, Dirs, SDs, Elev, Z = NULL, DrawTop = FALSE,
DrawFront = FALSE, DrawSide = FALSE)
}

\arguments{
	\item{x1}{
		must be a single value. Provide the \code{x} location of the first turbine.
	}
	\item{y1}{
		must be a single value. Provide the \code{y} location of the first turbine.
	}
	\item{x2}{
		must be a single value. Provide the \code{x} location of the second turbine.
	}
	\item{y2}{
		must be a single value. Provide the \code{y} location of the second turbine.
	}
	\item{Dirs}{
		a matrix containing average yearly wind directions. Usually, the third element of the list object \code{\link{FarmData}} will be used as this matrix.
	}
	\item{SDs}{
		a matrix containing average yearly wind direction standard deviations. Usually, the fourth element of the list object \code{\link{FarmData}} will be used as this matrix.
	}
	\item{Elev}{
		a matrix containing elevations. Usually, the fifth element of the list object \code{\link{FarmData}} will be used as this matrix.
	}
	\item{Z}{
		accepts a vector of two representing the heights of the turbines at the two points in meters. If \code{NULL} (the default), \code{PartialJensen} uses \code{e$FarmVars$z} as the heights information for both points. In all cases, the heights are added to the terrain topography information from \code{Elev}.
	}
	\item{DrawTop}{
		If \code{TRUE}, draws a top view of the partial wake situation. Defaults to \code{FALSE}. Should not be \code{TRUE} at the same time as \code{DrawFront} or \code{DrawSide}.
	}
	\item{DrawFront}{
		If \code{TRUE}, draws a front view of the partial wake situation. Defaults to \code{FALSE}. Should not be \code{TRUE} at the same time as \code{DrawTop} or \code{DrawSide}.
	}
	\item{DrawSide}{
		If \code{TRUE}, draws a side view of the partial wake situation. Defaults to \code{FALSE}. Should not be \code{TRUE} at the same time as \code{DrawTop} or \code{DrawFront}.
	}
}

\details{
The function first checks whether there is partial coverage. If so, it adjusts the penalty internally given by \code{\link{PairPenalty}} and if not, returns the full penalty. Therefore, the function is imposed over an existing Jensen model and refines it.\cr
Note that the penalty is the deduction to wind speed. It applies to wind power by its third power, so the user is responsible to take it to its cube himself. \code{\link{Profit}} does that automatically internally.
}

\value{
\code{PartialJensen} returns a penality between 0 and 1.
}

\references{
Frandsen, S. (1992). On the wind speed reduction in the center of large clusters of wind turbines. Journal of Wind Engineering and Industrial Aerodynamics, 39(1-3), pp. 251-265, https://doi.org/10.1016/0167-6105(92)90551-K.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{JensenTrapezoid}} to check whether there are wake effects present. \code{\link{FarmVars}} for the data object. \code{\link{PairPenalty}} for the non-partial wake penalty.
}

\examples{
P1 <- c(0.5868695, 0.9722714)
P2 <- c(0.4827957, 0.9552658)

if (exists("FarmData", envir = e, inherits = FALSE))
{
  Dirs <- e$FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
    e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  SDs <- e$FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
    e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  Elev <- e$FarmData[[5]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
    e$FarmVars$StartPoint:e$FarmVars$EndPoint]
} else
{
  Dirs <- FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
    e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  SDs <- FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
    e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  Elev <- FarmData[[5]][e$FarmVars$StartPoint:e$FarmVars$EndPoint,
    e$FarmVars$StartPoint:e$FarmVars$EndPoint]
}

## First, compute the non-partial penalty:
Penal <- PairPenalty(P2[1], P2[2], P1[1], P1[2], Dirs, SDs)

## Then, correct it for partial coverage:
Penal2 <- PartialJensen(P2[1], P2[2], P1[1], P1[2], Dirs, SDs, Elev)

## Now draw a top view:
Penal2 <- PartialJensen(P2[1], P2[2], P1[1], P1[2], Dirs, SDs, Elev, DrawTop = TRUE)

## Now draw a front view:
Penal2 <- PartialJensen(P2[1], P2[2], P1[1], P1[2], Dirs, SDs, Elev, DrawFront = TRUE)
## Note that the downwind point is somewhat elevated, seems 'right'
## from the upwind point of view, and far away (rotor disc seems smaller).

## Now draw a side view:
Penal2 <- PartialJensen(P2[1], P2[2], P1[1], P1[2], Dirs, SDs, Elev, DrawSide = TRUE)

## For elevation, see
Height(P1[1], P1[2], Elev) # (upwind point)
## and
Height(P2[1], P2[2], Elev) # (downwind point)

## Instead, for the next example, the downwind point is closer
## to the upwind point (larger impression of rotor disc), 'left'
## of it and lower in terms of elevation:
x1 <- 0.5
y1 <- 0.5
x2 <- 0.45
y2 <- 0.478
Penal <- PairPenalty(x2, y2, x1, y1, Dirs, SDs)
Penal2 <- PartialJensen(x2, y2, x1, y1, Dirs, SDs, Elev, DrawTop = TRUE)
Penal2 <- PartialJensen(x2, y2, x1, y1, Dirs, SDs, Elev, DrawFront = TRUE)
## For elevation, see
Height(x1, y1, Elev) # (upwind point)
## and
Height(x2, y2, Elev) # (downwind point)
}

\keyword{PartialJensen}
\keyword{Profit}
