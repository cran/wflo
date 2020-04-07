\name{PointInPolygon}
\alias{PointInPolygon}

\title{Checks whether a point is inside a polygon.}

\description{
Used the famous Jensen test to check whether a provided point is inside a polygon, the latter defined by a set of points.
}

\usage{
PointInPolygon(PointMat, TestPoint)
}

\arguments{
	\item{PointMat}{
		a 2 x 4 matrix, four corner points of each x and y coordinates.
	}
	\item{TestPoint}{
		a vector of two containing x and y coordinates of a point to test.
	}
}

\value{
\code{PointInPolygon} returns \code{TRUE} if the point is inside the polygon, or \code{FALSE}, else.
}

\references{
Jensen, N. O. (1983). A note on wind generator interaction. Roskilde: Risø National Laboratory. Risø-M, No. 2411
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{JensenTrapezoid}} to compute a Jensen trapezoid.
}

\examples{
set.seed(1357)
Angle <- runif(n = 1, min = 0, max = 360)
Point <- runif(n = 2, min = 0, max = 1)
Dist <- rnorm(n = 1, mean = 200, sd = 200)

MyTrapezoid = JensenTrapezoid(Angle, Point, Dist)

plot(x = MyTrapezoid[1, ], y = MyTrapezoid[2, ], xlab = "", ylab = "")
lines(x = c(MyTrapezoid[1, 1], MyTrapezoid[1, 2]), y = c(MyTrapezoid[2, 1], MyTrapezoid[2, 2]))
lines(x = c(MyTrapezoid[1, 2], MyTrapezoid[1, 3]), y = c(MyTrapezoid[2, 2], MyTrapezoid[2, 3]))
lines(x = c(MyTrapezoid[1, 3], MyTrapezoid[1, 4]), y = c(MyTrapezoid[2, 3], MyTrapezoid[2, 4]))
lines(x = c(MyTrapezoid[1, 4], MyTrapezoid[1, 1]), y = c(MyTrapezoid[2, 4], MyTrapezoid[2, 1]))

NumTest <- 50

xTest <- runif(n = NumTest, min = min(MyTrapezoid[1, ]), max = max(MyTrapezoid[1, ]))
yTest <- runif(n = NumTest, min = min(MyTrapezoid[2, ]), max = max(MyTrapezoid[2, ]))

for (i in 1:NumTest)
{
	ThisPoint <- c(xTest[i], yTest[i])
	if (PointInPolygon(MyTrapezoid, ThisPoint))
	{
		points(ThisPoint[1], ThisPoint[2], pch = 16, col = "green")
	} else
	{
		points(ThisPoint[1], ThisPoint[2], pch = 16, col = "red")
	}
}
}

\keyword{Jensen}
\keyword{Trapezoid}
