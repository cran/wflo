\name{WindspeedLog}
\alias{WindspeedLog}

\title{Scales wind speeds to certain heights.}

\description{
For a wind speed at given height, returns a scaled wind speed at some different height. Often used to obtain wind speed at hub height.
}

\usage{
WindspeedLog(v0, HH = 100, z0 = 0.1, refHeight = 10)
}

\arguments{
	\item{v0}{
		wind speed (in meters per second) at reference height.
	}
	\item{HH}{
		the height (in meters) at which wind speed is desired.
	}
	\item{z0}{
		roughness length (in meters). Usually set to 0.1 m for onshore sites.
	}
	\item{refHeight}{
		reference height (in meters). The height at which the actual wind speed (v0) was measured.
	}
}

\details{
This function simply implements \deqn{v = v_0 (\frac{\log{\frac{HH}{z_0}}}{\log{\frac{refHeight}{z_0}}})}\cr
Note that this way to scale wind speeds to certain heights is frequently considered deprecated in the literature. Use \code{\link{WindspeedHellmann}} instead.
}

\value{
\code{WindspeedLog} returns a wind speed (in meters per second) at the desired height.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{WindspeedHellmann}} for a different way to scale wind speeds to heights.
}

\examples{
WindspeedLog(v0 = 6, HH = 80, z0 = 0.1, refHeight = 20)
WindspeedHellmann(v0 = 6, HH = 80, refHeight = 20)
}

\keyword{WindspeedLog}
