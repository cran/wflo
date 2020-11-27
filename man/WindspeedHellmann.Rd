\name{WindspeedHellmann}
\alias{WindspeedHellmann}

\title{Scales wind speeds to certain heights.}

\description{
For a wind speed at given height, returns a scaled wind speed at some different height. Often used to obtain wind speed at hub height.
}

\usage{
WindspeedHellmann(v0, HH = 100, refHeight = 10, alpha = 1 / 7)
}

\arguments{
	\item{v0}{
		wind speed (in meters per second) at reference height.
	}
	\item{HH}{
		the height (in meters) at which wind speed is desired.
	}
	\item{refHeight}{
		reference height (in meters). The height at which the actual wind speed (v0) was measured.
	}
	\item{alpha}{
		power law parameter. Usually set to 1/7 for onshore sites.
	}
}

\details{
This function simply implements \deqn{v = v_0 (\frac{HH}{refHeight})^{\alpha}}
}

\value{
\code{WindspeedHellmann} returns a wind speed (in meters per second) at the desired height.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{WindspeedLog}} for a different way to scale wind speeds to heights.
}

\examples{
WindspeedLog(v0 = 6, HH = 80, z0 = 0.1, refHeight = 20)
WindspeedHellmann(v0 = 6, HH = 80, refHeight = 20)
}

\keyword{WindspeedHellmann}
