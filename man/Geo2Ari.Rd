\name{Geo2Ari}
\alias{Geo2Ari}

\title{Converts degrees between the arithmetic system and the azimuth system (and vice versa).}

\description{
Use this function to convert degrees from the arithmetic system (0° being east, ascending counterclockwise) into the azimuth system (0° being north, ascending clockwise) and vice versa.
}

\usage{
Geo2Ari(g)
}

\arguments{
	\item{g}{
		contains a single value degree (usually between 0° and 360°, decimal fractions allowed).
	}
}

\details{
\code{g} may contain degrees from both systems, the function turns the data into the respective other system.
}

\value{
\code{Geo2Ari} returns a single value of degrees.
}

\author{
Carsten Croonenbroeck
}

\seealso{
\code{\link{GetDirInfo}} for further degree information.
}

\examples{
Geo2Ari(0)
## In an arithmetic system, 0° means 'east', while 'east' in
## azimuth notation is 90°. This call returns 90.
Geo2Ari(90)
## In an azimuth system, 90° means 'east', while 'east' in
## arithmetic notation is 0°. This call returns 0.
Geo2Ari(Geo2Ari(123))
## Returns 123.
}

\keyword{GetDirInfo}
\keyword{Profit}
