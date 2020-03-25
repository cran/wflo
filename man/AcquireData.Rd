\name{AcquireData}
\alias{AcquireData}

\title{Downloads the larger data set for entire Germany}

\description{
Convenience function that downloads the larger data set covering the entire area of Germany and saves it to the specified directory. This data set can replace the smaller sub sample data set that is part of this package.
}

\usage{
AcquireData(Folder)
}

\arguments{
	\item{Folder}{
		must be a character string containing the folder location to where the file will be saved. For instance, use \code{\link{getwd}} to save to the current working directory.
	}
}

\details{
This function will interactively lead the user through the downloading process and also gives advice on how to use the larger data set. Make sure that R can write to the specified directory.
}

\value{
\code{AcquireData} returns nothing.
}

\author{
Carsten Croonenbroeck
}

\seealso{
Use \code{\link{FarmVars}} for the settings object and \code{\link{FarmData}} for the smaller, built-in data set that is part of this package.
}

\examples{
\dontrun{
AcquireData(tempdir())
# Will download the data file to the specified directory.
}
}

\keyword{Data set}
\keyword{FarmData}
