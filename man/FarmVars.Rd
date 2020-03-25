\name{FarmVars}
\docType{data}
\alias{FarmVars}

\title{Variables object for wind farm specifications.}

\description{
This list object collects all necessary variables for wind farm layout optimizations within this package. The object is loaded upon loading the package, embedded into a user changeable environment \code{e} and provides a set of default values.
}

\usage{
FarmVars
}

\format{
	A \code{list} object containing:\cr
	\cr
	\code{$UnitCost}: This contains yearly costs per turbine in actual currency. Mostly, this will be installation costs, maintenance cost and others. If, for example, total turbine installation costs are 2 mio. EUR and the projected life span of the turbine is 20 years, then the yearly cost is roughly 100,000 EUR, which is also the default value for this item. Refine the number at will to encompass annuity (taking interest effects into consideration), salvage, or the mentioned maintenance cost.\cr
	\cr
	\code{$Price}: This denotes the average sale price of electricity produced at the wind farm per megawatt hour. Defaults to 100, since 100 EUR per megawatt hour is a good (still rough) estimate in the European electricity market.\cr
	\cr
	\code{$StartPoint}: Denotes the point at which to start to 'cut out' a sample field from the data set. Since the full data set covers the entire region of Germany at a 200 x 200 meters resolution, the entire area is too big for a wind farm. For an actual wind farm layout optimization problem, a smaller region is indicated. This package will conveniently provide a square area of 5 km x 5 km in size. Note: Concerning the larger data set covering the entire area of Germany (conveniently downloadable using \code{\link{AcquireData}}), this area con be found at point 2000:2000 (full data covering a range of 4400 x 3250 points) and from there, spans a five kilometers square, i.e. covers the raster at 2000:2024, 2000:2024. Change the default value of 1 to 2000 to find that area or to something else for cross checking your algorithm with data covering a different area. However, for comparability, evaluate your algorithm result at the default of 2000 or using the built-in data set.\cr
	\cr
	\code{$Width}: Denotes the width (and, since the test area is a square, also the height) of the test area. Defaults to 25, which leads to a 5 km x 5 km square as a sample wind field, as the raster resolution is 200 x 200 m.\cr
	\cr
	\code{$MeterMinDist}: Denotes the minimum distance from one turbine to another in meters. Defaults to 500 meters, which is a rather conservative value.\cr
	\cr
	\code{$EndPoint}: The complement to '$StartPoint'. Is computed based on '$StartPoint' and '$Width' as follows: FarmVars$StartPoint + FarmVars$Width - 1. By default, $StartPoint = 1 and $Width = 25, so $EndPoint = 1 + 25 - 1 = 25.\cr
	\cr
	\code{$MeterWidth}: Contains the width (and, since the test area is a square, also the height) of the test area in meters. Since by default, FarmVars$Width = 25 and the raster resolution is 200 x 200 m, this defaults to $MeterWidth = 200 * 25 = 5000 meters (5 km).\cr
	\cr
	\code{$MinDist}: Contains the minimum distance in problem space. Since the problem space in a convenient unit square, this defaults to $MinDist = MeterMinDist / MeterWidth = 500 / 5000 = 0.1.\cr
	\cr
	\code{$z}: Contains the hub height of the turbine type under investigation. Default value is 100 meters.\cr
	\cr
	\code{$z0}: Contains the terrain's roughness length required for Jensen's wake model. In most studies, a default value of 0.1 is agreed to be a good guess for onshore wind farms, so this is the default value here.\cr
	\cr
	\code{$r0}: Contains the rotor radius for the turbine type under investigation. For instance, a Vestas V90 turbine has a rotor diameter of 90 meters (hence the name), so the radius is 45 mesters, which is also the default value here.
}

\examples{
# 'e$FarmVars' is generated using
e$FarmVars <- list()
e$FarmVars$UnitCost <- 100000
e$FarmVars$Price <- 100
e$FarmVars$StartPoint <- 1
e$FarmVars$Width <- 25
e$FarmVars$MeterMinDist <- 500
e$FarmVars$EndPoint <- e$FarmVars$StartPoint + e$FarmVars$Width - 1
e$FarmVars$MeterWidth <- 200 * e$FarmVars$Width
e$FarmVars$MinDist <- e$FarmVars$MeterMinDist / e$FarmVars$MeterWidth
e$FarmVars$z <- 100
e$FarmVars$z0 <- 0.1
e$FarmVars$r0 <- 45
}

\source{Carsten Croonenbroeck}

\keyword{datasets}
