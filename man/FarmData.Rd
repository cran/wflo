\name{FarmData}
\docType{data}
\alias{FarmData}

\title{Data set for wind farm layout optimization.}

\description{
This list object contains four matrices covering adjusted yield, wind speed, wind direction and standard deviations of wind directions in Germany.
}

\usage{
FarmData
}

\format{
	A \code{list} object containing seven matrices, each containing 25 x 25 values at a raster resolution of 200 x 200 m (note that for the larger data set downloadable using \code{\link{AcquireData}}, each matrix contains 4400 x 3250 values):\cr
	\cr
	\code{$AdjustedYield}: Yield is average annual energy production (AEP). According to FGW technical guidelines, AEP is adjusted due to different location qualities to obtain a better guess at the marketable energy output. Interpret these values directly as 'megawatt hours per year'.\cr
	\cr
	\code{$WindSpeed}: Average wind speeds in meters per second.\cr
	\cr
	\code{$WindDirection}: Average wind directions in degrees (azimuth system).\cr
	\cr
	\code{$SDDirection}: Standard deviations of wind directions in degrees (azimuth system).\cr
	\cr
	\code{$Elevation}: Terrain relief (elevation in meters).\cr
	\cr
	\code{$Slope}: Terrain slope in degrees.\cr
	\cr
	\code{$SlopeDirection}: Direction of hillside in degrees.\cr
}

\source{
DWD Climate Data Center (CDC): \url{ftp://opendata.dwd.de/climate_environment/CDC/grids_germany/multi_annual/wind_parameters/resol_200x200/}\cr
\cr
FGW Technical Guidelines: \url{https://wind-fgw.de/shop/technical-guidelines/?lang=en}\cr
\cr
Croonenbroeck, C. & Hennecke, D. Does the German renewable energy act provide a fair incentive system for onshore wind power? - A simulation analysis, Energy Policy, 2020, 114, 111663, 1-15, https://doi.org/10.1016/j.enpol.2020.111663
}

\author{
Carsten Croonenbroeck\cr
David Hennecke
}

\note{
If the full mode data set is present, it is loaded into the environment \code{e}. The built-in data set \code{FarmData} contains matrices of dimension 25 x 25, while the full data set \code{e$FarmData} consists of 4400 x 3250 matrices.
}

\examples{
# 'Profit' uses this data set internally:
NumTurbines <- 4
set.seed(1357)
Result <- optim(par = runif(NumTurbines * 2), fn = Profit, method = "L-BFGS-B",
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
Result
PlotResult(Result)
}

\keyword{datasets}
