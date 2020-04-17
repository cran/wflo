e <- new.env()

FarmVars <- list()
FarmVars$UnitCost <- 100000
FarmVars$Price <- 100
FarmVars$StartPoint <- 1
FarmVars$Width <- 25
FarmVars$MeterMinDist <- 500
FarmVars$EndPoint <- e$FarmVars$StartPoint + e$FarmVars$Width - 1
FarmVars$MeterWidth <- 200 * e$FarmVars$Width
FarmVars$MinDist <- e$FarmVars$MeterMinDist / e$FarmVars$MeterWidth
FarmVars$z <- 100
FarmVars$z0 <- 0.1
FarmVars$r0 <- 45
FarmVars$Bench <-
  c(0.6700295,  0.5224326, 0.4578065,   0.636406,     0.4339372,
    0.7446001,  0.7147827, 0.6249925,   0.04166674,   0.0854744,
    0.9016667,  0.105121,  0.0000000,   1.0000000,    0.8315406,
    0.5410363,  0.7901967, 0.4433062,   0.1416438,    0.10833333,
    0.8420252,  0.1918029, 0.5697075,   9.523577e-16, 0.552057,
    0.3994493,  0.4338577, 0.9155776,   0.6356543,    0.08333333,
    0.7954187,  0.338033,  0.8916667,   0.3003984,    0.3583333,
    0.2083333,  0.3939717, 0.3817109,   0.9586827,    1.080283e-05)

.onAttach <- function(libname, pkgname)
{
  Ver <- as.character(packageVersion("wflo")[1])
  packageStartupMessage(paste("This is wflo, version ", Ver, ". Looking for full data file...", sep = ""))

	if (file.exists("FarmData.RData"))
	{
		if (file.info("FarmData.RData")$size == 233636819)
		{
			packageStartupMessage("File found. Loading...")
			flush.console()

			load("FarmData.RData")
			e$FarmData <- FarmData

			e$FarmVars <- list()
			e$FarmVars$UnitCost <- 100000
			e$FarmVars$Price <- 100
			e$FarmVars$StartPoint <- 2000
			e$FarmVars$Width <- 25
			e$FarmVars$MeterMinDist <- 500
			e$FarmVars$EndPoint <- e$FarmVars$StartPoint + e$FarmVars$Width - 1
			e$FarmVars$MeterWidth <- 200 * e$FarmVars$Width
			e$FarmVars$MinDist <- e$FarmVars$MeterMinDist / e$FarmVars$MeterWidth
			e$FarmVars$z <- 100
			e$FarmVars$z0 <- 0.1
			e$FarmVars$r0 <- 45
			e$FarmVars$BenchmarkSolution <- FarmVars$Bench

			packageStartupMessage("Success. Operating in full mode now.")
		} else
		{
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
		  e$FarmVars$BenchmarkSolution <- FarmVars$Bench

			packageStartupMessage("File found, but seems to be corrupt. Operating in basic mode.")
		}
	} else
	{
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
	  e$FarmVars$BenchmarkSolution <- FarmVars$Bench

		packageStartupMessage("File not found, operating in basic mode.")
		packageStartupMessage("Type \"AcquireData(getwd())\" to get the full data file.")
		packageStartupMessage("If you already have the file, set the R working directory to the directory containing the file before loading wflo so that wflo can find the file.")
		packageStartupMessage("Alternatively, load the file manually using load(). Note that after that, you should make the necessary adjustments to find the benchmark area.")
		packageStartupMessage("Type 'vignette(\"wflo\")' for additional information.")
	}
}
