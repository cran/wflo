e <- new.env()

FarmVars <- list()
FarmVars$UnitCost <- 100000
FarmVars$Price <- 100
FarmVars$StartPoint <- 1
FarmVars$Width <- 25
FarmVars$MeterMinDist <- 500
FarmVars$EndPoint <- FarmVars$StartPoint + FarmVars$Width - 1
FarmVars$MeterWidth <- 200 * FarmVars$Width
FarmVars$MinDist <- FarmVars$MeterMinDist / FarmVars$MeterWidth
FarmVars$z <- 100
FarmVars$z0 <- 0.1
FarmVars$r0 <- 45

.onAttach <- function(libname, pkgname)
{
  packageStartupMessage("This is wflo. Looking for full data file...")

	if (file.exists("FarmData.RData"))
	{
		if (file.info("FarmData.RData")$size == 124606617)
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

			packageStartupMessage("Success. Operating in full mode now.")
		} else
		{
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

		packageStartupMessage("File not found, operating in basic mode.")
		packageStartupMessage("Type \"AcquireData(getwd())\" to get the full data file.")
		packageStartupMessage("If you already have the file, set the R working directory to the directory containing the file before loading wflo so that wflo can find the file.")
		packageStartupMessage("Alternatively, load the file manually using load(). Note that after that, you should make the following adjustments to find the benchmark area:")
		packageStartupMessage("\"e$FarmVars$StartPoint <- 2000\" and \"e$FarmVars$EndPoint <- 2024\"")
	}
}
