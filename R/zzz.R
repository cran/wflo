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
  c(0.66681197, 0.51920366, 0.40001173, 0.73883743, 0.39199855,
    0.29543266, 0.83333332, 0.36453646, 0.84200820, 0.17836044,
    0.58333332, 0.20833333, 0.77751433, 0.57246404, 0.13369234,
    0.08333335,          0,          1, 0.04166667, 0.12246495,
    0.53759769, 0.03688077, 0.57690332, 0.41578636, 0.75026990,
    0.42526027, 0.63112939, 0.07676504, 0.98465420, 0.02097486,
    0.44808668, 0.64066857, 0.43917715, 0.88938042, 0.37499991,
    0.19688800, 0.40500598, 0.40134872, 0.91059693, 0.08887331)

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
