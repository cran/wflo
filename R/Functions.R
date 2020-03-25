utils::globalVariables(c("FarmData"))

Cost <- function(x, y) #x, y \in R^n
{
	return(rep(e$FarmVars$UnitCost, length(x)))
}

Yield <- function(x, y, Adj) #x, y \in R
{
	return(Adj[as.integer((x * (e$FarmVars$Width - 1)) + 1), as.integer((y * (e$FarmVars$Width - 1)) + 1)])
}

GetDirInfo <- function(x, y, Dirs, SDs) #x, y \in R
{
	ThisDir <- Dirs[as.integer((x * (e$FarmVars$Width - 1)) + 1), as.integer((y * (e$FarmVars$Width - 1)) + 1)]
	ThisSD <- SDs[as.integer((x * (e$FarmVars$Width - 1)) + 1), as.integer((y * (e$FarmVars$Width - 1)) + 1)]

	return(c(ThisDir, ThisSD))
}

ValidSetup <- function(x, y) #x, y \in R^n
{
	n <- length(x)
	M <- matrix(ncol = n, nrow = n)

	for (i in 1:(n - 1))
	{
		for (j in (i + 1):n)
		{
			M[i, j] <- (x[i] - x[j]) ^ 2 + (y[i] - y[j]) ^ 2
		}
	}
	M[is.na(M)] <- 2 * (e$FarmVars$MinDist * e$FarmVars$MinDist)

	if (any(M < (e$FarmVars$MinDist * e$FarmVars$MinDist))) return(FALSE) else return(TRUE)
}

Geo2Ari <- function(g)
{
	if (g == -1) return(-1)
	return((450 - g) %% 360)
}

IsInWake <- function(Cone, Direction)
{
	if (Cone[1] < 0) Cone[1] <- 360 + Cone[1]
	Cone[2] <- Cone[2] %% 360

	u <- as.numeric()
	v <- as.numeric()

	u[1] <- sin(Cone[1] * pi / 180)
	u[2] <- sin(Cone[2] * pi / 180)

	v[1] <- cos(Cone[1] * pi / 180)
	v[2] <- cos(Cone[2] * pi / 180)

	u <- sort(u)
	v <- sort(v)

	Dir_u <- sin(Direction * pi / 180)
	Dir_v <- cos(Direction * pi / 180)

	if (Dir_u >= u[1] & Dir_u <= u[2]) return(TRUE)
	if (Dir_v >= v[1] & Dir_v <= v[2]) return(TRUE)

	return(FALSE)
}

JensenAngle <- function(x)
{
	alpha <- 0.5 / (log(e$FarmVars$z / e$FarmVars$z0))
	alpha_x <- alpha * x
	r <- e$FarmVars$r0 + alpha_x
	Norm <- sqrt((r ^ 2) + (x ^ 2))
	Cos_gamma <- x / Norm

	return(acos(Cos_gamma) / pi * 180)
}

JensenFactor <- function(x)
{
	alpha <- 0.5 / (log(e$FarmVars$z / e$FarmVars$z0))
	alpha_x <- alpha * x
	r <- e$FarmVars$r0 + alpha_x

	return((1 - ((2 / 3) * (e$FarmVars$r0 / r) ^ 2)))
}

GetAngle <- function(x1, y1, x2, y2) #As seen from point 2, where is point 1? Give target point first, location second.
{
	A <- c(x1, y1)
	B <- c(x2, y2)
	if (isTRUE(all.equal(A, B))) return(-1)

	APrime <- A - B
	lAPrime <- sqrt(APrime[1] ^ 2 + APrime[2] ^ 2)

	cosAlpha <- abs(APrime[1]) / (lAPrime)
	Degrees <- acos(cosAlpha) / pi * 180

	if (APrime[1] >= 0)
	{
		if (APrime[2] >= 0)
		{
			#Nothing to do here
		} else
		{
			Degrees <- 360 - Degrees
		}
	} else
	{
		if (APrime[2] >= 0)
		{
			Degrees <- 180 - Degrees
		} else
		{
			Degrees <- 180 + Degrees
		}
	}

	if (is.nan(Degrees) | is.na(Degrees)) Degrees <- 0

	return(Degrees)
}

PairPenalty <- function(x1, y1, x2, y2, Dirs, SDs) #As seen from point 2.
{
	PointsAngle <- Geo2Ari(GetAngle(x1, y1, x2, y2))
	if (PointsAngle == -1) return(0)

	PointInfo <- GetDirInfo(x2, y2, Dirs, SDs)

	Distance <- sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2) * e$FarmVars$MeterWidth
	JensenCone <- JensenAngle(Distance)
	JensenCone <- c(PointsAngle - JensenCone, PointsAngle + JensenCone)

	if (IsInWake(JensenCone, PointInfo[1]))
	{
		return(JensenFactor(Distance))
	} else
	{
		return(1)
	}
}

Profit <- function(X)
{
	xSel <- seq(from = 1, to = length(X) - 1, by = 2)
	x <- X[xSel]
	y <- X[xSel + 1]

	if (exists("FarmData", envir = e, inherits = FALSE))
	{
	  Adj <- e$FarmData[[1]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
	  Dirs <- e$FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
	  SDs <- e$FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
	} else
	{
	  Adj <- FarmData[[1]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
	  Dirs <- FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
	  SDs <- FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
	}

	n <- length(x)

	if (ValidSetup(x, y) == FALSE) return(sum(rep(e$FarmVars$UnitCost, n))) #Make it a minimization problem.

	M <- rep(1, n)
	M <- M %*% t(M)
	for (i in 1:n)
	{
		for (j in 1:n)
		{
			M[i, j] <- M[i, j] * PairPenalty(x[j], y[j], x[i], y[i], Dirs, SDs)
		}
	}
	M[M == 0] <- 1
	p <- rep(1, n)
	for (i in 1:n)
	{
		for (j in 1:n)
		{
			p[i] <- p[i] * M[j, i]
		}
	}

	PointYield <- as.numeric()
	for (i in 1:n) PointYield[i] <- Yield(x[i], y[i], Adj)

	Revenue <- (PointYield * e$FarmVars$Price) * p

	retVal <- sum(Revenue - Cost(x, y))

	return(-retVal) #Make it a minimization problem.
}

GetArrow <- function(BaseX, BaseY, Degrees, Frac = 25)
{
	XPrime <- BaseX + cos(Degrees * pi / 180)
	YPrime <- BaseY + sin(Degrees * pi / 180)

	dxfrac <- (XPrime - BaseX) / Frac
	dyfrac <- (YPrime - BaseY) / Frac

	x0 <- BaseX - dxfrac
	y0 <- BaseY - dyfrac
	x1 <- BaseX + dxfrac
	y1 <- BaseY + dyfrac

	return(c(x0, y0, x1, y1))
}

ImposeVectorField <- function(xNum = 5, yNum = 5, ColMain = "black", ColBand = "white", Frac = 25, DoSDs = TRUE)
{
  if (exists("FarmData", envir = e, inherits = FALSE))
  {
  	Dirs <- e$FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
	  SDs <- e$FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  } else
  {
    Dirs <- FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    SDs <- FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  }

	Counter <- 1

	MyArrows <- data.frame(x0 = rep(NA, xNum * yNum), y0 = NA, x1 = NA, y1 = NA)
	MyArrows_A <- MyArrows
	MyArrows_B <- MyArrows

	for (i in seq(from = 0, to = 1, length.out = xNum))
	{
		for (j in seq(from = 0, to = 1, length.out = yNum))
		{
			ThisArrowInfos <- GetDirInfo(i, j, Dirs, SDs)

			ThisArrow <- GetArrow(i, j, ThisArrowInfos[1], Frac)
			MyArrows$x0[Counter] <- ThisArrow[1]
			MyArrows$y0[Counter] <- ThisArrow[2]
			MyArrows$x1[Counter] <- ThisArrow[3]
			MyArrows$y1[Counter] <- ThisArrow[4]

			ThisArrow <- GetArrow(i, j, ThisArrowInfos[1] + (0.5 * ThisArrowInfos[2]), Frac)
			MyArrows_A$x0[Counter] <- ThisArrow[1]
			MyArrows_A$y0[Counter] <- ThisArrow[2]
			MyArrows_A$x1[Counter] <- ThisArrow[3]
			MyArrows_A$y1[Counter] <- ThisArrow[4]

			ThisArrow <- GetArrow(i, j, ThisArrowInfos[1] - (0.5 * ThisArrowInfos[2]), Frac)
			MyArrows_B$x0[Counter] <- ThisArrow[1]
			MyArrows_B$y0[Counter] <- ThisArrow[2]
			MyArrows_B$x1[Counter] <- ThisArrow[3]
			MyArrows_B$y1[Counter] <- ThisArrow[4]

			Counter <- Counter + 1
		}
	}

	if (DoSDs)
	{
		arrows(x0 = MyArrows_A$x0, y0 = MyArrows_A$y0, x1 = MyArrows_A$x1, y1 = MyArrows_A$y1, col = ColBand, length = 0.1)
		arrows(x0 = MyArrows_B$x0, y0 = MyArrows_B$y0, x1 = MyArrows_B$x1, y1 = MyArrows_B$y1, col = ColBand, length = 0.1)
	}
	arrows(x0 = MyArrows$x0, y0 = MyArrows$y0, x1 = MyArrows$x1, y1 = MyArrows$y1, col = ColMain, length = 0.1)
}

PlotResult <- function(Result)
{
  if (exists("FarmData", envir = e, inherits = FALSE))
  {
	  Adj <- e$FarmData[[1]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  } else
  {
    Adj <- FarmData[[1]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  }

	Z <- matrix(ncol = e$FarmVars$Width, nrow = e$FarmVars$Width)
	for (x in 0:e$FarmVars$Width)
	{
		for (y in 0:e$FarmVars$Width)
		{
			Z[x, y] <- Adj[x, y]
		}
	}

	xSel <- seq(from = 1, to = length(Result$par) - 1, by = 2)
	xOpt <- Result$par[xSel]
	yOpt <- Result$par[xSel + 1]

	Oldpar <- par(no.readonly = TRUE)
	on.exit(par(Oldpar))
	par(mar = c(0, 0, 0, 0), mai = c(0, 0, 0, 0))
	image(Adj, xaxt = "n", yaxt = "n", bty = "n", col = rgb(0, 0, seq(from = 0, to = 1, length.out = e$FarmVars$Width * e$FarmVars$Width)))
	contour(x = seq(from = 0, to = 1, length.out = e$FarmVars$Width), y = seq(from = 0, to = 1, length.out = e$FarmVars$Width), z = Z, add = TRUE)
	ImposeVectorField(xNum = 7, yNum = 7, ColMain = "gray", ColBand = "gray30", DoSDs = TRUE)
	xPlot <- xOpt
	yPlot <- yOpt

	a <- as.numeric()
	for (i in 0:e$FarmVars$Width) a[i] <- (i / (e$FarmVars$Width - 1)) - (1 / (e$FarmVars$Width - 1))

	for (i in 1:length(xOpt))
	{
		b <- xOpt[i]
		d <- a - b
		d[d > 0] <- NA
		d <- abs(d)
		e <- a[which(d == min(na.omit(d)))]
		xPlot[i] <- e

		b <- yOpt[i]
		d <- a - b
		d[d > 0] <- NA
		d <- abs(d)
		e <- a[which(d == min(na.omit(d)))]
		yPlot[i] <- e
	}

	points(xPlot, yPlot, col = "gold4", pch = 16)
}

AcquireData = function(Folder)
{
	if (!interactive())
	{
		message("Requires interactive mode.")
		return(invisible(NULL))
	}

	URL = "http://wflo.auf.uni-rostock.de/FarmData.RData"

	message(paste("Downloading data file from ", URL, sep = ""))
	message(paste("Writing to directory ", Folder, sep = ""))
	message("Proceed? Hit 'n' to cancel or any other key to continue.")
	retVal = readline()
	if (substr(retVal, 1, 1) == "n")
	{
		message("Cancelled.")
		return(invisible(NULL))
	}

	FileLocation = paste(Folder, "/FarmData.RData", sep = "")

	if (file.exists(FileLocation))
	{
		message("File exists. Really proceed? Hit 'n' to cancel or any other key to continue.")
		retVal = readline()
		if (substr(retVal, 1, 1) == "n")
		{
			message("Cancelled.")
			return(invisible(NULL))
		}

		message("Trying to delete the file...")
		file.remove(FileLocation)
		if (file.exists(FileLocation))
		{
			message("Could not delete file. Write protected? Aborting...")
			return(invisible(NULL))
		} else
		{
			message("Successfully deleted.")
		}
	}

	message("Starting download (118 MB), please wait...")
	flush.console()
	Ret = try(suppressWarnings(utils::download.file(URL, FileLocation, quiet = TRUE)), silent = TRUE)

	if (class(Ret) == "try-error")
	{
		message("File was not downloaded. Internet connection error? Please try again.")
		return(invisible(NULL))
	}

	if (Ret == 0 & file.exists(FileLocation))
	{
		if (file.info("FarmData.RData")$size != 124606617)
		{
			message("File seems to be corrupt. Download interrupted? Please try again.")
			return(invisible(NULL))
		} else
		{
			message("File downloaded successfully.\n")
			message("Type \"?FarmData\" for further information on the data.")
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
		}
	} else
	{
		message("File was not downloaded. Internet connection error? Please try again.")
		return(invisible(NULL))
	}
}
