utils::globalVariables(c("FarmData"))

Cost <- function(x, y) #x, y \in R^n
{
	return(rep(e$FarmVars$UnitCost, min(length(x), length(y))))
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

JensenTrapezoid <- function(WindDir, Point, x)
{
  x = x * 1.1 #Add a little margin.

  alpha <- 0.5 / (log(e$FarmVars$z / e$FarmVars$z0))
  alpha_x <- alpha * x
  r <- e$FarmVars$r0 + alpha_x

  r0U <- (1 / e$FarmVars$MeterWidth) * e$FarmVars$r0
  rU <- (1 / e$FarmVars$MeterWidth) * r
  xU <- (1 / e$FarmVars$MeterWidth) * x

  AAngle <- WindDir + 90
  AAngle <- AAngle %% 360

  BAngle <- WindDir + 270
  BAngle <- BAngle %% 360

  EAngle <- WindDir + 180
  EAngle <- EAngle %% 360

  ##

  Displace <- c(r0U, 0)

  Mat <- cbind(c(cos(AAngle * pi / 180), sin(AAngle * pi / 180)), c(-sin(AAngle * pi / 180), cos(AAngle * pi / 180)))
  A <- (Mat %*% Displace) + Point

  Mat <- cbind(c(cos(BAngle * pi / 180), sin(BAngle * pi / 180)), c(-sin(BAngle * pi / 180), cos(BAngle * pi / 180)))
  B <- (Mat %*% Displace) + Point

  ##

  Displace <- c(xU, 0)

  Mat <- cbind(c(cos(EAngle * pi / 180), sin(EAngle * pi / 180)), c(-sin(EAngle * pi / 180), cos(EAngle * pi / 180)))
  E <- (Mat %*% Displace) + Point

  ##

  Displace <- c(rU, 0)

  Mat <- cbind(c(cos(AAngle * pi / 180), sin(AAngle * pi / 180)), c(-sin(AAngle * pi / 180), cos(AAngle * pi / 180)))
  D <- (Mat %*% Displace) + E

  Mat <- cbind(c(cos(BAngle * pi / 180), sin(BAngle * pi / 180)), c(-sin(BAngle * pi / 180), cos(BAngle * pi / 180)))
  C <- (Mat %*% Displace) + E

  ##

  RetVal <- cbind(A, B, C, D)
  colnames(RetVal) <- c("A", "B", "C", "D")

  return(RetVal)
}

PointInPolygon <- function(PointMat, TestPoint) # Jordan test.
{
  CrossProdTest <- function(A, B, C)
  {
    if (B[2] > C[2])
    {
      Temp <- B
      B <- C
      C <- Temp
    }

    if (A[2] <= B[2] | A[2] > C[2]) return(1)

    Delta <- ((B[1] - A[1]) * (C[2] - A[2])) - ((B[2] - A[2]) * (C[1] - A[1]))

    if (Delta >= 0) return(1)

    return(-1)
  }

  P <- cbind(PointMat[, ncol(PointMat)], PointMat)

  t <- -1
  for (i in 1:(ncol(P) - 1))
  {
    t <- t * CrossProdTest(TestPoint, P[, i], P[, i + 1])
    if (t == 0) break
  }
  if (t < 0) return(FALSE) else return(TRUE)
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
			# Nothing to do here.
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
	PointInfo[1] <- Geo2Ari(PointInfo[1])
	PointInfo[1] <- PointInfo[1] + 180
	PointInfo[1] <- PointInfo[1] %% 360

	Distance <- sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2) * e$FarmVars$MeterWidth
	JensenCone <- JensenTrapezoid(PointInfo[1], c(x2, y2), Distance)

	if (PointInPolygon(JensenCone, c(x1, y1)))
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
			ThisArrowInfos[1] <- Geo2Ari(ThisArrowInfos[1])

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

ProfitContributors <- function(Result)
{
  if (inherits(Result, "list")) if (any(names(Result) == "par")) Result <- Result$par
  if (!inherits(Result, "numeric")) stop("Object is neither a vector nor optimizer result.")

  xSel <- seq(from = 1, to = length(Result) - 1, by = 2)
  x <- Result[xSel]
  y <- Result[xSel + 1]

  if (any(x < 0) | any(x > 1) | any(y < 0) | any(y > 1)) stop("At leat one point is out of bounds.")

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
  retVal <- Revenue - Cost(x, y)
  retVal <- cbind(1:length(x), retVal)
  colnames(retVal) <- c("Turbine", "Profit")
  return(retVal)
}

PlotResult <- function(Result, DoLabels = FALSE, Labels = "IDs")
{
  if (!any(names(Result) == "par")) stop("Not a valid optimizer result provided.")

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
	if (DoLabels)
	{
	  if (Labels[1] == "IDs") Labels <- 1:length(xPlot)
	  if (length(Labels) != length(xPlot))
	  {
	    warning("Length of labels object must match the number of points.")
	    return(invisible(NULL))
	  }
	  text(xPlot, yPlot, labels = Labels, pos = 3, col = "white")
	}
}

ShowWakePenalizers <- function(Result, Cones = TRUE, VectorField = TRUE)
{
  if (!any(names(Result) == "par")) stop("Not a valid optimizer result provided.")

  if (exists("FarmData", envir = e, inherits = FALSE))
  {
    Dirs <- e$FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    SDs <- e$FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  } else
  {
    Dirs <- FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    SDs <- FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  }

  xSel <- seq(from = 1, to = length(Result$par) - 1, by = 2)
  x <- Result$par[xSel]
  y <- Result$par[xSel + 1]

  n <- length(x)
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

  Causers <- 0
  Sufferers <- 0
  for (i in 1:n)
  {
    if (mean(M[, i]) != 1) Sufferers <- c(Sufferers, i)
    if (mean(M[i, ]) != 1) Causers <- c(Causers, i)
  }
  Causers <- Causers[2:length(Causers)]
  Sufferers <- Sufferers[2:length(Sufferers)]

  if (any(is.na(Causers)) | any(is.na(Sufferers)))
  {
    message("Nothing to show, no wake effects present.")
    return(invisible(NULL))
  }

  plot(x[c(Causers, Sufferers)], y[c(Causers, Sufferers)], ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "")
  for (i in 1:length(Sufferers))
  {
    points(x[Sufferers[i]], y[Sufferers[i]], col = "gold4", pch = 16)
    text(x[Sufferers[i]], y[Sufferers[i]], labels = Sufferers[i], pos = 3)
  }
  for (i in 1:length(Causers))
  {
    points(x[Causers[i]], y[Causers[i]], col = "grey", pch = 16)
    text(x[Causers[i]], y[Causers[i]], labels = Causers[i], pos = 3)

    if (Cones)
    {
      WindDir <- GetDirInfo(x[Causers[i]], y[Causers[i]], Dirs, SDs)
      WindDir <- Geo2Ari(WindDir[1])
      WindDir <- WindDir + 180
      WindDir <- WindDir %% 360
      JensenCone <- JensenTrapezoid(WindDir, c(x[Causers[i]], y[Causers[i]]), e$FarmVars$MeterWidth * 2)

      lines(x = c(JensenCone[1, 1], JensenCone[1, 2]), y = c(JensenCone[2, 1], JensenCone[2, 2]), col = "red")
      lines(x = c(JensenCone[1, 2], JensenCone[1, 3]), y = c(JensenCone[2, 2], JensenCone[2, 3]), col = "red")
      lines(x = c(JensenCone[1, 3], JensenCone[1, 4]), y = c(JensenCone[2, 3], JensenCone[2, 4]), col = "red")
      lines(x = c(JensenCone[1, 4], JensenCone[1, 1]), y = c(JensenCone[2, 4], JensenCone[2, 1]), col = "red")
    }
  }
  if (VectorField) ImposeVectorField(xNum = 7, yNum = 7, ColMain = "gray", ColBand = "gray30", DoSDs = TRUE)
  legend("topleft", c("Causers", "Sufferers"), bg = "white", pch = 16, col = c("grey", "gold4"))

  return(invisible(M))
}

AcquireData <- function(Folder)
{
	if (!interactive())
	{
		message("Requires interactive mode.")
		return(invisible(NULL))
	}

	URL <- "http://wflo.auf.uni-rostock.de/FarmData.RData"

	message(paste("Downloading data file from ", URL, sep = ""))
	message(paste("Writing to directory ", Folder, sep = ""))
	message("Proceed? Hit 'n' to cancel or any other key to continue.")
	retVal <- readline()
	if (substr(retVal, 1, 1) == "n")
	{
		message("Cancelled.")
		return(invisible(NULL))
	}

	FileLocation <- paste(Folder, "/FarmData.RData", sep = "")

	if (file.exists(FileLocation))
	{
		message("File exists. Really proceed? Hit 'n' to cancel or any other key to continue.")
		retVal <- readline()
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

	message("Starting download (222 MB), please wait...")
	flush.console()
	Ret <- try(suppressWarnings(utils::download.file(URL, FileLocation, quiet = TRUE)), silent = TRUE)

	if (inherits(Ret, "try-error"))
	{
		message("File was not downloaded. Internet connection error? Please try again.")
		return(invisible(NULL))
	}

	if (Ret == 0 & file.exists(FileLocation))
	{
		if (file.info("FarmData.RData")$size != 233636819)
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
