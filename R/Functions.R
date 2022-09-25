utils::globalVariables(c("FarmData", "ParallelProfile"))

#######################################################################################################################################

MosettiTurbineCost <- function(N) N * ((2 / 3) + ((1 / 3) * exp(-0.00174 * (N ^ 2))))

#######################################################################################################################################

WindspeedLog <- function(v0, HH = 100, z0 = 0.1, refHeight = 10) v0 * (log(HH / z0)) / (log(refHeight / z0))

#######################################################################################################################################

WindspeedHellmann <- function(v0, HH = 100, refHeight = 10, alpha = 1 / 7) v0 * (HH / refHeight) ^ (alpha)

#######################################################################################################################################

GK2LonLat <- function(X, Y)
{
  GK <- data.frame(cbind("X_GK" = X, "Y_GK" = Y))

  sp::coordinates(GK) <- c("X_GK", "Y_GK")
  sp::proj4string(GK) <- sp::CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +datum=potsdam +units=m +no_defs")
  LonLat <- sp::spTransform(GK, sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  LonLat <- drop(LonLat@coords)
  attr(LonLat, "names") <- NULL

  return(LonLat)
}

#######################################################################################################################################

LonLat2GK <- function(Lon, Lat)
{
  LonLat <- data.frame(cbind("Lon" = Lon, "Lat" = Lat))

  sp::coordinates(LonLat) <- c("Lon", "Lat")
  sp::proj4string(LonLat) <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  GK <- sp::spTransform(LonLat, sp::CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +datum=potsdam +units=m +no_defs"))
  GK <- drop(GK@coords)
  attr(GK, "names") <- NULL

  return(GK)
}

#######################################################################################################################################

Index2GK <- function(X, Y)
{
  Left <- 3280000
  Top <- 6109800

  X <- X - 1
  X <- X * 200
  X <- X + Left
  X <- round(X)

  Y <- Y - 1
  Y <- Y * 200
  Y <- Top - Y
  Y <- round(Y)

  return(c(X, Y))
}

#######################################################################################################################################

GK2Index <- function(X, Y)
{
  Left <- 3280000
  Top <- 6109800

  X <- X - Left
  X <- X / 200
  X <- X + 1
  X <- round(X)

  Y <- Top - Y
  Y <- Y / 200
  Y <- Y + 1
  Y <- round(Y)

  return(c(X, Y))
}

#######################################################################################################################################

GetFarmFromLonLat <- function(Top, Left, EdgeLen, DoPlot = FALSE)
{
#  Invert <- function(x)
#  {
#    y <- x
#
#    for (i in 1:nrow(x))
#    {
#      y[nrow(x) - i + 1, ] <- x[i, ]
#    }
#
#    return(y)
#  }
  Invert <- function(x) x[seq(from = nrow(x), to = 1, by = -1), ]

  Tiles <- EdgeLen / 200

  GK1 <- LonLat2GK(Left, Top)
  GK2 <- GK1
  GK2[1] <- GK2[1] + ((Tiles - 1) * 200)
  GK2[2] <- GK2[2] - ((Tiles - 1) * 200)
  I1 <- GK2Index(GK1[1], GK1[2])
  I2 <- GK2Index(GK2[1], GK2[2])

  if (!exists("FarmData", envir = e, inherits = FALSE)) stop("No full data set found.")
  if (I1[1] < 0 | I1[2] < 0 | I2[1] > ncol(e$FarmData$AdjustedYield) | I2[2] > nrow(e$FarmData$AdjustedYield)) stop("Selected farm out of boundaries.")

  Temp <- list()
  for (i in 1:length(e$FarmData))
  {
    Temp[[i]] <- e$FarmData[[i]][I1[1]:I2[1], I1[2]:I2[2]]
    Temp[[i]] <- t(Temp[[i]])
    Temp[[i]] <- Invert(Temp[[i]])
  }
  names(Temp) <- names(e$FarmData)

  if (DoPlot)
  {
    Test <- raster::raster(Temp$AdjustedYield)
    Proj <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"
    raster::crs(Test) <- Proj
    plot(Test)
  }

  return(Temp)
}

#######################################################################################################################################

AirDensity <- function(Altitude = 0, Temperature = 20, phi = 0.76)
{
  AirPressure <- function(H)
  {
    p0 <- 101325
    H0 <- 8435

    return(p0 * exp(-H / H0))
  }

  Humidity <- function(phi, p, Temperature)
  {
    Magnus <- function(Temperature)
    {
      return(6.112 * exp((17.62 * Temperature) / (243.12 + Temperature)))
    }

    RS <- 287.058
    Rd <- 461.523
    pd <- Magnus(Temperature)

    Rf <- RS / (1 - (phi * (pd / p) * (1 - (RS / Rd))))

    return(Rf)
  }

  if (Temperature < -45 | Temperature > 60) stop("Temperature outside valid range (-45 to +60 degrees Celsius).")
  if (phi < 0 | phi > 1) stop("phi outside valid range (0 to 1).")

  T <- 273.15 + Temperature
  p <- AirPressure(Altitude)
  R <- Humidity(phi, p, Temperature)

  return(p / (R * T))
}

#######################################################################################################################################

Cost <- function(x, y) #x, y \in R^n
{
	return(rep(e$FarmVars$UnitCost, min(length(x), length(y))))
}

#######################################################################################################################################

Yield <- function(x, y, Adj) #x, y \in R
{
	return(Adj[as.integer((x * (e$FarmVars$Width - 1)) + 1), as.integer((y * (e$FarmVars$Width - 1)) + 1)])
}

#######################################################################################################################################

Height <- function(x, y, Elev)
{
  return(Elev[as.integer((x * (e$FarmVars$Width - 1)) + 1), as.integer((y * (e$FarmVars$Width - 1)) + 1)])
}

#######################################################################################################################################

Area <- function(a, b, d) #http://walter.bislins.ch/blog/index.asp?page=Schnittfl%E4che+zweier+Kreise+berechnen
{
  if (d >= a + b) return(0)
  if (d <= abs(a - b))
  {
    if (a > b) r = b else r = a
    return(pi * (r ^ 2))
  }

  x <- (a * a - b * b + d * d) / (2 * d)
  if (suppressWarnings(is.nan(acos(x / a)))) return(0)
  if (suppressWarnings(is.nan(sqrt(a * a - x * x)))) return(0)

  As <- a * a * acos(x / a)
  Ad <- x * sqrt(a * a - x * x)
  Aa <- As - Ad

  y <- (b * b - a * a + d * d) / (2 * d)
  if (suppressWarnings(is.nan(acos(y / b)))) return(0)
  if (suppressWarnings(is.nan(sqrt(b * b - y * y)))) return(0)

  Bs <- b * b * acos(y / b)
  Bd <- y * sqrt(b * b - y * y)
  Ba <- Bs - Bd

  return (Aa + Ba)
}

#######################################################################################################################################

GetDirInfo <- function(x, y, Dirs, SDs) #x, y \in R
{
	ThisDir <- Dirs[as.integer((x * (e$FarmVars$Width - 1)) + 1), as.integer((y * (e$FarmVars$Width - 1)) + 1)]
	ThisSD <- SDs[as.integer((x * (e$FarmVars$Width - 1)) + 1), as.integer((y * (e$FarmVars$Width - 1)) + 1)]

	return(c(ThisDir, ThisSD))
}

#######################################################################################################################################

ValidSetup <- function(x, y) #x, y \in R^n
{
	n <- length(x)
	if (n < 2) return(TRUE)

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

#######################################################################################################################################

Geo2Ari <- function(g)
{
	if (g == -1) return(-1)

	return((450 - g) %% 360)
}

#######################################################################################################################################

JensenTrapezoid <- function(WindDir, Point, x, Margin = TRUE)
{
  if(Margin) x = x * 1.1 #Add a little margin.

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

#######################################################################################################################################

PartialJensen <- function(x2, y2, x1, y1, Dirs, SDs, Elev, Z = NULL, DrawTop = FALSE, DrawFront = FALSE, DrawSide = FALSE)
{
  Intersect <- function(ConeEdge1, ConeEdge2, Disc1, Disc2)
  {
    R <- ConeEdge2 - ConeEdge1
    S <- Disc2 - Disc1

    r1 <- R[1] / S[1]
    r2 <- R[2] / S[2]

    if (r1 == r2)
    {
      warning("Lines are collinear.")
      return(0)
    }

    M <- cbind(R, -S)
    b <- Disc1 - ConeEdge1
    p <- solve(M) %*% b

    return(ConeEdge1 + p[1] * R)
  }

  Norm <- function(x)
  {
    return(sqrt(sum(x ^ 2)))
  }

  if (length(Z) == 2)
  {
    z1 <- Z[2]
    z2 <- Z[1]
  } else
  {
    z1 <- e$FarmVars$z
    z2 <- e$FarmVars$z
  }

  PointInfo <- GetDirInfo(x1, y1, Dirs, SDs)
  PointInfo[1] <- Geo2Ari(PointInfo[1])
  PointInfo[1] <- PointInfo[1] + 180
  PointInfo[1] <- PointInfo[1] %% 360
  JensenCone1 <- JensenTrapezoid(PointInfo[1], c(x1, y1), max(sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2) * e$FarmVars$MeterWidth, e$FarmVars$MeterMinDist))

  PointInfo <- GetDirInfo(x2, y2, Dirs, SDs)
  PointInfo[1] <- Geo2Ari(PointInfo[1])
  PointInfo[1] <- PointInfo[1] + 180
  PointInfo[1] <- PointInfo[1] %% 360
  JensenCone2 <- JensenTrapezoid(PointInfo[1], c(x2, y2), 50)

  I1 <- Intersect(JensenCone1[, 3], JensenCone1[, 2], c(x2, y2), JensenCone2[, 2])
  I2 <- Intersect(JensenCone1[, 1], JensenCone1[, 4], c(x2, y2), JensenCone2[, 2])

  P1H <- Height(x1, y1, Elev)
  P2H <- Height(x2, y2, Elev)

  ConeCenter <- (I1 + I2) / 2
  ConeCenter <- c(ConeCenter[1] * e$FarmVars$MeterWidth, ConeCenter[2] * e$FarmVars$MeterWidth, P1H + z1)

  P2 <- c(x2 * e$FarmVars$MeterWidth, y2 * e$FarmVars$MeterWidth, P2H + z2)

  a <- (sqrt((I1[1] - I2[1]) ^ 2 + (I1[2] - I2[2]) ^ 2) / 2) * e$FarmVars$MeterWidth
  b <- e$FarmVars$r0
  d <- sqrt((ConeCenter[1] - P2[1]) ^ 2 + (ConeCenter[2] - P2[2]) ^ 2 + (ConeCenter[3] - P2[3]) ^ 2)

  Intersection <- Area(a, b, d)
  InWakeFrac <- Intersection / (pi * b ^ 2)

  Penal <- min(c(PairPenalty(x2, y2, x1, y1, Dirs, SDs), PairPenalty(JensenCone2[1, 1], JensenCone2[2, 1], x1, y1, Dirs, SDs), PairPenalty(JensenCone2[1, 2], JensenCone2[2, 2], x1, y1, Dirs, SDs))) # c(P2, A', B')

  Deduction <- 1 - Penal
  Deduction <- Deduction * InWakeFrac
  PartPenal <- 1 - Deduction

  if (DrawTop)
  {
    # Draw cone
    XLim <- range(JensenCone1[1, ])
    XLim[1] <- min(XLim[1], JensenCone2[1, 1], JensenCone2[1, 2])
    XLim[2] <- max(XLim[2], JensenCone2[1, 1], JensenCone2[1, 2])

    YLim <- range(JensenCone1[2, ])
    YLim[1] <- min(YLim[1], JensenCone2[2, 1], JensenCone2[2, 2])
    YLim[2] <- max(YLim[2], JensenCone2[2, 1], JensenCone2[2, 2])

    XLim[1] <- XLim[1] - (0.01 * XLim[1])
    XLim[2] <- XLim[2] + (0.01 * XLim[2])

    YLim[1] <- YLim[1] - (0.01 * YLim[1])
    YLim[2] <- YLim[2] + (0.01 * YLim[2])

    plot(x = JensenCone1[1, 1], y = JensenCone1[2, 1], xlim = XLim, ylim = YLim, xlab = "", ylab = "")
    text(JensenCone1[1, 1], JensenCone1[2, 1], labels = "A", pos = 4)
    points(x = JensenCone1[1, 2], y = JensenCone1[2, 2])
    text(JensenCone1[1, 2], JensenCone1[2, 2], labels = "B", pos = 1)
    points(x = JensenCone1[1, 3], y = JensenCone1[2, 3])
    text(JensenCone1[1, 3], JensenCone1[2, 3], labels = "C", pos = 2)
    points(x = JensenCone1[1, 4], y = JensenCone1[2, 4])
    text(JensenCone1[1, 4], JensenCone1[2, 4], labels = "D", pos = 3)
    points(x = x1, y = y1, pch = 16)
    text(x = x1, y = y1, labels = "P1", pos = 2)
    lines(x = JensenCone1[1, ], y = JensenCone1[2, ])
    lines(x = c(JensenCone1[1, 1], JensenCone1[1, 4]), y = c(JensenCone1[2, 1], JensenCone1[2, 4]))

    # Draw inner disc
    points(x = JensenCone2[1, 1], y = JensenCone2[2, 1])
    text(JensenCone2[1, 1], JensenCone2[2, 1], labels = "A'", pos = 4)
    points(x = JensenCone2[1, 2], y = JensenCone2[2, 2])
    text(JensenCone2[1, 2], JensenCone2[2, 2], labels = "B'", pos = 1)
    points(x = x2, y = y2, pch = 16)
    text(x = x2, y = y2, labels = "P2", pos = 2)
    lines(x = c(JensenCone2[1, 1], JensenCone2[1, 2]), y = c(JensenCone2[2, 1], JensenCone2[2, 2]))

    # Draw I1
    points(I1[1], I1[2], pch = 5)
    text(x = I1[1], y = I1[2], labels = "I1", pos = 2)
    abline(v = I1[1], col = "grey", lty = 2)
    abline(h = I1[2], col = "grey", lty = 2)

    # Draw I2
    points(I2[1], I2[2], pch = 5)
    text(x = I2[1], y = I2[2], labels = "I2", pos = 2)
    abline(v = I2[1], col = "grey", lty = 2)
    abline(h = I2[2], col = "grey", lty = 2)

    # Draw CC
    points(x = ConeCenter[1] / e$FarmVars$MeterWidth, y = ConeCenter[2] / e$FarmVars$MeterWidth, pch = 16, col = "blue")
    text(x = ConeCenter[1] / e$FarmVars$MeterWidth, y = ConeCenter[2] / e$FarmVars$MeterWidth, labels = "CC", pos = 2)
  }

  if (DrawFront)
  {
    CC <- c(d, 0)
    P2Draw <- c(d, (P2H + z2) - (P1H + z1))
    if (PointInPolygon(JensenCone1, JensenCone2[, 1])) P2Draw <- -P2Draw

    cos_theta <- (t(CC) %*% P2Draw) / (Norm(CC) * Norm(P2Draw))
    theta <- drop(acos(cos_theta) / pi * 180)

    theta <- (360 - theta) %% 360
    Mat <- cbind(c(cos(theta * pi / 180), sin(theta * pi / 180)), c(-sin(theta * pi / 180), cos(theta * pi / 180)))
    P2Draw <- Mat %*% CC

    CC <- c(0, 0)

    if (CC[1] < P2Draw[1])
    {
      Left <- CC[1]
      Right <- P2Draw[1]
    } else
    {
      Left <- P2Draw[1]
      Right <- CC[1]
    }

    if (CC[2] < P2Draw[2])
    {
      Bottom <- CC[2]
      Top <- P2Draw[2]
    } else
    {
      Bottom <- P2Draw[2]
      Top <- CC[2]
    }

    if (a < b) R <- b else R <- a
    R <- R * 1.4
    Left <- Left - R
    Right <- Right + R
    Bottom <- Bottom - R
    Top <- Top + R

    if ((P1H + z1) > (P2H + z2)) P2Draw[2] <- -P2Draw[2]
    P2Draw[2] <- -P2Draw[2]

    plot(x = c(CC[1], P2Draw[1]), y = c(CC[2], P2Draw[2]), xlim = c(Left, Right), ylim = c(Bottom, Top), xlab = "", ylab = "", xaxt = "n", yaxt = "n", pch = 4, bty = "n")

    plotrix::draw.circle(x = CC[1], y = CC[2], radius = a, col = rgb(red = 255, green = 127, blue = 80, alpha = 127, maxColorValue = 255))
    plotrix::draw.circle(x = P2Draw[1], y = P2Draw[2], radius = b, col = rgb(red = 0, green = 191, blue = 255, alpha = 127, maxColorValue = 255))
  }

  if (DrawSide)
  {
    YRange = c(min(c(P1H, P2H)), max(c(P1H, P2H)))
    YRange[2] = YRange[2] + z2 + e$FarmVars$r0
    YRange[2] = YRange[2] * 1.1

    Distance = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2) * e$FarmVars$MeterWidth
    XRange = c(0, Distance)

    MyCone = JensenTrapezoid(0, c(XRange[2] / e$FarmVars$MeterWidth, (P1H + z1) / e$FarmVars$MeterWidth), XRange[2]) * e$FarmVars$MeterWidth

    #plot.new plus terrain
    plot(x = c(0, max(XRange)), y = c(P2H, P1H), type = "l", lty = 2, xlim = XRange, ylim = YRange, xlab = "Meters", ylab = "Meters", main = "", sub = "(dashed line not actual terrain)")

    #towers
    lines(x = c(0, 0), y = c(P2H, P2H + z2), col = "grey", lwd = 2)
    lines(x = c(XRange[2], XRange[2]), y = c(P1H, P1H + z1), col = "grey", lwd = 2)

    #rotors
    lines(x = c(0, 0), y = c(P2H + z2 - e$FarmVars$r0, P2H + z2 + e$FarmVars$r0), lwd = 3)
    lines(x = c(XRange[2], XRange[2]), y = c(P1H + z1 - e$FarmVars$r0, P1H + z1 + e$FarmVars$r0), lwd = 3)

    #hubs
    points(0, P2H + z2, pch = 16)
    points(XRange[2], P1H + z1, pch = 16)

    #wake
    lines(x = c(MyCone[1, 1], MyCone[1, 4]), y = c(MyCone[2, 1], MyCone[2, 4]), col = "blue")
    lines(x = c(MyCone[1, 2], MyCone[1, 3]), y = c(MyCone[2, 2], MyCone[2, 3]), col = "blue")
  }

  return(PartPenal)
}

#######################################################################################################################################

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

#######################################################################################################################################

JensenAngle <- function(x)
{
	alpha <- 0.5 / (log(e$FarmVars$z / e$FarmVars$z0))
	alpha_x <- alpha * x
	r <- e$FarmVars$r0 + alpha_x
	Norm <- sqrt((r ^ 2) + (x ^ 2))
	Cos_gamma <- x / Norm

	return(acos(Cos_gamma) / pi * 180)
}

#######################################################################################################################################

JensenFactor <- function(x)
{
	alpha <- 0.5 / (log(e$FarmVars$z / e$FarmVars$z0))
	alpha_x <- alpha * x
	r <- e$FarmVars$r0 + alpha_x

	return((1 - ((2 / 3) * (e$FarmVars$r0 / r) ^ 2)))
}

#######################################################################################################################################

GenerateGauss <- function(u = 8, refHeight = 10, maxX = 500, resY = 100, resZ = 100, Verbose = TRUE)
{
  xMax <- maxX

  MultNorm <- function(x, mean, sigma)
  {
    if (nrow(sigma) != length(x)) stop("Non-matching arguments.")
    if (length(mean) != length(x)) stop("Non-matching arguments.")
    if (nrow(sigma) != ncol(sigma)) stop("Non-quadratic covariance matrix.")
    if (!isSymmetric(sigma)) stop("Non-symmetric covariance matrix.")
    if (any(eigen(sigma)$values <= 0)) stop("Covariance matrix is not positive definite.")

    p <- length(x)
    return(drop(1 / (sqrt(((2 * pi) ^ p) * det(sigma))) * exp(-0.5 * t(x - mean) %*% solve(sigma) %*% (x - mean))))
  }

  ## Begin prepare
  Alpha <- 0.5 / log(e$FarmVars$z / e$FarmVars$z0)
  r <- e$FarmVars$r0 + Alpha * xMax
  zMax <- e$FarmVars$z + r
  SeqZ <- seq(from = zMax, to = 0, length.out = resZ)
  SeqY <- seq(from = -0.5 * zMax, to = 0.5 * zMax, length.out = resY)
  u <- WindspeedHellmann(v0 = u, HH = SeqZ, refHeight = refHeight)
  Sigma2 <- cbind(c(e$FarmVars$r0 ^ 2, 0), c(0, e$FarmVars$r0 ^ 2))
  Mu2 <- c(e$FarmVars$z, 0)
  BaseGauss <- MultNorm(x = c(e$FarmVars$z, 0), mean = Mu2, sigma = Sigma2)
  ## End prepare

  SeqX <- 1:xMax

  Gauss <- array(dim = c(length(SeqZ), length(SeqY), length(SeqX))) #Tensor: Height, width, depth

  for (x in SeqX)
  {
    r <- e$FarmVars$r0 + Alpha * x
    Sigma <- cbind(c(r ^ 2, 0), c(0, r ^ 2))
    Mu <- c(e$FarmVars$z, 0)

    Dens <- matrix(ncol = resY, nrow = resZ)

    for (y in 1:resY)
    {
      for (z in 1:resZ)
      {
        Dens[z, y] <- MultNorm(x = c(SeqZ[z], SeqY[y]), mean = Mu, sigma = Sigma) / BaseGauss
      }
    }

    P <- JensenFactor(x)

    Dens <- Dens * (P - 1)
    Dens <- Dens + 1

    for (z in 1:resZ) Dens[z, ] <- Dens[z, ] * u[z]

    Gauss[, , x] <- Dens

    if (Verbose)
    {
      message(paste(x, " of ", xMax, " complete (", round(x / xMax * 100), " %).", sep = ""))
      flush.console()
    }
  }

  return(Gauss)
}

#######################################################################################################################################

GaussWS <- function(Gauss, x, y, z)
{
  xMax <- dim(Gauss)[3]

  Alpha <- 0.5 / log(e$FarmVars$z / e$FarmVars$z0)
  r <- e$FarmVars$r0 + Alpha * xMax
  zMax <- e$FarmVars$z + r

  SeqX <- 1:xMax
  SeqY <- seq(from = -zMax, to = zMax, length.out = dim(Gauss)[2])
  SeqZ <- seq(from = zMax, to = 0, length.out = dim(Gauss)[1])

  DiffX <- abs(x - SeqX)
  x <- which(DiffX == min(DiffX))

  DiffY <- abs(y - SeqY)
  y <- which(DiffY == min(DiffY))

  DiffZ <- abs(z - SeqZ)
  z <- which(DiffZ == min(DiffZ))

  return(Gauss[z, y, x])
}

#######################################################################################################################################

QuickGauss3D <- function(x, y, z, u = 8, refHeight = 10)
{
  MultNorm <- function(x, mean, sigma)
  {
    if (nrow(sigma) != length(x)) stop("Non-matching arguments.")
    if (length(mean) != length(x)) stop("Non-matching arguments.")
    if (nrow(sigma) != ncol(sigma)) stop("Non-quadratic covariance matrix.")
    if (!isSymmetric(sigma)) stop("Non-symmetric covariance matrix.")
    if (any(eigen(sigma)$values <= 0)) stop("Covariance matrix is not positive definite.")

    p <- length(x)
    return(drop(1 / (sqrt(((2 * pi) ^ p) * det(sigma))) * exp(-0.5 * t(x - mean) %*% solve(sigma) %*% (x - mean))))
  }

  Alpha <- 0.5 / log(e$FarmVars$z / e$FarmVars$z0)
  r <- e$FarmVars$r0 + Alpha * x
  u <- WindspeedHellmann(v0 = u, HH = z, refHeight = refHeight)

  Sigma <- cbind(c(r ^ 2, 0), c(0, r ^ 2))
  Mu <- c(e$FarmVars$z, 0)

  Sigma2 <- cbind(c(e$FarmVars$r0 ^ 2, 0), c(0, e$FarmVars$r0 ^ 2))
  Mu2 <- c(e$FarmVars$z, 0)

  Dens <- MultNorm(x = c(z, y), mean = Mu, sigma = Sigma) / MultNorm(x = c(e$FarmVars$z, 0), mean = Mu2, sigma = Sigma2)

  P <- JensenFactor(x)

  Dens <- Dens * (P - 1)
  Dens <- Dens + 1
  Dens <- Dens * u

  return(Dens)
}

#######################################################################################################################################

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

#######################################################################################################################################

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

#######################################################################################################################################

SwitchProfile <- function(On, Result)
{
  On[On < 0.5] <- 0
  On[On >= 0.5] <- 1

  if (exists("ParallelProfile"))
  {
    e$FarmVars <- ParallelProfile$FarmVars
    e$FarmData <- ParallelProfile$FarmData
  }

  xSel <- seq(from = 1, to = length(Result$par) - 1, by = 2)
  x <- Result$par[xSel]
  y <- Result$par[xSel + 1]

  n <- length(x)

  ThisSet <- as.numeric()
  for (i in seq(from = 1, to = (2 * n), by = 2))
  {
    Index <- (i + 1) / 2
    if (On[Index])
    {
      ThisSet[i] <- x[Index]
      ThisSet[i + 1] <- y[Index]
    }
  }
  ThisSet <- na.omit(ThisSet)

  if (length(ThisSet) >= 1) return(Profit(ThisSet)) else return(0)
}

#######################################################################################################################################

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
    Elev <- e$FarmData[[5]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  } else
  {
    Adj <- FarmData[[1]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    Dirs <- FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    SDs <- FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    Elev <- FarmData[[5]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
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
      if (e$FarmVars$Partial) if (i != j) M[i, j] <- PartialJensen(x[j], y[j], x[i], y[i], Dirs, SDs, Elev)
    }
  }
  M[M == 0] <- 1
  p <- rep(1, n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      p[i] <- p[i] * (M[j, i] ^ 3)
    }
  }

  PointYield <- as.numeric()

  if (exists("Yield", envir = e, inherits = FALSE))
  {
    for (i in 1:n) PointYield[i] <- e$Yield(x[i], y[i], Adj) * p[i]
    Revenue <- PointYield
  } else
  {
    for (i in 1:n) PointYield[i] <- Yield(x[i], y[i], Adj) * p[i]
    Revenue <- (PointYield * e$FarmVars$Price)
  }

  if (exists("Cost", envir = e, inherits = FALSE))
  {
    retVal <- sum(Revenue - e$Cost(x, y))
  } else
  {
    retVal <- sum(Revenue - Cost(x, y))
  }

  return(-retVal) #Make it a minimization problem.
}

#######################################################################################################################################

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

#######################################################################################################################################

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

#######################################################################################################################################

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
    Elev <- e$FarmData[[5]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  } else
  {
    Adj <- FarmData[[1]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    Dirs <- FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    SDs <- FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    Elev <- FarmData[[5]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  }

  n <- length(x)

  M <- rep(1, n)
  M <- M %*% t(M)
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      M[i, j] <- M[i, j] * PairPenalty(x[j], y[j], x[i], y[i], Dirs, SDs)
      if (e$FarmVars$Partial) if (i != j) M[i, j] <- PartialJensen(x[j], y[j], x[i], y[i], Dirs, SDs, Elev)
    }
  }
  M[M == 0] <- 1
  p <- rep(1, n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      p[i] <- p[i] * (M[j, i] ^ 3)
    }
  }

  PointYield <- as.numeric()

  if (exists("Yield", envir = e, inherits = FALSE))
  {
    for (i in 1:n) PointYield[i] <- e$Yield(x[i], y[i], Adj) * p[i]
    Revenue <- PointYield
  } else
  {
    for (i in 1:n) PointYield[i] <- Yield(x[i], y[i], Adj) * p[i]
    Revenue <- (PointYield * e$FarmVars$Price)
  }

  if (exists("Cost", envir = e, inherits = FALSE))
  {
    retVal <- Revenue - e$Cost(x, y)
  } else
  {
    retVal <- Revenue - Cost(x, y)
  }

  retVal <- cbind(1:n, retVal)
  colnames(retVal) <- c("Turbine", "Profit")
  return(retVal)
}

#######################################################################################################################################

PlotResult <- function(Result, ImageData = NULL, DoLabels = FALSE, Labels = "IDs")
{
  if (!any(names(Result) == "par")) stop("Not a valid optimizer result provided.")

  if (is.null(ImageData))
  {
    if (exists("FarmData", envir = e, inherits = FALSE))
    {
  	  Adj <- e$FarmData[[1]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    } else
    {
      Adj <- FarmData[[1]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    }
  } else
  {
    Adj <- ImageData
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

#######################################################################################################################################

ShowWakePenalizers <- function(Result, Cones = TRUE, All = FALSE, VectorField = TRUE, ImposeMST = FALSE, NS = FALSE, NSOnly = FALSE, Exaggerate = TRUE, DoubleLine = TRUE, Frames = 100, Alpha = 0.5, MaxContrast = TRUE, Soften = TRUE)
{
  if (!any(names(Result) == "par")) stop("Not a valid optimizer result provided.")

  if (ImposeMST) All <- TRUE
  if (ImposeMST == TRUE & Cones == FALSE)
  {
    #Nothing to do here yet.
  } else
  {
    if (All) Cones <- TRUE
  }
  if (NSOnly)
  {
    NS = TRUE
    VectorField = FALSE
  }

  if (exists("FarmData", envir = e, inherits = FALSE))
  {
    Dirs <- e$FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    SDs <- e$FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    Elev <- e$FarmData[[5]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
  } else
  {
    Dirs <- FarmData[[3]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    SDs <- FarmData[[4]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
    Elev <- FarmData[[5]][e$FarmVars$StartPoint:e$FarmVars$EndPoint, e$FarmVars$StartPoint:e$FarmVars$EndPoint]
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
      if (i != j)
      {
        M[i, j] <- M[i, j] * PairPenalty(x[j], y[j], x[i], y[i], Dirs, SDs)
        if (e$FarmVars$Partial) M[i, j] <- PartialJensen(x[j], y[j], x[i], y[i], Dirs, SDs, Elev)
      }
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

  if (NSOnly == FALSE)
  {
    if (All == FALSE)
    {
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
        points(x[Causers[i]], y[Causers[i]], col = "gray", pch = 16)
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
      legend("topleft", c("Causers", "Sufferers"), bg = "white", pch = 16, col = c("gray", "gold4"))

    } else
    {
      plot(x, y, ylim = c(0, 1), xlim = c(0, 1), xlab = "", ylab = "")
      points(x, y, col = "gold4", pch = 16)
      text(x, y, labels = 1:n, pos = 3)

      if (Cones)
      {
        for (i in 1:n)
        {
          WindDir <- GetDirInfo(x[i], y[i], Dirs, SDs)
          WindDir <- Geo2Ari(WindDir[1])
          WindDir <- WindDir + 180
          WindDir <- WindDir %% 360
          JensenCone <- JensenTrapezoid(WindDir, c(x[i], y[i]), e$FarmVars$MeterWidth * 2)

          lines(x = c(JensenCone[1, 1], JensenCone[1, 2]), y = c(JensenCone[2, 1], JensenCone[2, 2]), col = "red")
          lines(x = c(JensenCone[1, 2], JensenCone[1, 3]), y = c(JensenCone[2, 2], JensenCone[2, 3]), col = "red")
          lines(x = c(JensenCone[1, 3], JensenCone[1, 4]), y = c(JensenCone[2, 3], JensenCone[2, 4]), col = "red")
          lines(x = c(JensenCone[1, 4], JensenCone[1, 1]), y = c(JensenCone[2, 4], JensenCone[2, 1]), col = "red")
        }
      }
    }
  }

  if (VectorField) ImposeVectorField(xNum = 7, yNum = 7, ColMain = "gray", ColBand = "gray30", DoSDs = TRUE)

  if (ImposeMST)
  {
    d <- cbind(x, y)
    Tree <- invisible(emstreeR::ComputeMST(d, verbose = FALSE))
    for (i in 1:(nrow(Tree) - 1))
    {
      lines(x = c(Tree$x[Tree$from[i]], Tree$x[Tree$to[i]]), y = c(Tree$y[Tree$from[i]], Tree$y[Tree$to[i]]), lty = 2)
    }
  }

  if (NS)
  {
    Iter <- 4
    dT <- 0.1
    Diff <- 0.000001
    Visc <- 0.000001
    N <- 250
    StepConst <- N / ncol(Dirs)

    GodEnvir <- environment()

    s <- rep(0, times = N) %*% t(rep(0, times = N))
    Density <- s
    Vx <- s
    Vy <- s
    Vx0 <- s
    Vy0 <- s

    ObstaclesX <- x * N
    ObstaclesY <- y * N

    addDensity <- function(x, y, Amount)
    {
      Density[x, y] <- Amount

      assign(x = "Density", value = Density, envir = GodEnvir)
    }

    addVelocity <- function(x, y, AmountX, AmountY)
    {
      Vx[x, y] <- Vx[x, y] + AmountX
      Vy[x, y] <- Vy[x, y] + AmountY

      assign(x = "Vx", value = Vx, envir = GodEnvir)
      assign(x = "Vy", value = Vy, envir = GodEnvir)
    }

    BresLine <- function(Start, End)
    {
      x1 <- round(Start[1])
      y1 <- round(Start[2])
      x2 <- round(End[1])
      y2 <- round(End[2])

      RetValX <- as.numeric()
      RetValY <- as.numeric()

      Steep <- FALSE

      DeltaX <- abs(x2 - x1)
      if (x2 - x1 > 0) StepX <- 1 else StepX <- -1

      DeltaY <- abs(y2 - y1)
      if (y2 - y1 > 0) StepY <- 1 else StepY <- -1

      if (DeltaY > DeltaX)
      {
        Steep <- TRUE

        # Swap x1 and y1
        Temp <- x1
        x1 <- y1
        y1 <- Temp

        # Swap DeltaX and DeltaY
        Temp <- DeltaX
        DeltaX <- DeltaY
        DeltaY <- Temp

        # Swap StepX and StepY
        Temp <- StepX
        StepX <- StepY
        StepY <- Temp
      }

      Delta <- (DeltaY * 2) - DeltaX

      i <- 1
      for (Coord in 0:(DeltaX - 1))
      {
        if (Steep)
        {
          RetValX[i] <- y1
          RetValY[i] <- x1
        } else
        {
          RetValX[i] <- x1
          RetValY[i] <- y1
        }
        i <- i + 1

        while (Delta >= 0)
        {
          y1 <- y1 + StepY
          Delta <- Delta - (DeltaX * 2)
        }

        x1 <- x1 + StepX
        Delta <- Delta + (DeltaY * 2)
      }

      RetValX[i] <- x2
      RetValY[i] <- y2

      RetVal = cbind(RetValX, RetValY)
      return(RetVal)
    }

    AddParkData <- function()
    {
      X <- seq(from = 1, to = (N / StepConst), by = (N / StepConst) / 25)
      Y <- X

      Rad <- (450 - Dirs[X, Y]) %% 360
      Rad <- Rad * pi / 180
      WDx <- cos(Rad)
      WDy <- sin(Rad)
      addVelocity(round(X * StepConst), round(Y * StepConst), WDx, WDy)

      AddAmount <- 100
      i <- 1:N
      addDensity(N - 2, i, AddAmount) #right
      addDensity(2, i, AddAmount)     #left
      addDensity(i, N - 2, AddAmount) #top
      addDensity(i, 2, AddAmount)     #bottom
    }

    Diffuse1 <- function(b, x, x0, Diff, dT) #1, Vx0, Vx, Visc, dT
    {
      a <- dT * Diff * (N - 2) * (N - 2)
      x <- LinSolve(b, x, x0, a, 1 + 6 * a)

      assign(x = "Vx0", value = x, envir = GodEnvir)
    }

    Diffuse2 <- function(b, x, x0, Diff, dT) #2, Vy0, Vy, Visc, dT
    {
      a <- dT * Diff * (N - 2) * (N - 2)
      x <- LinSolve(b, x, x0, a, 1 + 6 * a)

      assign(x = "Vy0", value = x, envir = GodEnvir)
    }

    Diffuse3 <- function(b, x, x0, Diff, dT) #0, s, Density, Diff, dT
    {
      a <- dT * Diff * (N - 2) * (N - 2)
      x <- LinSolve(b, x, x0, a, 1 + 6 * a)

      assign(x = "s", value = x, envir = GodEnvir)
    }

    LinSolve <- function(b, x, x0, a, c)
    {
      cRecip <- 1 / c

      for (t in 1:Iter)
      {
        j <- 2:(N - 2)
        i <- 2:(N - 2)
        x[i, j] <- (x0[i, j] + a * (x[i + 1, j] + x[i - 1, j] + x[i, j + 1] + x[i, j - 1])) * cRecip

        x <- SetBound(b, x)
      }

      return(x)
    }

    Project1 <- function(velocX, velocY, p, div) #Vx0, Vy0, Vx, Vy
    {
      j <- 2:(N - 2)
      i <- 2:(N - 2)
      div[i, j] <- (-0.5 * (velocX[i + 1, j] - velocX[i - 1, j] + velocY[i, j + 1] - velocY[i, j - 1])) / N
      p[i, j] <- 0

      div <- SetBound(0, div)
      p <- SetBound(0, p)
      p <- LinSolve(0, p, div, 1, 6)

      assign(x = "Vx", value = p, envir = GodEnvir)
      assign(x = "Vy", value = div, envir = GodEnvir)

      j <- 2:(N - 2)
      i <- 2:(N - 2)
      velocX[i, j] <- velocX[i, j] - 0.5 * (p[i + 1, j] - p[i - 1, j]) * N
      velocY[i, j] <- velocY[i, j] - 0.5 * (p[i, j + 1] - p[i, j - 1]) * N

      velocX <- SetBound(1, velocX)
      velocY <- SetBound(2, velocY)

      assign(x = "Vx0", value = velocX, envir = GodEnvir)
      assign(x = "Vy0", value = velocY, envir = GodEnvir)
    }

    Project2 <- function(velocX, velocY, p, div) #Vx, Vy, Vx0, Vy0
    {
      j <- 2:(N - 2)
      i <- 2:(N - 2)
      div[i, j] <- (-0.5 * (velocX[i + 1, j] - velocX[i - 1, j] + velocY[i, j + 1] - velocY[i, j - 1])) / N
      p[i, j] <- 0

      div <- SetBound(0, div)
      p <- SetBound(0, p)
      p <- LinSolve(0, p, div, 1, 6)

      assign(x = "Vx0", value = p, envir = GodEnvir)
      assign(x = "Vy0", value = div, envir = GodEnvir)

      j <- 2:(N - 2)
      i <- 2:(N - 2)
      velocX[i, j] <- velocX[i, j] - 0.5 * (p[i + 1, j] - p[i - 1, j]) * N
      velocY[i, j] <- velocY[i, j] - 0.5 * (p[i, j + 1] - p[i, j - 1]) * N

      velocX <- SetBound(1, velocX)
      velocY <- SetBound(2, velocY)

      assign(x = "Vx", value = velocX, envir = GodEnvir)
      assign(x = "Vx", value = velocY, envir = GodEnvir)
    }

    Advect <- function(b, d, d0, velocX, velocY, dT)
    {
      dTX <- dT * (N - 2)
      dTY <- dT * (N - 2)

      NFloat <- N - 2

      jFloat <- 1
      for (j in 2:(N - 2))
      {
        jFloat <- jFloat + 1

        iFloat <- 1
        for (i in 2:(N - 2))
        {
          iFloat <- iFloat + 1

          Tmp1 <- dTX * velocX[i, j]
          Tmp2 <- dTY * velocY[i, j]
          X <- iFloat - Tmp1
          Y <- jFloat - Tmp2

          if (X < 0.5) X <- 0.5
          if (X > NFloat + 0.5) X <- NFloat + 0.5
          i0 <- trunc(X)
          i1 <- i0 + 1

          if (Y < 0.5) Y <- 0.5
          if (Y > NFloat + 0.5) Y <- NFloat + 0.5
          j0 <- trunc(Y)
          j1 <- j0 + 1

          s1 <- X - i0
          s0 <- 1 - s1
          t1 <- Y - j0
          t0 <- 1 - t1

          i0i <- i0
          i1i <- i1
          j0i <- j0
          j1i <- j1

          Temp <- s0 * (t0 * d0[i0i, j0i] + t1 * d0[i0i, j1i]) + s1 * (t0 * d0[i1i, j0i] + t1 * d0[i1i, j1i])
          if (length(Temp) > 0) d[i, j] <- Temp else d[i, j] <- 0
        }
      }

      d <- SetBound(b, d)

      return(d)
    }

    SetBound <- function(b, x)
    {
      i <- 2:(N - 2)
      x[i, 0] <- ifelse(b == 2, -x[i, 1], x[i, 1])
      x[i, N - 1] <- ifelse(b == 2, -x[i, N - 2], x[i, N - 2])

      j <- 2:(N - 2)
      x[0, j] <- ifelse(b == 1, -x[1, j], x[1, j])
      x[N - 1, j] <- ifelse(b == 1, -x[N - 2, j], x[N - 2, j])

      x[0, 0] <- 0.5 * (x[1, 0] + x[0, 1])
      x[0, N - 1] <- 0.5 * (x[1, N - 1] + x[0, N - 2])
      x[N - 1, 0] <- 0.5 * (x[N - 2, 0] + x[N - 1, 1])
      x[N - 1, N - 1] <- 0.5 * (x[N - 2, N - 1] + x[N - 1, N - 2])

      return(x)
    }

    Step <- function()
    {
      applyObstacles(Density, Vx, Vy, Vx0, Vy0)

      Diffuse1(1, Vx0, Vx, Visc, dT)
      Diffuse2(2, Vy0, Vy, Visc, dT)

      Project1(Vx0, Vy0, Vx, Vy)

      Vx <- Advect(1, Vx, Vx0, Vx0, Vy0, dT)
      Vy <- Advect(2, Vy, Vy0, Vx0, Vy0, dT)

      assign(x = "Vx", value = Vx, envir = GodEnvir)
      assign(x = "Vy", value = Vy, envir = GodEnvir)

      Project2(Vx, Vy, Vx0, Vy0)
      Diffuse3(0, s, Density, Diff, dT)
      Density <- Advect(0, Density, s, Vx, Vy, dT)

      assign(x = "Density", value = Density, envir = GodEnvir)

      applyObstacles(Density, Vx, Vy, Vx0, Vy0)
    }

    applyObstacles <- function(DensityO, VxO, VyO, Vx0O, Vy0O)
    {
      if (length(ObstaclesX) < 1) return(0)

      if (Exaggerate)
      {
        for (i in 1:length(ObstaclesX))
        {
          for (x in -2:2)
          {
            for (y in -2:2)
            {
              if (ObstaclesX[i] + x >= 1 & ObstaclesX[i] + x <= N & ObstaclesY[i] + y >= 1 & ObstaclesY[i] + y <= N)
              {
                DensityO[ObstaclesX[i] + x, ObstaclesY[i] + y] <- 0
                VxO[ObstaclesX[i] + x, ObstaclesY[i] + y] <- 0
                VyO[ObstaclesX[i] + x, ObstaclesY[i] + y] <- 0
                Vx0O[ObstaclesX[i] + x, ObstaclesY[i] + y] <- 0
                Vy0O[ObstaclesX[i] + x, ObstaclesY[i] + y] <- 0
              }
            }
          }
        }
      } else
      {
        for (i in 1:length(ObstaclesX))
        {
          MyDir <- GetDirInfo(x[i], y[i], Dirs, SDs)[1]
          MyDir <- Geo2Ari(MyDir)
          MyDir <- MyDir + 180
          MyDir <- MyDir %% 360

          Trapez <- JensenTrapezoid(MyDir, c(x[i], y[i]), 500)
          A <- Trapez[, 1] * N
          B <- Trapez[, 2] * N

          MyLine <- BresLine(A, B)

          for (j in 1:length(MyLine[, 1]))
          {
            if (MyLine[j, 1] >= 1 & MyLine[j, 1] <= N & MyLine[j, 2] >= 1 & MyLine[j, 2] <= N)
            {
              DensityO[MyLine[j, 1], MyLine[j, 2]] <- 0
              VxO[MyLine[j, 1], MyLine[j, 2]] <- 0
              VyO[MyLine[j, 1], MyLine[j, 2]] <- 0
              Vx0O[MyLine[j, 1], MyLine[j, 2]] <- 0
              Vy0O[MyLine[j, 1], MyLine[j, 2]] <- 0
            }
          }

          if (DoubleLine)
          {
            if (abs(B[2] - A[2]) > abs(B[1] - A[1]))
            {
              A[1] <- A[1] + 1
              B[1] <- B[1] + 1
            } else
            {
              A[2] <- A[2] + 1
              B[2] <- B[2] + 1
            }
            MyLine <- BresLine(A, B)

            for (j in 1:length(MyLine[, 1]))
            {
              if (MyLine[j, 1] >= 1 & MyLine[j, 1] <= N & MyLine[j, 2] >= 1 & MyLine[j, 2] <= N)
              {
                DensityO[MyLine[j, 1], MyLine[j, 2]] <- 0
                VxO[MyLine[j, 1], MyLine[j, 2]] <- 0
                VyO[MyLine[j, 1], MyLine[j, 2]] <- 0
                Vx0O[MyLine[j, 1], MyLine[j, 2]] <- 0
                Vy0O[MyLine[j, 1], MyLine[j, 2]] <- 0
              }
            }
          }
        }
      }

      assign(x = "Density", value = DensityO, envir = GodEnvir)
      assign(x = "Vx", value = VxO, envir = GodEnvir)
      assign(x = "Vy", value = VyO, envir = GodEnvir)
      assign(x = "Vx0", value = Vx0O, envir = GodEnvir)
      assign(x = "Vy0", value = Vy0O, envir = GodEnvir)
    }

    MainLoop <- function(Frames)
    {
      pb <- progress::progress_bar$new(total = Frames)

      for (i in 1:Frames)
      {
        AddParkData()
        Step()
        pb$tick()
      }
    }

    MainLoop(Frames)

    if (MaxContrast) Density <- ((Density - min(Density)) / (max(Density) - min(Density))) * 255

    if (Soften)
    {
      a <- Density
      f = 1 / 9
      for (x in 2:(N - 1))
      {
        for (y in 2:(N - 1))
        {
          a[x, y] = f * Density[x - 1, y - 1] + f * Density[x, y - 1] + f * Density[x + 1, y - 1] + f * Density[x - 1, y] + f * Density[x, y] + f * Density[x + 1, y] + f * Density[x - 1, y + 1] + f * Density[x, y + 1] + f * Density[x + 1, y + 1]
        }
      }

      Density <- a
    }

    if (NSOnly)
    {
      image(Density, col = rainbow(255, alpha = 1))
    } else
    {
      image(Density, col = rainbow(255, alpha = Alpha), add = TRUE)
    }
  }

  return(invisible(M))
}

#######################################################################################################################################

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

			Bench <- e$FarmVars$BenchmarkSolution
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
			e$FarmVars$BenchmarkSolution <- Bench
		}
	} else
	{
		message("File was not downloaded. Internet connection error? Please try again.")
		return(invisible(NULL))
	}
}
