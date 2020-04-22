### R code from vignette source 'wflo.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("MASS")


###################################################
### code chunk number 2: Newfield (eval = FALSE)
###################################################
## e$FarmVars$StartPoint <- 1800
## e$FarmVars$Width <- 100
## e$FarmVars$EndPoint <- e$FarmVars$StartPoint + e$FarmVars$Width - 1
## e$FarmVars$MeterWidth <- 200 * e$FarmVars$Width
## e$FarmVars$MinDist <- e$FarmVars$MeterMinDist / e$FarmVars$MeterWidth


###################################################
### code chunk number 3: Result
###################################################
library(wflo)
NumTurbines <- 4
set.seed(2763)
Result <- optim(par = runif(NumTurbines * 2), fn = Profit,
  method = "L-BFGS-B", lower = rep(0, NumTurbines * 2),
  upper = rep(1, NumTurbines * 2))
Result


###################################################
### code chunk number 4: Slots
###################################################
Oldpar <- par(no.readonly = TRUE)

par(mfrow = c(4, 2), mar = c(0, 0, 0, 0), mai = c(0.2, 0.2, 0.2, 0.2))

MyNames <- c("AEP", "Wind speed", "Wind direction", "SD of wind direction", "Elevation", "Slope", "Hillside")

for (i in 1:7)
{
	X <- FarmData[[i]][1:25, 1:25]
	image(X, main = MyNames[i], xaxt = "n", yaxt = "n", bty = "n", col = rgb(0, 0, seq(from = 0, to = 1, length.out = e$FarmVars$Width * e$FarmVars$Width)))
}

par(Oldpar)


###################################################
### code chunk number 5: PlotResult
###################################################
PlotResult(Result)


###################################################
### code chunk number 6: ShowWake
###################################################
ShowWakePenalizers(Result)


###################################################
### code chunk number 7: CRS
###################################################
library(nloptr)
set.seed(1357)
Result <- crs2lm(x0 = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
Result


###################################################
### code chunk number 8: GA
###################################################
library(rgenoud)
set.seed(1357)
Dom = cbind(rep(0, 2 * NumTurbines), rep(1, 2 * NumTurbines))
Result <- genoud(fn = Profit, nvars = 2 * NumTurbines,
  starting.values = runif(NumTurbines * 2), Domains = Dom,
  boundary.enforcement = 2, print.level = 0)
Result


###################################################
### code chunk number 9: PSO
###################################################
library(pso)
set.seed(1357)
Result <- psoptim(par = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
Result


###################################################
### code chunk number 10: SANN
###################################################
lower <- rep(0, NumTurbines * 2)
upper <- rep(1, NumTurbines * 2)
Wrapper <- function(X)
{
	xSel <- seq(from = 1, to = length(X) - 1, by = 2)
	x <- X[xSel]
	y <- X[xSel + 1]

	if (any(x < lower) | any(x > upper) | any(y < lower) | any(y > upper))
	{
		return(sum(rep(e$FarmVars$UnitCost, length(x))))
	}

	return(Profit(X))
}
set.seed(1357)
Result <- optim(par = runif(NumTurbines * 2), fn = Wrapper,
  method = "SANN")
Result


###################################################
### code chunk number 11: Cost
###################################################
e$Cost <- function(x, y) #x, y \in R^n
{
	retVal <- rep(e$FarmVars$UnitCost, min(length(x), length(y)))
	retVal[x > 0.5] <- retVal[x > 0.5] * 2
	return(retVal)
}


###################################################
### code chunk number 12: CostOptim
###################################################
set.seed(1357)
Result <- psoptim(par = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
Result
rm(Cost, envir = e)


###################################################
### code chunk number 13: Cost2
###################################################
e$Cost <- function(x, y)
{
	n <- min(length(x), length(y))

	retVal <- rep(e$FarmVars$UnitCost, n)

	DistMat <- matrix(ncol = n, nrow = n)
	for (i in 1:n)
	{
		for (j in 1:n)
		{
			DistMat[i, j] <- sqrt((x[i] - x[j]) ^ 2 + (y[i] - y[j]) ^ 2)
		}
	}
	SumDist <- as.numeric()
	for (i in 1:n) SumDist[i] <- sum(DistMat[i, ])

	retVal <- retVal * SumDist
	return(retVal)
}
rm(Cost, envir = e)


###################################################
### code chunk number 14: Yield
###################################################
e$Yield <- function(x, y, AEP) #x, y \in R
{
	return(x + y)
}


###################################################
### code chunk number 15: YieldOptim
###################################################
set.seed(1357)
Result <- psoptim(par = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
Result
rm(Yield, envir = e)


###################################################
### code chunk number 16: Contribs
###################################################
NumTurbines <- 4
set.seed(1235)
Result <- optim(par = runif(NumTurbines * 2), fn = Profit,
  method = "L-BFGS-B", lower = rep(0, NumTurbines * 2),
  upper = rep(1, NumTurbines * 2))
Contribs <- ProfitContributors(Result)
Contribs
PlotResult(Result, DoLabels = TRUE, Labels = Contribs[, 2])


###################################################
### code chunk number 17: TurbsOutput
###################################################
Result <- list(par = e$FarmVars$BenchmarkSolution)
Result$value <- Profit(Result$par)
Result
PlotResult(Result)
ShowWakePenalizers(Result)


###################################################
### code chunk number 18: PSO-Fig
###################################################
set.seed(1357)
Result <- psoptim(par = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
PlotResult(Result)


###################################################
### code chunk number 19: TurbsContribs
###################################################
NumTurbines <- 4
set.seed(1235)
Result <- optim(par = runif(NumTurbines * 2), fn = Profit,
  method = "L-BFGS-B", lower = rep(0, NumTurbines * 2),
  upper = rep(1, NumTurbines * 2))
MyLabels <- ProfitContributors(Result)
PlotResult(Result, DoLabels = TRUE, Labels = MyLabels[, 2])


###################################################
### code chunk number 20: TurbsPlot
###################################################
Result <- list(par = e$FarmVars$BenchmarkSolution)
PlotResult(Result)


###################################################
### code chunk number 21: TurbsWake
###################################################
Result <- list(par = e$FarmVars$BenchmarkSolution)
ShowWakePenalizers(Result)


