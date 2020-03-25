### R code from vignette source 'wflo.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("MASS")


###################################################
### code chunk number 2: Result
###################################################
library(wflo)
NumTurbines <- 4
set.seed(1357)
Result <- optim(par = runif(NumTurbines * 2), fn = Profit,
  method = "L-BFGS-B", lower = rep(0, NumTurbines * 2),
  upper = rep(1, NumTurbines * 2))
Result


###################################################
### code chunk number 3: PlotResult
###################################################
PlotResult(Result)


###################################################
### code chunk number 4: CRS
###################################################
library(nloptr)
set.seed(1357)
Result <- crs2lm(x0 = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
Result


###################################################
### code chunk number 5: GA
###################################################
library(rgenoud)
set.seed(1357)
Dom = cbind(rep(0, 2 * NumTurbines), rep(1, 2 * NumTurbines))
Result <- genoud(fn = Profit, nvars = 2 * NumTurbines,
  starting.values = runif(NumTurbines * 2), Domains = Dom,
  boundary.enforcement = 2, print.level = 0)
Result


###################################################
### code chunk number 6: PSO
###################################################
library(pso)
set.seed(1357)
Result <- psoptim(par = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
Result


###################################################
### code chunk number 7: SANN
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
### code chunk number 8: PSO-Fig
###################################################
set.seed(1357)
Result <- psoptim(par = runif(NumTurbines * 2), fn = Profit,
  lower = rep(0, NumTurbines * 2), upper = rep(1, NumTurbines * 2))
PlotResult(Result)


