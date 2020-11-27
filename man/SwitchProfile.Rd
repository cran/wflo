\name{SwitchProfile}
\alias{SwitchProfile}

\title{Helper function for ex-post wind direction dependence analysis}

\description{
This function helps with analzing an already computed turbine setup. As \code{\link{Profit}} (together with an optimizer) finds an optimal setup given actual wind speeds in that area, the question arises whether that setup is 'robust' against changing wind directions. If for a given setup, there is no or only little wake dependence (analyze e.g. with \code{\link{ShowWakePenalizers}}), this may change if wind direction changes. Given a different wind direction, a situation may come up at which there is an upwind turbine and a local set of, say, two, three, or even more downwind turbines. As these downwind turbines are in the wake of the upwind turbines ('receiving' only reduced wind speeds), it may be beneficial to shut down the upwind turbine in order to have the downwind turbine receive full wind, thus possibly generating greater profit.\cr
\cr
\code{SwitchProfile} expects a 'bit' pattern of turbines (so-called 'profiles') at which '1' means that a turbine is to be used, while '0' means, the turbine is to be shut down. For a total of four tubines in a setup, profile \code{c(1, 1, 1, 1)} means that all turbines are used, while profile \code{c(1, 0, 0, 1)} means that turbines 1 and 4 are used, turbines 2 and 3 are to be shut down. If a wind farm is 'robust' against wind direction changes, it should return profiles of all ones for all wind directions.\cr
\cr
In the examples we show a short code sample that uses \code{SwitchProfile} to loop through 360 integer wind direction degrees and optimizes actual profiles. It is up to the researcher to use alternative optimizers in order to find (possibly) better profiles.
}

\usage{
SwitchProfile(On, Result)
}

\arguments{
	\item{On}{
		must be a numeric vector of length \code{n}, where \code{n} is the number of turbines in a setup.
	}
	\item{Result}{
		the actual setup as an optimizer result. Must be a list that contains the turbine locations is a slot "$par". For example, use \code{Result = list(par = e$FarmVars$BenchmarkSolution)}, see \code{\link{FarmVars}} for details.
	}
}

\details{
For compatibility with non-integer solution optimizers, \code{SwitchProfile} expects not only values equal to 0 or 1, but will assume that values below 0.5 are meant to be 0, and values >= 0.5 are coerced to 1.
}

\value{
\code{SwitchProfile} returns a single number. The result is the negative profit for the profile given. The result can be optimized in order to find the maximum profit profile, which should be identical to a profile of all ones if the setup is robust against wind direction changes.
}

\author{
Carsten Croonenbroeck
}

\seealso{
Use \code{\link{Profit}} for the target function used internally.
}

\examples{
Result <- list(par = e$FarmVars$BenchmarkSolution)

Profile1 <- rep(1, times = 20)

#Will be identical to Profit(Result$par):
SwitchProfile(Profile1, Result)

Profile2 <- Profile1
Profile2[8] <- 0 #Disable turbine 8.

#Returns farm profit if turbine 8 is off:
SwitchProfile(Profile2, Result)

##########################################

#Warning, this will overwrite e$FarmData. If not dispensable,
#backup before use.
#Computation may be slow.

\donttest{
Result <- list(par = e$FarmVars$BenchmarkSolution)
e$FarmData <- FarmData

n <- length(e$FarmVars$BenchmarkSolution) / 2

Profiles <- matrix(ncol = n, nrow = 360)
Dom <- cbind(rep(0, n), rep(1, n))

for (j in 1:360)
{
	e$FarmData$WindDirection <- matrix(data = j, ncol = e$FarmVars$Width,
	  nrow = e$FarmVars$Width)
	message(paste("Processing wind direction: ", j, " degrees.", sep = ""))
	flush.console()

	Result <- rgenoud::genoud(fn = SwitchProfile, nvars = n, starting.values = runif(n),
	  Domains = Dom, boundary.enforcement = 2, MemoryMatrix = FALSE, print.level = 0,
	  Result = Result)
	On <- Result$par
	On[On < 0.5] <- 0
	On[On >= 0.5] <- 1
	Profiles[j, ] <- On
}

#Now 'Profiles' is a matrix of 360 rows (degrees) and n columns.
#If the setup is robust against wind direction changes, all values
#in the columns should be ones. If there are zeros left, double
#check these profiles (compare to full turbines setup) manually,
#e.g. using

for (i in 1:nrow(Profiles))
{
	if (any(Profiles[i, ] == 0))
	{
		e$FarmData$WindDirection <- matrix(data = i,
    ncol = e$FarmVars$Width, nrow = e$FarmVars$Width)

		print(paste("Is unrestricted profit (",
    abs(Profit(e$FarmVars$BenchmarkSolution)), ") still
    greater than 'profile profit' (",
    abs(SwitchProfile(Profiles[i, ], Result)), ")? ",
    abs(Profit(e$FarmVars$BenchmarkSolution)) >=
    abs(SwitchProfile(Profiles[i, ], Result)), ".", sep = ""))
	}
}

#Is absolute profile profit actually greater than the absolute value
#function Profit() returns? If yes, then it is beneficial to turn off
#turbines (use profiles) if wind comes from the given angle.
}

##########################################

#As an alternative to using an optimizer to find the optimal profile,
#all possible profiles can be 'looped through' deterministically using
#the following few lines of code.

#Warning, this will overwrite e$FarmData. If not dispensable,
#backup before use.

#Be advised, computation may be very slow, even if run on modern
#machines and even though using parallelization.

\donttest{
Result <- list(par = e$FarmVars$BenchmarkSolution)
e$FarmData <- FarmData

ComputeProfile <- function(i)
{
	e$FarmData$WindDirection <- matrix(data = i, ncol = e$FarmVars$Width,
	  nrow = e$FarmVars$Width)

	BestProfit <- n * e$FarmVars$UnitCost
	BestProfile <- rep(0, n)
	for (j in 0:MaxNum)
	{
		ThisProfile <- as.numeric(rev(intToBits(j)))[(32 - (n - 1)):32]
		if (length(which(ThisProfile == 1)) > 1)
		{
			ThisProfit <- SwitchProfile(ThisProfile, Result)
		} else
		{
			ThisProfit <- n * e$FarmVars$UnitCost
		}

		if (ThisProfit < BestProfit)
		{
			BestProfit <- ThisProfit
			BestProfile <- ThisProfile
		}
	}
	return(BestProfile)
}

n <- length(e$FarmVars$BenchmarkSolution) / 2
MaxNum <- (2 ^ n) - 1

NumCores <- parallel::detectCores() - 1
cl <- snow::makeCluster(NumCores)

ParallelProfile <- new.env()
ParallelProfile$FarmVars <- e$FarmVars
ParallelProfile$FarmData <- e$FarmData

snow::clusterExport(cl, list = ls())
snow::clusterExport(cl, list = ls(envir = as.environment("package:wflo")))

doSNOW::registerDoSNOW(cl)

iterations <- 360
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
`\%dopar\%` <- foreach::`\%dopar\%`
AllProfiles <- foreach::foreach(i = 1:iterations,
  .options.snow = opts) \%dopar\% ComputeProfile(i)

close(pb)
snow::stopCluster(cl)
}

#Afterward, check the resulting list of profiles for remaining zeros as above.

##########################################

#Thinking the thought of wind direction dependent layouts further,
#Profit() can be used to optimize a farm setup that takes all wind
#direction wake patterns into account during the optimization
#stage already.

#Warning, this will overwrite \code{e$FarmData}. If not dispensable,
#backup before use.
#Computation may be slow.

\donttest{
e$FarmData <- FarmData

SumProfit <- function(X)
{
	MySum <- as.numeric()
	for (i in 1:360)
	{
		e$FarmData$WindDirection <- matrix(data = i, ncol = e$FarmVars$Width,
		  nrow = e$FarmVars$Width)
		MySum[i] <- Profit(X)
	}
	MySum <- sum(MySum)
	return(MySum)
}

n <- 20
Dom <- cbind(rep(0, 2 * n), rep(1, 2 * n))
Result <- rgenoud::genoud(fn = SumProfit, nvars = 2 * n, starting.values = runif(2 * n),
  Domains = Dom, boundary.enforcement = 2, MemoryMatrix = FALSE, print.level = 1)
}
}

\keyword{SwitchProfile}
