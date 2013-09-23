#Description: This is a script that performs the maximum likelihood analysis to calculate the slope and intercepts.
require("maxLik") # library required for maximum likelihood estimation

# definition of lambda as a function of r
lambda <- function( r, slope, int ) { int + slope*r }

# calculates the likelihood function. Assumes data are stored in global variable "x",
# in the right format (as used by the function "simulate.data")
logLikFun <- function( param )
{
	slope <- param[1]
	int <- param[2]
	logLik <- 0
	for ( i in 1:nrow(x) )
	{
		prob <- exp(-lambda(x[i,2], slope, int)*(0:19))
		logLik <- logLik + dmultinom(unlist(x[i,3:22]), prob=prob, log = TRUE)
	}
	#cat( "log-Likelihood:", logLik, "\n" )
	logLik
}

fit.data <- function( filename )
{
	x <<- read.table( filename ) # assign to global variable x
	#print( x )
	mle <- maxLik( logLik = logLikFun, start = c(slope=0, int=1) )
	print( summary(mle) )

	return(mle)
}

filename = commandArgs(trailingOnly=TRUE)

#print(filename)

# fit model to data using maximum likelihood
mle = fit.data( filename )

print(mle$estimate)
