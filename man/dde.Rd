\name{dde}
\alias{dde}
\title{Solve Delay Differential Equations}

\description{
  A solver for systems of delay differential equations based off numerical routines from 
  Simon Wood's solv95 program. This solver is also capable of solving systems of ordinary 
  differential equations. 
  
  Please see the included demos for examples of how to use dde. To view available demos
  run \code{demo(package="ddesolve")}. The supplied demos require that the package
  PBSmodelling be installed.
}

\usage{
dde(y, func, parms=NULL, from=0, to=10, by=0.01, tol=1e-8, dt=0.1, hbsize=10000)
}

\arguments{
\item{y}{vector of initial values of the DDE system. The size of the supplied vector
    determines the number of variables in the system.}
\item{func}{a user supplied function that computes the gradients in the DDE system at time t. 
	The function must be defined either as: yprime <- func(t,y) or yprime <- func(t,y,parms). 
	where t is the current time point in the integration;
	y is a vector of the current estimated variables of the DDE system;
	and parms is any R object representing additional parameters (optional). 

	func must return one of the two following return types:
	1) a vector containing the calculated gradients for	each variable; 
	or 2) a list whose first element is a vector of calculated gradients, 
	and whose second element is a (possibly named) vector of values required 
	at each point in the integration.
},
\item{parms}{any parameters to pass to func}
\item{from}{start time}
\item{to}{stop time}
\item{by}{save results approximately every \'by\' timestep}
\item{tol}{maximum error tolerated at each timestep (as a proportion of the state variable concerned)}
\item{dt}{maximum initial timestep}
\item{hbsize}{history buffer size required for solving DDEs)}
}

\details{
  The user supplied function \code{func} can access past values (lags) of \code{y}
  by calling the \code{\link{pastvalue}} function. Past gradients are accessible
  by the \code{\link{pastgradient}} function. These functions can only be called from \code{func}
  and can only be passed values of \code{t} greater or equal to the start time, but less than the
  current time of the integration point. For example, calling \code{pastvalue(t)} is not allowed,
  since these values are the current values which are passed in as \code{"y"}.
}

\value{
  a data.frame with one column for t, a colomn for every variable in the system,
  and a column for every additional value that may (or may not) have been returned
  by \code{func} in the second element of the list.
  
  If the initial \"y\" values parameter was named, then the solved values column 
  will use the same names. Otherwise "y1", "y2", ... will be used.
  
  If \code{func} returned a list, with a named vector as the second element, then
  those names will be used as the column names. If the vector was not named, then
  "extra1", "extra2", ... will be used.
}

\seealso{ \code{\link{pastvalue}} }

\examples{
##################################################
# This is just a single example of using dde.
# For more examples see demo(package="ddesolve")
# the demos require the package PBSmodelling
##################################################

#create a func to return dde gradient
yprime <- function(t,y,parms) {
	if (t < parms$tau)
		lag <- parms$initial
	else
		lag <- pastvalue(t - parms$tau)
	y1 <- parms$a * y[1] - (y[1]^3/3) + parms$m * (lag[1] - y[1])
	y2 <- y[1] - y[2]
	return(c(y1,y2))
}

#define initial values and parameters
yinit <- c(1,1)
parms <- list(tau=3, a=2, m=-10, initial=yinit)

# solve the dde system
yout <- dde(y=yinit, func=yprime, parms=parms, from=0, to=30)

# and display the results
plot(yout$t, yout$y1, type="l", col="red", xlab="t", ylab="y", 
     ylim=c(min(yout$y1, yout$y2), max(yout$y1, yout$y2)))
lines(yout$t, yout$y2, col="blue")
legend("topleft", legend = c("y1", "y2"),lwd=2, lty = 1, 
       xjust = 1, yjust = 1, col = c("red","blue"))
}

\keyword{math}
