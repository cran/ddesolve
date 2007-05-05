# Newton's Law of cooling
# models the cooling of a cup of coffee left in a room
# see http://en.wikipedia.org/wiki/Heat_conduction#Newton.27s_law_of_cooling
require(ddesolve)
if (!require(PBSmodelling)) stop("The package PBSmodelling must be installed for this demo")

#close any existing windows
closeWin("window")

#store working directory
if (!exists("oldwd") || getwd() != system.file(package = "ddesolve")) {
	oldwd <- getwd()
	setwd(system.file(package = "ddesolve"))
}

runPlot <- function()
{
	#extract variables from window
	getWinVal(scope="L")
	
	# `y' is the estimated value of the variable at time `t'
	# y[1] represents the tempurature of a cup of coffee at a given time `t'.
	# dy/dt ~ y - room tempurature. 
	# so we can write dy/dt = -rho(y - roomTemp)
	yprime <- function(t, y) {
	    y1 <- -rho*(y[1]-Tenv)
	    return(y1)
	}
	#we must also specify the initial value or conditions.
	yinit <- c(Tcup)
	
	#solve the ODE from t0..t1 - (ignore hbsize, setting to zero may crash)
	x <- dde(y=yinit, func=yprime, from=t0, to=t1, hbsize=0) 
	resetGraph()
	plot(x, type="l", main="Cooling of a cup of coffee", ylab="Tempurature", xlab="Time")
}

#restore working directory once demo is done
onClose <- function() { setwd(oldwd); }

createWin("demo_files/cooling_win.txt")