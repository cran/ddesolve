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
	#extract variables from GUI
	getWinVal(scope="L")
	
	#function to calculate gradient at a given time
	yprime <- function(t, y, parms=NULL) {
    	y1 <- sigma*(y[2]-y[1]) 
    	y2 <- y[1]*(tau-y[3]) - y[2]
    	y3 <- y[1]*y[2] - rho*y[3]
    	if (derivative=="yes")
    		return(list(c(y1,y2,y3), c(dy1=y1,dy2=y2,dy3=y3)))
    	else
    		return(list(c(y1,y2,y3), NULL))
	}
	#initial values
	yinit <- c(y1=y1,y2=y2,y3=y3)
	#solve ODE
	if (solver=="ddesolve") {
		x <- dde(y=yinit, func=yprime, from=t0, to=t1, by=timestep, hbsize=0)
	}
	else if (solver=="odesolve") {
		require(odesolve) || stop("odesolve is required")
		x <- lsoda(y=yinit, times=seq(from=t0, to=t1, by=timestep), func=yprime, parms=NULL, rtol=1e-6, atol=1e-4)
		x <- as.data.frame(x)
		if (derivative=="yes") #something weird happens to the labels with odesolve
			colnames(x)<-c("time", "y1", "y2", "y3", "dy1", "dy2", "dy3")
	}
	resetGraph()
	pairs(x, pch=183)
}

#function to restore working directory once demo is done
onClose <- function() { setwd(oldwd); }
createWin("demo_files/lorenz_win.txt")