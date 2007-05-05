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
	getWinVal(scope="L")
	
	yprime <- function(t, y) {
	    if (t-t0 > tau)
	    	lag <- pastvalue(t-tau)
	    else
	    	lag <- 0
	    y1 <- P*lag[1]*exp(-lag[1]/theta)-delta*y[1]
	    return(list(y1, c(dy=y1, exp=exp(-lag[1]/theta))))
	}
	#we must also specify the initial value or conditions.
	yinit <- c(y=initPop)
	
	#solve the ODE from t0..t1
	x <- dde(y=yinit, func=yprime, from=t0, to=t1, hbsize=1000) 
	
	if (ptype=="t") {
		resetGraph()
		par(mfrow=c(3,1))
		plot(x=x$t, y=x$y, type="l", main="Adult Blowfly Population", xlab="Time", ylab="Population (y)")
		plot(x=x$t, y=x$dy, type="l", main="Rate of Change of Adult Population", xlab="Time", ylab="delta Population (dy)")
		plot(x=x$t, y=x$exp, type="l", main="exp(-lag[1]/theta)", xlab="Time", ylab="exp value")
	}
	else if (ptype=="p") {
		resetGraph()
		plot(x, main="Adult Blowfly Population")
	}
}

#restore working directory once demo is done
onClose <- function() { setwd(oldwd); }

createWin("demo_files/blowflies_win.txt")