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
	
	myGrad <- function(t, y) {
		if (t-t0 >= tau) ylag <- pastvalue(t-tau)
		else ylag <- 0
		yexp <- exp(-theta*ylag[1]/A0)
		yp <- P*ylag[1]*yexp-delta*y[1]
		return( list(yp, c(dy=yp, yexp=yexp)) ) }
	
	#solve the ODE from t0..t1
	x <- dde(y=A0, func=myGrad, times=seq(t0,t1), hbsize=1000) 
	
	if (ptype=="t") {
		frame(); resetGraph();
		expandGraph(mfrow=c(3,1),mar=c(4,4,2,1),mgp=c(2.75,.75,0),cex.main=1.5,cex.lab=1.5)
		plot(x=x$t, y=x$y1, type="l", main="Adult Blowfly Population", xlab="Time", ylab="Population (y)")
		plot(x=x$t, y=x$dy, type="l", main="Rate of Change of Adult Population", xlab="Time", ylab="delta Population (dy)")
		plot(x=x$t, y=x$yexp, type="l", main="exp(-theta*ylag[1] / A0)", xlab="Time", ylab="exp value")
	}
	else if (ptype=="p") {
		frame(); resetGraph();
		plot(x, main="Adult Blowfly Population") }
}

#restore working directory once demo is done
onClose <- function() { setwd(oldwd); }

createWin("demo_files/blowflies_win.txt")