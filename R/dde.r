
dde <- function (y, func, parms=NULL, from=0, to=10, by=0.01, tol=1e-8, dt=0.1, hbsize=10000) 
{
	if (hbsize<=0) hbsize=1 #if 0 or less, C code crashes
    out <- .Call("startDDE", 
                 gradFunc=func, 
                 env=new.env(), 
                 yinit=y,
                 parms=parms,
                 settings=c(tol, from, to, dt, by, hbsize),
                 PACKAGE = "ddesolve")

    return(data.frame(out))
}

# t - at what time
# markno - used for optimization
pastvalue <- function(t)
{
	markno=0 #used for optimization when more than one lag time is used
	         #but will need code to update `data.nlag' in me95.c:setupglobaldata
    out <- .Call("getPastValue", 
                 t=t, 
                 markno=as.integer(markno),
                 PACKAGE = "ddesolve")
	return(out)
}

# t - at what time
# markno - used for optimization
pastgradient <- function(t)
{
	markno=0 #used for optimization when more than one lag time is used
	         #but will need code to update `data.nlag' in me95.c:setupglobaldata
    out <- .Call("getPastGradient", 
                 t=t, 
                 markno=as.integer(markno),
                 PACKAGE = "ddesolve")
	return(out)
}
