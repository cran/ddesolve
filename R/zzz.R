
.First.lib <- function(lib, pkg)
{
	library.dynam("ddesolve", pkg, lib)
	cat("
ddesolve 1.00 -- based on solv95 by Simon Wood

A complete User Guide appears as ddesolve-UG.pdf in the root library
directory of ddesolve. To use this package effectively, please consult 
the guide.

Last built on Tue Jul 17, 2007
")
}
