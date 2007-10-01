
.First.lib <- function(lib, pkg)
{
	library.dynam("ddesolve", pkg, lib)
	cat("
ddesolve 1.01 -- based on solv95 by Simon Wood

A complete User Guide appears as ddesolve-UG.pdf in the root library
directory of ddesolve. To use this package effectively, please consult 
the guide.

Last built on Mon Oct 1, 2007 (23:15:27)
")
}
