
.First.lib <- function(lib, pkg)
{
	library.dynam("ddesolve", pkg, lib)
	cat("
ddesolve 1.02 -- based on solv95 by Simon Wood

A complete User Guide appears as ddesolve-UG.pdf in the root/doc 
directory of ddesolve. To use this package effectively, please consult 
the guide.

Last built on Mon April 16, 2007 (19:49:42)
")
}
