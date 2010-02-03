.First.lib <- function(libname, pkgname)
{
  library.dynam("ExtremalProc", package = pkgname, lib.loc = libname)
  invisible()
}


