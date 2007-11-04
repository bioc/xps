#===========================================================================#
# xps package initialization
#===========================================================================#

.onLoad <- function(libname, pkgname) {
  require(methods);
  cat(paste("\nWelcome to", pkgname, "version 0.3.10", "\n"));
  cat("    an R wrapper for XPS - eXpression Profiling System\n");
  cat("    (c) Copyright 2001-2007 by Christian Stratowa\n");
  cat("    \n");
}

.onUnload <- function( libpath ) {
  library.dynam.unload( "xps", libpath )
}
