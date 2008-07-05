#===========================================================================#
# xps package initialization
#===========================================================================#

.onLoad <- function(libname, pkgname) {
  require(methods);
  cat(paste("\nWelcome to", pkgname, "version 1.1.3", "\n"));
  cat("    an R wrapper for XPS - eXpression Profiling System\n");
  cat("    (c) Copyright 2001-2008 by Christian Stratowa\n");
  cat("    \n");
}

.onUnload <- function( libpath ) {
  library.dynam.unload( "xps", libpath )
}
