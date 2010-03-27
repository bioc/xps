#===========================================================================#
# xps package initialization
#===========================================================================#

.onLoad <- function(libname, pkgname) {
   require(methods);
   require(utils);
   cat(paste("\nWelcome to", pkgname, "version", packageDescription(pkgname, field="Version"), "\n"));
   cat("    an R wrapper for XPS - eXpression Profiling System\n");
   cat("    (c) Copyright 2001-2010 by Christian Stratowa\n");
   cat("    \n");
}

.onUnload <- function(libpath) {
   library.dynam.unload( "xps", libpath )
}
