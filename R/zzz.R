#===========================================================================#
# xps package initialization
#===========================================================================#

.onLoad <- function(libname, pkgname) {
   require(methods);
   require(utils);
   packageStartupMessage(paste("\nWelcome to", pkgname, "version", packageDescription(pkgname, field="Version")));
   packageStartupMessage("    an R wrapper for XPS - eXpression Profiling System");
   packageStartupMessage("    (c) Copyright 2001-2012 by Christian Stratowa");
   packageStartupMessage("    ");
}

.onUnload <- function(libpath) {
   library.dynam.unload( "xps", libpath )
}
