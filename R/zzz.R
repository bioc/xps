#===========================================================================#
# xps package initialization
#===========================================================================#

#.onLoad <- function(libname, pkgname) {
#   require(methods);
#   require(utils);
#}

.onAttach <- function(libname, pkgname) {
   packageStartupMessage(paste("\nWelcome to", pkgname, "version", packageDescription(pkgname, fields="Version")));
   packageStartupMessage("    an R wrapper for XPS - eXpression Profiling System");
   packageStartupMessage("    (c) Copyright 2001-2014 by Christian Stratowa");
   packageStartupMessage("    ");
}

.onUnload <- function(libpath) {
   library.dynam.unload( "xps", libpath )
}
