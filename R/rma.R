#==============================================================================#
# rma.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# rma: RMA preprocessing
#==============================================================================#

"rma" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         background = "pmonly",
         normalize  = TRUE,
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         verbose    = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {
      set <- xpsRMA(xps.data,
                    filename   = filename,
                    filedir    = filedir,
                    tmpdir     = tmpdir,
                    background = background,
                    normalize  = normalize,
                    option     = option,
                    exonlevel  = exonlevel,
                    xps.scheme = xps.scheme,
                    verbose    = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#rma
