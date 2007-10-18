#==============================================================================#
# normalize.R: normalization functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# normalize:           methods: mean, median, quantile, lowess, supsmu
# normalize.constant:  mean or median normalization
# normalize.quantiles: RMA-tpye quantile normalization
# normxpress.lowess:   lowess normalization
# normxpress.supsmu:   supsmu normalization
#==============================================================================#

"normalize" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "all",
         method    = "mean",
         option    = "all",
         logbase   = "0",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(0.02, 0),
         verbose   = TRUE)
{
   if (is(xps.data, "DataTreeSet") || is(xps.data, "ExprTreeSet")) {
      set <- xpsNormalize(xps.data,
                          filename  = filename,
                          filedir   = filedir,
                          tmpdir    = tmpdir,
                          update    = update,
                          select    = select,
                          method    = method,
                          option    = option,
                          logbase   = logbase,
                          exonlevel = exonlevel,
                          refindex  = refindex,
                          refmethod = refmethod,
                          params    = params,
                          verbose   = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class <DataTreeSet,ExprTreeSet>"));
   }#if
}#normalize

#------------------------------------------------------------------------------#

"normalize.constant" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         method    = "mean",
         logbase   = "0",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(0.02, 0),
         verbose   = TRUE)
{
   set <- normalize(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = "all",
                    method    = method,
                    option    = "all",
                    logbase   = logbase,
                    exonlevel = exonlevel,
                    refindex  = refindex,
                    refmethod = refmethod,
                    params    = params,
                    verbose   = verbose);
   return(set);
}#normalize.constant

#------------------------------------------------------------------------------#

"normalize.quantiles" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         exonlevel = "",
         verbose   = TRUE)
{
   set <- normalize(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = "pmonly",
                    method    = "quantile",
                    option    = "together:none",
                    logbase   = "0",
                    exonlevel = exonlevel,
                    refindex  = 1,
                    refmethod = "mean",
                    params    = list(0.0),
                    verbose   = verbose);
   return(set);
}#normalize.quantiles

#------------------------------------------------------------------------------#

"normalize.lowess" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         logbase   = "log2",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(0.67, 3),
         verbose   = TRUE)
{
   set <- normalize(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = "all",
                    method    = "lowess",
                    option    = "all",
                    logbase   = logbase,
                    exonlevel = exonlevel,
                    refindex  = refindex,
                    refmethod = refmethod,
                    params    = params,
                    verbose   = verbose);
   return(set);
}#normalize.lowess

#------------------------------------------------------------------------------#

"normalize.supsmu" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         logbase   = "log2",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(0.67, 3),
         verbose   = TRUE)
{
   set <- normalize(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = "all",
                    method    = "supsmu",
                    option    = "all",
                    logbase   = logbase,
                    exonlevel = exonlevel,
                    refindex  = refindex,
                    refmethod = refmethod,
                    params    = params,
                    verbose   = verbose);
   return(set);
}#normalize.supsmu

