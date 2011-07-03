#------------------------------------------------------------------------------#
# plotDensity: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotDensity" <-
function(x,
         which      = "",
         size       = 0,
         transfo    = log2,
         ylab       = "density",
         xlab       = "log intensity",
         names      = "namepart",
         type       = "l",
         col        = 1:6,
         lty        = 1:5,
         add.legend = FALSE,
         dev        = "screen",
         outfile    = "DensityPlot",
         w          = 540,
         h          = 540,
         verbose    = TRUE,
         ...) 
{
   if (debug.xps()) print("------plotDensity------")

   ## check for correct class
   if (!extends(class(x), "ProcesSet")) {
      stop(paste(sQuote("x"), "is not derived from class", sQuote("ProcesSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(file=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(file=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   hist(x,
        which      = which,
        size       = size,
        transfo    = transfo,
        ylab       = ylab,
        xlab       = xlab,
        names      = names,
        type       = type,
        col        = col,
        lty        = lty,
        add.legend = add.legend,
        verbose    = verbose,
         ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotDensity

#------------------------------------------------------------------------------#
