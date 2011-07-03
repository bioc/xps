#------------------------------------------------------------------------------#
# plotPM: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotPM" <-
function(x,
         which   = "",
         size    = 0,
         transfo = NULL,
         method  = mean,
         names   = "namepart",
         beside  = TRUE,
         col     = c("red", "blue"),
         legend  = c("PM","MM"),
         las     = 2,
         ylab    = "mean intensities",
         dev     = "screen",
         outfile = "PMPlot",
         w       = 540,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------plotPM------")

   ## check for correct class
   if (!is(x, "DataTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("DataTreeSet")));
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
   pmplot(x,
          which   = which,
          size    = size,
          transfo = transfo,
          method  = method,
          names   = names,
          beside  = beside,
          col     = col,
          legend  = legend,
          las     = las,
          ylab    = ylab,
          ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotPM

#------------------------------------------------------------------------------#
