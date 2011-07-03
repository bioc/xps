#------------------------------------------------------------------------------#
# plotCall: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotCall" <-
function(x,
         beside  = TRUE,
         names   = "namepart",
         col     = c("red","green","blue"),
         legend  = c("P","M","A"),
         ylim    = c(0,100),
         ylab    = "detection call [%]",
         las     = 2,
         dev     = "screen",
         outfile = "CallPlot",
         w       = 800,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------plotCall------")

   ## check for correct class
   if (!is(x, "CallTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("CallTreeSet")));
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
   callplot(x,
            beside = beside,
            names  = names,
            col    = col,
            legend = legend,
            ylim   = ylim,
            ylab   = ylab,
            las    = las,
            ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotCall

#------------------------------------------------------------------------------#
