#------------------------------------------------------------------------------#
# plotMAD: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotMAD" <-
function(x,
         which      = "UnitName",
         transfo    = log2,
         col        = NULL,
         names      = "namepart",
         sort       = FALSE,
         bmar       = NULL,
         add.legend = FALSE,
         dev        = "screen",
         outfile    = "MADPlot",
         w          = 540,
         h          = 540,
         ...) 
{
   if (debug.xps()) print("------plotMAD------")

   ## check for correct class
   if (!is(x, "ExprTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("ExprTreeSet")));
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
   madplot(x,
           which      = which,
           transfo    = transfo,
           col        = col,
           names      = names,
           sort       = sort,
           bmar       = bmar,
           add.legend = add.legend,
           ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotMAD

#------------------------------------------------------------------------------#
