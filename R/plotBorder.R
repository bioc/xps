#------------------------------------------------------------------------------#
# plotBorder: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotBorder" <-
function(x,
         type    = c("pos", "neg"),
         qualopt = "raw",
         transfo = log2,
         range   = 0,
         names   = "namepart",
         ylim    = NULL,
         bmar    = NULL,
         las     = 2,
         dev     = "screen",
         outfile = "BorderPlot",
         w       = 800,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------plotBorder------")

   ## check for correct class
   if (!is(x, "QualTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("QualTreeSet")));
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
   borderplot(x,
              type    = type,
              qualopt = qualopt,
              transfo = transfo,
              range   = range,
              names   = names,
              ylim    = ylim,
              bmar    = bmar,
              las     = las,
              ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotBorder

#------------------------------------------------------------------------------#
