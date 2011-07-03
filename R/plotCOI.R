#------------------------------------------------------------------------------#
# plotCOI: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotCOI" <-
function(x,
         type    = c("pos", "neg"),
         qualopt = "raw",
         radius  = 0.5,
         linecol = "gray70",
         visible = TRUE,
         dev     = "screen",
         outfile = "CenterOfIntensityPlot",
         w       = 540,
         h       = 540,
         ...)
{
   if (debug.xps()) print("------plotCOI------")

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
   coiplot(x,
           type    = type,
           qualopt = qualopt,
           radius  = radius,
           linecol = linecol,
           visible = visible,
           ...) 

   if (dev != "screen") {
      dev.off();
   }#if
}#plotCOI

#------------------------------------------------------------------------------#
