#------------------------------------------------------------------------------#
# plotNUSE: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotNUSE" <-
function(x,
         which    = "UnitName",
         size     = 0,
         range    = 0,
         names    = "namepart",
         main     = "NUSE Plot",
         ylim     = c(0.8,1.2),
         las      = 2,
         add.line = TRUE,
         outline  = FALSE,
         dev      = "screen",
         outfile  = "NUSEPlot",
         w        = 800,
         h        = 540,
         ...) 
{
   if (debug.xps()) print("------plotNUSE------")

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
   if (is (x, "ExprTreeSet")) {
      nuseplot(x,
               which    = which,
               size     = size,
               range    = range,
               names    = names,
               main     = main,
               ylim     = ylim,
               las      = las,
               add.line = add.line,
               outline  = outline,
               ...);
   } else if (is (x, "QualTreeSet")) {
      nuseplot(x,
               range    = range,
               names    = names,
               main     = main,
               ylim     = ylim,
               las      = las,
               add.line = add.line,
               ...);
   }#if

   if (dev != "screen") {
      dev.off();
   }#if
}#plotNUSE

#------------------------------------------------------------------------------#
