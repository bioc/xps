#------------------------------------------------------------------------------#
# plotCorr: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotCorr" <-
function(x,
         which      = "UnitName",
         transfo    = log2,
         method     = "spearman",
         col        = NULL,
         names      = "namepart",
         sort       = FALSE,
         reverse    = TRUE,
         bmar       = NULL,
         add.legend = FALSE,
         dev        = "screen",
         outfile    = "CorrelationPlot",
         w          = 540,
         h          = 540,
         ...) 
{
   if (debug.xps()) print("------plotCorr------")

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
   corplot(x,
           which      = which,
           transfo    = transfo,
           method     = method,
           col        = col,
           names      = names,
           sort       = sort,
           reverse    = reverse,
           bmar       = bmar,
           add.legend = add.legend,
            ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotCorr

#------------------------------------------------------------------------------#
