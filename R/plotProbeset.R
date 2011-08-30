#------------------------------------------------------------------------------#
# plotProbeset: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotProbeset" <-
function(x,
         unitID,
         unittype   = "transcript",
         which      = "pm",
         transfo    = log2,
         names      = "namepart",
         ylim       = NULL,
         col        = 1:6,
         lty        = 1:5,
         add.legend = FALSE,
         dev        = "screen",
         outfile    = "ProbesetPlot",
         w          = 540,
         h          = 540,
         ...) 
{
   if (debug.xps()) print("------plotProbeset------")

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

   ## plot probeset
   probesetplot(x,
        unitID     = unitID,
        unittype   = unittype,
        which      = which,
        transfo    = transfo,
        names      = names,
        ylim       = ylim,
        col        = col,
        lty        = lty,
        add.legend = add.legend,
         ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotProbeset

#------------------------------------------------------------------------------#
