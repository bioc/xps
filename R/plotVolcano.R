#------------------------------------------------------------------------------#
# plotVolcano: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotVolcano" <-
function(x,
         labels      = "",
         p.value     = "pval",
         mask        = FALSE,
         show.cutoff = TRUE,
         cex.text    = 0.7,
         col.text    = "blue",
         col.cutoff  = "grey",
         xlim        = NULL,
         xlab        = "Log2(Fold-Change)",
         ylab        = "-Log10(P-Value)",
         pch         = '.',
         dev         = "screen",
         outfile     = "VolcanoPlot",
         w           = 540,
         h           = 540,
         ...) 
{
   if (debug.xps()) print("------plotVolcano------")

   ## check for correct class
   if (!is(x, "AnalysisTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("AnalysisTreeSet")));
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
   volcanoplot(x,
               labels      = labels,
               p.value     = p.value,
               mask        = mask,
               show.cutoff = show.cutoff,
               cex.text    = cex.text,
               col.text    = col.text,
               col.cutoff  = col.cutoff,
               xlim        = xlim,
               xlab        = xlab,
               ylab        = ylab,
               pch         = pch,
               ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotVolcano

#------------------------------------------------------------------------------#
