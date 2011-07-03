#------------------------------------------------------------------------------#
# plotPCA: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotPCA" <-
function(x,
         which      = "UnitName",
         transfo    = log2,
         method     = "none",
         groups     = NULL,
         screeplot  = FALSE,
         squarepca  = FALSE,
         pcs        = c(1,2),
         add.labels = FALSE,
         add.legend = FALSE,
         col        = NULL,
         names      = "namepart",
         as.list    = FALSE,
         dev        = "screen",
         outfile    = "PCAPlot",
         w          = 540,
         h          = 540,
         ...) 
{
   if (debug.xps()) print("------plotPCA------")

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
   pcaplot(x,
           which      = which,
           transfo    = transfo,
           method     = method,
           groups     = groups,
           screeplot  = screeplot,
           squarepca  = squarepca,
           pcs        = pcs,
           add.labels = add.labels,
           add.legend = add.legend,
           col        = col,
           names      = names,
           as.list    = as.list,
           ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotPCA

#------------------------------------------------------------------------------#
