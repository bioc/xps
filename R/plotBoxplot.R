#------------------------------------------------------------------------------#
# plotBoxplot: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotBoxplot" <-
function(x,
         which   = "",
         size    = 0,
         transfo = log2,
         range   = 0,
         names   = "namepart",
         mar     = c(10,5,2,1),
         las     = 2,
         cex     = 1.0,
         dev     = "screen",
         outfile = "BoxPlot",
         w       = 800,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------plotBoxplot------")

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

   if (is.null(mar)) {
      mar  <- c(5, 4, 4, 2) + 0.1;
      bmar <- NULL;
   } else {
      bmar <- list(b=mar[1], cex=cex, w=w);
   }#if

   oldpar <- par(no.readonly=TRUE, mar=mar);

   ## plot data
   boxplot(x,
           which   = which,
           size    = size,
           transfo = transfo,
           range   = range,
           names   = names,
           bmar    = bmar,
           las     = las,
           ...)

   par(oldpar);
   if (dev != "screen") {
      dev.off();
   }#if
}#plotBoxplot

#------------------------------------------------------------------------------#
