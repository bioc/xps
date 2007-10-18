"mvaplot.dev" <-
function(x,
         transfo = log2,
         method  = "median",
         names   = "namepart",
         ylim    = c(-6,6),
         xlab    = "A",
         ylab    = "M",
         pch     = '.',
         mar     = c(3,3,2,1),
         dev     = "screen",
         outfile = "MvAPlot",
         w       = 540,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------mvaplot.dev------")

   ## check for correct class
   if (!is (x, "ExprTreeSet")) {
      stop(paste(sQuote("x"), "is not  class", sQuote("ExprTreeSet")));
   }#if

   ## get expression levels from data
   ds <- validData(x);

   if (is.function(transfo)) ds <- transfo(ds);
   if (method == "median")   mn <- apply(ds, 1, median)
   else                      mn <- rowMeans(ds);

   if (is.null(names))              names <- colnames(ds)
   else if (names[1] == "namepart") names <- namePart(colnames(ds))
   else                             ds    <- ds[, names, drop=F];

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
   oldpar <- par(no.readonly=TRUE);

   size <- function(x){c(ceiling(sqrt(x)), round(sqrt(x)))}
   nc <-size(ncol(ds))[1];
   nr <-size(ncol(ds))[2];
#   nc <- nr <- ceiling(sqrt(ncol(ds)));
   nfig   <- layout(matrix(c(1:(nr*nc)), nr, nc, byrow=TRUE), respect=TRUE);

   ## images
   layout.show(nfig);
   for (i in 1:ncol(ds)) {
      plot(x    = mn,
           y    = (ds[,i] - mn),
           main = names[i],
           ylim = ylim,
           xlab = xlab,
           ylab = ylab,
           pch  = pch,
           ...);
      abline(h=0, col=gray(0.5));
   }#for

   par(oldpar);
   if (dev != "screen") {
      dev.off();
   }#if
}#mvaplot.dev
