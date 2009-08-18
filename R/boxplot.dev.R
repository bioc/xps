"boxplot.dev" <-
function(x,
         which   = "",
         size    = 0,
         transfo = log2,
         range   = 0,
         names   = "namepart",
         mar     = c(10,5,4,1),
         las     = 2,
         dev     = "screen",
         outfile = "BoxPlot",
         w       = 800,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------boxplot.dev------")

   ## check for correct class
   if (!extends(class(x), "ProcesSet")) {
      stop(paste(sQuote("x"), "is not derived from class", sQuote("ProcesSet")));
   }#if

   ## get expression levels from data
   m <- validData(x, which=which);
   if (size > 1)             m <- m[seq(1,nrow(m),len=size),];
   if (is.function(transfo)) m <- transfo(m);

   if (is.null(names))              names <- colnames(m)
   else if (names[1] == "namepart") names <- namePart(colnames(m))
   else                             m     <- m[, names, drop=F];

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
   oldpar <- par(no.readonly=TRUE, mar=mar);

   ## plot data
   boxplot(data.frame(m),
           range = range,
           names = names,
           las   = las,
           ...)

   par(oldpar);
   if (dev != "screen") {
      dev.off();
   }#if
}#boxplot.dev
