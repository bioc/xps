"image.dev" <-
function(x,
         bg      = FALSE,
         transfo = log2,
         col     = gray((0:64)/64),
         names   = "namepart",
         xlab    = "",
         ylab    = "",
         mar     = c(1,1,2,1),
         dev     = "screen",
         outfile = "Image",
         w       = 540,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------image.dev------")

   ## check for correct class
   if (!is (x, "DataTreeSet")) {
      stop(paste(sQuote("x"), "is not  class", sQuote("DataTreeSet")));
   }#if

   ## get expression levels from data or bgrd
   if (bg == FALSE)          ds <- validData(x, which="")
   else                      ds <- validBgrd(x, which="");
   if (is.function(transfo)) ds <- transfo(ds);

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
      m <- as.numeric(ds[,i]);
      m <- matrix(m, ncol=ncols(x), nrow=nrows(x));
      m <- m[,nrows(x):1];

      image(m,
            col  = col,
            main = names[i],
            xlab = xlab,
            ylab = ylab,
            xaxt = 'n',
            yaxt = 'n',
            ...);
   }#for

   par(oldpar);
   if (dev != "screen") {
      dev.off();
   }#if
}#image.dev
