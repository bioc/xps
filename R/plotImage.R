#------------------------------------------------------------------------------#
# plotImage: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"plotImage" <-
function(x,
         type    = character(),
         qualopt = c("raw", "adjusted", "normalized"),
         transfo = log2,
         col     = NULL,
         names   = character(),
         dev     = "screen",
         outfile = "Image",
         w       = 800,
         h       = 800,
         ...) 
{
   if (debug.xps()) print("------plotImage------");

   qualopt <- match.arg(qualopt);
   if (names[1] == "namepart") {
      stop("<names='namepart'> is not possible for this function");
   }#if

   ## check for correct class
   if (is (x, "DataTreeSet")) {
      if (is.na(match(type, c("intensity", "background")))) {
         stop(paste(sQuote(type), "is not a valid residual option"));
      }#if

      if (names[1] == "*") names <- unlist(treeNames(x));

      bg <- FALSE;
   } else if (is (x, "ExprTreeSet")) {
      if (is.na(match(type, c("intensity", "background")))) {
         stop(paste(sQuote(type), "is not a valid data option"));
      }#if

      if (type == "intensity") {
         if (names[1] == "*") names <- getTreeNames(rootFile(x), treetype=ADJTYPE[1]);

         bg <- FALSE;
      } else {
         if (names[1] == "*") {
            exten <- unique(extenPart(getTreeNames(rootFile(x))));
            id    <- match(exten, BGDTYPE);
            id    <- id[!is.na(id)];
            names <- getTreeNames(rootFile(x), treetype=BGDTYPE[id]);
         }#if

         bg <- TRUE;
      }#if
   } else if (is (x, "QualTreeSet")) {
      if (is.na(match(type, RESIOPT))) {
         stop(paste(sQuote(type), "is not a valid residual option"));
      }#if

      if (names[1] == "*") names <- getTreeNames(rootFile(x), treetype=QUATYPE[4]);

      names <- names[grep(qualopt, names)];
   } else {
      stop(paste(sQuote("x"), "must be one of class <DataTreeSet, QualTreeSet>"));
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
   oldpar <- par(no.readonly=TRUE);

   ## layout for images
   size <- function(x){c(ceiling(sqrt(x)), round(sqrt(x)))}
   nr   <- size(length(names))[1];
   nc   <- size(length(names))[2];
   nfig <- layout(matrix(c(1:(nr*nc)), nr, nc, byrow=TRUE), respect=TRUE);

   ## images
   layout.show(nfig);
   for (i in 1:length(names)) {
      if (is (x, "DataTreeSet") | is(x, "ExprTreeSet")) {
         image(x,
               bg         = bg,
               transfo    = transfo,
               col        = col,
               names      = names[i],
               xlab       = "",
               ylab       = "",
               add.legend = FALSE,
               ...);
      } else if (is (x, "QualTreeSet")) {
         image(x,
               type       = type,
               qualopt    = qualopt,
               transfo    = transfo,
               col        = col,
               names      = names[i],
               xlab       = "",
               ylab       = "",
               add.legend = FALSE,
               ...);
      }#if
   }#for

   par(oldpar);
   if (dev != "screen") {
      dev.off();
   }#if
}#plotImage

#------------------------------------------------------------------------------#
