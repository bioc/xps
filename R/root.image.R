"root.image" <-
function(x,
         treename   = character(0),
         leafname   = "fInten",
         logbase    = "log2",
         option     = "COLZ",
         canvasname = "Image",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.image------")

   ## check for correct class
   if (!is (x, "DataTreeSet")) {
      stop(paste(sQuote("x"), "is not  class", sQuote("DataTreeSet")));
   }#if

   rootfile <- rootFile(x);
   varlist  <- paste("fX", "fY", leafname, sep=":");

   ## test for valid treename
   treenames <- getTreeNames(rootfile);
   if (treename%in%treenames == FALSE) {
      stop(paste("treename", sQuote(treename), "is not present in", sQuote(rootfile)));
   }#if
   ## need to add ROOT directory in rootfile
   treename <- paste("DataSet", treename, sep="/");

   ## get path to library xps.so
   xpsso <- system.file("libs", .Platform$r_arch,
                  paste("xps", .Platform$dynlib.ext, sep=''), package="xps");
   xpsso <- gsub("//", "/", xpsso);

   ## need to put character variables in parenthesis
   temp       <- as.character("\""); #"
   xpsso      <- paste(temp, xpsso, temp, sep="");
   rootfile   <- paste(temp, rootfile, temp, sep="");
   canvasname <- paste(temp, canvasname, temp, sep="");
   treename   <- paste(temp, treename, temp, sep="");
   varlist    <- paste(temp, varlist, temp, sep="");
   logbase    <- paste(temp, logbase, temp, sep="");
   option     <- paste(temp, option, temp, sep="");
   w          <- as.integer(w);
   h          <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDrawImage.C", sep="/");
   macro <- paste("'", macro, "(", xpsso,",", rootfile,",", canvasname,",",
                         treename,",", varlist,",", logbase,",", option,",",
                         w,",", h, ")", "'", sep="");
   system(paste("root -l", macro));

}#root.image
