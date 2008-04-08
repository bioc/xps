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

   ## unix vs windows settings
   is.win <- (.Platform$OS.type == "windows");
   if (is.win) {
      sq  <- as.character("");  #
      dq  <- as.character("\\\""); #\"
   } else {
      sq  <- as.character("'");  #'
      dq  <- as.character("\""); #"
   }#if

   ## get path to library xps.so
   xpsso <- system.file("libs", .Platform$r_arch,
                  paste("xps", .Platform$dynlib.ext, sep=''), package="xps");
   xpsso <- gsub("//", "/", xpsso);

   ## need to put character variables in parenthesis
#   temp       <- as.character("\""); #"
   xpsso      <- paste(dq, xpsso, dq, sep="");
   rootfile   <- paste(dq, rootfile, dq, sep="");
   canvasname <- paste(dq, canvasname, dq, sep="");
   treename   <- paste(dq, treename, dq, sep="");
   varlist    <- paste(dq, varlist, dq, sep="");
   logbase    <- paste(dq, logbase, dq, sep="");
   option     <- paste(dq, option, dq, sep="");
   w          <- as.integer(w);
   h          <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDrawImage.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",",
                         treename,",", varlist,",", logbase,",", option,",",
                         w,",", h, ")", sq, sep="");
   system(paste("root -l", macro));

}#root.image
