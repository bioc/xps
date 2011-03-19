"root.image" <-
function(x,
         treename   = character(0),
         leafname   = "fInten",
         logbase    = "log2",
         option     = "COLZ",
         zlim       = NULL,
         canvasname = "Image",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.image------")

   ## check for correct class
   if (!(is(x, "DataTreeSet") | is(x, "QualTreeSet"))) {
      stop(paste(sQuote("x"), "must be one of class <DataTreeSet, QualTreeSet>"));
   }#if

   rootfile <- rootFile(x);
   varlist  <- paste("fX", "fY", leafname, sep=":");
   savename <- rootDrawName(canvasname, save.as);

   ## test for valid treename
   treenames <- getTreeNames(rootfile);
   if (treename%in%treenames == FALSE) {
      stop(paste("treename", sQuote(treename), "is not present in", sQuote(rootfile)));
   }#if

### to do ###
# test for extenPart suitable for images, ev validImageExtension()
#############

   ## need to add ROOT directory in rootfile
   treename <- paste(setName(x), treename, sep="/");

   ## test for zlim
   if (is.null(zlim)) {
      minz <- 0.0;
      maxz <- 0.0;
   } else {
      minz <- zlim[1];
      maxz <- zlim[2];
   }#if

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
   xpsso      <- paste(dq, xpsso, dq, sep="");
   rootfile   <- paste(dq, rootfile, dq, sep="");
   canvasname <- paste(dq, canvasname, dq, sep="");
   treename   <- paste(dq, treename, dq, sep="");
   varlist    <- paste(dq, varlist, dq, sep="");
   logbase    <- paste(dq, logbase, dq, sep="");
   option     <- paste(dq, option, dq, sep="");
   savename   <- paste(dq, savename, dq, sep="");
   minz       <- as.double(minz);
   maxz       <- as.double(maxz);
   w          <- as.integer(w);
   h          <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDrawImage.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",",
                  treename,",", varlist,",", logbase,",", option,",",
                  savename,",", minz,",", maxz,",", w,",", h, ")", sq, sep="");

   if (save.as == "") {
      system(paste("root -l", macro));
   } else {
      system(paste("root -l -q", macro));
   }#if
}#root.image
