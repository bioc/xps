"root.density" <-
function(x,
         treename   = "*",
         logbase    = "log2",
         canvasname = "DensityPlot",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.density------")

   rootfile  <- rootFile(x);
   treenames <- treeNames(x);

   ## get tree extension
   exten <- unlist(treenames)[1];
   exten <- unlist(strsplit(exten, "\\."))[2];

   ## test for valid treename
   namepart <- unlist(strsplit(treename, "\\."))[1];
   if (namepart == "*") {
      treename <- paste(namepart, exten, sep=".");
   } else if (treename%in%treenames == FALSE) {
      stop(paste("treename", sQuote(treename), "is not present in", sQuote(rootfile)));
   }#if

   ## add ROOT directory in rootfile and select leafname
   treename <- paste(setName(x), treename, sep="/");
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else {
      stop(paste(sQuote("x"), "must be  class <DataTreeSet,ExprTreeSet>"));
   }#if

   varlist <- leafname;
   option  <- "L";

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
   macro <- paste(system.file(package="xps"), "rootsrc/macroDrawDensity.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",",
                  treename,",", varlist,",", logbase,",", option,",",
                  w,",", h, ")", sq, sep="");
   system(paste("root -l", macro));
}#root.density
