"root.profile" <-
function(x,
         treename    = "*",
         varlist     = NULL,
         as.log      = TRUE,
         globalscale = TRUE,
         boxes       = TRUE,
         ylim        = NULL,
         canvasname  = "ProfilePlot",
         save.as     = "",
         w           = 800,
         h           = 600) 
{
   if (debug.xps()) print("------root.profile------")

   rootfile  <- rootFile(x);
   treenames <- treeNames(x);
   savename  <- rootDrawName(canvasname, save.as);

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

   ## add ROOT directory in rootfile
   treename <- paste(setName(x), treename, sep="/");

   ## select leafname as varlist
   if (is.null(varlist)) {
      if (is (x, "DataTreeSet")) {
         varlist <- "fInten";
      } else if (is (x, "ExprTreeSet")) {
         varlist <- "fLevel";
      } else {
         stop(paste("for <varlist=NULL>", sQuote("x"), 
                    "must be  class <DataTreeSet,ExprTreeSet>"));
      }#if
   }#if

   ## test for ylim
   if (is.null(ylim)) {
      miny <- 0.0;
      maxy <- 0.0;
   } else {
      miny <- ylim[1];
      maxy <- ylim[2];
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
   xpsso       <- paste(dq, xpsso, dq, sep="");
   rootfile    <- paste(dq, rootfile, dq, sep="");
   canvasname  <- paste(dq, canvasname, dq, sep="");
   treename    <- paste(dq, treename, dq, sep="");
   varlist     <- paste(dq, varlist, dq, sep="");
   savename    <- paste(dq, savename, dq, sep="");
   miny        <- as.double(miny);
   maxy        <- as.double(maxy);
   aslog       <- as.integer(as.log);
   globalscale <- as.integer(globalscale);
   boxes       <- as.integer(boxes);
   w           <- as.integer(w);
   h           <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDrawProfilePlot.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",",
                  treename,",", varlist,",", savename,",", miny,",", maxy,",", aslog,",",
                  globalscale,",", boxes,",", w,",", h, ")", sq, sep="");

   if (save.as == "") {
      system(paste("root -l", macro));
   } else {
      system(paste("root -l -q", macro));
   }#if
}#root.profile
