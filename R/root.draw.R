#==============================================================================#
# root.draw.R: plotting functions based on ROOT TCanvas:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# root.graph1D: 
# root.graph2D: 
# root.mvaplot: 
# root.hist1D: 
# root.hist2D: 
# root.hist3D: 
# root.draw: 
#==============================================================================#

"root.graph1D" <-
function(x,
         treename   = character(0),
         logbase    = "log2",
         option     = "P",
         canvasname = "Graph1D",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.graph1D------")

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
      leafname <- "fPValue";
   }#if

   root.draw(x,
             treename   = treename,
             leafname   = leafname,
             logbase    = logbase,
             type       = "graph",
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.graph1D

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.graph2D" <-
function(x,
         treename1  = character(0),
         treename2  = character(0),
         logbase    = "log2",
         option     = "P",
         canvasname = "Graph2D",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.graph2D------")

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
      leafname <- "fPValue";
   }#if

   root.draw(x,
             treename   = c(treename1, treename2),
             leafname   = leafname,
             logbase    = logbase,
             type       = "graph",
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.graph2D

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.mvaplot" <-
function(x,
         treename1  = character(0),
         treename2  = character(0),
         logbase    = "log2",
         option     = "P",
         canvasname = "MvAPlot",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.mvaplot------")

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
#      leafname <- "fPValue";
      stop(paste(sQuote("mvaplot"), "does not make sense for class", sQuote("CallTreeSet")));
   }#if

   root.draw(x,
             treename   = c(treename1, treename2),
             leafname   = leafname,
             logbase    = logbase,
             type       = "mvaplot",
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.mvaplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.hist1D" <-
function(x,
         treename   = character(0),
         logbase    = "log2",
         type       = "hist",
         option     = "HIST",
         canvasname = "Histogram1D",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.hist1D------")

   if (!(identical(type, "hist") || identical(type, "density"))) {
      stop(paste(sQuote(type), "is not a valid type option"));
   }#if

   ## density requires special option
   if (identical(type, "density")) {
      option = "AL";
   }#if

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
      leafname <- "fPValue";
   }#if

   root.draw(x,
             treename   = treename,
             leafname   = leafname,
             logbase    = logbase,
             type       = type,
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.hist1D

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.hist2D" <-
function(x,
         treename1  = character(0),
         treename2  = character(0),
         logbase    = "log2",
         option     = "COLZ",
         canvasname = "Histogram2D",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.hist2D------")

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
      leafname <- "fPValue";
   }#if

   root.draw(x,
             treename   = c(treename1, treename2),
             leafname   = leafname,
             logbase    = logbase,
             type       = "hist",
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.hist2D

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.hist3D" <-
function(x,
         treename1  = character(0),
         treename2  = character(0),
         treename3  = character(0),
         logbase    = "log2",
         option     = "HIST",
         canvasname = "Histogram3D",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.hist3D------")

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
      leafname <- "fPValue";
   }#if

   root.draw(x,
             treename   = c(treename1, treename2, treename3),
             leafname   = leafname,
             logbase    = logbase,
             type       = "hist",
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.hist3D

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.draw" <-
function(x,
         treename   = "",
         leafname   = "fInten",
         logbase    = "log2",
         type       = "graph",
         option     = "L",
         canvasname = "DrawCanvas",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.draw------")

   ## check for correct class
   if (!extends(class(x), "ProcesSet")) {
      stop(paste(sQuote("x"), "is not derived from class", sQuote("ProcesSet")));
   }#if

   rootfile <- rootFile(x);
   varlist  <- leafname;
   savename <- rootDrawName(canvasname, save.as);

   ## test for valid treenames
   ntree <- length(treename);
   if (ntree > 3) {
      stop(paste("function", dQuote("root.draw"), "can only be applied to at most 3 trees."));
   }#if
   treenames <- getTreeNames(rootfile);
   validtree <- treename%in%treenames;
   if (length(validtree[validtree==TRUE]) < ntree) {
      stop(paste("at least one tree", dQuote("treename"),
                 "is missing or not present in", sQuote(rootfile)));
   }#if

   ## need to add ROOT directory in rootfile
   if (is (x, "DataTreeSet")) {
      treename <- paste("DataSet", treename, sep="/");
   } else if (is (x, "ExprTreeSet")) {
      treename <- paste("PreprocesSet", treename, sep="/");
   } else if (is (x, "CallTreeSet")) {
      treename <- paste("CallSet", treename, sep="/");
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
   ntree      <- as.integer(ntree);
   treename1  <- paste(dq, treename[1], dq, sep="");
   treename2  <- paste(dq, treename[2], dq, sep="");
   treename3  <- paste(dq, treename[3], dq, sep="");
   varlist    <- paste(dq, varlist, dq, sep="");
   logbase    <- paste(dq, logbase, dq, sep="");
   type       <- paste(dq, type, dq, sep="");
   option     <- paste(dq, option, dq, sep="");
   savename   <- paste(dq, savename, dq, sep="");
   w          <- as.integer(w);
   h          <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDraw.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",",
                  ntree,",", treename1,",", treename2,",", treename3,",",
                  varlist,",", logbase,",", type,",", option,",", savename,",", 
                  w,",", h, ")", sq, sep="");

   if (save.as == "") {
      system(paste("root -l", macro));
   } else {
      system(paste("root -l -q", macro));
   }#if
}#root.draw
