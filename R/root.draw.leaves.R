#==============================================================================#
# root.draw.leaves.R: plotting functions based on ROOT TCanvas:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# root.density: 
# root.draw.leaves: 
#==============================================================================#

#"root.density" <-
#function(x,
#         treename   = "*",
#         logbase    = "log2",
#         canvasname = "DensityPlot",
#         w          = 540,
#         h          = 540) 
#{
#   if (debug.xps()) print("------root.density------")
#
#   ## select leafname for class
#   if (is (x, "DataTreeSet")) {
#      leafname <- "fInten";
#   } else if (is (x, "ExprTreeSet")) {
#      leafname <- "fLevel";
#   } else {
#      stop(paste(sQuote("x"), "must be  class <DataTreeSet,ExprTreeSet>"));
#   }#if
#
#   root.draw.leaves(x,
#                    treename   = treename,
#                    leafname   = leafname,
#                    logbase    = logbase,
#                    type       = "density",
#                    option     = "AL",
#                    sortopt    = 0,
#                    canvasname = canvasname,
#                    w          = w,
#                    h          = h);
#}#root.density

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.draw.leaves" <-
function(x,
         treename   = "*",
         leafname   = "fInten",
         logbase    = "log2",
         type       = "graph",
         option     = "L",
         sortopt    = 0,
         canvasname = "DrawCanvas",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.draw.leaves------")

   ## check for correct class
   if (!extends(class(x), "ProcesSet")) {
      stop(paste(sQuote("x"), "is not derived from class", sQuote("ProcesSet")));
   }#if

   rootfile  <- rootFile(x);
   treenames <- treeNames(x);
   varlist   <- leafname;

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
#   temp       <- as.character("\""); #"
   xpsso      <- paste(dq, xpsso, dq, sep="");
   rootfile   <- paste(dq, rootfile, dq, sep="");
   canvasname <- paste(dq, canvasname, dq, sep="");
   treename   <- paste(dq, treename, dq, sep="");
   varlist    <- paste(dq, varlist, dq, sep="");
   logbase    <- paste(dq, logbase, dq, sep="");
   type       <- paste(dq, type, dq, sep="");
   option     <- paste(dq, option, dq, sep="");
   sortopt    <- as.integer(sortopt);
   w          <- as.integer(w);
   h          <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDrawLeaves.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",",
                  treename,",", varlist,",", logbase,",", type,",", option,",",
                  sortopt,",", w,",", h, ")", sq, sep="");
   system(paste("root -l", macro));
}#root.draw.leaves
