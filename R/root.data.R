#==============================================================================#
# root.data.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# root.data: create new DataTreeSet for existing root data file
#==============================================================================#

"root.data" <-
function(xps.scheme,
         rootfile = character(0),
         celnames = "*")
 {
   ## check for presence of class xps.scheme
   xps.scheme <- validSchemeTreeSet(xps.scheme);

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## check for correct chip name
   chipname  <- getChipName(rootFile(xps.scheme));
   treetitle <- unique(getTreeNames(rootfile, "cel", gettitle=TRUE));
   if (chipname != treetitle) {
      stop(paste("wrong", sQuote("xps.scheme"), "for chip", sQuote(treetitle)));
   }#if

   ## get treenames
   treenames <- getTreeNames(rootfile, "cel");

   if (celnames[1] != "*") {
      celnames  <- paste(namePart(celnames), "cel", sep=".");
      dupnames  <- celnames[duplicated(celnames)];
      if (length(dupnames) > 0) {
         stop(paste("duplicated celnames: ", dupnames, "\n"));
      }#if

      treenames <- treenames[match(celnames, treenames)];
      if (length(treenames[is.na(treenames)]) > 0) {
         stop(paste("some", sQuote("celnames"), " are not found in",
                    sQuote(rootfile)));
      }#if
   }#if

   ## create new class
   set <- new("DataTreeSet",
              setname   = "DataSet",
              settype   = "rawdata",
              rootfile  = rootfile,
              filedir   = dirname(rootfile),
              numtrees  = length(treenames),
              treenames = as.list(treenames),
              scheme    = xps.scheme);
   return(set);
}#root.data


