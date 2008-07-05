#==============================================================================#
# root.call.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# root.call: create new CallTreeSet for existing root present call file
#==============================================================================#

"root.call" <-
function(xps.scheme,
         rootfile  = character(0),
         treetype  = character(0),
         treenames = "*")
 {
   ## get schemefile and chiptype
   xps.scheme <- validSchemeTreeSet(xps.scheme);
   schemefile <- validROOTFile(rootFile(xps.scheme), "SchemeTreeSet");
   chiptype   <- chipType(xps.scheme);

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "CallTreeSet");

   ## check for correct treetype
   if (!is.na(match(treetype, CALTYPE))) {
      setname <- "CallSet";
      settype <- "preprocess";
   } else {
      stop(paste("invalid treetype", sQuote(treetype)));
   }#if

   ## get calltype
   if (treetype == "dc5") {
      calltype <- "mas5";
   } else if (treetype == "dab") {
      calltype <- "dabg";
   } else {
      calltype <- "none";
   }#if

   ## get trees with extension "treetype"
   nametypes <- getTreeNames(rootfile, treetype);
   if (length(nametypes) == 0) {
      stop(paste("root file", sQuote(rootfile), "has no trees of type", sQuote(treetype)));
   }#if

   ## get treenames
   if (treenames[1] == "*") {
      treenames <- nametypes;
   } else {
      treenames <- paste(namePart(treenames), treetype, sep=".");
      nametypes <- nametypes[match(treenames, nametypes)];
      if (length(nametypes[is.na(nametypes)]) > 0) {
         stop(paste("some", sQuote("treenames"), " are not found in",
                    sQuote(rootfile)));
      }#if
   }#if

   ## get tree names "treeset.treename.treetype"
   numtrees <- length(treenames);
   if (numtrees == 1) {
      treeset <- paste(setname, treenames, sep=".");
   } else if (numtrees > 1) {
      treeset <- paste(rootfile, setname, treenames, sep="/");
   }#if

   ## export pvalue data from root file
   outfile  <- "tmp_pvalue.txt";

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treeset),
           as.integer(numtrees),
           as.character(treetype),
           as.character("fUnitName:fPValue"),
           as.character(outfile),
           as.character("\t"),
           as.integer(FALSE),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## import outfile as dataframe
   data <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);

   ## export present call from root file
   outfile  <- "tmp_call.txt";

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treeset),
           as.integer(numtrees),
           as.character(treetype),
           as.character("fUnitName:fCall"),
           as.character(outfile),
           as.character("\t"),
           as.integer(FALSE),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## import outfile as dataframe
   dc <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);

   ## create new class
   set <- new("CallTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = dirname(rootfile),
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = xps.scheme,
              data      = data,
              params    = list(),
              calltype  = calltype,
              detcall  = dc);

   return(set);
}#root.call

