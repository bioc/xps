#==============================================================================#
# root.expr.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# root.expr: create new ExprTreeSet for existing root expression file
#==============================================================================#

"root.expr" <-
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
   rootfile <- validROOTFile(rootfile, "ExprTreeSet");

   ## check for correct treetype
   if (!is.na(match(treetype, NRMTYPE))) {
      setname <- "NormSet";
      settype <- "normation";
   } else if (!is.na(match(treetype, PRETYPE))) {
      setname <- "PreprocesSet";
      settype <- "preprocess";
   } else {
      stop(paste("invalid treetype", sQuote(treetype)));
   }#if

   ## get exprtype
   if (treetype == "adf") {
      exprtype <- "mas4";
   } else if (treetype == "tbw") {
      exprtype <- "mas5";
   } else if (treetype == "mdp") {
      setname  <- "PreprocesSet";
      settype  <- "preprocess";
      exprtype <- "rma";
   } else {
      exprtype <- "none";
   }#if

   ## get normtype
   if (treetype == "tmn") {
      normtype <- "mean";
   } else if (treetype == "med") {
      normtype <- "median";
   } else if (treetype == "low") {
      normtype <- "lowess";
   } else if (treetype == "sup") {
      normtype <- "supsmu";
   } else {
      normtype <- "none";
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
      dupnames  <- treenames[duplicated(treenames)];
      if (length(dupnames) > 0) {
         stop(paste("duplicated treenames: ", dupnames, "\n"));
      }#if

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

   ## export data from root file
   outfile  <- "tmp_expr.txt";

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treeset),
           as.integer(numtrees),
           as.character(treetype),
           as.character("fUnitName:fLevel"),
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

   ## create new class
   set <- new("ExprTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = dirname(rootfile),
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = xps.scheme,
              data      = data,
              params    = list(),
              exprtype  = exprtype,
              normtype  = normtype);

   return(set);
}#root.expr


