#==============================================================================#
# import.expr.scheme.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# import.expr.scheme: import CDF-file and annotation files for exression array
#==============================================================================#

"import.expr.scheme" <-
function(filename   = character(0),
         filedir    = getwd(),
         schemefile = character(0),
         probefile  = character(0),
         annotfile  = character(0),
         add.mask   = FALSE,
         verbose    = TRUE)
{
   ## check for presence of filename
   if (filename == "") {
      stop(paste(sQuote("filename"), "is missing"))
   }#if

   ## root scheme file
   rootfile <- rootDirFile(filename, filedir);

   ## check for presence of parameters
   if (!(is.character(schemefile) &&
         is.character(probefile) &&
         is.character(annotfile))) {
      stop("all arguments must be of type character");
   }#if

   ## check for presence of scheme file
   if (!file.exists(schemefile)) {
      stop("missing chip definition file *.CDF");
   }#if

   ## get chipname and extension
   scheme <- unlist(strsplit(schemefile, "/"));
   scheme <- scheme[length(scheme)];
   scheme <- unlist(strsplit(scheme, "\\."));

   chipname <- scheme[1];
   exten    <- scheme[length(scheme)];

   if (!((exten == "CDF") || (exten == "cdf"))) {
      stop(paste(sQuote("schemefile"), "is not a CDF-file"))
   }#if

   ## check for presence of probe file
   if (!file.exists(probefile)) {
      stop("missing probe file");
   }#if

   probe <- unlist(strsplit(probefile, "/"));
   probe <- probe[length(probe)];
   exten <- substr(probe, nchar(probe)-8, nchar(probe));
   if (exten != "probe.tab") {
      stop(paste(sQuote("probefile"),
                 "does not end with", sQuote("probe.tab")));
   }#if

   ## check for presence of annotation file
   if (annotfile == "") {
      warning("no annotation file is given");
   } else if (!file.exists(annotfile)) {
      stop("missing chip annotation file");
   } else {
      ann   <- unlist(strsplit(annotfile, "/"));
      ann   <- ann[length(ann)];
      exten <- substr(ann, nchar(ann)-8, nchar(ann));

      if (exten != "annot.csv") {
         stop(paste(sQuote("annotfile"),
                    "does not end with", sQuote("annot.csv")));
      }#if
   }#if

   ## create root scheme file for expression array
   r <- .C("ImportExprSchemes",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(schemefile),
           as.character(probefile),
           as.character(annotfile),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ImportExprSchemes")));
      return(NULL);
   }#if

   ## get treenames and probe information
   treenames <- as.list(getTreeNames(rootfile));
   probeinfo <- as.list(getProbeInfo(rootfile));

   ## create new class
   set <- new("SchemeTreeSet",
              setname   = chipname,
              settype   = "scheme",
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = length(treenames),
              treenames = as.list(treenames),
              chipname  = chipname,
              chiptype  = "GeneChip",
              probeinfo = probeinfo);

   ## get mask for scheme
   if (add.mask) {
      chipMask(set) <- export(set,
                              treetype     = "scm",
                              varlist      = "fMask",
                              as.dataframe = TRUE,
                              verbose      = TRUE);
   }#if

   return(set);
}#import.expr.scheme


