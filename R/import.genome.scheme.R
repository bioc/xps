#==============================================================================#
# import.genome.scheme.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# import.genome.scheme: import CLF, PGF and annotation files for whole genome array
#==============================================================================#

"import.genome.scheme" <-
function(filename   = character(0),
         filedir    = getwd(),
         layoutfile = character(0),
         schemefile = character(0),
         transcript = character(0),
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
   if (!(is.character(layoutfile) && is.character(schemefile) &&
         is.character(transcript))) {
      stop("all arguments must be of type character");
   }#if

   ## check for presence of layout file
   if (!file.exists(layoutfile)) {
      stop("missing CEL Layout File *.CLF");
   }#if

   scheme <- unlist(strsplit(layoutfile, "/"));
   scheme <- scheme[length(scheme)];
   scheme <- unlist(strsplit(scheme, "\\."));
   exten  <- scheme[length(scheme)];

   if (!((exten == "CLF") || (exten == "clf"))) {
      stop(paste("argument", sQuote("layoutfile"), "is not a CLF-file"))
   }#if

   ## check for presence of scheme file
   if (!file.exists(schemefile)) {
      stop("missing Probe Group File *.PGF");
   }#if

   ## get chipname and extension
   scheme <- unlist(strsplit(schemefile, "/"));
   scheme <- scheme[length(scheme)];
   scheme <- unlist(strsplit(scheme, "\\."));

   chipname <- scheme[1];
   exten    <- scheme[length(scheme)];

   if (!((exten == "PGF") || (exten == "pgf"))) {
      stop(paste(sQuote("schemefile"), "is not a PGF-file"))
   }#if

   ## check for presence of transcript annotation file
   if (!file.exists(transcript)) {
      stop("missing transcript annotation file");
   }#if

   annot <- unlist(strsplit(transcript, "/"));
   annot <- annot[length(annot)];
   exten <- substr(annot, nchar(annot)-2, nchar(annot));
   if (!((length(grep("transcript",annot)) > 0) && (exten == "csv"))) {
      stop(paste(sQuote("transcript"), "is not a transcript annotation file"));
   }#if

   ## create root scheme file for whole genome array
   r <- .C("ImportGenomeSchemes",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(layoutfile),
           as.character(schemefile),
           as.character(transcript),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ImportGenomeSchemes")));
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
              chiptype  = "GenomeChip",
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
}#import.genome.scheme
