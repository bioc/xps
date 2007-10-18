#==============================================================================#
# root.scheme.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# root.scheme: create new SchemeTreeSet for existing root scheme file
#==============================================================================#

"root.scheme" <-
function(rootfile = character(0),
         add.mask = FALSE)
{
   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## check for presence of chip name
   chipinfo <- as.list(getNameType(rootfile));
   chipname <- chipinfo$chipname;
   if (length(chipname) == 0 || nchar(chipname) < 1) {
      stop("missing chip name");
   }#if

   ## check for presence of correct chip type
   chiptype <- chipinfo$chiptype;
   TYPE <- c("GeneChip", "GenomeChip", "ExonChip");
   if (is.na(match(chiptype, TYPE))) {
      stop(paste(sQuote("chiptype"), "must be <GeneChip,GenomeChip,ExonChip>"));
   }#if

   ## get treenames and probe information
   treenames <- as.list(getTreeNames(rootfile));
   probeinfo <- as.list(getProbeInfo(rootfile));

   ## create new class
   set <- new("SchemeTreeSet",
              setname   = chipname,
              settype   = "scheme",
              rootfile  = rootfile,
              filedir   = dirname(rootfile),
              numtrees  = length(treenames),
              treenames = as.list(treenames),
              chipname  = chipname,
              chiptype  = chiptype,
              probeinfo = probeinfo);

   ## get mask for scheme
   if (add.mask) {
      chipMask(set) <- export(set,
                              treetype     = "scm",
                              varlist      = "fMask",
                              as.dataframe = TRUE,
                              verbose      = FALSE);
   }#if

   return(set);
}#root.scheme


