#==============================================================================#
# class SchemeTreeSet: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# chipName:
# chipType:
# chipType<-:
# chipMask:
# chipMask<-:
# probeInfo:
# nrows:
# ncols:
# attachMask:
# removeMask:
# export:
#==============================================================================#


#------------------------------------------------------------------------------#
# SchemeTreeSet initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "SchemeTreeSet", 
   function(.Object,
            chipname  = character(),
            chiptype  = character(),
            probeinfo = list(),
            ...) 
   {
      if (debug.xps()) print("------initialize:SchemeTreeSet------")

      ## set default chiptype
      if (missing(chiptype) || chiptype == "") {
         chiptype <- "GeneChip";
      }#if

      .Object <- callNextMethod(.Object,
                                chipname  = chipname,
                                chiptype  = chiptype,
                                probeinfo = probeinfo,
                                ...);
      .Object@chipname  <- chipname;
      .Object@chiptype  <- chiptype;
      .Object@probeinfo <- probeinfo;
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("SchemeTreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:SchemeTreeSet------")

      msg <- NULL;

      ## check for correct chiptype
      if(!isChipType(object@chiptype)) {
         msg <- validMsg(msg, paste(sQuote("chiptype"), "must be",
                         sQuote("<GeneChip,GenomeChip,ExonChip>")));
      }#if

      ## check for correct settype
      if (object@settype != "scheme") {
         msg <- validMsg(msg,
                         paste(sQuote("settype"), "is not", sQuote("scheme")));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# SchemeTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("chipName", signature(object="SchemeTreeSet"),
   function(object) object@chipname
)#chipName

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("chipType", signature(object="SchemeTreeSet"),
   function(object) object@chiptype
)#chipType

setReplaceMethod("chipType", signature(object="SchemeTreeSet", value="character"),
   function(object, value) {
      object@chiptype <- value;
      return(object);
   }
)#chipType<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("chipMask", signature(object="SchemeTreeSet"),
   function(object) object@mask
)#chipMask

setReplaceMethod("chipMask", signature(object="SchemeTreeSet", value="data.frame"),
   function(object, value) {
      object@mask <- value;
      return(object);
   }
)#chipMask<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("probeInfo", signature(object="SchemeTreeSet"),
   function(object) object@probeinfo
)#probeInfo

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("nrows", signature(object="SchemeTreeSet"),
   function(object) object@probeinfo$nrows
)#nrows

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("ncols", signature(object="SchemeTreeSet"),
   function(object) object@probeinfo$ncols
)#ncols


#------------------------------------------------------------------------------#
# SchemeTreeSet methods:
#------------------------------------------------------------------------------#

setMethod("attachMask", signature(object="SchemeTreeSet"),
   function(object) {
      if (debug.xps()) print("------attachMask.SchemeTreeSet------")

      chipMask(object) <- export(object,
                                 treetype     = "scm",
                                 varlist      = "fMask",
                                 as.dataframe = TRUE,
                                 verbose      = FALSE);
      return(object);
   }
)#attachMask

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeMask", signature(object="SchemeTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeMask.SchemeTreeSet------")

      chipMask(object) <- data.frame(matrix(nr=0,nc=0));
      gc(); #????
      return(object);
   }
)#removeMask

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"exportSchemeTreeSet" <-
function(object,
         treetype     = character(0),
         varlist      = "*",
         outfile      = character(0),
         sep          = "\t",
         as.dataframe = FALSE,
         verbose      = TRUE,
         ...) 
{
   if (debug.xps()) print("------export.SchemeTreeSet------")

   callNextMethod();

   ## get scheme file
   scheme <- rootFile(object);

   ## check for presence of parameters
   if (!(is.character(treetype) && is.character(varlist))) {
      stop("arguments <treetype,varlist> must be of type character");
   }#if

   ## check for presence of valid tree type (see utils.R)
   TYPE <- SCMTYPE;
   if (chipType(object) == "ExonChip") {
      TYPE <- c(SCMTYPE, EXNTYPE);
   }#if
   type <- match(treetype, TYPE);
   if (length(type) == 0 || is.na(type)) {
      stop(paste("invalid parameter", sQuote("treetype")));
   }#if

   if (nchar(varlist) < 1) {
      varlist <- "*";
   }#if

   ## get tree name "treeset.treename.treetype"
   treename <- paste(chipName(object), chipName(object), treetype, sep=".");

   ## check for presence of outfile=/path/outname
   outname <- paste(chipName(object), treetype, sep="_");
   outfile <- validOutfile(outname, outfile);

   ## check for presence of valid separator
   validSeparator(sep);

   ## export scheme from root file
   r <- .C("ExportScheme",
           as.character(scheme),
           as.character(treename),
           as.character(varlist),
           as.character(outfile),
           as.character(sep),
           as.integer(verbose),
           PACKAGE="xps");

   ## import outfile as dataframe
   ds <- NULL;
   if (as.dataframe) {
      ## use quote = "\"", since annotation has single quotes e.g. "2'-PDE"
      ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep=sep, row.names=NULL, quote = "\""); #"
   }#if

   return(ds);
}#exportSchemeTreeSet

setMethod("export", "SchemeTreeSet", exportSchemeTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
