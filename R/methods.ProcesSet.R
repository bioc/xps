#==============================================================================#
# methods.ProcesSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# schemeFile:
# schemeFile<-:
# schemeSet:
# schemeSet<-:
# chipName:
# chipType:
# getTreeData:
# validData:
# export:
# boxplot:
# mboxplot:
# hist:
#==============================================================================#


#------------------------------------------------------------------------------#
# ProcesSet initialization:
#------------------------------------------------------------------------------#

"initialize.ProcesSet" <-
function(.Object,
         scheme = NULL,
         data   = data.frame(),
         params = list(),
         ...) 
{
   if (debug.xps()) print("------initialize:ProcesSet------")

   ## prevent effects of multiple initialization
   if (is.null(scheme))  return(.Object);
   if (is.null(data))    data   <- data.frame(matrix(nr=0,nc=0));
   if (!is.list(params)) params <- as.list(params);

   .Object@scheme = scheme;
   .Object@data   = data;
   .Object@params = params;

   .Object <- callNextMethod(.Object,
                             scheme = scheme,
                             data   = data,
                             params = params,
                             ...);
   .Object@scheme = scheme;
   .Object@data   = data;
   .Object@params = params;
   .Object;
}#initialize.ProcesSet

setMethod("initialize", "ProcesSet", initialize.ProcesSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("ProcesSet",
   function(object) 
   {
      if (debug.xps()) print("------setValidity:ProcesSet------")

      msg <- NULL;

      ## check setname
      if (!(is(object@setname, "character") && nchar(object@setname) > 0)) {
         msg <- validMsg(msg, paste(sQuote("setname"), "is missing"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# ProcesSet accessors:
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("schemeFile", signature(object="ProcesSet"),
   function(object) {
      scheme <- object@scheme;
      return(scheme@rootfile);
   }
)#schemeFile

setReplaceMethod("schemeFile", signature(object="ProcesSet", value="character"),
   function(object, value) {
      scheme <- object@scheme;
      scheme@rootfile <- value;
      object@scheme <- scheme;
      return(object);
   }
)#schemeFile<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("schemeSet", signature(object="ProcesSet"),
   function(object) object@scheme
)#schemeSet

setReplaceMethod("schemeSet", signature(object="ProcesSet", value="SchemeTreeSet"),
   function(object, value) {
      object@scheme <- value;
      return(object);
   }
)#schemeSet<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("chipName", signature(object="ProcesSet"),
   function(object) {
      scheme <- object@scheme;
      return(scheme@chipname);
   }
)#chipName

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("chipType", signature(object="ProcesSet"),
   function(object) {
      scheme <- object@scheme;
      return(scheme@chiptype);
   }
)#chipType


#------------------------------------------------------------------------------#
# ProcesSet methods:
#------------------------------------------------------------------------------#

setMethod("getTreeData", signature(object="ProcesSet"),
   function(object,
            treetype = "cel",
            varlist  = "fInten")
   {
      if (debug.xps()) print("------getTreeData.ProcesSet------")

      ds <- export(object,
                   treetype     = treetype,
                   varlist      = varlist,
                   as.dataframe = TRUE);
      return(ds);
   }
)#getTreeData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"dataProcesSet" <-
function(object,
         which = "UnitName") 
{
   if (debug.xps()) print("------dataProcesSet------")

   ## check for presence of data
   data  <- object@data;
   if (min(dim(data)) == 0) {
      stop(paste("slot", sQuote("data"), "has no data"));
   }#if

   ## use names from column "which" as rownames
   if (!is.na(match(which, colnames(data)))) {
      rownames(data) <- data[,which];
   }#if

   treenames <- namePart(object@treenames);
   datanames <- namePart(colnames(data));

   return(data[,!is.na(match(datanames, treenames))]);
}#dataProcesSet

setMethod("validData", "ProcesSet", dataProcesSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"exportProcesSet" <-
function(object,
         treenames    = "*",
         treetype     = character(0),
         varlist      = "*",
         outfile      = character(0),
         sep          = "\t",
         as.dataframe = FALSE,
         verbose      = TRUE,
         ...) 
{
   if (debug.xps()) print("------export.ProcesSet------")

   callNextMethod();

   ## get name of root file
   filename <- rootFile(object);

   ## get name of root scheme file
   schemefile <- validROOTFile(schemeFile(object), class(object));

   ## get chiptype
   chiptype <- chipType(object);

   ## get treeset
   treeset <- setName(object);

   ## check for treenames
   numtrees <- length(treenames);
   if (treenames[1] != "*") {
      treenames <- namePart(treenames);
   }#if

   ## check for presence of correct tree type
   datatype <- setType(object);
   treetype <- validTreetype(treetype, datatype);

   ## check for varlist
   if (nchar(varlist) < 1) {
      varlist <- "*";
   }#if

   ## check for presence of valid separator and outfile
   exten <- validSeparator(sep);
   if (length(outfile) && outfile != "") {
      if (regexpr(".txt",outfile) > 0 || regexpr(".csv",outfile) > 0) {
         outfile <- substr(outfile, 1, nchar(outfile)-4);
      }#if
      outfile <- paste(outfile, exten, sep=".");
   }#if

   ## check for presence of outfile=/path/outname
   outname <- "";
   if (numtrees > 1 || treenames[1] == "*") {
      outname <- "PivotTable";
   } else {
      outname <- treenames;
   }#if
   outfile <- validOutfile(outname, outfile);

   ## get tree names "treeset.treename.treetype"
   if (numtrees == 1) {
      treenames <- paste(treeset, treenames, treetype, sep=".");
   } else if (numtrees > 1) {
      treenames <- paste(treenames, treetype, sep=".");
      treenames <- paste(filename, treeset, treenames, sep="/");
   }#if

   ## export data from root file
   r <- .C("ExportData",
           as.character(filename),
           as.character(schemefile),
           as.character(chiptype),
           as.character(datatype),
           as.character(treenames),
           as.integer(numtrees),
           as.character(treetype),
           as.character(varlist),
           as.character(outfile),
           as.character(sep),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## import outfile as dataframe
   ds <- NULL;
   if (as.dataframe) {
      ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep=sep, row.names=NULL);
   }#if

   return(ds);
}#exportProcesSet

setMethod("export", "ProcesSet", exportProcesSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("boxplot", signature(x="ProcesSet"),
   function(x,
            which   = "",
            size    = 0,
            transfo = log2,
            range   = 0,
            names   = "namepart",
            ...) 
   {
      if (debug.xps()) print("------boxplot.ProcesSet------")

      ds <- validData(x, which=which);
      if (size > 1)             ds <- ds[seq(1,nrow(ds),len=size),];
      if (is.function(transfo)) ds <- transfo(ds);

      if (is.null(names))              names <- colnames(ds)
      else if (names[1] == "namepart") names <- namePart(colnames(ds))
      else                             ds    <- ds[, names, drop=F];

      boxplot(ds,
              range = range,
              names = names,
              ...);
   }
)#boxplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("mboxplot", signature(x="ProcesSet"),
   function(x,
            which   = "",
            size    = 0,
            transfo = log2,
            method  = "mean",
            range   = 0,
            ylim    = c(-1,1),
            outline = FALSE,
            names   = "namepart",
            ...) 
   {
      if (debug.xps()) print("------mboxplot.ProcesSet------")

      m <- validData(x, which=which);
      if (size > 1)             m <- m[seq(1,nrow(m),len=size),];
      if (is.function(transfo)) m <- transfo(m);

      if (is.null(names))              names <- colnames(m)
      else if (names[1] == "namepart") names <- namePart(colnames(m))
      else                             m     <- m[, names, drop=F];

      if (method == "median") {
         mn  <- apply(m, 1, median);
      } else {
         mn  <- rowMeans(m);
      }#if
      rle <- data.frame(sweep(m, MARGIN=1, STATS=mn, FUN='-'));

      boxplot(rle,
              range   = range,
              xlab    = '',
              xaxt    = 'n',
              ylim    = ylim,
              outline = outline);
      abline(h=0);
      axis(side=1, outer=F, at=1:length(names), labels=names, ...);
   }
)#mboxplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("hist", signature(x="ProcesSet"),
   function(x,
            which   = "",
            size    = 0,
            transfo = log2,
            ylab    = "density",
            xlab    = "log intensity",
            type    = "l",
            col     = 1:6,
            ...) 
   {
      if (debug.xps()) print("------hist.ProcesSet------")

      m <- validData(x, which=which);
      if (size > 1) m <- m[seq(1,nrow(m),len=size),];
      if (is.function(transfo)) m <- transfo(m);

      plotDensity(m,
                  ylab = ylab,
                  xlab = xlab,
                  type = type,
                  col  = col,
                  ...);
   }
)#hist

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
