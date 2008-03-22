#==============================================================================#
# methods.TreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# rootFile:
# rootFile<-:
# fileDir:
# fileDir<-:
# setName:
# setName<-:
# setType:
# setType<-:
# treeNames:
# export:
# root.browser: 
#==============================================================================#


#------------------------------------------------------------------------------#
# TreeSet initialization:
#------------------------------------------------------------------------------#

"initialize.TreeSet" <-
function(.Object,
         setname   = character(),
         settype   = character(),
         rootfile  = NULL,
         filedir   = getwd(),
         numtrees  = 0,
         treenames = list(),
         ...) 
{
   if (debug.xps()) print("------initialize:TreeSet------")

   ## prevent effects of multiple initialization
   if (is.null(rootfile)) return(.Object);
   .Object@setname  <- setname;
   .Object@settype  <- settype;
   .Object@rootfile <- rootfile;
   .Object@filedir  <- filedir;

   rootfile <- rootDirFile(rootfile, filedir);

   .Object <- callNextMethod(.Object,
                             setname   = setname,
                             settype   = settype,
                             rootfile  = rootfile,
                             filedir   = filedir,
                             numtrees  = numtrees,
                             treenames = treenames,
                             ...);
   .Object@setname   <- setname;
   .Object@settype   <- settype;
   .Object@rootfile  <- rootfile;
   .Object@filedir   <- filedir;
   .Object@numtrees  <- numtrees;
   .Object@treenames <- treenames;
   .Object;
}#initialize.TreeSet

setMethod("initialize", "TreeSet", initialize.TreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("TreeSet",
   function(object) 
   {
      if (debug.xps()) print("------setValidity:TreeSet------")

      msg <- NULL;

      ## check rootfile
      if (!isROOTFile(object@rootfile)) {
         msg <- validMsg(msg, paste(sQuote("rootfile"),
                                   "is missing or is not a ROOT file"));
      }#if

      ## check for presence of directory for file
      if (!file.exists(object@filedir)) {
         msg <- validMsg(msg, paste(sQuote("filedir"),
                                   "is not a system directory"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# TreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("rootFile", signature(object="TreeSet"),
   function(object) object@rootfile
)#rootFile

setReplaceMethod("rootFile", signature(object="TreeSet", value="character"),
   function(object, value) {
      object@rootfile <- value;
      return(object);
   }
)#rootFile<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("fileDir", signature(object="TreeSet"),
   function(object) object@filedir
)#fileDir

setReplaceMethod("fileDir", signature(object="TreeSet", value="character"),
   function(object, value) {
      object@filedir <- value;
      return(object);
   }
)#fileDir<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("setName", signature(object="TreeSet"),
   function(object) object@setname
)#setName

setReplaceMethod("setName", signature(object="TreeSet", value="character"),
   function(object, value) {
      object@setname <- value;
      return(object);
   }
)#setName<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("setType", signature(object="TreeSet"),
   function(object) object@settype
)#setType

setReplaceMethod("setType", signature(object="TreeSet", value="character"),
   function(object, value) {
      object@settype <- value;
      return(object);
   }
)#setType<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("treeNames", signature(object="TreeSet"),
   function(object) object@treenames
)#treeNames


#------------------------------------------------------------------------------#
# TreeSet methods:
#------------------------------------------------------------------------------#

"exportTreeSet" <-
function(object, ...) 
{
   if (debug.xps()) print("------export.TreeSet------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootFile(object), class(object));

   ## check for number of trees
   if (nchar(object@numtrees) < 1) {
      stop("number of trees is null");
   }#if
}#exportTreeSet

setMethod("export", "TreeSet", exportTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"rootBrowser" <-
function(object) 
{
   if (debug.xps()) print("------rootBrowser.TreeSet------")

   savedir <- getwd();
   setwd(object@filedir);

   macro <- paste(system.file(package="xps"), "rootsrc/macroOpenBrowser.C", sep="/");
   xpsso <- system.file("libs", .Platform$r_arch,
                  paste("xps", .Platform$dynlib.ext, sep=''), package="xps");
   xpsso <- gsub("//", "/", xpsso);
   temp  <- as.character("\""); #"
   xpsso <- paste(temp, xpsso, temp, sep="");
   ## sQuote is ok on Mac ('') but false on Linux (´`)!!
#   macro <- sQuote(paste(macro, "(", xpsso, ")", sep=""));
   macro <- paste("'", macro, "(", xpsso, ")", "'", sep="");
   system(paste("root -l", macro));

   setwd(savedir);
}#rootBrowser

setMethod("root.browser", "TreeSet", rootBrowser);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# other potential methods:
#setMethod("root.canvas", "TreeSet", rootCanvas);
#setMethod("root.image", "TreeSet", rootImage);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




