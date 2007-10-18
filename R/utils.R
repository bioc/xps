#==============================================================================#
# utils.R: contains utilty functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# xpsOptions:
# isChipType:
# isROOTFile:
# rootDirFile:
# validLogbase:
# validMsg:
# validOutfile:
# validROOTFile:
# validSchemeTreeSet:
# validSeparator:
# validTempDir: 
# validTranscriptOption:
# validTreenames:
# validTreetype:
# getDatatype:
# exonLevel:
# extenPart:
# namePart:
# type2Exten:
# plotDensity:
# getNameType:
# getChipName:
# getChipType:
# getProbeInfo:
# getNumberTrees:
# getTreeNames:
#==============================================================================#


#------------------------------------------------------------------------------#
# treetype: parameter used to get tree with corresponding tree extension
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ## scheme tree extensions are defined in src/XPSSchemes.cxx
   SCMTYPE <- c("scm", "idx", "prb", "ann");
   EXNTYPE <- c("cxy", "exn", "pbs", "anx", "anp");

   ## data tree extensions are defined in src/XPSProcessing.cxx
   RAWTYPE <- c("cel");
   PRETYPE <- c("int",
                "sbg", "wbg", "rbg", "gbg",
                "cmn", "cmd", "clw", "css", "cqu",
                "dc5", "dab",
                "amn", "gmn", "wmn", "wdf", "adf", "tbw", "mdp");
   NRMTYPE <- c("tmn", "med", "ksm", "low", "sup", "qua", "mdp");


#------------------------------------------------------------------------------#
# xpsOptions: set options for xps
# currently only used to enable/disable debug messages
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
OPTIONS <- new.env(hash=TRUE, parent=emptyenv());
OPTIONS$debug <- FALSE;

xpsOptions <- function(debug=FALSE) {
   if (debug == FALSE) {
      OPTIONS$debug <- FALSE;
   } else {
      OPTIONS$debug <- TRUE;
   }#if
}#xpsOptions

debug.xps <- function(x) {OPTIONS$debug};

#------------------------------------------------------------------------------#
# isChipType: check if chiptype is "GeneChip", "GenomeChip" or "ExonChip"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
isChipType <- function(chiptype) {
   if (debug.xps()) print("------isChipType------")

   TYPE <- c("GeneChip", "GenomeChip", "ExonChip");
   if (is.na(match(chiptype, TYPE))) return(FALSE);
   return(TRUE);
}#isChipType

#------------------------------------------------------------------------------#
# isROOTFile: check if filename is a ROOT file, i.e. "filename.root"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
isROOTFile <- function(filename) {
   if (debug.xps()) print("------isROOTFile------")

   if (!file.exists(filename)) {
      return(FALSE);
   }#if

   rtf <- unlist(strsplit(filename, "/"));
   rtf <- rtf[length(rtf)];
   xtn <- substr(rtf, nchar(rtf)-3, nchar(rtf));
   if (length(xtn) == 0 || xtn != "root") {
      return(FALSE);
   }#if

   return(TRUE);
}#isROOTFile

#------------------------------------------------------------------------------#
# rootDirFile: return root file as /filedir/filename.root
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootDirFile <- function(filename=character(0), filedir=character(0)) {
   if (debug.xps()) print("------rootDirFile------")

   ## for case that filename=/path/filename
   if (length(filename) == 0) {
      filename <- "ERROR_MISSING_FILENAME";
   } else {
      filename <- basename(filename);
   }#if

   ## set default directory
   if (dirname(filename) != ".") {
      filedir <- dirname(filename);
   }#if
   if (filedir == "" || filedir == ".") {
      filedir <- as.character(getwd());
   }#if
   ## remove trailing "/"
   if (substr(filedir, nchar(filedir), nchar(filedir)) == "/") {
      filedir <- substr(filedir, 0, nchar(filedir)-1);
   }#if

   ## set rootfile to /fullpath/filename.root
   exten <- unlist(strsplit(filename, "\\."));
   if (exten[length(exten)] == "root") {
      rootfile <- paste(filedir, "/", filename, sep="");
   } else {
      rootfile <- paste(filedir, "/", filename, ".root", sep="");
   }#if

   return(rootfile);
}#rootDirFile

#------------------------------------------------------------------------------#
# validLogbase: check for presence of valid logbase
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validLogbase <- function(logbase) {
   if (debug.xps()) print("------validLogbase------")

   if (!(identical(logbase, "0")    || identical(logbase, "log") ||
         identical(logbase, "log2") || identical(logbase, "log10"))) {
      stop(paste(sQuote(logbase), "is not a valid logbase option"));
   }#if

   return(as.character(logbase));
}#validLogbase

#------------------------------------------------------------------------------#
# validMsg: utility function for setValidity
# taken from package BioBase: tools.R
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validMsg <- function(msg, result) {
   if (is.character(result)) {
      append(msg, result);
   } else {
      msg;
   }#if
}#validMsg

#------------------------------------------------------------------------------#
# validOutfile: check for presence of valid outfile = /filepath/filename
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validOutfile <- function(outname, outfile=character(0)) {
   if (debug.xps()) print("------validOutfile------")

   filename <- "";
   filepath <- paste(getwd(), "/", sep="");

   if (!(length(outfile) <= 0 || outfile == "")) {
      nsep <- gregexpr("/", outfile)[[1]];

      if (!(length(nsep) ==1 && nsep < 0)) {
         nsep <- nsep[length(nsep)];
         filepath <- substr(outfile, 1, nsep);
      }#if

      if (!(filepath != "" && file.exists(filepath))) {
         stop(paste("path of", sQuote("outfile"), "is not a valid path"));
      }#if

      filename <- substring(outfile, nsep+1);
   }#if

   if (debug.xps()) print(paste("filepath = ",filepath));
   if (debug.xps()) print(paste("filename = ",filename));

   if (filename == "") {
      filename <- paste(outname, "txt", sep=".");
      outfile  <- paste(filepath, filename, sep="");
   }#if

   return(outfile);
}#validOutfile

#------------------------------------------------------------------------------#
# validROOTFile: check if rootfile is a valid ROOT TFile, i.e. "myfile.root"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validROOTFile <- function(rootfile, treeset) {
   if (debug.xps()) print("------validROOTFile------")

   if (!isROOTFile(rootfile)) {
      stop(paste(sQuote("rootfile"), "of class", sQuote(treeset),
          "is missing or not a ROOT file *.root"));
   }#if

   filedir <- dirname(rootfile);
   if (!(is(filedir, "character") && file.exists(filedir))) {
      stop(paste(sQuote("filedir"), "is not a system directory"));
   }#if

   return(rootfile);
}#validROOTFile

#------------------------------------------------------------------------------#
# validSchemeTreeSet: check if "object" is a valid class with correct parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validSchemeTreeSet <- function(object) {
   if (debug.xps()) print("------validSchemeTreeSet------")

   if (!is (object, "SchemeTreeSet")) {
      stop(paste(sQuote("object"), "is not class", sQuote("SchemeTreeSet")));
   }#if

   ## check for correct settype
   settype  <- as.character(object@settype);
   if (settype != "scheme") {
      stop(paste(sQuote("object"), "is not of settype", sQuote("scheme")));
   }#if

   ## check for presence of root scheme file
   scheme   <- as.character(object@rootfile);
   scheme <- validROOTFile(scheme);

   ## check for presence of chip name
   chipname <- as.character(object@chipname);
   if (length(chipname) == 0 || nchar(chipname) < 1) {
      stop("missing chip name");
   }#if

   ## check for presence of correct chip type
   if(!isChipType(object@chiptype)) {
      stop(paste(sQuote("chiptype"), "of class", sQuote("SchemeTreeSet"),
                 "must be", sQuote("<GeneChip,GenomeChip,ExonChip>")));
   }#if

   return(object);
}#validSchemeTreeSet

#------------------------------------------------------------------------------#
# validSeparator: check for presence of valid separator
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validSeparator <- function(sep) {
   if (debug.xps()) print("------validSeparator------")

   if(!(identical(sep, " ") || identical(sep, "\t") ||
        identical(sep, ",") || identical(sep, ";"))) {
      stop(paste(sQuote(sep), "is not a valid separator"));
   }#if
}#validSeparator

#------------------------------------------------------------------------------#
# validTempDir: check for presence of valid temporary directory
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validTempDir <- function(tmpdir) {
   if (debug.xps()) print("------validTempDir------")

   if (tmpdir != "") {
      if (!file.exists(tmpdir)) {
         stop(paste(sQuote("tmpdir"), "is not a system directory"));
      }#if
      if (substr(tmpdir,nchar(tmpdir), nchar(tmpdir)) == "/") {
         tmpdir <- substr(tmpdir, 0, nchar(tmpdir)-1);
      }#if
   }#if

   return(tmpdir);
}#validTempDir

#------------------------------------------------------------------------------#
# validTranscriptOption: check for presence of valid transcript option
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validTranscriptOption <- function(transcript) {
   if (debug.xps()) print("------validTranscriptOption------")

   if (!(identical(transcript, "transcript") ||
         identical(transcript, "exon") ||
         identical(transcript, "probeset"))) {
      stop(paste(sQuote(transcript), "is not a valid transcript option"));
   }#if

   return(as.character(transcript));
}#validTranscriptOption

#------------------------------------------------------------------------------#
# validTreenames: remove potential reference tree
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validTreenames <- function(treenames, refname) {
   if (debug.xps()) print("------validTreenames------")

   treexten  <- extenPart(treenames);
   treenames <- namePart(treenames);
   reference <- match(refname, treenames);

   if (!is.na(reference)) treenames <- treenames[-reference];

   return(paste(treenames, treexten, sep="."));
}#validTreenames

#------------------------------------------------------------------------------#
# validTreetype: check for presence of valid tree extension
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validTreetype <- function(treetype, datatype) {
   if (debug.xps()) print("------validTreetype------")

   if (datatype == "scheme") {
      TYPE <- c(SCMTYPE, EXNTYPE);
   } else if (datatype == "rawdata") {
      TYPE <- RAWTYPE;
   } else if (datatype == "preprocess") {
      TYPE <- PRETYPE;
   } else if (datatype == "normation") {
      TYPE <- NRMTYPE;
   } else if (datatype == "all") {
      TYPE <- c(RAWTYPE, PRETYPE, NRMTYPE);
   } else {
      stop(paste("invalid parameter", sQuote("datatype")));
   }#if

   if (is.na(match(treetype, TYPE))) {
      stop(paste("invalid treetype", sQuote(treetype),
                 "for datatype", sQuote(datatype)));
   }#if

   return(treetype);
}#validTreetype

#------------------------------------------------------------------------------#
# getDatatype: get data type for tree extension treetype
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getDatatype <- function(treetype) {
   if (debug.xps()) print("------getDatatype------")

   if (!is.na(match(treetype, RAWTYPE))) {
      datatype <- "rawdata";
   } else if (!is.na(match(treetype, PRETYPE))) {
      datatype <- "preprocess";
   } else if (!is.na(match(treetype, NRMTYPE))) {
      datatype <- "normation";
   } else {
      stop(paste("invalid parameter", sQuote("treetype")));
   }#if

   return(datatype);
}#getDatatype

#------------------------------------------------------------------------------#
# exonLevel: utility function for setValidity
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exonLevel <- function(exonlevel="", chiptype, as.sum=TRUE) {
   if (debug.xps()) print("------exonLevel------")

   exlevel <- 0;
   if (chiptype == "GenomeChip") {
      if (exonlevel == "") {
         stop(paste("invalid argument", sQuote("exonlevel")));
      }#if

      ## levels are defined in src/XPSSchemes.h
      level <- unlist(strsplit(exonlevel, "\\+"));
      LEVEL <- c("core", "metacore", "affx", "all");
      CODE  <- c( 1024,   8192,       60,     9276);
      level <- match(level, LEVEL);

      if (is.na(all(level))) {
         stop(paste("invalid argument", sQuote("exonlevel")));
      }#if

      ## replace "core" with "core+metacore"
      if (any(level == 1)) level <- c(level, 1, 2);
      if (any(level == 4)) level <- c(1:3); 
      level <- unique(level);

      if (as.sum) exlevel <- sum(CODE[level])
      else        exlevel <- CODE[level];
   } else if (chiptype == "ExonChip") {
      if (exonlevel == "") {
         stop(paste("invalid argument", sQuote("exonlevel")));
      }#if

      ## levels are defined in src/XPSSchemes.h
      level <- unlist(strsplit(exonlevel, "\\+"));
      LEVEL <- c("core", "metacore", "extended", "metaextended", "full", "metafull", "ambiguous", "affx", "all");
      CODE  <- c( 1024,   8192,       512,        4096,           256,    2048,       128,         60,    16256);
      level <- match(level, LEVEL);

      if (is.na(all(level))) {
         stop(paste("invalid argument", sQuote("exonlevel")));
      }#if

      ## replace e.g. "core" with "core+metacore"
      if (any(level == 1)) level <- c(level, 1, 2);
      if (any(level == 3)) level <- c(level, 3, 4);
      if (any(level == 5)) level <- c(level, 5, 6);
      if (any(level == 9)) level <- c(1:8);
      level <- unique(level);

      if (as.sum) exlevel <- sum(CODE[level])
      else        exlevel <- CODE[level];
   }#if

   return(exlevel);
}#exonLevel

#------------------------------------------------------------------------------#
# extenPart: utility function to extract extension part from "name.exten"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
extenPart <- function(names, as.unique=TRUE) {
   if (debug.xps()) print("------extenPart------")

   strg <- unlist(names);
   len  <- length(strg);
   strg <- sapply(strsplit(strg,"\\."), function(x) x[2]);
   if (as.unique == TRUE) {
      strg <- unique(strg);
   }#if
 
   return(strg[!is.na(strg)]);
}#extenPart

#------------------------------------------------------------------------------#
# namePart: utility function to extract name part from "name.exten"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
namePart <- function(names) {
   if (debug.xps()) print("------namePart------")

   strg <- unlist(names);
   len  <- length(strg);
   strg <- sapply(strsplit(strg,"\\."), function(x) x[1]);
 
   return(strg[!is.na(strg)]);
}#namePart

#------------------------------------------------------------------------------#
# type2Exten: utility function to get extension for method type
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type2Exten <- function(type, datatype) {
   if (debug.xps()) print("------type2Exten------")

   ## type and extensions are defined in src/XPSProcessing.cxx
   if (datatype == "preprocess") {
      TYPE <- c("mean", "median", "lowess", "supsmu", "quantile",
                "avgdiff", "tukeybiweight", "medianpolish");
      XTEN <- c("cmn", "cmd", "clw", "css", "cqu", "adf", "tbw", "mdp");
   } else if (datatype == "normation") {
      TYPE <- c("mean", "median", "ksmooth", "lowess", "supsmu",
                "quantile", "medianpolish");
      XTEN <- c("tmn", "med", "ksm", "low", "sup", "qua", "mdp");
   } else {
      stop(paste("invalid parameter", sQuote("datatype")));
   }#if
 
   return(XTEN[match(type, TYPE)]);
}#type2Exten

#------------------------------------------------------------------------------#
# plotDensity: function to plot density
# taken from package affy: plot.density.R
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plotDensity <- function(mat,
                        ylab="density", xlab="x", type="l", col=1:6,
                        ...) {
  
   x.density <- apply(mat, 2, density);

   all.x <- do.call("cbind", lapply(x.density, function(x) x$x));
   all.y <- do.call("cbind", lapply(x.density, function(x) x$y));
  
   matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, ...);

   invisible(list(all.x=all.x, all.y=all.y));
}#plotDensity


################################################################################
# utility functions based on rwrapper.h C-functions
################################################################################

#------------------------------------------------------------------------------#
# getNameType: get chip name and type from root scheme file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getNameType <- function(rootfile) {
   if (debug.xps()) print("------getNameType------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## get chip name and type
   info <- .C("ChipNameType",
              as.character(rootfile),
              nametype=character(2),
              PACKAGE="xps")$nametype;

   return(list(chipname=info[1], chiptype=info[2]));
}#getNameType

#------------------------------------------------------------------------------#
# getChipName: get chip name from root scheme file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getChipName <- function(rootfile) {
   return(getNameType(rootfile)$chipname);
}#getChipName

#------------------------------------------------------------------------------#
# getChipType: get chip type from root scheme file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getChipType <- function(rootfile) {
   return(getNameType(rootfile)$chiptype);
}#getChipType

#------------------------------------------------------------------------------#
# getProbeInfo: get GeneChip probe information from root scheme file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getProbeInfo <- function(rootfile) {
   if (debug.xps()) print("------probeInfo------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## get probe info
   info <- .C("GeneChipProbeInfo",
              as.character(rootfile),
              value=integer(8),
              PACKAGE="xps")$value;

   return(list(nrows=info[1], ncols=info[2], nprobes=info[3], ncontrols=info[4],
            ngenes=info[5], nunits=info[6], nprobesets=info[7], naffx=info[8]));
}#getProbeInfo

#------------------------------------------------------------------------------#
# getNumberTrees: get number of trees saved in root file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getNumberTrees <- function(rootfile, treetype="*", setname=NULL) {
   if (debug.xps()) print("------getNumberTrees------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## get number of trees
   if (is.null(setname)) {
      ntrees <- .C("GetNumberOfTrees4Exten",
                   as.character(rootfile),
                   as.character(treetype),
                   numtrees=integer(1),
                   PACKAGE="xps")$numtrees;
   } else {
      ntrees <- .C("GetNumberOfTrees",
                   as.character(rootfile),
                   as.character(setname),
                   as.character(treetype),
                   numtrees=integer(1),
                   PACKAGE="xps")$numtrees;
   }#if
   return(ntrees);
}#getNumberTrees

#------------------------------------------------------------------------------#
# getTreeNames: get names of trees saved in root file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getTreeNames <- function(rootfile, treetype="*", setname=NULL, gettitle=FALSE) {
   if (debug.xps()) print("------getTreeNames------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## get tree names
   if (is.null(setname)) {
      ntrees <- .C("GetNumberOfTrees4Exten",
                   as.character(rootfile),
                   as.character(treetype),
                   numtrees=integer(1),
                   PACKAGE="xps")$numtrees;

      tnames <- .C("GetTreeNames4Exten",
                   as.character(rootfile),
                   as.character(treetype),
                   as.integer(gettitle),
                   treenames=character(ntrees),
                   PACKAGE="xps")$treenames;
   } else {
      ntrees <- .C("GetNumberOfTrees",
                   as.character(rootfile),
                   as.character(setname),
                   as.character(treetype),
                   numtrees=integer(1),
                   PACKAGE="xps")$numtrees;

      tnames <- .C("GetTreeNames",
                   as.character(rootfile),
                   as.character(setname),
                   as.character(treetype),
                   as.integer(gettitle),
                   treenames=character(ntrees),
                   PACKAGE="xps")$treenames;
   }#if

#   return(tnames);
   return(as.vector(tnames, mode="character"));
}#getTreeNames

#------------------------------------------------------------------------------#

