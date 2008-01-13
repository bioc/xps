#==============================================================================#
# methods.ExprTreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# exprType:
# exprType<-:
# normType:
# normType<-:
# exprs:
# exprs<-:
# se.exprs:
# xpsNormalize:
# xpsPreFilter:
# xpsUniFilter:
# mvaplot:
#==============================================================================#


#------------------------------------------------------------------------------#
# ExprTreeSet initialization:
#------------------------------------------------------------------------------#

"initialize.ExprTreeSet" <-
function(.Object,
         exprtype = "none",
         normtype = "none",
         ...) 
{
   if (debug.xps()) print("------initialize:ExprTreeSet------")

   ## set default expression/normation type
   if (exprtype == "") exprtype <- "none";
   if (normtype == "") normtype <- "none";

   .Object <- callNextMethod(.Object,
                             exprtype = exprtype,
                             normtype = normtype,
                             ...);
   .Object@exprtype = exprtype;
   .Object@normtype = normtype;
   .Object;
}#initialize.ExprTreeSet

setMethod("initialize", "ExprTreeSet", initialize.ExprTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("ExprTreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:ExprTreeSet------")

      msg <- NULL;

      ## check for correct settype
      if (!(object@settype == "preprocess" ||
            object@settype == "normation")) {
         msg <- validMsg(msg, paste(sQuote("settype"),
                        "must be of type <preprocess, normation>"));
      }#if

      ## check expression method
      TYPE <- c("none", "rma", "mas4", "mas5", "custom");
      if (is.na(match(object@exprtype, TYPE))) {
         msg <- validMsg(msg, paste(sQuote("exprtype"), "must be one of",
                         "<none,rma,mas4,mas5,custom>"));
      }#if

      ## check normation method
      TYPE <- c("none", "mean", "median", "lowess", "supsmu");
      if (is.na(match(object@normtype, TYPE))) {
         msg <- validMsg(msg, paste(sQuote("normtype"), "must be one of",
                         "<none,mean,median,lowess,supsmu>"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# ExprTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("exprType", signature(object="ExprTreeSet"),
   function(object) object@exprtype
)#exprType

setReplaceMethod("exprType", signature(object="ExprTreeSet", value="character"),
   function(object, value) {
      object@exprtype <- value;
      return(object);
   }
)#exprType<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("normType", signature(object="ExprTreeSet"),
   function(object) object@normtype
)#normType

setReplaceMethod("normType", signature(object="ExprTreeSet", value="character"),
   function(object, value) {
      object@normtype <- value;
      return(object);
   }
)#normType<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("exprs", signature(object="ExprTreeSet"),
   function(object) object@data
)#exprs

setReplaceMethod("exprs", signature(object="ExprTreeSet", value="data.frame"),
   function(object, value) {
      object@data <- value;
      return(object);
   }
)#exprs<-


#------------------------------------------------------------------------------#
# ExprTreeSet methods:
#------------------------------------------------------------------------------#

setMethod("se.exprs", signature(object="ExprTreeSet"),
   function(object) {
      if (debug.xps()) print("------se.exprs.ExprTreeSet------")

      se <- export(object,
                   treetype     = extenPart(object@treenames),
                   varlist      = "fUnitName:fStdev",
                   as.dataframe = TRUE,
                   verbose      = FALSE);
      return(se);
   }
)#se.exprs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"normalize.ExprTreeSet" <-
function(object,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "separate",
         method    = "mean",
         option    = "all",
         logbase   = "0",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(), 
         verbose   = TRUE)
{
   if (debug.xps()) print("------normalize.ExprTreeSet------")

   ## get schemefile and chiptype for object
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chiptype   <- chipType(object);

   ## treenames to normalize as fullnames=/datadir/treenames
   datafile  <- object@rootfile;
   setname   <- object@setname;
   treenames <- as.character(object@treenames);
   numtrees  <- length(treenames);
   fullnames <- paste(datafile, setname, treenames, sep="/");

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## update: need to set filename to rootfile
   if (!is.logical(update)) {
      stop(paste(sQuote("update"), "must be TRUE or FALSE"));
   } else if (update == TRUE) {
      filename <- rootfile;
   }#if

   ## check for presence of valid selector option
   TYPE <- c("separate", "together");
   if (is.na(match(select, TYPE))) {
      stop(paste(sQuote(select), "is not a valid selector option"));
   }#if

   ## check for valid normalization method
   TYPE <- c("mean", "median", "lowess", "supsmu");
   if (is.na(match(method, TYPE))) {
      stop(paste(sQuote("method"), "must be <mean, median,lowess,supsmu>"));
   }#if

   ## percent units to be selected
   if (option == "all") {
      pc <- 1.0;
   } else if (option == "sel") {
      pc <- 0.3; #default
   } else {
      pc <- as.integer(as.character(option));
   }#if
   if (!(pc > 0 && pc <= 1)) {
      stop(paste(sQuote("option"), "must be <all,sel> or number in range (0,1]"));
   }#if

   ## check parameter list
   if (length(params) == 0) {
      stop(paste("empty parameter list", sQuote("params")));
   }#if

   ## option, tree extension and parameters for normalization method
   logbase <- validLogbase(logbase);
   switch(method,
      mean = {
         if (!(length(params) == 2 &&
               params[[1]] >= 0 && params[[1]] <= 0.5 &&
               params[[2]] >= 0)) {
            warning(paste(sQuote("params"),
                    "are not correct, using default parameters <0.02,0>"));
            params[[1]] <- 0.02;
            params[[2]] <- 0;
         }#if
         params <- list(trim=params[[1]], sc=params[[2]]);
         option <- paste(ifelse(pc<1, "sel", "all"), logbase, sep=":");
         exten  <- "tmn";
      }, median = {
         if (length(params) == 1) {
            if (params[[1]] < 0) {
               warning(paste(sQuote("params"),
                       "is not correct, using default setting sc=0"));
               params[[1]] <- 0;
            }#if
            params <- list(trim=0.5, sc=params[[1]]);
         } else if (!(length(params) == 2 &&
                    params[[1]] == 0.5 && params[[2]] >= 0)) {
            warning(paste(sQuote("params"),
                    "are not correct, using default parameters <0.5,0>"));
            params[[1]] <- 0.5;
            params[[2]] <- 0;
         }#if
         params <- list(trim=params[[1]], sc=params[[2]]);
         method <- "mean";
         option <- paste(ifelse(pc<1, "sel", "all"), logbase, sep=":");
         exten  <- "med";
      }, lowess = {
         if (!(length(params) == 2 &&
               params[[1]] > 0 && params[[1]] <= 1 &&
               params[[2]] > 0)) {
            warning(paste(sQuote("params"),
                          "are not correct, using default parameters <0.67,3>"));
            params[[1]] <- 0.67;
            params[[2]] <- 3;
         }#if
         params <- list(f=params[[1]], iter=params[[2]]);
         option <- logbase;
         exten  <- "low";
      }, supsmu = {
         if (!(length(params) == 2 &&
               params[[1]] >= 0 && params[[1]] <= 1 &&
               params[[2]] >= 0)) {
            warning(paste(sQuote("params"),
                          "are not correct, using default parameters <0.,0.>"));
            params[[1]] <- 0.0;
            params[[2]] <- 0.0;
         }#if
         params <- list(span=params[[1]], bass=params[[2]]);
         option <- logbase;
         exten  <- "sup";
      }
   )#switch

   ## exon level
   if (exonlevel != "") {
      cat("note:", sQuote("exonlevel"), "has no effect on ExprTreeSet\n");
   }#if
   exonlevel <- 0;

   ## reference tree
   reftree <- "";
   if (!(refindex >= 0 && refindex <= numtrees)) {
      stop(paste(sQuote("refindex"), "must be 0 or less than number of trees"));
   }#if
   if (refindex == 0) {
      xten    <- unlist(strsplit(treenames[[1]], "\\."));
      reftree <- paste(setname, "/*.", xten[2], sep="");
   } else {
      reftree <- paste(setname, treenames[[refindex]], sep="/");
   }#if

   ## reference method
   TYPE <- c("mean", "median");
   if (is.na(match(refmethod, TYPE))) {
      stop(paste(sQuote("refmethod"), "must be <mean, median>"));
   }#if

   ## define setname and settype for new treeset
   setname <- "NormSet";
   settype <- "normation";

   ## normalize
   r <- .C("Normxpress",
           as.character(filename),
           as.character(filedir),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.character(select),
           as.double(pc),
           as.character(method),
           as.character(option),
           as.integer(length(params)),
           as.double(params),
           as.integer(exonlevel),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.character(reftree),
           as.character(refmethod),
           as.integer(update),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("Normxpress")));
      return(NULL);
   }#if

   ## export/import result as dataframe
   outfile  <- sub(".root", ".txt", rootfile);
   # get treename "treeset.treename.treetype"
   treename <- paste(setname, "*", exten, sep=".");
   numtrees <- 1; # must be one for treename="*"

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treename),
           as.integer(numtrees),
           as.character(exten),
           as.character("fUnitName:fLevel"),
           as.character(outfile),
           as.character("\t"),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## get treenames after normalization
   treenames <- getTreeNames(rootfile, exten);
   numtrees  <- length(treenames);

   ds <- data.frame(matrix(nr=0,nc=0));
   if (file.exists(outfile)) {
      ds <- read.table(outfile, header=TRUE, sep="\t", row.names=NULL);
   } else {
      warning(paste("could not export results as", sQuote(outfile)));
   }#if

   ## create new class ExprTreeSet
   set <- new("ExprTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              data      = ds,
              params    = params,
              exprtype  = "none",
              normtype  = method);

   return(set);
}#normalize.ExprTreeSet

setMethod("xpsNormalize", "ExprTreeSet", normalize.ExprTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"prefilter.ExprTreeSet" <-
function(object,
         filename   = character(0),
         filedir    = getwd(),
         filter     = NULL,
         minfilters = 999,
         logbase    = "log2",
         treename   = "PreFilter",
         xps.call   = NULL,
         verbose    = TRUE)
{
   if (debug.xps()) print("------prefilter.ExprTreeSet------")

   ## get schemefile and chiptype for object
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chiptype   <- chipType(object);
   chipname   <- chipName(scheme);

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check for prefilter
   if (is.null(filter) || !is(filter, "PreFilter")) {
      stop(paste("missing filter of class", sQuote("PreFilter")));
   }#if

   ## check for number of filters
   if (filter@numfilters <= 0) {
      stop(paste(sQuote("PreFilter"), "must have at least one initialized filter"));
   }#if

   ## check for minfilters
   if (minfilters <= 0) {
      stop(paste(sQuote("minfilters"), "must be at least one"));
   }#if

   ## check for valid logbase
   logbase <- validLogbase(logbase);

   ## get exprtree names to add to filter set
   datafile  <- object@rootfile;
   setname   <- object@setname;
   treenames <- as.character(object@treenames);
   numexprs  <- length(treenames);
   exprnames <- paste(datafile, setname, treenames, sep="/");

   ## check for presence of optional root detection-call file
   numcalls  <- 0;
   callnames <- "";
   if (!is.null(xps.call)) {
      if (!(is(xps.call, "CallTreeSet") &&
            (schemeFile(xps.call) == schemefile) &&
            (chipType(xps.call)   == chiptype) &&
            (file.exists(rootFile(xps.call))))) {
         stop(paste("wrong detection call file", sQuote("xps.call")));
      }#if

      ## get calltree names to add to filter set
      datafile  <- rootFile(xps.call);
      setname   <- setName(xps.call);
      treenames <- treeNames(xps.call);
      numcalls  <- length(treenames);
      callnames <- paste(datafile, setname, treenames, sep="/");
   }#if

   ## prefilter data
   nvarparams  <- 0; varparams  <- 0; varoption  = "";
   nlowparams  <- 0; lowparams  <- 0; lowoption  = "";
   nhiparams   <- 0; hiparams   <- 0; hioption   = "";
   nquanparams <- 0; quanparams <- 0; quanoption = "";
   ncallparams <- 0; callparams <- 0; calloption = "";
   trim    <- NULL;
   epsilon <- NULL;

   ## variance filter
   madev      <- madFilter(filter);
   cv         <- cvFilter(filter);
   variance   <- varFilter(filter);
   difference <- diffFilter(filter);
   ratio      <- ratioFilter(filter);
   gap        <- gapFilter(filter);
   if (length(madev) == 2) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- madev[1];
      epsilon    <- madev[2];
      varoption  <- paste(varoption, "mad", sep=":", collapse=":");
   }#if
   if (length(cv) == 3) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- cv[1];
      trim       <- cv[2];
      epsilon    <- cv[3];
      varoption  <- paste(varoption, "cov2mn", sep=":", collapse=":");
   }#if
   if (length(variance) == 3) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- variance[1];
      trim       <- variance[2];
      epsilon    <- variance[3];
      varoption  <- paste(varoption, "var2mn", sep=":", collapse=":");
   }#if
   if (length(difference) == 3) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- difference[1];
      trim       <- difference[2];
      epsilon    <- difference[3];
      varoption  <- paste(varoption, "dif2mn", sep=":", collapse=":");
   }#if
   if (length(ratio) == 1) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- ratio[1];
      varoption  <- paste(varoption, "max2min", sep=":", collapse=":");
   }#if
   if (length(gap) == 4) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- gap[1];
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- gap[2];
      trim       <- gap[3];
      epsilon    <- gap[4];
      varoption  <- paste(varoption, "gap2mn", sep=":", collapse=":");
   }#if
   if (!is.null(trim)) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams] <- trim;
   }#if
   if (!is.null(epsilon)) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams] <- epsilon;
   }#if
   varoption <- substring(varoption, 2);

   ## lowFilter
   lothreshold <- lowFilter(filter);
   if (length(lothreshold) == 3) {
      nlowparams <- 2;
      lowparams  <- lothreshold[1:2];
      lowoption  <- lothreshold[3];
   }#if

   ## highFilter
   hithreshold <- highFilter(filter);
   if (length(hithreshold) == 3) {
      nhiparams <- 2;
      hiparams  <- hithreshold[1:2];
      hioption  <- hithreshold[3];
   }#if

   ## quantileFilter
   quantil <- quantileFilter(filter);
   if (length(quantil) == 3) {
      nquanparams <- 3;
      quanparams  <- quantil;
   }#if

   ## callFilter
   prescall <- callFilter(filter);
   if (length(prescall) == 3) {
      ncallparams <- 2;
      callparams  <- prescall[1:2];
      calloption  <- prescall[3];
   }#if

   ## define setname and settype for new treeset
   setname <- "PreFilterSet";
   settype <- "prefilter";
   exten   <- "pfr";

   ## prefilter
   r <- .C("PreFilter",
           as.character(filename),
           as.character(filedir),
           as.character(chiptype),
           as.character(chipname),
           as.character(schemefile),
           as.character(setname),
           as.character(treename),
           as.integer(minfilters),
           as.character(logbase),
           as.integer(nvarparams),
           as.double(varparams),
           as.character(varoption),
           as.integer(nlowparams),
           as.double(lowparams),
           as.character(lowoption),
           as.integer(nhiparams),
           as.double(hiparams),
           as.character(hioption),
           as.integer(nquanparams),
           as.double(quanparams),
           as.character(quanoption),
           as.integer(ncallparams),
           as.double(callparams),
           as.character(calloption),
           as.character(exprnames),
           as.integer(numexprs),
           as.character(callnames),
           as.integer(numcalls),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("PreFilter")));
      return(NULL);
   }#if

   ## export/import filter mask as dataframe
   outfile  <- sub(".root", ".txt", rootfile);
   # get treename "treeset.treename.treetype"
   treename <- paste(setname, treename, exten, sep=".");
   numtrees <- 1; 

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treename),
           as.integer(numtrees),
           as.character(exten),
           as.character("fUnitName:fFlag"),
           as.character(outfile),
           as.character("\t"),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## get treenames after normalization
   treenames <- getTreeNames(rootfile, exten);
   numtrees  <- length(treenames);

   ds <- data.frame(matrix(nr=0,nc=0));
   if (file.exists(outfile)) {
      ds <- read.table(outfile, header=TRUE, sep="\t", row.names=NULL);
   } else {
      warning(paste("could not export results as", sQuote(outfile)));
   }#if

   ## create new class FilterTreeSet
   set <- new("FilterTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              data      = ds,
#?              params    = params,
              filter    = filter,
              exprset   = object);

   ## add callset
   if (!is.null(xps.call) && is(xps.call, "CallTreeSet")) {
      set@callset <- xps.call;
   }#if

   return(set);
}#prefilter.ExprTreeSet

setMethod("xpsPreFilter", "ExprTreeSet", prefilter.ExprTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"unifilter.ExprTreeSet" <-
function(object,
         filename   = character(0),
         filedir    = getwd(),
         update     = FALSE,
         filter     = NULL,
         minfilters = 999,
         logbase    = "log2",
         group      = character(0),
         treename   = "UniTest",
         xps.fltr   = NULL,
         xps.call   = NULL,
         verbose    = TRUE)
{
   if (debug.xps()) print("------unifilter.ExprTreeSet------")

   ## get schemefile and chiptype for object
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chiptype   <- chipType(object);
   chipname   <- chipName(scheme);

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## update: need to set filename to rootfile
   if (!is.logical(update)) {
      stop(paste(sQuote("update"), "must be TRUE or FALSE"));
   } else if (update == TRUE) {
      filename <- rootfile;
   }#if

   ## check for unifilter
   if (is.null(filter) || !is(filter, "UniFilter")) {
      stop(paste("missing filter of class", sQuote("UniFilter")));
   }#if

   ## check for minfilters
   if (minfilters <= 0) {
      stop(paste(sQuote("minfilters"), "must be at least one"));
   }#if

   ## check for valid logbase
   logbase <- validLogbase(logbase);

   ## check for presence of optional root prefilter file
   fltrname <- "";
   numfltr  <- 0;
   if (!is.null(xps.fltr)) {
      if (setType(xps.fltr) != "prefilter") {
         stop(paste(sQuote("xps.fltr"), "must be of type", sQuote("PreFilter")));
      }#if
      if (chipName(xps.fltr@scheme) != chipname) {
         stop(paste(sQuote("xps.fltr"), "was not created for chip", sQuote(chipname)));
      }#if

      # get treename "rootfile.treeset.treename"
      rootname <- rootFile(xps.fltr);
      setname  <- setName(xps.fltr);
      fltrname <- treeNames(xps.fltr);
      fltrname <- paste(rootname, setname, fltrname, sep="/");
      numfltr  <- 1;
   }#if

   ## get exprtree names to add to filter set
   datafile  <- object@rootfile;
   setname   <- object@setname;
   treenames <- as.character(object@treenames);
   numexprs  <- length(treenames);
   exprnames <- paste(datafile, setname, treenames, sep="/");

   ## check for presence of optional root detection call file
   numcalls  <- 0;
   callnames <- "";
   if (!is.null(xps.call)) {
      if (!(is(xps.call, "CallTreeSet") &&
            (schemeFile(xps.call) == schemefile) &&
            (chipType(xps.call)   == chiptype) &&
            (file.exists(rootFile(xps.call))))) {
         stop(paste("wrong detection call file", sQuote("xps.call")));
      }#if

      ## get calltree names to add to filter set
      datafile  <- rootFile(xps.call);
      setname   <- setName(xps.call);
      treenames <- treeNames(xps.call);
      numcalls  <- length(treenames);
      if (numcalls != numexprs) {
         stop("number of call trees must be equal to number of expression trees");
      }#if
      callnames <- paste(datafile, setname, treenames, sep="/");
   }#if

   ## unifilter data
   nuniparams  <- 0; uniparams  <- 0; unioption  = "";
   nfcparams   <- 0; fcparams   <- 0; fcoption   = "";
   nufparams   <- 0; ufparams   <- 0; ufoption   = "";
   ncallparams <- 0; callparams <- 0; calloption = "";

   ## uniTest
   unitest <- uniTest(filter);
   if (length(unitest) == 8) {
      nuniparams <- 5;
      uniparams <- as.double(unitest[4:8]);
      unioption <- sub("\\.","",paste(unitest[2], unitest[3], sep=":"));
      unitest   <- sub("\\.","",unitest[1]);
   } else {
         stop(paste("missing or incorrect parameter", sQuote("UniFilter@unitest")));
   }#if

   ## fcFilter
   foldchange <- fcFilter(filter);
   if (length(foldchange) == 2) {
      nfcparams <- 2;
      fcparams[1] <- foldchange[1];
      fcparams[2] <- ifelse(foldchange[2]=="up",1, ifelse(foldchange[2]=="down",-1,0));
   }#if

   ## unitestFilter
   unifltr <- unitestFilter(filter);
   if (length(foldchange) == 2) {
      nufparams <- 1;
      ufparams <- unifltr[1];
      ufoption <- unifltr[2];
   }#if

   ## callFilter
   prescall <- callFilter(filter);
   if (length(prescall) == 3) {
      ncallparams <- 3;
      callparams[1] <- prescall[1];
      callparams[2] <- prescall[2];
      callparams[3] <- prescall[2];
      calloption    <- paste(prescall[3], prescall[3], sep=":");
   }#if

   ## check group and get group index
   if (length(unique(group)) != 2) {
         stop("number of groups must be two");
   }#if
   grpidx <- as.integer(as.factor(group));
   
   ## define setname and settype for unifilter
   setname <- "UniFilterSet";
   settype <- "unifilter";
   exten   <- "ufr";

   ## unifilter
   r <- .C("UniFilter",
           as.character(filename),
           as.character(filedir),
           as.character(chiptype),
           as.character(chipname),
           as.character(schemefile),
           as.character(setname),
           as.character(treename),
           as.integer(minfilters),
           as.character(logbase),
           as.character(unitest),
           as.integer(nuniparams),
           as.double(uniparams),
           as.character(unioption),
           as.integer(nfcparams),
           as.double(fcparams),
           as.character(fcoption),
           as.integer(nufparams),
           as.double(ufparams),
           as.character(ufoption),
           as.integer(ncallparams),
           as.double(callparams),
           as.character(calloption),
           as.character(exprnames),
           as.integer(numexprs),
           as.character(callnames),
           as.integer(numcalls),
           as.character(group),
           as.integer(grpidx),
           as.character(fltrname),
           as.integer(numfltr),
           as.character("*"),
           as.integer(update),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("UniFilter")));
      return(NULL);
   }#if

   ## export optional unifilter
   ds <- data.frame(matrix(nr=0,nc=0));
   if (numberFilters(filter) > 0) {
      ## export/import unifilter as dataframe
      outfile  <- sub(".root", ".txt", rootfile);
      # get treename "treeset.treename.treetype"
      ufltrname <- paste(setname, treename, exten, sep=".");
      numtrees  <- 1;

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(ufltrname),
              as.integer(numtrees),
              as.character(exten),
              as.character("fUnitName:fFlag"),
              as.character(outfile),
              as.character("\t"),
              as.integer(verbose),
              err=integer(1),
              PACKAGE="xps")$err;

      if (r != 0) {
         stop(paste("error in function", sQuote("ExportData")));
         return(NULL);
      }#if

      if (file.exists(outfile)) {
         ds <- read.table(outfile, header=TRUE, sep="\t", row.names=NULL);
      } else {
         warning(paste("could not export results as", sQuote(outfile)));
      }#if
   }#if

   ## get treenames after normalization
   treenames <- getTreeNames(rootfile, exten);
   numtrees  <- length(treenames);

   ## create new class FilterTreeSet
   fltrset <- new("FilterTreeSet",
                  setname   = setname,
                  settype   = settype,
                  rootfile  = rootfile,
                  filedir   = filedir,
                  numtrees  = numtrees,
                  treenames = as.list(treenames),
                  scheme    = scheme,
                  data      = ds,
#?                  params    = params,
                  filter    = filter,
                  exprset   = object);

   ## add callset to class FilterTreeSet
   if (!is.null(xps.call) && is(xps.call, "CallTreeSet")) {
      fltrset@callset <- xps.call;
   }#if

   ## varlist for export
#ev depending on uniparams, unioption
   varlist <- "fUnitName:df:mn1:mn2:fc:pval:padj:pcha:se:stat";

   ## define setname and settype for new treeset
   setname <- "UniFilterSet";
   settype <- "UnivariateAnalysis";
#??   settype <- "unifilter";
   exten   <- "stt";

   ## export/import unitest result as dataframe
   outfile  <- sub(".root", ".txt", rootfile);
   # get treename "treeset.treename.treetype"
   treename <- paste(setname, treename, exten, sep=".");
   numtrees <- 1;

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treename),
           as.integer(numtrees),
           as.character(exten),
           as.character(varlist),
           as.character(outfile),
           as.character("\t"),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## get treenames
   treenames <- getTreeNames(rootfile, exten);
   numtrees  <- length(treenames);

   ds <- data.frame(matrix(nr=0,nc=0));
   if (file.exists(outfile)) {
      ds <- read.table(outfile, header=TRUE, sep="\t", row.names=NULL);
   } else {
      warning(paste("could not export results as", sQuote(outfile)));
   }#if

   ## create new class AnalysisTreeSet
   set <- new("AnalysisTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              data      = ds,
#?              params    = params,
              fltrset   = fltrset);

   return(set);
}#unifilter.ExprTreeSet

setMethod("xpsUniFilter", "ExprTreeSet", unifilter.ExprTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("mvaplot", signature(x="ExprTreeSet"),
   function(x,
            transfo = log2,
            method  = "median",
            names   = "namepart",
            ylim    = c(-6,6),
            xlab    = "A",
            ylab    = "M",
            pch     = '.',
            ...) 
   {
      if (debug.xps()) print("------mvaplot.ExprTreeSet------")

      m <- validData(x);
      if (is.function(transfo)) m <- transfo(m);

      if (method == "median") mn  <- apply(m, 1, median)
      else                    mn  <- rowMeans(m);

      if (is.null(names))              names <- colnames(m)
      else if (names[1] == "namepart") names <- namePart(colnames(m))
      else                             m     <- m[, names, drop=F];

      ## plot images
      if (ncol(m) > 1 && interactive()) par(ask=TRUE) else par(ask=FALSE);
      for (i in 1:ncol(m)) {
         plot(x    = mn,
              y    = (m[,i] - mn),
              main = names[i],
              ylim = ylim,
              xlab = xlab,
              ylab = ylab,
              pch  = pch,
              ...);
         abline(h=0);
      }#for
      par(ask=FALSE);
   }
)#mvaplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

