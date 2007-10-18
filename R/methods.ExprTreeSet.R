#==============================================================================#
# methods.ExprTreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

