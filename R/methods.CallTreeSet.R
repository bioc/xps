#==============================================================================#
# methods.CallTreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# pvalData:
# pvalData<-:
# presCall:
# presCall<-:
# validCall:
# callplot:
#==============================================================================#


#------------------------------------------------------------------------------#
# CallTreeSet initialization:
#------------------------------------------------------------------------------#

"initialize.CallTreeSet" <-
function(.Object,
         calltype = "none",
         ...) 
{
   if (debug.xps()) print("------initialize:CallTreeSet------")

   ## set default calltype
   if (calltype == "") {
      calltype <- "mas5";
   }#if

   .Object <- callNextMethod(.Object, calltype=calltype, ...);
   .Object@calltype = calltype;
   .Object;
}#initialize.CallTreeSet

setMethod("initialize", "CallTreeSet", initialize.CallTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("CallTreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:CallTreeSet------")

      msg <- NULL;

      ## check for correct settype
      if (object@settype != "preprocess") {
         msg <- validMsg(msg, paste(sQuote("settype"),
                        "must be of type <preprocess>"));
      }#if

      ## check calltype
      TYPE <- c("mas5", "dabg");
      if (is.na(match(object@calltype, TYPE))) {
         msg <- validMsg(msg, paste(sQuote("normation"), "must be one of",
                         "<mas5,dabg>"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# CallTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("pvalData", signature(object="CallTreeSet"),
   function(object) object@data
)#pvalData

setReplaceMethod("pvalData", signature(object="CallTreeSet", value="data.frame"),
   function(object, treenames = NULL, value) {
      ## keep columns "UNIT_ID" and "UnitName"
      ds <- value[, c("UNIT_ID", "UnitName")];

      ## remove columns "UNIT_ID" and "UnitName"
      value <- value[, is.na(match(colnames(value), "UNIT_ID"))];
      value <- value[, is.na(match(colnames(value), "UnitName"))];

      ## result dataframe ds
      if (is.null(treenames)) {
         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      } else {
         treenames <- namePart(treenames);
         datanames <- namePart(colnames(value));

         pos <- match(datanames, treenames);
         if (length(pos[!is.na(pos)]) != length(treenames)) {
            stop("at least one treename is not present in value");
         }#if

         value <- value[,!is.na(match(datanames, treenames))];

         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      }#if

      ## correct extension
      treenames <- colnames(value);
      treenames <- sub("_PVALUE", "", treenames);

      object@data      <- ds;
      object@treenames <- as.list(treenames);
      object@numtrees  <- length(treenames);
      return(object);
   }
)#pvalData<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("presCall", signature(object="CallTreeSet"),
   function(object) object@detcall
)#presCall

setReplaceMethod("presCall", signature(object="CallTreeSet", value="data.frame"),
   function(object, treenames = NULL, value) {
      ## keep columns "UNIT_ID" and "UnitName"
      ds <- value[, c("UNIT_ID", "UnitName")];

      ## remove columns "UNIT_ID" and "UnitName"
      value <- value[, is.na(match(colnames(value), "UNIT_ID"))];
      value <- value[, is.na(match(colnames(value), "UnitName"))];

      ## result dataframe ds
      if (is.null(treenames)) {
         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      } else {
         treenames <- namePart(treenames);
         datanames <- namePart(colnames(value));

         pos <- match(datanames, treenames);
         if (length(pos[!is.na(pos)]) != length(treenames)) {
            stop("at least one treename is not present in value");
         }#if

         value <- value[,!is.na(match(datanames, treenames))];

         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      }#if

      object@detcall <- ds;
      return(object);
   }
)#presCall<-


#------------------------------------------------------------------------------#
# CallTreeSet methods:
#------------------------------------------------------------------------------#

"callCallTreeSet" <-
function(object) {
   if (debug.xps()) print("------callCallTreeSet------")

   ## check for presence of detcall
   dcall <- object@detcall;
   ntree <- object@numtrees;
   if (min(dim(dcall)) == 0 || ncol(dcall) < ntree) {
      stop(paste("slot", sQuote("detcall"), "has no data"));
   }#if

   treenames <- namePart(object@treenames);
#x   treenames <- make.names(treenames);  #to compare names with colnames of data.frame
   callnames <- namePart(colnames(dcall));

   return(dcall[,!is.na(match(callnames, treenames))]);
}#callCallTreeSet

setMethod("validCall", "CallTreeSet", callCallTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("callplot", signature(x="CallTreeSet"),
   function(x,
            beside = TRUE,
            names  = "namepart",
            col    = c("red","green","blue"),
            legend = c("P","M","A"),
            las    = 2,
            ylim   = c(0,100),
            ylab   = "detection call [%]",
            ...) 
   {
      if (debug.xps()) print("------callplot.CallTreeSet------")

      dc <- validCall(x);

      if (is.null(names))              names <- colnames(dc)
      else if (names[1] == "namepart") names <- namePart(colnames(dc))
      else                             dc    <- dc[, names, drop=F];

      P  <- 100*apply(dc, 2, function(x)length(which(x=="P")))/nrow(dc);
      M  <- 100*apply(dc, 2, function(x)length(which(x=="M")))/nrow(dc);
      A  <- 100*apply(dc, 2, function(x)length(which(x=="A")))/nrow(dc);

      barplot(t(cbind(P, M, A)),
              beside    = beside,
              col       = col,
              names.arg = names,
              legend    = legend,
              las       = las,
              ylim      = ylim,
              ylab      = ylab,
              ...);
   }
)#callplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
