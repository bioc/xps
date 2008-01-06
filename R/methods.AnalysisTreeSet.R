#==============================================================================#
# methods.AnalysisTreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# filterTreeset:
#==============================================================================#


#------------------------------------------------------------------------------#
# AnalysisTreeSet initialization:
#------------------------------------------------------------------------------#

"initialize.AnalysisTreeSet" <-
function(.Object, ...) 
{
   if (debug.xps()) print("------initialize:AnalysisTreeSet------")

   .Object <- callNextMethod(.Object, ...);
   .Object;
}#initialize.AnalysisTreeSet

setMethod("initialize", "AnalysisTreeSet", initialize.AnalysisTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("AnalysisTreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:AnalysisTreeSet------")
      msg <- NULL;

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# AnalysisTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("filterTreeset", signature(object="AnalysisTreeSet"),
   function(object) object@fltrset
)#filterTreeset

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("validFilter", signature(object="AnalysisTreeSet"),
   function(object) validData(object@fltrset)
)#validFilter


#------------------------------------------------------------------------------#
# AnalysisTreeSet methods:
#------------------------------------------------------------------------------#

setMethod("getTreeData", signature(object="AnalysisTreeSet"),
   function(object,
            treename = "UniTest",
            treetype = "stt",
            varlist  = "fUnitName:mn1:mn2:fc:pval")
   {
      if (debug.xps()) print("------getTreeData.AnalysisTreeSet------")

      ds <- export(object,
                   treenames    = treename,
                   treetype     = treetype,
                   varlist      = varlist,
                   as.dataframe = TRUE);
      return(ds);
   }
)#getTreeData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"dataAnalysisTreeSet" <-
function(object,
         which = "UnitName:mask") 
{
   if (debug.xps()) print("------dataAnalysisTreeSet------")

   ## check for presence of data
   data  <- object@data;
   if (min(dim(data)) == 0) {
      stop(paste("slot", sQuote("data"), "has no data"));
   }#if

   strg <- unlist(strsplit(which,":"));

   ## for which="mask" export only data with mask=1
   if (strg[1] == "mask" || (length(strg) == 2 && strg[2] == "mask")) {
      mask <- validFilter(object);
      data <- cbind(data, mask[,"FLAG",drop=F]);
      data <- data[data[,"FLAG"]==1,];
      data <- data[,-which(colnames(data) == "FLAG")];
   }#if

   ## use UnitName from column "which" as rownames
   if (!is.na(match(strg[1], colnames(data)))) {
      rownames(data) <- data[,strg[1]];
      return(data[,-which(colnames(data) == strg[1])]);
   }#if

   return(data);
}#dataAnalysisTreeSet

setMethod("validData", "AnalysisTreeSet", dataAnalysisTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
