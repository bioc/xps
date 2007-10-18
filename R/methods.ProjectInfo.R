#==============================================================================#
# methods.ProjectInfo.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# projectInfo:
# projectInfo<-:
# authorInfo:
# authorInfo<-:
# datasetInfo:
# datasetInfo<-:
# sourceInfo:
# sourceInfo<-:
# arrayInfo:
# arrayInfo<-:
# show:
#==============================================================================#


#------------------------------------------------------------------------------#
# ProjectInfo initialization:
#------------------------------------------------------------------------------#

"initialize.ProjectInfo" <-
function(.Object, ...) 
{
   if (debug.xps()) print("------initialize:ProjectInfo------")

#   callNextMethod();
   .Object <- callNextMethod(.Object, ...);
   .Object;
}#initialize.ProjectInfo

setMethod("initialize", "ProjectInfo", initialize.ProjectInfo);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("ProjectInfo",
   function(object) {
      if (debug.xps()) print("------setValidity:ProjectInfo------")
      msg <- NULL;

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# ProjectInfo accessors:
#------------------------------------------------------------------------------#

setMethod("projectInfo", signature(object="ProjectInfo"),
   function(object) object@project
)#projectInfo

setReplaceMethod("projectInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 3) {
         stop(paste(sQuote("project"),
              "must have <name,date,description>"));
      }#if

      object@project <- list(name        = value[1],
                             date        = value[2],
                             description = value[3]);
      return(object);
   }
)#projectInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("authorInfo", signature(object="ProjectInfo"),
   function(object) object@author
)#authorInfo

setReplaceMethod("authorInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 6) {
         stop(paste(sQuote("author"),
              "must have <lastname,firstname,company,lab,email,phone>"));
      }#if

      object@author <- list(lastname  = value[1],
                            firstname = value[2],
                            company   = value[3],
                            lab       = value[4],
                            email     = value[5],
                            phone     = value[6]);
      return(object);
   }
)#authorInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("datasetInfo", signature(object="ProjectInfo"),
   function(object) object@dataset
)#datasetInfo

setReplaceMethod("datasetInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 6) {
         stop(paste(sQuote("dataset"),
              "must have <name,type,sample,submitter,date,description>"));
      }#if

      TYPE <- c("UD", "TS", "DR", "MC");
      if (is.na(match(value[2], TYPE))) {
         stop(paste(sQuote("type"), "must be <UD,TS,DR,MC>"));
      }#if

      object@dataset <- list(name        = value[1],
                             type        = value[2],
                             sample      = value[3],
                             submitter   = value[4],
                             date        = value[5],
                             description = value[6]);
      return(object);
   }
)#datasetInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("sourceInfo", signature(object="ProjectInfo"),
   function(object) object@source
)#sourceInfo

setReplaceMethod("sourceInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 4) {
         stop(paste(sQuote("source"),
              "must have <name,species,subspecies,description>"));
      }#if

      object@source <- list(name        = value[1],
                            species     = value[2],
                            subspecies  = value[3],
                            description = value[4]);
      return(object);
   }
)#sourceInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("arrayInfo", signature(object="ProjectInfo"),
   function(object) object@arraytype
)#arrayInfo

setReplaceMethod("arrayInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 3) {
         stop(paste(sQuote("arraytype"),
              "must have <chipname,chiptype,description>"));
      }#if

      if (!isChipType(value[2])) {
         stop(paste(sQuote("arraytype"), "must be <GeneChip,ExonChip>"));
      }#if

      object@arraytype <- list(chipname    = value[1],
                               chiptype    = value[2],
                               description = value[3]);
      return(object);
   }
)#arrayInfo<-


#------------------------------------------------------------------------------#
# ProjectInfo methods:
#------------------------------------------------------------------------------#

setMethod("show", signature(object="ProjectInfo"),
   function(object) {
      if (debug.xps()) print("------show.ProjectInfo------")

      type  <- c("project", "author", "dataset", "source", "arraytype");
      index <- c(length(object@project)   > 0,
                 length(object@author)    > 0,
                 length(object@dataset)   > 0,
                 length(object@source)    > 0,
                 length(object@arraytype) > 0);

      cat("Project information:\n");
      cat("   Submitter: ", object@submitter,  "\n");
      cat("   Laboratory:", object@laboratory, "\n");
      cat("   Contact:   ", object@contact,    "\n");

      if (any(index)) {
         cat("   Information is available on:",
             paste(type[index],collapse=", "),"\n");
      }#if
   }
)#show

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

