#==============================================================================#
# methods.Filter.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# filters:
# filters<-:
#==============================================================================#


#------------------------------------------------------------------------------#
# Filter initialization:
#------------------------------------------------------------------------------#

"initialize.Filter" <-
function(.Object, ...) 
{
   if (debug.xps()) print("------initialize:Filter------")

   .Object <- callNextMethod(.Object, ...);
   .Object;
}#initialize.Filter

setMethod("initialize", "Filter", initialize.Filter);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("Filter",
   function(object) {
      if (debug.xps()) print("------setValidity:Filter------")
      msg <- NULL;

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# Filter accessors:
#------------------------------------------------------------------------------#

setMethod("numberFilters", signature(object="Filter"),
   function(object) object@numfilters
)#numberFilters

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

