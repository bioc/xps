#==============================================================================#
# TreeSetClasses.R: contains all class definitions and method definitions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# methods.ProjectInfo.R:   contains methods for ProjectInfo
# methods.TreeSet.R:       contains methods for TreeSet
# methods.SchemeTreeSet.R: contains methods for SchemeTreeSet
# methods.ProcesSet.R:     contains methods for ProcesSet
# methods.DataTreeSet.R:   contains methods for DataTreeSet
# methods.ExprTreeSet.R:   contains methods for ExprTreeSet
# methods.CallTreeSet.R:   contains methods for CallTreeSet
#==============================================================================#


#------------------------------------------------------------------------------#
# ProjectInfo:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# class ProjectInfo
setClass("ProjectInfo",
   representation(submitter      = "character",
                  laboratory     = "character",
                  contact        = "character",
                  project        = "list",
                  author         = "list",
                  dataset        = "list",
                  source         = "list",
                  sample         = "list",
                  celline        = "list",
                  primarycell    = "list",
                  tissue         = "list",
                  biopsy         = "list",
                  arraytype      = "list",
                  hybridizations = "data.frame",
                  treatments     = "data.frame"
   ),
   prototype(submitter      = "",
             laboratory     = "",
             contact        = "",
             project        = list(),
             author         = list(),
             dataset        = list(),
             source         = list(),
             sample         = list(),
             celline        = list(),
             primarycell    = list(),
             tissue         = list(),
             biopsy         = list(),
             arraytype      = list(),
             hybridizations = data.frame(matrix(nr=0,nc=0)),
             treatments     = data.frame(matrix(nr=0,nc=0))
   )
)#ProjectInfo

# generic methods for class ProjectInfo
setGeneric("projectInfo",     function(object) standardGeneric("projectInfo"));
setGeneric("projectInfo<-",   function(object, value) standardGeneric("projectInfo<-"));
setGeneric("authorInfo",      function(object) standardGeneric("authorInfo"));
setGeneric("authorInfo<-",    function(object, value) standardGeneric("authorInfo<-"));
setGeneric("datasetInfo",     function(object) standardGeneric("datasetInfo"));
setGeneric("datasetInfo<-",   function(object, value) standardGeneric("datasetInfo<-"));
setGeneric("sourceInfo",      function(object) standardGeneric("sourceInfo"));
setGeneric("sourceInfo<-",    function(object, value) standardGeneric("sourceInfo<-"));
setGeneric("sampleInfo",      function(object) standardGeneric("sampleInfo"));
setGeneric("sampleInfo<-",    function(object, value) standardGeneric("sampleInfo<-"));
setGeneric("cellineInfo",     function(object) standardGeneric("cellineInfo"));
setGeneric("cellineInfo<-",   function(object, value) standardGeneric("cellineInfo<-"));
setGeneric("primcellInfo",    function(object) standardGeneric("primcellInfo"));
setGeneric("primcellInfo<-",  function(object, value) standardGeneric("primcellInfo<-"));
setGeneric("tissueInfo",      function(object) standardGeneric("tissueInfo"));
setGeneric("tissueInfo<-",    function(object, value) standardGeneric("tissueInfo<-"));
setGeneric("biopsyInfo",      function(object) standardGeneric("biopsyInfo"));
setGeneric("biopsyInfo<-",    function(object, value) standardGeneric("biopsyInfo<-"));
setGeneric("arrayInfo",       function(object) standardGeneric("arrayInfo"));
setGeneric("arrayInfo<-",     function(object, value) standardGeneric("arrayInfo<-"));
setGeneric("hybridizInfo",    function(object) standardGeneric("hybridizInfo"));
setGeneric("hybridizInfo<-",  function(object, value) standardGeneric("hybridizInfo<-"));
setGeneric("treatmentInfo",   function(object) standardGeneric("treatmentInfo"));
setGeneric("treatmentInfo<-", function(object, value) standardGeneric("treatmentInfo<-"));


#------------------------------------------------------------------------------#
# TreeSet: virtual superset for 'SchemeTreeSet', 'DataTreeSet' etc
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# class TreeSet
setClass("TreeSet",
   representation(setname   = "character",
                  settype   = "character",
                  rootfile  = "character",
                  filedir   = "character",
                  numtrees  = "numeric",
                  treenames = "list",
                  "VIRTUAL"
   ),
   prototype(setname   = "",
             settype   = "",
             rootfile  = "ROOTFile",
             filedir   = getwd(),
             numtrees  = 0,
             treenames = list()
   )
)#TreeSet

# generic methods for class TreeSet
setGeneric("rootFile",     function(object) standardGeneric("rootFile"));
setGeneric("rootFile<-",   function(object, value) standardGeneric("rootFile<-"));
setGeneric("fileDir",      function(object) standardGeneric("fileDir"));
setGeneric("fileDir<-",    function(object, value) standardGeneric("fileDir<-"));
setGeneric("setName",      function(object) standardGeneric("setName"));
setGeneric("setName<-",    function(object, value) standardGeneric("setName<-"));
setGeneric("setType",      function(object) standardGeneric("setType"));
setGeneric("setType<-",    function(object, value) standardGeneric("setType<-"));
setGeneric("treeNames",    function(object) standardGeneric("treeNames"));
setGeneric("export",       function(object,...) standardGeneric("export"));
setGeneric("root.browser", function(object) standardGeneric("root.browser"));


#------------------------------------------------------------------------------#
# SchemeTreeSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("SchemeTreeSet",
   representation(chipname  = "character",
                  chiptype  = "character",
                  probeinfo = "list",
                  mask      = "data.frame"
   ),
   contains=c("TreeSet"),
   prototype(chipname  = "",
             chiptype  = "GeneChip",
             probeinfo = list(),
             mask      = data.frame(matrix(nr=0,nc=0))
   )
)#SchemeTreeSet

# generic methods for class SchemeTreeSet
setGeneric("chipName",   function(object) standardGeneric("chipName"));
setGeneric("chipType",   function(object) standardGeneric("chipType"));
setGeneric("chipType<-", function(object, value) standardGeneric("chipType<-"));
setGeneric("chipMask",   function(object) standardGeneric("chipMask"));
setGeneric("chipMask<-", function(object, value) standardGeneric("chipMask<-"));
setGeneric("probeInfo",  function(object) standardGeneric("probeInfo"));
setGeneric("nrows",      function(object) standardGeneric("nrows"));
setGeneric("ncols",      function(object) standardGeneric("ncols"));
setGeneric("attachMask", function(object) standardGeneric("attachMask"));
setGeneric("removeMask", function(object) standardGeneric("removeMask"));


#------------------------------------------------------------------------------#
# ProcesSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("ProcesSet",
   representation(scheme = "SchemeTreeSet",
                  data   = "data.frame",
                  params = "list"
   ),
   contains=c("TreeSet"),
   prototype(scheme = new("SchemeTreeSet"),
             data   = data.frame(matrix(nr=0,nc=0)),
             params = list()
   )
)#ProcesSet

# generic methods for class ProcesSet
setGeneric("schemeFile",  function(object) standardGeneric("schemeFile"));
setGeneric("schemeFile<-",function(object, value) standardGeneric("schemeFile<-"));
setGeneric("schemeSet",   function(object) standardGeneric("schemeSet"));
setGeneric("schemeSet<-", function(object, value) standardGeneric("schemeSet<-"));
setGeneric("getTreeData", function(object, ...) standardGeneric("getTreeData"));
setGeneric("validData",   function(object, ...) standardGeneric("validData"));
setGeneric("boxplot",     function(x, ...) standardGeneric("boxplot"));
setGeneric("mboxplot",    function(x, ...) standardGeneric("mboxplot"));
setGeneric("hist",        function(x, ...) standardGeneric("hist"));


#------------------------------------------------------------------------------#
# DataTreeSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("DataTreeSet",
   representation(bgtreenames = "list",
                  bgrd        = "data.frame",
                  projectinfo = "ProjectInfo"
   ),
   contains=c("ProcesSet"),
   prototype(bgtreenames = list(),
             bgrd        = data.frame(matrix(nr=0,nc=0)),
             projectinfo = new("ProjectInfo")
   )
)#DataTreeSet

# generic methods for class DataTreeSet
setGeneric("intensity",     function(object) standardGeneric("intensity"));
setGeneric("intensity<-",   function(object, value) standardGeneric("intensity<-"));
setGeneric("background",    function(object) standardGeneric("background"));
setGeneric("background<-",  function(object, value) standardGeneric("background<-"));
setGeneric("bgtreeNames",   function(object) standardGeneric("bgtreeNames"));
setGeneric("attachInten",   function(object, ...) standardGeneric("attachInten"));
setGeneric("removeInten",   function(object) standardGeneric("removeInten"));
setGeneric("attachBgrd",    function(object, ...) standardGeneric("attachBgrd"));
setGeneric("removeBgrd",    function(object) standardGeneric("removeBgrd"));
setGeneric("validBgrd",     function(object, ...) standardGeneric("validBgrd"));
setGeneric("addData",       function(object, ...) standardGeneric("addData"));
setGeneric("pm",            function(object, ...) standardGeneric("pm"));
setGeneric("mm",            function(object, ...) standardGeneric("mm"));
setGeneric("xpsRMA",        function(object, ...) standardGeneric("xpsRMA"));
setGeneric("xpsMAS4",       function(object, ...) standardGeneric("xpsMAS4"));
setGeneric("xpsMAS5",       function(object, ...) standardGeneric("xpsMAS5"));
setGeneric("xpsMAS5Call",   function(object, ...) standardGeneric("xpsMAS5Call"));
setGeneric("xpsDABGCall",   function(object, ...) standardGeneric("xpsDABGCall"));
setGeneric("xpsPreprocess", function(object, ...) standardGeneric("xpsPreprocess"));
setGeneric("xpsBgCorrect",  function(object, ...) standardGeneric("xpsBgCorrect"));
setGeneric("xpsNormalize",  function(object, ...) standardGeneric("xpsNormalize"));
setGeneric("xpsSummarize",  function(object, ...) standardGeneric("xpsSummarize"));
setGeneric("pmplot",        function(x, ...) standardGeneric("pmplot"));
setGeneric("image",         function(x, ...) standardGeneric("image"));


#------------------------------------------------------------------------------#
# ExprTreeSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("ExprTreeSet",
   representation(exprtype = "character",
                  normtype = "character"
   ),
   contains=c("ProcesSet"),
   prototype(exprtype = "none",
             normtype = "none"
   )
)#ExprTreeSet

# generic methods for class ExprTreeSet
setGeneric("exprType",   function(object) standardGeneric("exprType"));
setGeneric("exprType<-", function(object, value) standardGeneric("exprType<-"));
setGeneric("normType",   function(object) standardGeneric("normType"));
setGeneric("normType<-", function(object, value) standardGeneric("normType<-"));
setGeneric("exprs",      function(object) standardGeneric("exprs"));
setGeneric("exprs<-",    function(object, value) standardGeneric("exprs<-"));
setGeneric("se.exprs",   function(object) standardGeneric("se.exprs"));
setGeneric("mvaplot",    function(x, ...) standardGeneric("mvaplot"));


#------------------------------------------------------------------------------#
# CallTreeSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("CallTreeSet",
   representation(calltype = "character",
                  detcall  = "data.frame"
   ),
   contains=c("ProcesSet"),
   prototype(calltype = "mas5",
             detcall  = data.frame(matrix(nr=0,nc=0))
   )
)#CallTreeSet

# generic methods for class CallTreeSet
setGeneric("pvalData",   function(object) standardGeneric("pvalData"));
setGeneric("pvalData<-", function(object, value) standardGeneric("pvalData<-"));
setGeneric("presCall",   function(object) standardGeneric("presCall"));
setGeneric("presCall<-", function(object, value) standardGeneric("presCall<-"));
setGeneric("validCall",  function(object) standardGeneric("validCall"));
setGeneric("callplot",   function(x, ...) standardGeneric("callplot"));



#------------------------------------------------------------------------------#


