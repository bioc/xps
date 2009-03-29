#==============================================================================#
# methods.DataTreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# nrows:
# ncols:
# intensity:
# intensity<-:
# background:
# background<-:
# bgtreeNames:
# attachInten:
# removeInten:
# attachBgrd:
# removeBgrd:
# attachMask:
# removeMask:
# validData:
# validBgrd:
# addData:
# pm:
# mm:
# rawCELName:
# xpsRMA:
# xpsMAS4:
# xpsMAS5:
# xpsMAS5Call:
# xpsDABGCall:
# xpsINICall:
# xpsPreprocess:
# xpsBgCorrect:
# xpsNormalize:
# xpsSummarize:
# pmplot:
# image:
#==============================================================================#


#------------------------------------------------------------------------------#
# DataTreeSet initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "DataTreeSet", 
   function(.Object, ...) {
      if (debug.xps()) print("------initialize:DataTreeSet------")

      .Object <- callNextMethod(.Object, ...);
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("DataTreeSet",
   function(object) 
   {
      if (debug.xps()) print("------setValidity:DataTreeSet------")

      msg <- NULL;

      ## check for correct settype
      TYPE <- c("rawdata", "preprocess");
      if (is.na(match(object@settype, TYPE))) {
         msg <- validMsg(msg,
                     paste(sQuote("settype"), "must be <rawdata,preprocess>"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# DataTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("nrows", signature(object="DataTreeSet"),
   function(object) nrows(object@scheme)
)#nrows

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("ncols", signature(object="DataTreeSet"),
   function(object) ncols(object@scheme)
)#ncols

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("intensity", signature(object="DataTreeSet"),
   function(object) object@data
)#intensity

setReplaceMethod("intensity", signature(object="DataTreeSet", value="data.frame"),
   function(object, filename = NULL, verbose = FALSE, value)  {
      ## check if dataframe value has same dimension as data
      if (!(length(rownames(value)) == 0 || length(rownames(object@data)) == 0) &&
          (length(rownames(value)) != length(rownames(object@data)))) {
         stop("replacement data.frame value must have same number of rows as data.");
      }#if

      ## check if dataframe value has columns X,Y
      if (dim(value)[2] != 0 && colnames(value)[1] != "X") {
         stop(paste("first column of value must be", sQuote("X")));
      }#if
      if (dim(value)[2] != 0 && colnames(value)[2] != "Y") {
         stop(paste("second column of value must be", sQuote("Y")));
      }#if

      ## sort for Y then X to ensure correct order
      if (dim(value)[2] > 0) {
         value <- value[order(value[,"Y"],value[,"X"]),];
      }#if

      ## replace data with value
      if (is.null(filename)) {
         object@data <- value;
         return(object);
      } else {
         if (dim(value)[1] == 0 || dim(value)[2] == 0) {
            warning(paste(sQuote("value"), "has dimension zero, object will not be changed."));
            return(object);
         }#if

         ## root data file
         rootdata <- paste(filename, "_cel.root", sep="");
         filedir  <- getwd();
         rootfile <- rootDirFile(rootdata, filedir);

         ## check if root file exists (necessary for WinXP to test already here)
         if (existsROOTFile(rootfile)) {
            stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
         }#if

         ## get names of CEL-files
         celnames <- namePart(colnames(value)[3:dim(value)[2]]);
         celfiles <- paste(getwd(), "/", celnames, ".CEL", sep="");

         ## create CEL-files (version 3)
         for (i in 1:length(celnames)) {
            celname <- paste(celnames[i], ".CEL", sep="");
            if (verbose) print(paste("exporting value as CEL-file:", celname))

            ## create CEL header
            cel.header <- CELHeader(celnames[i], object@scheme);

            ## create CEL data: STDV=3*sqrt(MEAN), NPIXELS=20
#x            cel.data <- value[,c(1,2,i+2)];
            cel.data <- cbind(value[,1:2], formatC(value[,i+2], digits=1, format="f"));
            cel.data <- cbind(cel.data, formatC(3*sqrt(as.numeric(cel.data[,3])), digits=1, format="f"), 20);
            colnames(cel.data) <- c("CellHeader=X", "Y", "MEAN", "STDV", "NPIXELS");

            write.table(cel.header, celname, quote=FALSE, col.names=FALSE, row.names=FALSE);
            suppressWarnings(write.table(cel.data, celname, quote=FALSE, sep="\t", append=TRUE, col.names=TRUE, row.names=FALSE));
         }#for

         ## create root data file
         r <- .C("ImportData",
                 as.character(filename),
                 as.character(filedir),
                 as.character(chipType(object@scheme)),
                 as.character(chipName(object@scheme)),
                 as.character(rootFile(object@scheme)),
                 as.character("DataSet"),
                 as.character(celfiles),
                 as.character(celnames),
                 as.integer(length(celfiles)),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.character(""),
                 as.integer(0),
                 as.integer(0),
                 as.integer(0),
                 as.integer(verbose),
                 result=character(2),
                 PACKAGE="xps")$result;

         ## returned result: saved rootfile and error
         rootfile <- r[1];
         error    <- as.integer(r[2]);

         if (error != 0) {
            stop(paste("error in function", sQuote("ImportData")));
            return(NULL);
         }#if

         ## replace object parameters
         object@data        <- value;
         object@rootfile    <- rootfile;
         object@filedir     <- filedir;
         object@treenames   <- as.list(paste(celnames, ".cel", sep=""));
         object@bgtreenames <- list();
         object@bgrd        <- data.frame(matrix(nr=0,nc=0));

         return(object);
      }#if
   }
)#intensity<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("background", signature(object="DataTreeSet"),
   function(object) object@bgrd
)#background

setReplaceMethod("background", signature(object="DataTreeSet", value="data.frame"),
   function(object, value) {
      object@bgrd <- value;
      return(object);
   }
)#background<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("bgtreeNames", signature(object="DataTreeSet"),
   function(object) object@bgtreenames
)#bgtreeNames

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("projectInfo", signature(object="DataTreeSet"),
   function(object) object@projectinfo
)#projectInfo

setReplaceMethod("projectInfo", signature(object="DataTreeSet", value="ProjectInfo"),
   function(object, value) {
      object@projectinfo <- value;
      return(object);
   }
)#projectInfo<-

#------------------------------------------------------------------------------#
# DataTreeSet methods:
#------------------------------------------------------------------------------#

setMethod("attachInten", signature(object="DataTreeSet"),
   function(object, treenames="*") {
      if (debug.xps()) print("------attachInten.DataTreeSet------")

      treetype <- extenPart(object@treenames);
      if (treenames[1] == "*") treenames <- object@treenames;
      if (length(treenames) > 0) {
         intensity(object, NULL, FALSE) <- export(object,
                                                  treenames    = treenames,
                                                  treetype     = treetype,
                                                  varlist      = "fInten",
                                                  as.dataframe = TRUE,
                                                  verbose      = FALSE);
      } else {
         warning("missing data tree names, data will not be added.");
      }#if
      return(object);
   }
)#attachInten

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeInten", signature(object="DataTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeInten.DataTreeSet------")

      intensity(object, NULL, FALSE) <- data.frame(matrix(nr=0,nc=0));
      gc(); #????
      return(object);
   }
)#removeInten

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("attachBgrd", signature(object="DataTreeSet"),
   function(object, treenames="*") {
      if (debug.xps()) print("------attachBgrd.DataTreeSet------")

      treetype <- extenPart(object@bgtreenames);
      if (treenames[1] == "*") treenames <- object@bgtreenames;
      if (length(treenames) > 0) {
         background(object) <- export(object,
                                      treenames    = treenames,
                                      treetype     = treetype,
                                      varlist      = "fBg",
                                      as.dataframe = TRUE,
                                      verbose      = FALSE);
      } else {
         warning("missing background tree names, data will not be added.");
      }#if
      return(object);
   }
)#attachBgrd

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeBgrd", signature(object="DataTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeBgrd.DataTreeSet------")

      object@bgrd <- data.frame(matrix(nr=0,nc=0));
      gc(); #????
      return(object);
   }
)#removeBgrd

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("attachMask", signature(object="DataTreeSet"),
   function(object) {
      if (debug.xps()) print("------attachMask.DataTreeSet------")

      chipMask(object@scheme) <- export(object@scheme,
                                        treetype     = "scm",
                                        varlist      = "fMask",
                                        as.dataframe = TRUE,
                                        verbose      = FALSE);
      return(object);
   }
)#attachMask

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeMask", signature(object="DataTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeMask.DataTreeSet------")

      chipMask(object@scheme) <- data.frame(matrix(nr=0,nc=0));
      gc(); #????
      return(object);
   }
)#removeMask

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"dataDataTreeSet" <-
function(object,
         which="") 
{
   if (debug.xps()) print("------dataDataTreeSet------")

   ## check for presence of data
   data  <- object@data;
   if (min(dim(data)) == 0) {
      stop(paste("slot", sQuote("data"), "has no data"));
   }#if

   ## check for presence of scheme mask
   msk <- chipMask(object@scheme);
   if (nrow(msk) == 0 || ncol(msk) == 0) {
      cat("slot", sQuote("mask"), "of", sQuote("scheme"),
          "is empty, importing mask from scheme.root...\n");
      msk <- export(object@scheme,
                    treetype     = "scm",
                    varlist      = "fMask",
                    as.dataframe = TRUE,
                    verbose      = FALSE);
   }#if
   ## need to sort msk to Y then X
   msk <- msk[order(msk[,"Y"], msk[,"X"]),];

   ## get data indices from scheme mask
   id <- which(msk[,"Mask"] < 999999); ##all
   if (chipType(object) == "GeneChip") {
      if (which == "pm") {
         id <- which(msk[,"Mask"] == 1);
      } else if (which == "mm") {
         id <- which(msk[,"Mask"] == 0);
      } else if (which == "both") {
         id <- c(which(msk[,"Mask"] == 1), which(msk[,"Mask"] == 0));
         id <- id[order(id)];
      }#if
   } else if (chipType(object) == "GenomeChip") {
      if (which == "") {
         id <- 1:nrow(data); ##all
      } else if (which == "antigenomic") {
         id <- which(msk[,"Mask"] == -2);
      } else {
         level <- exonLevel(which, "GenomeChip", as.sum=FALSE);
         id    <- unlist(sapply(level, function(x) {which(msk[,"Mask"] == x)}));
         id    <- id[order(id)];
      }#if
   } else if (chipType(object) == "ExonChip") {
      if (which == "") {
         id <- 1:nrow(data); ##all
      } else if (which == "genomic") {
         id <- which(msk[,"Mask"] == -1);
      } else if (which == "antigenomic") {
         id <- which(msk[,"Mask"] == -2);
      } else {
         level <- exonLevel(which, "ExonChip", as.sum=FALSE);
         id    <- unlist(sapply(level, function(x) {which(msk[,"Mask"] == x)}));
         id    <- id[order(id)];
      }#if
   }#if

   treenames <- namePart(object@treenames);
   datanames <- namePart(colnames(data));

   ds <- data[,!is.na(match(datanames, treenames))];
   return(ds[id,]);
}#dataDataTreeSet

setMethod("validData", "DataTreeSet", dataDataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"bgrdDataTreeSet" <-
function(object,
         which="") 
{
   if (debug.xps()) print("------bgrdDataTreeSet------")

   ## check for presence of bgrd
   bgrd  <- object@bgrd;
   if (min(dim(bgrd)) == 0) {
      stop(paste("slot", sQuote("bgrd"), "has no data"));
   }#if

   treenames <- namePart(object@bgtreenames);
   bgrdnames <- namePart(colnames(bgrd));

   return(bgrd[,!is.na(match(bgrdnames, treenames))]);
}#bgrdDataTreeSet

setMethod("validBgrd", "DataTreeSet", bgrdDataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"addDataDataTreeSet" <-
function(object,
         celdir   = NULL,
         celfiles = "",
         celnames = NULL,
         project  = NULL,
         verbose  = TRUE)
{
   if (debug.xps()) print("------addDataDataTreeSet------")

   ## check for presence of cel-file directory
   if (!is.null(celdir)) {
      if (celdir == "") {
         celdir <- getwd();
      }#if
      if (!file.exists(celdir)) {
         stop(paste("argument", sQuote("celdir"), "is not a system directory"));
      }#if
      if (substr(celdir,nchar(celdir),nchar(celdir)) != "/") {
         celdir <- paste(celdir, "/", sep="");
      }#if
   }#if

   ## get cel-files
   if (celfiles[1] == "*") {
      celfiles <- list.files(celdir, pattern = "\\.[cC][eE][lL]");
   } else if (is.null(celdir)) {
      existing <- file.exists(celfiles);
      if (length(existing[existing==TRUE]) != length(celfiles)) {
         stop(paste("some", sQuote("celfiles"), " are not found."));
      }#if
   } else {
      tmpfiles <- list.files(celdir, pattern = "\\.[cC][eE][lL]");
      tmpfiles <- match(celfiles, tmpfiles);
      if (length(tmpfiles[is.na(tmpfiles)]) > 0) {
         stop(paste("some", sQuote("celfiles"), " are not found in directory",
                    sQuote(celdir)));
      }#if
   }#if

   ## get number of cel-files
   ncel <- length(celfiles);
   if (ncel == 0) {
      stop(paste("directory", sQuote(celdir), "does not contain CEL-files"));
   }#if

   ## check for presence of potential celnames
   if (!is.null(celnames) && (length(unique(celnames)) != ncel)) {
      stop(paste(sQuote("celnames"), "must have same number of entries as",
                 sQuote("celfiles")));
   }#if
   if (is.null(celnames)) {
      celnames <- CELNames(celfiles);
   } else {
      celnames <- unlist(strsplit(celnames, "\\.[cC][eE][lL]"));
   }#if

   ## concatenate celfiles with celdir
   tmpfiles <- unlist(celfiles);
   celfiles <- paste(celdir, celfiles, sep="");

   ## tree names from celfiles
   if (length(celnames) == 1 && celnames == "") {
      celnames <- unique(unlist(strsplit(tmpfiles, "\\.[cC][eE][lL]")));
   }#if

   ## get root files and setname
   rootfile <- object@rootfile;
   scheme   <- object@scheme;
   setname  <- object@setname;
   chipname <- chipName(scheme);
   chiptype <- chipType(scheme);

   ## project information
   replace <- 0;
   projct  <- ""; nproject <- 0;
   author  <- ""; nauthor  <- 0;
   dataset <- ""; ndataset <- 0;
   source  <- ""; nsource  <- 0;
   sample  <- ""; nsample  <- 0;
   cell    <- ""; ncell    <- 0;
   pcell   <- ""; npcell   <- 0;
   tissue  <- ""; ntissue  <- 0;
   biopsy  <- ""; nbiopsy  <- 0;
   array   <- ""; narray   <- 0;
   hyb     <- ""; nhyb     <- 0;
   treat   <- ""; ntreat   <- 0;

   if (!is.null(project) && is(project, "ProjectInfo")) {
      replace  <- 1;
      nproject <- length(projectInfo(project));
      nauthor  <- length(authorInfo(project));
      ndataset <- length(datasetInfo(project));
      nsource  <- length(sourceInfo(project));
      nsample  <- length(sampleInfo(project));
      ncell    <- length(cellineInfo(project));
      npcell   <- length(primcellInfo(project));
      ntissue  <- length(tissueInfo(project));
      nbiopsy  <- length(biopsyInfo(project));
      narray   <- length(arrayInfo(project));
      nhyb     <- length(unlist(hybridizInfo(project)));
      ntreat   <- length(unlist(treatmentInfo(project)));

      if (nproject == 5)  projct  <- projectInfo(project)  else projct  <- "";
      if (nauthor  == 8)  author  <- authorInfo(project)   else author  <- "";
      if (ndataset == 7)  dataset <- datasetInfo(project)  else dataset <- "";
      if (nsource  == 6)  source  <- sourceInfo(project)   else source  <- "";
      if (nsample  == 12) sample  <- sampleInfo(project)   else sample  <- "";
      if (ncell    == 15) cell    <- cellineInfo(project)  else cell    <- "";
      if (npcell   == 14) pcell   <- primcellInfo(project) else pcell   <- "";
      if (ntissue  == 19) tissue  <- tissueInfo(project)   else tissue  <- "";
      if (nbiopsy  == 19) biopsy  <- biopsyInfo(project)   else biopsy  <- "";
      if (narray   == 4)  array   <- arrayInfo(project)    else array   <- "";

      if ((narray == 4) & (array$chipname != chipname)) {
         stop(paste(sQuote("arrayInfo$chipname"), " must be", sQuote(chipname)));
      }#if
      if ((narray == 4) & (array$chiptype != chiptype)) {
         stop(paste(sQuote("arrayInfo$chiptype"), " must be", sQuote(chiptype)));
      }#if

###########
# to do: check if hybnames and treatnames are identical to celnames
###########
      if (nhyb >= 9) {
         hyb <- as.vector(unlist(t(hybridizInfo(project))));
      } else {
         hyb <- "";
      }#if
      if (ntreat >= 8) {
         treat <- as.vector(unlist(t(treatmentInfo(project))));
      } else {
         treat <- "";
      }#if
   }#if

   ## update root data file
   r <- .C("ImportData",
           as.character(rootfile),
           as.character(""),
           as.character(chipType(scheme)),
           as.character(chipName(scheme)),
           as.character(rootFile(scheme)),
           as.character(setname),
           as.character(celfiles),
           as.character(celnames),
           as.integer(ncel),
           as.character(projct),
           as.integer(nproject),
           as.character(author),
           as.integer(nauthor),
           as.character(dataset),
           as.integer(ndataset),
           as.character(source),
           as.integer(nsource),
           as.character(sample),
           as.integer(nsample),
           as.character(cell),
           as.integer(ncell),
           as.character(pcell),
           as.integer(npcell),
           as.character(tissue),
           as.integer(ntissue),
           as.character(biopsy),
           as.integer(nbiopsy),
           as.character(array),
           as.integer(narray),
           as.character(hyb),
           as.integer(nhyb),
           as.character(treat),
           as.integer(ntreat),
           as.integer(replace),
           as.integer(1),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("ImportData")));
      return(NULL);
   }#if

   ## get all treenames (original and newly added treenames)
   celnames <- getTreeNames(rootfile, "cel", setname);

   ## create new class DataTreeSet
   set <- new("DataTreeSet",
              setname   = setname,
              settype   = object@settype,
              rootfile  = rootfile,
              filedir   = object@filedir,
              numtrees  = as.numeric(length(celnames)),
              treenames = as.list(celnames),
              scheme    = scheme);

   ## add project information
   if (!is.null(project) && is(project, "ProjectInfo")) {
      projectInfo(set) <- project;
   }#if

   return(set);
}#addDataDataTreeSet

setMethod("addData", "DataTreeSet", addDataDataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"pm.DataTreeSet" <-
function(object,
         which = "pm") 
{
   if (debug.xps()) print("------pm.DataTreeSet------")

   level <- unlist(strsplit(which, "\\+"));
   LEVEL <- c("pm",
              "core",     "metacore",
              "extended", "metaextended",
              "full",     "metafull",
              "ambiguous",
              "affx");
   if (is.na(all(match(level, LEVEL)))) {
      stop(paste("invalid argument", sQuote("which")));
   }#if

   return(validData(object, which=which));
}#pm.DataTreeSet

setMethod("pm", "DataTreeSet", pm.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"mm.DataTreeSet" <-
function(object,
         which = "mm") 
{
   if (debug.xps()) print("------mm.DataTreeSet------")

   LEVEL <- c("mm", "genomic", "antigenomic");
   if (is.na(all(match(which, LEVEL)))) {
      stop(paste("invalid argument", sQuote("which")));
   }#if

   return(validData(object, which=which));
}#mm.DataTreeSet

setMethod("mm", "DataTreeSet", mm.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"rawCELName.DataTreeSet" <-
function(object,
         treename = "*",
         fullpath = TRUE) 
{
   if (debug.xps()) print("------rawCELName.DataTreeSet------")

   ## get valid tree names
   treenames <- treeNames(object);
   if (treename == "*") {
      treename <- treenames;
   } else if (is.na(match(treename, treenames))) {
      stop(paste("invalid tree name", sQuote(treename)));
   }#if

   ## setname.treename
   treename <- paste(setName(object), treename, sep=".");

   rawnames <- .C("GetRawCELNames",
                as.character(rootFile(object)),
                as.integer(length(treename)),
                as.character(treename),
                celnames=character(length(treename)),
                PACKAGE="xps")$celnames;

   ## names only
   if (fullpath == FALSE) {
#      raw <- unlist(strsplit(rawnames, "/"));
      rawnames <- sapply(strsplit(rawnames, "/"), function(x)x[length(x)]);
   }#if

#   return(rawnames);
   return(as.vector(rawnames, mode="character"));
}#mm.DataTreeSet

setMethod("rawCELName", "DataTreeSet", rawCELName.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"rma.DataTreeSet" <-
function(object,
         filename   = character(),
         filedir    = getwd(),
         tmpdir     = "",
         background = "pmonly",
         normalize  = TRUE,
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE) 
{
   if (debug.xps()) print("------rma.DataTreeSet------")

   ## get schemefile and chipname
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chipname   <- chipName(object);
   chiptype   <- chipType(object);

   ## check for presence of alternative root scheme file
   if ((!is.null(xps.scheme)) &&
       is(xps.scheme, "SchemeTreeSet") &&
       (chipType(xps.scheme) == chiptype) &&
       (file.exists(rootFile(xps.scheme)))) {
      scheme     <- xps.scheme;
      schemefile <- rootFile(xps.scheme);
      chipname   <- chipName(xps.scheme);
      chiptype   <- chipType(xps.scheme);
   }#if

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check if root file exists (necessary for WinXP to test already here)
   if (existsROOTFile(rootfile)) {
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## check for presence of valid background option
   if (!(identical(background, "pmonly") ||
         identical(background, "mmonly") ||
         identical(background, "genomic") ||
         identical(background, "antigenomic") ||
         identical(background, "none"))) {
      stop(paste(sQuote(background), "is not a valid background option"));
   }#if

   ## check for valid normalize
   if (!is.logical(normalize)) {
      stop(paste(sQuote("normalize"), "must be TRUE or FALSE"));
   }#if

   ## check for presence of valid transcript option
   transcript <- validTranscriptOption(option);

   ## check for correct exonlevel
   exlevel <- exonLevel(exonlevel, chiptype);

   ## get treenames to normalize as fullnames=/datadir/treenames
   listnames <- listTreeNames(object);
   treenames <- listnames$treenames;
   fullnames <- listnames$fullnames;
   numtrees  <- length(treenames);

   ## define setname and settype for new treeset
   setname <- "PreprocesSet";
   settype <- "preprocess";

   ## preprocess RMA
   r <- .C("PreprocessRMA",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.character(background),
           as.character(transcript),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.integer(normalize),
           as.integer(exlevel[1]),
           as.integer(exlevel[2]),
           as.integer(exlevel[3]),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("PreprocessRMA")));
      return(NULL);
   }#if

   ## export result to outfile and import as dataframe ds
   ds <- data.frame(matrix(nr=0,nc=0));
   if (add.data) {
      outfile  <- sub("\\.root", ".txt", rootfile);
      ## get treename "treeset.treename.treetype"
      treetype <- "mdp";
      treename <- paste(setname, "*", treetype, sep=".");
      numtrees <- 1; # must be one for treename="*"

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(treename),
              as.integer(numtrees),
              as.character(treetype),
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

      if (file.exists(outfile)) {
         ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("could not export results as", sQuote(outfile)));
      }#if
   }#if

   ## get treenames after preprocessing
   treenames <- as.list(getTreeNames(rootfile, "mdp"));
   numtrees  <- length(treenames);

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
              params    = list(),
              exprtype  = "rma",
              normtype  = "none");

   return(set);
}#rma.DataTreeSet

setMethod("xpsRMA", "DataTreeSet", rma.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"mas4.DataTreeSet" <-
function(object,
         filename   = character(),
         filedir    = getwd(),
         tmpdir     = "",
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   if (debug.xps()) print("------mas4.DataTreeSet------")

   ## get schemefile and chipname
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chipname   <- chipName(object);
   chiptype   <- chipType(object);

   ## check for presence of alternative root scheme file
   if ((!is.null(xps.scheme)) &&
       is(xps.scheme, "SchemeTreeSet") &&
       (chipType(xps.scheme) == chiptype) &&
       (file.exists(rootFile(xps.scheme)))) {
      scheme     <- xps.scheme;
      schemefile <- rootFile(xps.scheme);
      chipname   <- chipName(xps.scheme);
      chiptype   <- chipType(xps.scheme);
   }#if

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check if root file exists (necessary for WinXP to test already here)
   if (existsROOTFile(rootfile)) {
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## check for presence of valid transcript option
   transcript <- validTranscriptOption(option);

   ## check for correct exonlevel
   exlevel <- exonLevel(exonlevel, chiptype);

   ## get treenames to normalize as fullnames=/datadir/treenames
   listnames <- listTreeNames(object);
   treenames <- listnames$treenames;
   fullnames <- listnames$fullnames;
   numtrees  <- length(treenames);

   ## define setname and settype for new treeset
   setname <- "PreprocesSet";
   settype <- "preprocess";

   ## preprocess MAS4
   r <- .C("PreprocessMAS4",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.character(transcript),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.integer(exlevel[1]),
           as.integer(exlevel[3]),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("PreprocessMAS4")));
      return(NULL);
   }#if

   ## export result to outfile and import as dataframe ds
   ds <- data.frame(matrix(nr=0,nc=0));
   if (add.data) {
      outfile  <- sub("\\.root", ".txt", rootfile);
      ## get treename "treeset.treename.treetype"
      treetype <- "adf";
      treename <- paste(setname, "*", treetype, sep=".");
      numtrees <- 1; # must be one for treename="*"

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(treename),
              as.integer(numtrees),
              as.character(treetype),
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

      if (file.exists(outfile)) {
         ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("error: could not export results as", sQuote(outfile)));
      }#if
   }#if

   ## normalized treenames 
   treenames <- as.list(getTreeNames(rootfile, "adf"));
   numtrees  <- length(treenames);

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
              params    = list(),
              exprtype  = "mas4",
              normtype  = "none");

   return(set);
}#mas4.DataTreeSet

setMethod("xpsMAS4", "DataTreeSet", mas4.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"mas5.DataTreeSet" <-
function(object,
         filename   = character(),
         filedir    = getwd(),
         tmpdir     = "",
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   if (debug.xps()) print("------mas5.DataTreeSet------")

   ## get schemefile and chipname
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chipname   <- chipName(object);
   chiptype   <- chipType(object);

   ## check for presence of alternative root scheme file
   if ((!is.null(xps.scheme)) &&
       is(xps.scheme, "SchemeTreeSet") &&
       (chipType(xps.scheme) == chiptype) &&
       (file.exists(rootFile(xps.scheme)))) {
      scheme     <- xps.scheme;
      schemefile <- rootFile(xps.scheme);
      chipname   <- chipName(xps.scheme);
      chiptype   <- chipType(xps.scheme);
   }#if

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check if root file exists (necessary for WinXP to test already here)
   if (existsROOTFile(rootfile)) {
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## check for presence of valid transcript option
   transcript <- validTranscriptOption(option);

   ## check for correct exonlevel
   exlevel <- exonLevel(exonlevel, chiptype);

   ## get treenames to normalize as fullnames=/datadir/treenames
   listnames <- listTreeNames(object);
   treenames <- listnames$treenames;
   fullnames <- listnames$fullnames;
   numtrees  <- length(treenames);

   ## define setname and settype for new treeset
   setname <- "PreprocesSet";
   settype <- "preprocess";

   ## preprocess MAS5
   r <- .C("PreprocessMAS5",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.character(transcript),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.integer(exlevel[1]),
           as.integer(exlevel[3]),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("PreprocessMAS5")));
      return(NULL);
   }#if

   ## export result to outfile and import as dataframe ds
   ds <- data.frame(matrix(nr=0,nc=0));
   if (add.data) {
      outfile  <- sub("\\.root", ".txt", rootfile);
      ## get treename "treeset.treename.treetype"
      treetype <- "tbw";
      treename <- paste(setname, "*", treetype, sep=".");
      numtrees <- 1; # must be one for treename="*"

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(treename),
              as.integer(numtrees),
              as.character(treetype),
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

      if (file.exists(outfile)) {
         ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("error: could not export results as", sQuote(outfile)));
      }#if
   }#if

   ## normalized treenames 
   treenames <- as.list(getTreeNames(rootfile, "tbw"));
   numtrees  <- length(treenames);

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
              params    = list(),
              exprtype  = "mas5",
              normtype  = "none");

   return(set);
}#mas5.DataTreeSet

setMethod("xpsMAS5", "DataTreeSet", mas5.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"mas5Call.DataTreeSet" <-
function(object,
         filename         = character(),
         filedir          = getwd(),
         tmpdir           = "",
         tau              = 0.015,
         alpha1           = 0.04,
         alpha2           = 0.06,
         ignore.saturated = TRUE,
         option           = "transcript",
         exonlevel        = "",
         xps.scheme       = NULL,
         add.data         = TRUE,
         verbose          = TRUE)
{
   if (debug.xps()) print("------mas5Call.DataTreeSet------")

   ## get schemefile and chipname
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chipname   <- chipName(object);
   chiptype   <- chipType(object);

   ## check for presence of alternative root scheme file
   if ((!is.null(xps.scheme)) &&
       is(xps.scheme, "SchemeTreeSet") &&
       (chipType(xps.scheme) == chiptype) &&
       (file.exists(rootFile(xps.scheme)))) {
      scheme     <- xps.scheme;
      schemefile <- rootFile(xps.scheme);
      chipname   <- chipName(xps.scheme);
      chiptype   <- chipType(xps.scheme);
   }#if

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check if root file exists (necessary for WinXP to test already here)
   if (existsROOTFile(rootfile)) {
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## check for valid tau, alpha1, alpha2, ignore.saturated
   if (!(is.numeric(tau) && tau > 0 && tau < 1)) {
      stop(paste(sQuote("tau"), "must be numeric and in range (0,1)"));
   }#if
   if (!(is.numeric(alpha1) && alpha1 > 0 && alpha1 < 1)) {
      stop(paste(sQuote("alpha1"), "must be numeric and in range (0,1)"));
   }#if
   if (!(is.numeric(alpha2) && alpha2 > 0 && alpha2 < 1)) {
      stop(paste(sQuote("alpha2"), "must be numeric and in range (0,1)"));
   }#if
   if (!is.logical(ignore.saturated)) {
      stop(paste(sQuote("ignore.saturated"), "must be TRUE or FALSE"));
   }#if

   ## check for presence of valid transcript option
   transcript <- validTranscriptOption(option);

   ## check for correct exonlevel
   exlevel <- exonLevel(exonlevel, chiptype);

   ## get treenames as fullnames=/datadir/treenames
   listnames <- listTreeNames(object);
   treenames <- listnames$treenames;
   fullnames <- listnames$fullnames;
   numtrees  <- length(treenames);

   ## define setname and settype for new treeset
   setname <- "CallSet";
   settype <- "preprocess";

   ## MAS5 detection call
   r <- .C("PreprocessMAS5Call",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.character(transcript),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.double(tau),
           as.double(alpha1),
           as.double(alpha2),
           as.integer(ignore.saturated),
           as.integer(exlevel[1]),
           as.integer(exlevel[3]),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("PreprocessMAS5Call")));
      return(NULL);
   }#if

   ## export p-value to outfile and import as dataframe ds
   ds <- data.frame(matrix(nr=0,nc=0));
   dc <- data.frame(matrix(nr=0,nc=0));
   if (add.data) {
      outfile  <- sub("\\.root", "_pval.txt", rootfile);
      ## get treename "treeset.treename.treetype"
      treetype <- "dc5";
      treename <- paste(setname, "*", treetype, sep=".");
      numtrees <- 1; # must be one for treename="*"

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(treename),
              as.integer(numtrees),
              as.character(treetype),
              as.character("fUnitName:fPValue"),
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
         ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("error: could not export results as", sQuote(outfile)));
      }#if

      ## export detection call to outfile and import as dataframe dc
      outfile  <- sub("\\.root", "_call.txt", rootfile);
      # get treename "treeset.treename.treetype"
      treename <- paste(setname, "*", treetype, sep=".");

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(treename),
              as.integer(numtrees),
              as.character(treetype),
              as.character("fUnitName:fCall"),
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
         dc <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("error: could not export results as", sQuote(outfile)));
      }#if
   }#if

   ## mas5 call treenames 
   treenames <- as.list(getTreeNames(rootfile, "dc5"));
   numtrees  <- length(treenames);

   ## named parameter list
   params <- list(tau    = tau,
                  alpha1 = alpha1,
                  alpha2 = alpha2,
                  ignore = ignore.saturated);

   ## create new class CallTreeSet
   set <- new("CallTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              data      = ds,
              params    = params,
              calltype  = "mas5",
              detcall   = dc);
   return(set);
}#mas5Call.DataTreeSet

setMethod("xpsMAS5Call", "DataTreeSet", mas5Call.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"dabgCall.DataTreeSet" <-
function(object,
         filename   = character(),
         filedir    = getwd(),
         alpha1     = 0.04,
         alpha2     = 0.06,
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   if (debug.xps()) print("------dabgCall.DataTreeSet------")

   ## get schemefile and chipname
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chipname   <- chipName(object);
   chiptype   <- chipType(object);

   ## check for presence of alternative root scheme file
   if ((!is.null(xps.scheme)) &&
       is(xps.scheme, "SchemeTreeSet") &&
       (chipType(xps.scheme) == chiptype) &&
       (file.exists(rootFile(xps.scheme)))) {
      scheme     <- xps.scheme;
      schemefile <- rootFile(xps.scheme);
      chipname   <- chipName(xps.scheme);
      chiptype   <- chipType(xps.scheme);
   }#if

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check if root file exists (necessary for WinXP to test already here)
   if (existsROOTFile(rootfile)) {
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for valid alpha1, alpha2
   if (!(is.numeric(alpha1) && alpha1 > 0 && alpha1 < 1)) {
      stop(paste(sQuote("alpha1"), "must be numeric and in range (0,1)"));
   }#if
   if (!(is.numeric(alpha2) && alpha2 > 0 && alpha2 < 1)) {
      stop(paste(sQuote("alpha2"), "must be numeric and in range (0,1)"));
   }#if

   ## check for presence of valid transcript option
   transcript <- validTranscriptOption(option);

   ## check for correct exonlevel
   exlevel <- exonLevel(exonlevel, chiptype);

   ## get treenames as fullnames=/datadir/treenames
   listnames <- listTreeNames(object);
   treenames <- listnames$treenames;
   fullnames <- listnames$fullnames;
   numtrees  <- length(treenames);

   ## define setname and settype for new treeset
   setname <- "CallSet";
   settype <- "preprocess";

   ## DABG detection call
   r <- .C("PreprocessDABGCall",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(chiptype),
           as.character(schemefile),
           as.character(transcript),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.double(alpha1),
           as.double(alpha2),
           as.integer(exlevel[3]),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("PreprocessDABGCall")));
      return(NULL);
   }#if

   ## export p-value to outfile and import as dataframe ds
   ds <- data.frame(matrix(nr=0,nc=0));
   dc <- data.frame(matrix(nr=0,nc=0));
   if (add.data) {
      outfile  <- sub("\\.root", "_pval.txt", rootfile);
      ## get treename "treeset.treename.treetype"
      treetype <- "dab";
      treename <- paste(setname, "*", treetype, sep=".");
      numtrees <- 1; # must be one for treename="*"

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(treename),
              as.integer(numtrees),
              as.character(treetype),
              as.character("fUnitName:fPValue"),
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
         ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("error: could not export results as", sQuote(outfile)));
      }#if

      ## export detection call to outfile and import as dataframe dc
      outfile  <- sub("\\.root", "_call.txt", rootfile);
      ## get treename "treeset.treename.treetype"
      treename <- paste(setname, "*", treetype, sep=".");

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(treename),
              as.integer(numtrees),
              as.character(treetype),
              as.character("fUnitName:fCall"),
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
         dc <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("error: could not export results as", sQuote(outfile)));
      }#if
   }#if

   ## mas5 call treenames 
   treenames <- as.list(getTreeNames(rootfile, "dab"));
   numtrees  <- length(treenames);

   ## named parameter list
   params <- list(alpha1 = alpha1,
                  alpha2 = alpha2);

   ## create new class CallTreeSet
   set <- new("CallTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              data      = ds,
              params    = params,
              calltype  = "dabg",
              detcall   = dc);
   return(set);
}#dabgCall.DataTreeSet

setMethod("xpsDABGCall", "DataTreeSet", dabgCall.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"iniCall.DataTreeSet" <-
function(object,
         filename   = character(),
         filedir    = getwd(),
         tmpdir     = "",
         weight     = 0.5,
         mu         = 0.0,
         scale      = 1.0,
         tol        = 0.00001,
         cyc        = 100,
         alpha1     = 0.4,
         alpha2     = 0.6,
         version    = "1.3.1",
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   if (debug.xps()) print("------iniCall.DataTreeSet------")

   ## get schemefile and chipname
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chipname   <- chipName(object);
   chiptype   <- chipType(object);

   ## check for presence of alternative root scheme file
   if ((!is.null(xps.scheme)) &&
       is(xps.scheme, "SchemeTreeSet") &&
       (chipType(xps.scheme) == chiptype) &&
       (file.exists(rootFile(xps.scheme)))) {
      scheme     <- xps.scheme;
      schemefile <- rootFile(xps.scheme);
      chipname   <- chipName(xps.scheme);
      chiptype   <- chipType(xps.scheme);
   }#if

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check if root file exists (necessary for WinXP to test already here)
   if (existsROOTFile(rootfile)) {
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## check for valid version
   version <- as.integer(gsub("\\.","",version));
   if (!(version == 131 || version == 130)) {
      stop(paste("wrong version, currently only", dQuote("1.3.1"), "or", dQuote("1.3.0"),
                 "are allowed"));
   }#if

   ## check for valid parameters
   if (!(is.numeric(weight) && is.numeric(mu) && is.numeric(scale) && is.numeric(tol))) {
      stop("parameters <weight,mu,scale,tol> must be numeric");
   }#if
   if (!(is.numeric(cyc) && cyc >= 0)) {
      stop(paste(sQuote("cyc"), "must be integer >= 0"));
   }#if
   cyc <- as.integer(cyc);

   if (!(is.numeric(alpha1) && alpha1 > 0 && alpha1 < 1)) {
      stop(paste(sQuote("alpha1"), "must be numeric and in range (0,1)"));
   }#if
   if (!(is.numeric(alpha2) && alpha2 > 0 && alpha2 < 1)) {
      stop(paste(sQuote("alpha2"), "must be numeric and in range (0,1)"));
   }#if

   ## check for presence of valid transcript option
   transcript <- validTranscriptOption(option);

   ## check for correct exonlevel
   exlevel <- exonLevel(exonlevel, chiptype);

   ## get treenames as fullnames=/datadir/treenames
   listnames <- listTreeNames(object);
   treenames <- listnames$treenames;
   fullnames <- listnames$fullnames;
   numtrees  <- length(treenames);

   ## define setname and settype for new treeset
   setname <- "CallSet";
   settype <- "preprocess";

   ## Informative call
   r <- .C("PreprocessINICall",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.character(transcript),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.integer(version),
           as.double(weight),
           as.double(mu),
           as.double(scale),
           as.double(tol),
           as.integer(cyc),
           as.double(alpha1),
           as.double(alpha2),
           as.integer(exlevel[2]),
           as.integer(exlevel[3]),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("PreprocessINICall")));
      return(NULL);
   }#if

   ## export p-value to outfile and import as dataframe ds
   ds <- data.frame(matrix(nr=0,nc=0));
   dc <- data.frame(matrix(nr=0,nc=0));
   if (add.data) {
      outfile  <- sub("\\.root", "_pval.txt", rootfile);
      ## get treename "treeset.treename.treetype"
      treetype <- "ini";
      treename <- paste(setname, "*", treetype, sep=".");
      numtrees <- 1; # must be one for treename="*"

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(treename),
              as.integer(numtrees),
              as.character(treetype),
              as.character("fUnitName:fPValue"),
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
         ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("error: could not export results as", sQuote(outfile)));
      }#if

      ## export detection call to outfile and import as dataframe dc
      outfile  <- sub("\\.root", "_call.txt", rootfile);
      # get treename "treeset.treename.treetype"
      treename <- paste(setname, "*", treetype, sep=".");

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(treename),
              as.integer(numtrees),
              as.character(treetype),
              as.character("fUnitName:fCall"),
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
         dc <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("error: could not export results as", sQuote(outfile)));
      }#if
   }#if

   ## informative call treenames 
   treenames <- as.list(getTreeNames(rootfile, "ini"));
   numtrees  <- length(treenames);

   ## named parameter list
   params <- list(version = version,
                  weight  = weight,
                  mu      = mu,
                  scale   = scale,
                  tol     = tol,
                  cyc     = cyc,
                  alpha1  = alpha1,
                  alpha2  = alpha2);

   ## create new class CallTreeSet
   set <- new("CallTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              data      = ds,
              params    = params,
              calltype  = "ini",
              detcall   = dc);
   return(set);
}#iniCall.DataTreeSet

setMethod("xpsINICall", "DataTreeSet", iniCall.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"preprocess.DataTreeSet" <-
function(object,
         filename          = character(),
         filedir           = getwd(),
         tmpdir            = "",
         update            = FALSE,
         bgcorrect.method  = NULL,
         bgcorrect.select  = character(),
         bgcorrect.option  = character(),
         bgcorrect.params  = list(),
         normalize.method  = NULL,
         normalize.select  = character(),
         normalize.option  = character(),
         normalize.logbase = character(),
         normalize.params  = list(),
         summarize.method  = NULL,
         summarize.select  = character(),
         summarize.option  = character(),
         summarize.logbase = character(),
         summarize.params  = list(),
         reference.index   = 0,
         reference.method  = "mean",
         reference.params  = list(0.0),
         exonlevel         = "",
         xps.scheme        = NULL,
         add.data          = TRUE,
         verbose           = TRUE)
{
   if (debug.xps()) print("------preprocess.DataTreeSet------")

   if (is.null(bgcorrect.method) && 
       is.null(normalize.method) &&
       is.null(summarize.method)) {
      stop(paste("at least one", sQuote("xps.method"), "must not be NULL"));
   }#if

   ## get schemefile and chiptype for object
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chipname   <- chipName(object);
   chiptype   <- chipType(object);

   ## check for presence of alternative root scheme file
   if ((!is.null(xps.scheme)) &&
       is(xps.scheme, "SchemeTreeSet") &&
       (chipType(xps.scheme) == chiptype) &&
       (file.exists(rootFile(xps.scheme)))) {
      scheme     <- xps.scheme;
      schemefile <- rootFile(xps.scheme);
      chipname   <- chipName(xps.scheme);
      chiptype   <- chipType(xps.scheme);
   }#if

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## update: need to set filename to rootfile
   if (!is.logical(update)) {
      stop(paste(sQuote("update"), "must be TRUE or FALSE"));
   } else if (update == TRUE) {
      filename <- rootfile;
   } else if (update == FALSE && existsROOTFile(rootfile)) {
      ## check if root file exists (necessary for WinXP to test already here)
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## bgcorrect method
   if (is.null(bgcorrect.method)) {
      bgcorrect.method <- "none";
   } else {
      TYPE <- c("sector", "weightedsector", "rma", "gccontent");
      if (is.na(match(bgcorrect.method, TYPE))) {
         stop(paste(sQuote("bgcorrect.method"),
                   "must be one of <sector,weightedsector,rma,gccontent>"));
      }#if

      TYPE <- c("pmonly", "mmonly", "both", "genomic", "antigenomic", "all", "none");
      if (is.na(match(bgcorrect.select, TYPE))) {
         stop(paste(sQuote(bgcorrect.select), "is not a valid selector option"));
      }#if

      TYPE <- c("subtractbg", "correctbg", "attenuatebg", "pmonly:epanechnikov");
      if (is.na(match(bgcorrect.option, TYPE))) {
         warning(paste(sQuote("bgcorrect.option"),
                      "is different from <pmonly:epanechnikov> for rma"));
      }#if
      if (length(bgcorrect.params) == 0) {
         stop(paste("empty parameter list", sQuote("bgcorrect.params")));
      }#if
   }#if

   ## normalize method
   if (is.null(normalize.method)) {
      normalize.method <- "none";
   } else {
      TYPE <- c("mean", "median", "quantile");
      if (is.na(match(normalize.method, TYPE))) {
         cat("methods <lowess,supsmu> have not been tested yet for DataTreeSet\n");
         stop(paste(sQuote("normalize.method"), "must be one of <mean,median,quantile>"));
      }#if

      TYPE <- c("pmonly", "mmonly", "both", "all");
      if (is.na(match(normalize.select, TYPE))) {
         stop(paste(sQuote("normalize.select"), "is not a valid selector option"));
      }#if

      ## check for valid normalization "option:logbase"
      normalize.option  <- validOption(normalize.option);
      normalize.logbase <- validLogbase(normalize.logbase);
      normalize.option  <- paste(normalize.option, normalize.logbase, sep=":");

      if (length(normalize.params) == 0) {
         stop(paste("empty parameter list", sQuote("normalize.params")));
      }#if
   }#if

   ## summarize method
   if (is.null(summarize.method)) {
      summarize.method <- "none";
   } else {
      TYPE <- c("avgdiff", "tukeybiweight", "medianpolish", "farms", "dfw");
      if (is.na(match(summarize.method, TYPE))) {
         stop(paste(sQuote("summarize.method"),
              "must be one of <avgdiff,tukeybiweight,medianpolish,farms,dfw>"));
      }#if

      TYPE <- c("pmonly", "mmonly", "both", "all", "none");
      if (is.na(match(summarize.select, TYPE))) {
         stop(paste(sQuote(summarize.select), "is not a valid selector option"));
      }#if

      ## check for valid summarization "option:logbase"
      summarize.option  <- validTranscriptOption(summarize.option);
      summarize.logbase <- validLogbase(summarize.logbase);
      summarize.option  <- paste(summarize.option, summarize.logbase, sep=":");

      if (length(summarize.params) == 0) {
         stop(paste("empty parameter list", sQuote("summarize.params")));
      }#if
   }#if

   ## check for correct exonlevel
   exlevel <- exonLevel(exonlevel, chiptype);

   ## get treenames to normalize as fullnames=/datadir/treenames
   listnames <- listTreeNames(object);
   treenames <- listnames$treenames;
   fullnames <- listnames$fullnames;
   numtrees  <- length(treenames);

   ## reference tree
   reference.tree <- "";
   if (!(reference.index >= 0 & reference.index <= numtrees)) {
      stop(paste(sQuote("reference.index"),
                        "must be 0 or less than number of trees"));
   }#if

   setname <- object@setname;
   if (reference.index == 0) {
      reference.tree <- paste(setname, "/*.", extenPart(treenames), sep="");
   } else {
      reference.tree <- paste(setname, treenames[[reference.index]], sep="/");
   }#if

   ## reference method
   TYPE <- c("mean", "median");
   if (is.na(match(reference.method, TYPE))) {
      stop(paste(sQuote("reference.method"), "must be <mean, median>"));
   }#if

   ## check reference parameter list
   if (length(reference.params) == 0) {
      stop(paste("empty parameter list", sQuote("reference.params")));
   }#if

   ## define setname and settype for new treeset
   setname <- "PreprocesSet";
   settype <- "preprocess";

   ## preprocess
   r <- .C("Preprocess",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.integer(update),
           as.character(bgcorrect.method),
           as.character(bgcorrect.select),
           as.character(bgcorrect.option),
           as.integer(length(bgcorrect.params)),
           as.double(bgcorrect.params),
           as.character(normalize.method),
           as.character(normalize.select),
           as.character(normalize.option),
           as.integer(length(normalize.params)),
           as.double(normalize.params),
           as.character(summarize.method),
           as.character(summarize.select),
           as.character(summarize.option),
           as.integer(length(summarize.params)),
           as.double(summarize.params),
           as.character(reference.tree),
           as.character(reference.method),
           as.double(reference.params),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.integer(exlevel[1]),
           as.integer(exlevel[2]),
           as.integer(exlevel[3]),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("Preprocess")));
      return(NULL);
   }#if

   if (summarize.method != "none") {
      exten <- type2Exten(summarize.method, settype);

      ds <- data.frame(matrix(nr=0,nc=0));
      if (add.data) {
         ## export result to outfile and import as dataframe
         outfile  <- sub("\\.root", ".txt", rootfile);
         ## get treename "treeset.treename.treetype"
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

         if (file.exists(outfile)) {
            ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
         } else {
            warning(paste("could not export results as", sQuote(outfile)));
         }#if
      }#if

      ## expression treenames 
      treenames <- as.list(getTreeNames(rootfile, exten));
      numtrees  <- length(treenames);

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
                 params    = as.list(summarize.params),
                 exprtype  = "custom",
                 normtype  = "none");
   } else if (normalize.method != "none") {
      ## get treenames with extension for method
      exten     <- type2Exten(normalize.method, settype);
      treenames <- getTreeNames(rootfile, exten);
      numtrees  <- length(treenames);

      ## create new class DataTreeSet
      set <- new("DataTreeSet",
                 setname   = setname,
                 settype   = settype,
                 rootfile  = rootfile,
                 filedir   = filedir,
                 numtrees  = numtrees,
                 treenames = as.list(treenames),
                 scheme    = scheme,
                 params    = as.list(normalize.params));
   } else if (bgcorrect.method != "none") {
      ## get treenames after background correction
      alltrees  <- getTreeNames(rootfile, "*");
      treenames <- alltrees[grep("\\.int",alltrees)];
      bgrdtrees <- alltrees[setdiff(1:length(alltrees), grep("\\.int",alltrees))];
      numtrees  <- length(treenames);
      numbgrds  <- length(bgrdtrees);
      if (numbgrds != numtrees) {
         stop(paste("number of bgrdtrees <", numbgrds,
                    "> is not equal to <", numtrees,">",sep=""));
      }#if

      ## create new class DataTreeSet
      set <- new("DataTreeSet",
                 setname   = setname,
                 settype   = settype,
                 rootfile  = rootfile,
                 filedir   = filedir,
                 numtrees  = numtrees,
                 treenames = as.list(treenames),
                 scheme    = scheme,
                 params    = as.list(bgcorrect.params));

      ## add background trees
      set@bgtreenames <- as.list(bgrdtrees);
   }#if

   return(set);
}#preprocess.DataTreeSet

setMethod("xpsPreprocess", "DataTreeSet", preprocess.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"bgCorrect.DataTreeSet" <-
function(object,
         filename  = character(),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "none",
         method    = character(),
         option    = character(),
         exonlevel = "",
         params    = list(),
         verbose   = TRUE)
{
   if (debug.xps()) print("------bgCorrect.DataTreeSet------")

   ## get schemefile and chipname
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chiptype   <- chipType(object);

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## update: need to set filename to rootfile
   if (!is.logical(update)) {
      stop(paste(sQuote("update"), "must be TRUE or FALSE"));
   } else if (update == TRUE) {
      filename <- rootfile;
   } else if (update == FALSE && existsROOTFile(rootfile)) {
      ## check if root file exists (necessary for WinXP to test already here)
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for presence of valid selector option
   TYPE <- c("pmonly", "mmonly", "both", "genomic", "antigenomic", "all", "none");
   if (is.na(match(select, TYPE))) {
      stop(paste(sQuote(select), "is not a valid selector option"));
   }#if

   ## check for valid background method
   TYPE <- c("sector", "weightedsector", "rma", "gccontent");
   if (is.na(match(method, TYPE))) {
      stop(paste(sQuote("method"),
                "must be one of <sector,weightedsector,rma,gccontent>"));
   }#if

   ## check for valid background correction option
   TYPE <- c("subtractbg", "correctbg", "attenuatebg", "pmonly:epanechnikov");
   if (is.na(match(option, TYPE))) {
      warning(paste(sQuote("option"),
                   "is different from <pmonly:epanechnikov> for rma"));
   }#if

   ## check for correct exonlevel
   exlevel <- exonLevel(exonlevel, chiptype);

   ## check parameter list
   if (length(params) == 0) {
      stop(paste("empty parameter list", sQuote("params")));
   }#if

   ## get treenames to normalize as fullnames=/datadir/treenames
   listnames <- listTreeNames(object);
   treenames <- listnames$treenames;
   fullnames <- listnames$fullnames;
   numtrees  <- length(treenames);

   ## define setname and settype for new treeset
   setname <- "BackgroundSet";
   settype <- "preprocess";

   ## background correct
   r <- .C("BgCorrect",
           as.character(filename),
           as.character(filedir),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.character(select),
           as.character(method),
           as.character(option),
           as.integer(length(params)),
           as.double(params),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.integer(update),
           as.integer(exlevel[1]),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("BgCorrect")));
      return(NULL);
   }#if

   ## get treenames after background correction
   alltrees  <- getTreeNames(rootfile, "*", setname);
   treenames <- alltrees[grep("\\.int",alltrees)];
   bgrdtrees <- alltrees[setdiff(1:length(alltrees), grep("\\.int",alltrees))];
   numtrees  <- length(treenames);
   numbgrds  <- length(bgrdtrees);
   if (numbgrds != numtrees) {
      stop(paste("number of bgrdtrees <", numbgrds,
                 "> is not equal to <", numtrees,">",sep=""));
   }#if

   ## create new class DataTreeSet
   set <- new("DataTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              params    = as.list(params));

   ## add background trees
   set@bgtreenames <- as.list(bgrdtrees);

   return(set);
}#bgCorrect.DataTreeSet

setMethod("xpsBgCorrect", "DataTreeSet", bgCorrect.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"normalize.DataTreeSet" <-
function(object,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "all",
         method    = "mean",
         option    = "all",
         logbase   = "0",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(0.02, 0),
         add.data  = TRUE,
         verbose   = TRUE)
{
   if (debug.xps()) print("------normalize.DataTreeSet------")

   ## get schemefile and chiptype for object
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chiptype   <- chipType(object);

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## update: need to set filename to rootfile
   if (!is.logical(update)) {
      stop(paste(sQuote("update"), "must be TRUE or FALSE"));
   } else if (update == TRUE) {
      filename <- rootfile;
   } else if (update == FALSE && existsROOTFile(rootfile)) {
      ## check if root file exists (necessary for WinXP to test already here)
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for presence of valid selector option
   TYPE <- c("pmonly", "mmonly", "both", "all");
   if (is.na(match(select, TYPE))) {
      stop(paste(sQuote("select"), "is not a valid selector option"));
   }#if

   ## check for valid normalization method
   TYPE <- c("mean", "median", "quantile");
   if (is.na(match(method, TYPE))) {
      cat("methods <lowess,supsmu> have not been tested yet for DataTreeSet\n");
      stop(paste(sQuote("method"), "must be one of <mean,median,quantile>"));
   }#if

   ## check for valid normalization "option:logbase"
   option  <- validOption(option);
   logbase <- validLogbase(logbase);
   option  <- paste(option, logbase, sep=":");

   ## check exon level
   exonlevel <- exonLevel(exonlevel, chiptype);

   ## get treenames to normalize as fullnames=/datadir/treenames
   listnames <- listTreeNames(object);
   treenames <- listnames$treenames;
   fullnames <- listnames$fullnames;
   numtrees  <- length(treenames);

   ## reference tree
   reftree <- "";
   if (!(refindex >= 0 && refindex <= numtrees)) {
      stop(paste(sQuote("refindex"), "must be 0 or less than number of trees"));
   }#if

   setname <- object@setname;
   if (refindex == 0) {
      reftree <- paste(setname, "/*.", extenPart(treenames), sep="");
   } else {
      reftree <- paste(setname, treenames[[refindex]], sep="/");
   }#if

   ## reference method
   TYPE <- c("mean", "median");
   if (is.na(match(refmethod, TYPE))) {
      stop(paste(sQuote("refmethod"), "must be <mean, median>"));
   }#if

   ## check parameter list
   if (length(params) == 0) {
      stop(paste("empty parameter list", sQuote("params")));
   }#if

   ## define setname and settype for new treeset
   setname <- "NormationSet";
   settype <- "preprocess";

   ## normalize
   r <- .C("Normalize",
           as.character(filename),
           as.character(filedir),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.character(select),
           as.character(method),
           as.character(option),
           as.integer(length(params)),
           as.double(params),
           as.integer(exonlevel[2]),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.character(reftree),
           as.character(refmethod),
           as.integer(update),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("Normalize")));
      return(NULL);
   }#if

   if (add.data) tmp <- NULL; #dummy only

   ## get treenames with extension for method
   treenames <- getTreeNames(rootfile, type2Exten(method, settype), setname);
   numtrees  <- length(treenames);

   ## create new class DataTreeSet
   set <- new("DataTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              params    = as.list(params));

   return(set);
}#normalize.DataTreeSet

setMethod("xpsNormalize", "DataTreeSet", normalize.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"summarize.DataTreeSet" <-
function(object,
         filename   = character(),
         filedir    = getwd(),
         tmpdir     = "",
         update     = FALSE,
         select     = "none",
         method     = character(),
         option     = "transcript",
         logbase    = "0",
         exonlevel  = "",
         params     = list(),
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   if (debug.xps()) print("------summarize.DataTreeSet------")

   ## get schemefile and chipname
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chipname   <- chipName(object);
   chiptype   <- chipType(object);

   ## check for presence of alternative root scheme file
   if ((!is.null(xps.scheme)) &&
       is(xps.scheme, "SchemeTreeSet") &&
       (chipType(xps.scheme) == chiptype) &&
       (file.exists(rootFile(xps.scheme)))) {
      scheme     <- xps.scheme;
      schemefile <- rootFile(xps.scheme);
      chipname   <- chipName(xps.scheme);
      chiptype   <- chipType(xps.scheme);
   }#if

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## update: need to set filename to rootfile
   if (!is.logical(update)) {
      stop(paste(sQuote("update"), "must be TRUE or FALSE"));
   } else if (update == TRUE) {
      filename <- rootfile;
   } else if (update == FALSE && existsROOTFile(rootfile)) {
      ## check if root file exists (necessary for WinXP to test already here)
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for presence of valid selector option
   TYPE <- c("pmonly", "mmonly", "both", "all", "none");
   if (is.na(match(select, TYPE))) {
      stop(paste(sQuote(select), "is not a valid selector option"));
   }#if

   ## check for valid summarization method
   TYPE <- c("avgdiff", "tukeybiweight", "medianpolish", "farms", "dfw");
   if (is.na(match(method, TYPE))) {
      stop(paste(sQuote("method"),
           "must be one of <avgdiff,tukeybiweight,medianpolish,farms,dfw>"));
   }#if

   ## check for valid option "transcript:logbase"
   transcript <- validTranscriptOption(option);
   logbase    <- validLogbase(logbase);
   sumoption  <- paste(transcript, logbase, sep=":");

   ## check parameter list
   if (length(params) == 0) {
      stop(paste("empty parameter list", sQuote("params")));
   }#if

   ## get treenames to summarize as fullnames=/datadir/treenames
   listnames <- listTreeNames(object);
   treenames <- listnames$treenames;
   fullnames <- listnames$fullnames;
   numtrees  <- length(treenames);

   ## check exon level
   exonlevel <- exonLevel(exonlevel, chiptype);

   ## define setname and settype for new treeset
   setname <- "SummarySet";
   settype <- "preprocess";

   ## summarize
   r <- .C("Summarize",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.character(select),
           as.character(method),
           as.character(sumoption),
           as.integer(length(params)),
           as.double(params),
           as.integer(exonlevel[3]),
           as.character(setname),
           as.character(fullnames),
           as.integer(numtrees),
           as.integer(update),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("Summarize")));
      return(NULL);
   }#if

   ## get extension
   exten <- type2Exten(method, settype);

   ## export result to outfile and import as dataframe
   ds <- data.frame(matrix(nr=0,nc=0));
   if (add.data) {
      outfile  <- sub("\\.root", ".txt", rootfile);
      ## get treename "treeset.treename.treetype"
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

      if (file.exists(outfile)) {
         ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("error: could not export results as", sQuote(outfile)));
      }#if
   }#if

   ## expression treenames 
   treenames <- getTreeNames(rootfile, exten, setname);
   numtrees  <- length(treenames);

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
              params    = list(),
              exprtype  = "custom",
              normtype  = "none");

   return(set);
}#summarize.DataTreeSet

setMethod("xpsSummarize", "DataTreeSet", summarize.DataTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("pmplot", signature(x="DataTreeSet"),
   function(x,
           which   = "",
           size    = 0,
           transfo = NULL,
           method  = mean,
           names   = "namepart",
           beside  = TRUE,
           col     = c("red", "blue"),
           legend  = c("PM","MM"),
           las     = 2,
           ylab    = "mean intensities",
           ...) 
{
      if (debug.xps()) print("------pmplot.DataTreeSet------")

      if (chipType(x) == "GeneChip") {
         pm <- validData(x, which="pm");
         mm <- validData(x, which="mm");
      } else if (chipType(x) == "GenomeChip") {
         pm <- validData(x, which=which);
         mm <- validData(x, which="antigenomic");
      } else if (chipType(x) == "ExonChip") {
         pm <- validData(x, which=which);
         mm <- validData(x, which="antigenomic");
      }#if

      if (is.null(names)) {
         names <- colnames(pm);
      } else if (names[1] == "namepart") {
         names <- namePart(colnames(pm));
      } else {
         pm <- pm[, names, drop=F];
         mm <- mm[, names, drop=F];
      }#if

      if (size > 1) {
         pm <- pm[seq(1,nrow(pm),len=size), , drop=F];
         mm <- mm[seq(1,nrow(mm),len=size), , drop=F];
      }#if

      if (is.function(transfo)) {
         pm <- transfo(pm);
         mm <- transfo(mm);
      }#if

      if (is.function(method)) {
         pm <- apply(pm, 2, method);
         mm <- apply(mm, 2, method);
      }#if

      barplot(t(cbind(pm, mm)),
              beside    = beside,
              col       = col,
              names.arg = names,
              legend    = legend,
              las       = las,
              ylab      = ylab,
              ...);
   }
)#pmplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("image", signature(x="DataTreeSet"),
   function(x,
            bg      = FALSE,
            transfo = log,
            col     = gray((0:64)/64),
            names   = "namepart",
            xlab    = "",
            ylab    = "",
            ...) 
{
      if (debug.xps()) print("------image.DataTreeSet------")

      ## get intensities from data or bgrd
      if (bg == FALSE)          ds <- validData(x, which="")
      else                      ds <- validBgrd(x, which="");
      if (is.function(transfo)) ds <- transfo(ds);

      if (is.null(names))              names <- colnames(ds)
      else if (names[1] == "namepart") names <- namePart(colnames(ds))
      else                             ds    <- ds[, names, drop=F];

      ## plot images
      if (ncol(ds) > 1 && interactive()) par(ask=TRUE) else par(ask=FALSE);
      for (i in 1:ncol(ds)) {
         m <- as.numeric(ds[,i]);
         m <- matrix(m, ncol=ncols(x), nrow=nrows(x));
         m <- m[,nrows(x):1];

         image(m,
               col  = col,
               main = names[i],
               xlab = xlab,
               ylab = ylab,
               xaxt = 'n',
               yaxt = 'n',
                ...);
      }#for
      par(ask=FALSE);
   }
)#image

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

