#==============================================================================#
# import.data.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# import.data: import CEL-files
#==============================================================================#

"import.data" <-
function(xps.scheme,
         filename = character(0),
         filedir  = getwd(),
         celdir   = NULL,
         celfiles = "*",
         celnames = NULL,
         project  = NULL,
         verbose  = TRUE)
{
   ## check for presence of class xps.scheme
   xps.scheme <- validSchemeTreeSet(xps.scheme);

   ## extract scheme, chipname and chiptype from xps.scheme
   settype  <- setType(xps.scheme);
   scheme   <- rootFile(xps.scheme);
   chipname <- chipName(xps.scheme);
   chiptype <- chipType(xps.scheme);

   ## check for presence of parameters
   if (!(is.character(filename) &&
         is.character(filedir))) {
      stop("arguments <filename,filedir> must be of type character");
   }#if

   ## root data file
   rootdata <- paste(filename, "_cel.root", sep="");
   rootfile <- rootDirFile(rootdata, filedir);

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
#      celnames <- "";
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

   ## project information
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

   ## define treeset
   setname <- "DataSet";

   ## create root data file
   r <- .C("ImportData",
           as.character(filename),
           as.character(filedir),
           as.character(chiptype),
           as.character(chipname),
           as.character(scheme),
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
           as.integer(0),
           as.integer(0),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ImportData")));
      return(NULL);
   }#if

   celnames <- paste(celnames, ".cel", sep="");

   ## create new class
   set <- new("DataTreeSet",
              setname   = setname,
              settype   = "rawdata",
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = as.numeric(ncel),
              treenames = as.list(celnames),
              scheme    = xps.scheme);

   ## add project information
   if (!is.null(project) && is(project, "ProjectInfo")) {
      projectInfo(set) <- project;
   }#if

   return(set);
}#import.data
