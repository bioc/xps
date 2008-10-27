#==============================================================================#
# root.merge.data.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# root.merge.data: merge ROOT files containing trees from CEL-files
#==============================================================================#

"root.merge.data" <- 
function(xps.scheme,
         rootfiles = list(),
         celnames  = "*") 
{
   if (debug.xps()) print("------root.merge.data------")

   ## check for presence of class xps.scheme
   xps.scheme <- validSchemeTreeSet(xps.scheme);

   ## extract chipname from xps.scheme
   chipname <- chipName(xps.scheme);
   setname  <- "DataSet";

   ## check if root files contain trees for chipname
   for (rootfile in rootfiles) {
      treetitle <- getTreeNames(rootfile, treetype = "cel", setname = setname, gettitle = TRUE)[1];
      if (treetitle != chipname) {
         stop(paste("ROOT file", sQuote(rootfile), "has no trees with chipname", sQuote(chipname)));
      }#if
   }#for

   ## get correct treenames from rootfiles
   celnames  <- unique(celnames);
   tmpnames  <- vector();
   fullnames <- vector();

   for (rootfile in rootfiles) {
      treenames <- getTreeNames(rootfile, treetype = "cel", setname = setname, gettitle = FALSE);

      if (length(celnames) == 0) next;
      if (celnames[1] == "*") {
         id <- match(tmpnames, treenames);
         id <- id[!is.na(id)];
         tmpnames <- c(tmpnames, treenames);
         if (length(id) > 0) {
            tmpnames  <- tmpnames[-id];
            treenames <- treenames[-id];
         }#if

         if (length(treenames) > 0) {
            fullnames <- c(fullnames, paste(rootfile, setname, treenames, sep="/"));
         }#if
      } else {
         id <- match(treenames, celnames);
         id <- id[!is.na(id)];

         if (length(id) > 0) {
            treenames <- celnames[id];
            celnames  <- celnames[-id];
            tmpnames  <- c(tmpnames, treenames);
            fullnames <- c(fullnames, paste(rootfile, setname, treenames, sep="/"));
         }#if
      }#if
   }#for

   ## check for presence of celnames
   if (celnames[1] != "*" && length(celnames) > 0) {
      stop(paste(sQuote(celnames), "not found in", sQuote("rootfiles"), "\n"));
   }#if

   ## create new class
   set <- new("DataTreeSet",
              setname   = setname,
              settype   = "rawdata",
              rootfile  = basename(rootfiles[1]),
              filedir   = dirname(rootfiles[1]),
              numtrees  = length(fullnames),
              treenames = as.list(fullnames),
              scheme    = xps.scheme);

   return(set);
}#root.merge.data
