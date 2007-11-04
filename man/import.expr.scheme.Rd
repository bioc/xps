\name{import.expr.scheme}
\alias{import.expr.scheme}
\title{Import CDF, probe and annotation files into a SchemeTreeSet}
\description{
Import the Affymetrix CDF, probe and annotation files into a ROOT file and create S4 class SchemeTreeSet
}
\usage{
import.expr.scheme(filename   = character(0),
                   filedir    = getwd(),
                   schemefile = character(0),
                   probefile  = character(0),
                   annotfile  = character(0),
                   chipname   = NULL,
                   add.mask   = FALSE,
                   verbose    = TRUE)
}
\arguments{
  \item{filename}{file name of ROOT scheme file.}
  \item{filedir}{system directory where ROOT scheme file should be stored.}
  \item{schemefile}{name of CDF-file, including full path.}
  \item{probefile}{name of probe-file, including full path.}
  \item{annotfile}{name of annotation-file, including full path.}
  \item{chipname}{optional chip name when using an alternative CDF-file.}
  \item{add.mask}{logical. If \code{TRUE} mask information will be included as slot \code{mask}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
\code{import.expr.scheme} is used to import all information for an Affymetrix 
expression array into a \code{\link{ROOT}} scheme file, including CDF-file, the corresponding
probe file, and the current Afymetrix annotation file. 

Usually, \code{chipname} is extracted from the name of the CDF-file, however, when using
an alternative CDF-file, e.g. from BrainArray or AffyProbeMiner, a \code{chipname} must be
supplied which starts with (or contains) the exact Affymetrix chip name.

An S4 class \code{\link{SchemeTreeSet}} will be created, serving as R wrapper to the 
\code{\link{ROOT}} scheme file \code{filename}.

Since a new \code{\link{ROOT}} scheme file needs only to be created when a new annotation file
is available from the Affymetrix website, it is recommended to store all \code{\link{ROOT}}
scheme files in a commonly accessible system directory \code{filedir}.

Use function \code{\link{root.scheme}} to access the \code{\link{ROOT}} scheme file from new
R sessions to avoid creating a new \code{\link{ROOT}} scheme file for every session.
}
\value{
  A \code{SchemeTreeSet} object.
}
\author{Christian Stratowa}
\note{As mentioned above, use function \code{\link{root.scheme}} to access the 
 \code{\link{ROOT}} scheme file from new R sessions to avoid creating a new \code{\link{ROOT}} 
 scheme file for every R session.

 Do not separate \code{filename} of ROOT files with dots, use underscores, e.g. do not use 
 \code{filename="Scheme.Test3.na21"} but use \code{filename="Scheme_Test3_na21"} or simply
 \code{filename="SchemeTest3na21"} instead. Extension \dQuote{root} is added automatically,
 so that ROOT is able to recognize the file as ROOT file.
}
\seealso{\code{\link{import.exon.scheme}}, \code{\link{import.genome.scheme}},
 \code{\link{root.scheme}}, \code{\link{SchemeTreeSet}}}
\examples{
\dontrun{
## define paths
scmdir <- "/common/path/schemes"
libdir <- "/my/path/Affy/libraryfiles"
anndir <- "/my/path/Affy/Annotation"

## create scheme for Test3 GeneChip
scheme.test3.na22 <- import.expr.scheme("Scheme_Test3_na22",filedir=scmdir,
                     schemefile=paste(libdir,"Test3.CDF",sep="/"),
                     probefile=paste(libdir,"Test3_probe.tab",sep="/"),
                     annotfile=paste(anndir,"Test3.na22.annot.csv",sep="/"))

## access ROOT scheme file from new R session
scheme.test3 <- root.scheme(paste(scmdir,"Scheme_Test3_na22.root",sep="/"))

## create scheme for HG-U133_Plus_2 GeneChip
scheme.hgu133p2.na22 <- import.expr.scheme("Scheme_HGU133p2_na22",filedir=scmdir,
                        schemefile=paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),
                        probefile=paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),
                        annotfile=paste(anndir,"HG-U133_Plus_2.na22.annot.csv",sep="/"))
## access ROOT scheme file from new R session
scheme.hgu133p2 <- root.scheme(paste(scmdir,"Scheme_HGU133p2_na22.root",sep="/"))
}
}
\keyword{manip}
