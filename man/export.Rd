\name{export}
\alias{export.scheme}
\alias{export.data}
\alias{export.expr}
\alias{export.call}
\alias{export-methods}
\alias{export}
\title{Export data as text files}
\description{
  Export data from classes \code{\link{SchemeTreeSet}}, \code{\link{DataTreeSet}},
 \code{\link{ExprTreeSet}}, or \code{\link{CallTreeSet}} to \code{outfile}.
}
\usage{
export.scheme(xps.scheme, treetype = character(0), varlist = "*", outfile = character(0), sep = "\t", as.dataframe = FALSE, verbose = TRUE)

export.data(xps.data, treename = "*", treetype = "cel", varlist = "*", outfile = character(0), sep = "\t", as.dataframe = FALSE, verbose = TRUE)

export.expr(xps.expr, treename = "*", treetype = character(0), varlist = "*", outfile = character(0), sep = "\t", as.dataframe = FALSE, verbose = TRUE)

export.call(xps.call, treename = "*", treetype = character(0), varlist = "*", outfile = character(0), sep = "\t", as.dataframe = FALSE, verbose = TRUE)

export(object, ...)
}
\arguments{
  \item{xps.scheme}{an object of type \code{\link{SchemeTreeSet}}.}
  \item{xps.data}{an object of type \code{\link{DataTreeSet}}.}
  \item{xps.expr}{an object of type \code{\link{ExprTreeSet}}.}
  \item{xps.call}{an object of type \code{\link{CallTreeSet}}.}
  \item{treename}{vector of tree names to export.}
  \item{treetype}{type of tree(s) to export, see \code{\link{validTreetype}}}
  \item{varlist}{names of tree leaves to export}
  \item{outfile}{name of output file.}
  \item{sep}{column separator}
  \item{as.dataframe}{if \code{TRUE} a data.frame will be returned.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{arguments \code{treenames},\code{treetype},\code{varlist},\code{outfile},\code{sep},\code{as.dataframe}.}
}
\details{
  Export data from classes \code{\link{SchemeTreeSet}}, \code{\link{DataTreeSet}},
 \code{\link{ExprTreeSet}}, or \code{\link{CallTreeSet}} to \code{outfile}.

\code{export} is a generic method to export data from \code{\link{ROOT}} trees as text file.
}
\value{
If \code{as.dataframe} is \code{TRUE}, the data will be imported into the
current R session as \code{data.frame}. Otherwise, \code{NULL} will be returned.
}
\author{Christian Stratowa}
\seealso{\code{\link{export-methods}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## export as table only
export(scheme.test3, treetype="idx", outfile="Test3_idx.txt",verbose=FALSE)

## export as table and import as data.frame
ann <- export.scheme(scheme.test3, treetype="ann", outfile="Test3_ann.txt",as.dataframe=TRUE,verbose=FALSE)
head(ann)
data <- export.data(data.test3, outfile="Test3_cel.txt",as.dataframe=TRUE,verbose=FALSE)
head(data)
}
\keyword{manip}
