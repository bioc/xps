\name{indexUnits-methods}
\docType{methods}
\alias{indexUnits-methods}
\alias{indexUnits}
\alias{pmindex-methods}
\alias{pmindex}
\alias{mmindex-methods}
\alias{mmindex}
\title{Unit Locations}
\description{
 Returns a data.frame or list with locations of the probes in each probe set. 

\emph{Usage}

\code{indexUnits(object, which = "", unitid = NULL, as.list = TRUE, data = NULL)}

\code{pmindex(object, unitid = NULL, as.list = TRUE)}

\code{mmindex(object, unitid = NULL, as.list = TRUE)}
}
\arguments{
  \item{object}{Object of class \code{"DataTreeSet"}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{unitid}{optional vector of UNIT_IDs.}
  \item{as.list}{if TRUE a \code{list} will be returned (default is \code{data.frame}).}
  \item{data}{optional \code{data.frame} containing (x,y)-coordinates.}
}
\details{
 Function \code{indexUnits} returns a \code{data.frame} or \code{list} with locations of the probes in
 each probe set. 

 By default a \code{data.frame} for selected \code{unitid}s or all unitIDs (\code{unitid="*"}) will be
 returned with columns <UNIT_ID, X, Y, XY>. Here "XY" are the selected rows of slot \code{data}.
 
 For \code{as.list=TRUE} a \code{list} of unitIDs will be returned containing the selected rows "XY". The
 names of the elements in the list returned are the \code{UNIT_ID}s.

 For \code{unitid=NULL} a vector of \code{data} rows "XY" will be returned. 

 For expression arrays which can be one of "pm", "mm", or "both". Alternatively, functions \code{pmindex}
 and \code{mmindex} can be used for PM probes or MM probes, respectively. 

 For exon arrays \code{which} is described in \code{\link{validData}}. However, in this case slot \code{data}
 must contain the (x,y)-coordinates of the probesetIDs.
}
\value{
  A \code{list} or \code{data.frame}.
}
\author{Christian Stratowa}
\seealso{\code{\link{unitID2transcriptID}}, \code{\link{unitID2probesetID}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## dataXY not attached
id <- indexUnits(data.test3, which="pm", unitid=c(34,36,122))
id

## dataXY attached
data.test3 <- attachDataXY(data.test3)
xy <- treeData(data.test3)
id <- indexUnits(data.test3, which="pm", unitid=c(34,36,122), data=xy)
id
id <- indexUnits(data.test3, which="", unitid=c(34,36,122), data=xy)
id
id <- indexUnits(data.test3, which="", unitid=34, as.list=FALSE, data=xy)
id
data.test3 <- removeDataXY(data.test3)

rm(scheme.test3, data.test3)
gc()
}
\keyword{methods}
