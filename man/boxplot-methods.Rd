\name{boxplot-methods}
\docType{methods}
\alias{boxplot-methods}
\alias{boxplot}
\title{Box Plots}
\description{
  Produce box-and-whisker plot(s) of the samples.

\emph{Usage}

\code{
boxplot(x, which = "", size = 0, transfo = log2, range = 0, names = "namepart", ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{names}{optional vector of sample names.}
  \item{\dots}{optional arguments to be passed to \code{boxplot}.}
}
\details{
 Creates a boxplot for slot \code{data} for an object of class \code{\link{DataTreeSet}}
 or \code{\link{ExprTreeSet}}.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as boxplot.
}
\author{Christian Stratowa}
\note{
 For a \code{\link{DataTreeSet}} object, data must first be attached using method 
 \code{\link{attachInten}}.
}
\seealso{\code{\link{boxplot.dev}}, \code{\link[graphics]{boxplot}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## need to attach scheme mask and probe intensities
data.test3 <- attachMask(data.test3)
data.test3 <- attachInten(data.test3)

if (interactive()) {
boxplot(data.test3)
}

## optionally remove mask and data to free memory
data.test3 <- removeInten(data.test3)
data.test3 <- removeMask(data.test3)

rm(scheme.test3, data.test3)
gc()
}
\keyword{methods}
