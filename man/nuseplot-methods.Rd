\name{nuseplot-methods}
\docType{methods}
\alias{nuseplot-methods}
\alias{nuseplot}
\title{Box Plots of Normalized Unscaled Standard Errors (NUSE)}
\description{
  Produce boxplot of Normalized Unscaled Standard Errors (NUSE) for the set of arrays.

\emph{Usage}

\code{
nuseplot(x, which = "UnitName", size = 0, range = 0, ylim = c(0.9,1.2), outline = FALSE, names = "namepart", las = 2, ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{ylim}{range for the plotted y values.}
  \item{outline}{if \code{outline} is not true, the outliers are not drawn.}
  \item{names}{optional vector of sample names.}
  \item{las}{the style of axis labels.}
  \item{\dots}{optional arguments to be passed to \code{boxplot}.}
}
\details{
 Create boxplots of Normalized Unscaled Standard Errors (NUSE) for the set of arrays.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as boxplot.
}
\author{Christian Stratowa}
\seealso{\code{\link{rleplot}}}
\examples{
# load existing ROOT scheme file and ROOT expression file for rma
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.rma <- root.expr(scheme.test3, paste(.path.package("xps"),"rootdata/tmp_Test3RMA.root",sep="/"), "mdp")

if (interactive()) {
nuseplot(data.rma)
}
}
\keyword{methods}
