\name{hist-methods}
\docType{methods}
\alias{hist-methods}
\alias{hist}
\title{Density Estimate}
\description{
  Plot the density estimates for each sample.

\emph{Usage}

\code{
hist(x, which = "", size = 0, transfo = log2, ylab = "density", xlab = "log intensity", type = "l", col = 1:6, ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{ylab}{a title for the y axis.}
  \item{xlab}{a title for the x axis.}
  \item{type}{type for the plot.}
  \item{col}{colors to use for the different arrays.}
  \item{\dots}{optional arguments to be passed to \code{plot}.}
}
\details{
 Plots the non-parametric density estimates using values contained in the columns of 
 slot \code{data}.

 For a \code{\link{DataTreeSet}} object, data must first be attached using method 
 \code{\link{attachInten}}.
}
\author{Christian Stratowa}
\note{
 For exon array raw data only a limited number of samples can be displayed as density plot due to
 memory limitations. To display all samples it is proposed to use function \code{\link{root.density}}
 instead.
}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## need to attach scheme mask and probe intensities
data.test3 <- attachMask(data.test3)
data.test3 <- attachInten(data.test3)

if (interactive()) {
hist(data.test3)
}

## optionally remove mask and data to free memory
data.test3 <- removeInten(data.test3)
data.test3 <- removeMask(data.test3)
}
\keyword{methods}
