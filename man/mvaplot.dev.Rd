\name{mvaplot.dev}
\alias{mvaplot.dev}
\title{MvA Scatter Plot for Device}
\description{
  Produce scatter plots of M values vs A values of the samples for the selected device.
}
\usage{
mvaplot.dev(x, transfo = log2, method = "median", names = "namepart", ylim = c(-6, 6), xlab = "A", ylab = "M", pch = ".", mar = c(3, 3, 2, 1), dev = "screen", outfile = "MvAPlot", w = 540, h = 540, ...)
}
\arguments{
  \item{x}{object of class \code{\link{ExprTreeSet}}.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{method}{method to compute M, \dQuote{mean} or \dQuote{median}.}
  \item{names}{optional vector of sample names.}
  \item{ylim}{range for the plotted M values.}
  \item{xlab}{a title for the x axis.}
  \item{ylab}{a title for the y axis.}
  \item{pch}{either an integer specifying a symbol or a single character to be used in plotting points.}
  \item{mar}{plot margin.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{\dots}{optional arguments to be passed to \code{plot}.}
}
\details{
 Produces mvaplots for slot \code{data} for an object of class \code{\link{ExprTreeSet}}
 for the selected graphics device.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as mvaplot.
}
\author{Christian Stratowa}
\seealso{\code{\link{mvaplot}}}
\keyword{device}
