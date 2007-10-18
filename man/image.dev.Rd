\name{image.dev}
\alias{image.dev}
\title{Image for Device}
\description{
  Creates an image for each sample for the selected device.
}
\usage{
image.dev(x, bg = FALSE, transfo = log2, col = gray((0:64)/64), names = "namepart", xlab = "", ylab = "", mar = c(1, 1, 2, 1), dev = "screen", outfile = "Image", w = 540, h = 540, ...)
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}}.}
  \item{bg}{logical. If \code{FALSE}, intensities from slot \code{data} will be used;
   if \code{TRUE}, background intensities from slot \code{bgrd} will be used.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{col}{color range for intensities.}
  \item{names}{optional vector of sample names.}
  \item{xlab}{a title for the x axis.}
  \item{ylab}{a title for the y axis.}
  \item{mar}{plot margins.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{\dots}{optional arguments to be passed to \code{image}.}
}
\details{
 Creates an image for each array for the selected graphics device.

 For \code{bgrd=TRUE} the distribution of the computed background intensities will be shown;
 this can be useful to see potential density gradients caused by hybridization conditions.
 For the computation of background intensities see function \code{\link{bgcorrect}};
 it is suggested to use \code{\link{bgcorrect.mas4}} to identify density gradients.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as image.
}
\author{Christian Stratowa}
\note{
 Data must first be attached to class \code{\link{DataTreeSet}} using method \code{\link{attachInten}}.
}
\seealso{\code{\link{image-methods}}, \code{\link{image}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## need to attach scheme mask and data
data.test3 <- attachMask(data.test3)
data.test3 <- attachInten(data.test3)

if (interactive()) {
image.dev(data.test3)
}

## to avoid memory comsumption of R remove data:
data.test3 <- removeInten(data.test3)
data.test3 <- removeMask(data.test3)
}
\keyword{device}
