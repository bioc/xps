\name{image-methods}
\docType{methods}
\alias{image-methods}
\alias{image}
\title{Display an Image}
\description{
 Creates an image for each sample.

\emph{Usage}

\code{
image(x, bg = FALSE, transfo = log2, col = gray((0:64)/64), names = "namepart", xlab = "", ylab = "", ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}}.}
  \item{bg}{logical. If \code{FALSE}, intensities from slot \code{data} will be used;
   if \code{TRUE}, background intensities from slot \code{bgrd} will be used.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{col}{color range for intensities.}
  \item{names}{optional vector of sample names.}
  \item{\dots}{optional arguments to be passed to \code{image}.}
}
\details{
 Creates an image for each array. For \code{bgrd=TRUE} the distribution of the computed
 background intensities will be shown; this can be useful to see potential density gradients
 caused by hybridization conditions. For the computation of background intensities see function
 \code{\link{bgcorrect}}; it is suggested to use \code{\link{bgcorrect.mas4}} to identify
 density gradients.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as image.
}
\author{Christian Stratowa}
\note{
 Data must first be attached to class \code{\link{DataTreeSet}} using method \code{\link{attachInten}}.
}
\seealso{\code{\link{image.dev}}, \code{\link{image}}}
\keyword{methods}
