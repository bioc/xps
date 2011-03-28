\name{AffyRNAdeg}
\alias{AffyRNAdeg}
\alias{summaryAffyRNAdeg}
\alias{plotAffyRNAdeg}
\alias{xpsRNAdeg-methods}
\alias{xpsRNAdeg}
\title{Functions to assess RNA Degradation.}
\description{
 Functions to detect possible RNA degradation.
}
\usage{
AffyRNAdeg(xps.data, treename = "*", qualopt = "raw", log.it = TRUE)

summaryAffyRNAdeg(rna.deg, signif.digits=3)

plotAffyRNAdeg(rna.deg, transform = "shift.scale", cols = NULL, ...)

xpsRNAdeg(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{\link{QualTreeSet}}.}
  \item{treename}{vector of tree names to export.}
  \item{qualopt}{option determining the data to which to apply qualification,
   one of \sQuote{raw}, \sQuote{adjusted}, \sQuote{normalized}.}
  \item{log.it}{logical, if TRUE, then probe data is log2 transformed.}
  \item{rna.deg}{\code{list}, output from \code{AffyRNAdeg}.}
  \item{signif.digits}{number of significant digits to show.}
  \item{transform}{transform data before plotting, one of "shift.scale",
   "shift.only", "none". }
  \item{cols}{vector of colors for plot, length is number of trees.}
  \item{object}{object of class \code{QualTreeSet}.}
  \item{\dots}{optional arguments to be passed to \code{plotAffyRNAdeg}.}
}
\details{
 Since probes within a probeset are ordered directionally from the 5' end to
 the 3' end, it is possible to estimate the quality (degradation status) of
 the RNA.

 Function \code{AffyRNAdeg} averages the probe intensities by location in the
 probeset, with the average taken over all probesets with identical number of
 probes. 

 Function \code{summaryAffyRNAdeg} produces a single summary statistic for 
 each array.

 Function \code{plotAffyRNAdeg} produces a side-by-side plot of the averaged
 intensities. Option \code{transform = "none"} shows the averaged intensities
 for each array while option "shift" staggers the plots for individual arrays
 vertically to make the display easier to read, and option "scale" normalizes
 the averaged intensities so that the standard deviation is equal to one.
}
\value{
  \code{AffyRNAdeg} returns a \code{list} with following components:
  \item{N}{number of probesets with identical number of probes}
  \item{sample.names }{names of samples, derived from affy batch object}
  \item{mns}{average intensity by probe position}
  \item{ses}{standard errors for probe position averages}
  \item{slope}{from linear regression of means.by.number}
  \item{pvalue}{from linear regression of means.by.number}
}
\author{Christian Stratowa, adapted from package affy}
\examples{
\dontrun{
rnadeg <- xpsRNAdeg(rlm.all, treename="*", qualopt="raw")
plotAffyRNAdeg(rnadeg)

rnadeg <- AffyRNAdeg(rlm.all)
result <- summaryAffyRNAdeg(rnadeg)
plotAffyRNAdeg(rnadeg)
}
}
\keyword{manip}
