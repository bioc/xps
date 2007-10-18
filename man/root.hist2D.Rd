\name{root.hist2D}
\alias{root.hist2D}
\title{ROOT 2D-Histogram}
\description{
  Creates a ROOT 2D-histogram for a ROOT tree.
}
\usage{
root.hist2D(x, treename1 = character(0), treename2 = character(0),logbase = "log2", option = "COLZ", canvasname = "Histogram2D", w = 540, h = 540)
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{treename1}{name of first tree, must be present in \code{rootfile} of object \code{x}.}
  \item{treename2}{name of second tree, must be present in \code{rootfile} of object \code{x}.}
  \item{logbase}{usually \dQuote{log2}, or \dQuote{0}, determines if leaf data should be converted to log.}
  \item{option}{ROOT hist TH1 option, usually one of \dQuote{SCAT}, \dQuote{COLZ}, \dQuote{BOX}.}
  \item{canvasname}{name of ROOT canvas}
  \item{w}{the width of the canvas in pixels.}
  \item{h}{the height of the canvas in pixels.}
}
\details{
 Creates a ROOT 2D-histogram for trees \code{treename1} and \code{treename2} present in \code{rootfile}.

 By selecting menu \dQuote{File->Save->canvasname.xxx} you can save the figure as e.g. *gif, *.jpg, *.pdf, *.ps or even as C++ macro.
}
\author{Christian Stratowa}
\note{
Always select menu item \dQuote{Quit ROOT} from menu \dQuote{File} to close the ROOT browser,
otherwise you are in the CINT C/C++ interpreter from \code{\link{ROOT}}. To exit CINT,
you need to type \dQuote{.q}.
}
\seealso{\code{\link{root.hist1D}}, \code{\link{root.hist3D}}}
\examples{
\dontrun{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

root.hist2D(data.test3, "TestA1.cel", "TestB1.cel", option="COLZ")
}
}
\keyword{device}
