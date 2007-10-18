\name{mvaplot-methods}
\docType{methods}
\alias{mvaplot-methods}
\alias{mvaplot}
\title{M vs A Plot}
\description{
  Produce scatter plots of M values vs A values of the samples.

\emph{Usage}

\code{
mvaplot(x, transfo = log2, method  = "median", names = "namepart", ylim = c(-6,6), ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{ExprTreeSet}}.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{method}{method to compute M, \dQuote{mean} or \dQuote{median}.}
  \item{names}{optional vector of sample names.}
  \item{ylim}{range for the plotted M values.}
  \item{\dots}{optional arguments to be passed to \code{plot}.}
}
\details{
 Produces mvaplots for slot \code{data} for an object of class \code{\link{ExprTreeSet}}.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as mvaplot.
}
\author{Christian Stratowa}
\seealso{\code{\link{mvaplot.dev}}}
\keyword{methods}
