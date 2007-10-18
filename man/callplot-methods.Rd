\name{callplot-methods}
\docType{methods}
\alias{callplot-methods}
\alias{callplot}
\title{Barplot of Percent Present and Absent Calls.}
\description{
Creates a barplot of percent Present/Marginal/Absent calls.

\emph{Usage}

\code{
callplot(x, beside = TRUE, names = "namepart", col = c("red","green","blue"), legend = c("P","M","A"), ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{CallTreeSet}}.}
  \item{beside}{logical. If \code{FALSE}, the columns of height are portrayed as stacked bars,
   and if \code{TRUE} the columns are portrayed as juxtaposed bars.}
  \item{names}{optional vector of sample names.}
  \item{col}{color for P/M/A bars}
  \item{legend}{legend for the plot, defaults to P/M/A.}
  \item{\dots}{optional arguments to be passed to \code{barplot}.}
}
\details{
Creates a barplot of percent Present/Marginal/Absent calls.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as callplot.
}
\author{Christian Stratowa}
\seealso{\code{\link{pmplot}}}
\keyword{methods}
