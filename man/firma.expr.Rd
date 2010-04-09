\name{firma.expr}
\alias{firma.expr}
\title{Get Expression Levels from FIRMA}
\description{
 Extracts FIRMA expression levels from data.frame \code{data}.
}
\usage{
firma.expr(xps.data,
           probeset = NULL,
           option   = "probeset")
}
\arguments{
  \item{xps.data}{object of class \code{\link{ExprTreeSet}}.}
  \item{probeset}{transcriptID or probesetID or NULL.}
  \item{option}{option determining the \code{probeset} type for which to extract expression levels, 
   one of \sQuote{transcript}, \sQuote{probeset}, \sQuote{exon}.}
}
\details{
Function \code{firma.expr} returns the expression levels from slot \code{data} for a given \code{probeset},
or for all probesets or transcripts in case of \code{probeset=NULL}. Row names will be the Affymetrix 
transcriptIDs, probesetIDs or exonIDs, respectively, dependent on the selected \code{option}.
}
\value{
  A \code{\link{data.frame}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{firma}}}
\keyword{manip}
