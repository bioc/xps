\name{firma.score}
\alias{firma.score}
\title{Get Splice Score from FIRMA}
\description{
 Extracts the FIRMA splice score from data.frame \code{data}.
}
\usage{
firma.score(xps.data,
            probeset = NULL,
            option   = "probeset")
}
\arguments{
  \item{xps.data}{object of class \code{\link{ExprTreeSet}}.}
  \item{probeset}{probesetID or NULL.}
  \item{option}{option determining the \code{probeset} type for which to extract expression levels, 
   one of \sQuote{probeset}, \sQuote{exon}.}
}
\details{
Function \code{firma.score} returns the FIRMA splice score described in Purdom et al. from slot \code{data} 
for a given \code{probeset}, or for all probesets in case of \code{probeset=NULL}. Row names will be the Affymetrix 
probesetIDs or exonIDs, respectively, dependent on the selected \code{option}.
}
\value{
  A \code{\link{data.frame}}.
}
\author{Christian Stratowa}
\references{
 Purdom, E., Simpson K.M., Robinson M.D., Conboy J.G., Lapuk A.V. and Speed, T.P. (2008),
 FIRMA: a method for detection of alternative splicing from exon array data.
 Bioinformatics 24(15):1707-1714
}
\seealso{\code{\link{firma}}}
\keyword{manip}
