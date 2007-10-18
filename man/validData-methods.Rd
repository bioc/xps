\name{validData-methods}
\docType{methods}
\alias{validData}
\title{Extract Subset of Data}
\description{
 Extracts a subset of valid data from data.frame \code{data}.

\emph{Usage}

\code{
validData(object, which = "")
}
}
\arguments{
  \item{object}{object of class \code{DataTreeSet}.}
  \item{which}{type of probes to be returned.}
}
\details{
For expression arrays, \code{validData} returns all the perfect match or mismatch  
probes on the arrays the object represents as data.frame, i.e. \code{which} can have 
the following values:
\tabular{lll}{
  \tab \code{pm}:\tab perfect match probes. \cr
  \tab \code{mm}:\tab mismatch probes. \cr
  \tab \code{both}:\tab both perfect match and mismatch probes. 
}

For exon arrays, \code{validData} returns the probes of the different exon levels as 
data.frame, i.e. \code{which} can have one of the following values:
\tabular{lll}{
  \tab \code{core}:\tab probesets supported by RefSeq and full-length GenBank transcripts. \cr
  \tab \code{metacore}:\tab core meta-probesets. \cr
  \tab \code{extended}:\tab probesets with other cDNA support. \cr
  \tab \code{metaextended}:\tab extended meta-probesets. \cr
  \tab \code{full}:\tab probesets supported by gene predictions only. \cr
  \tab \code{metafull}:\tab full meta-probesets. \cr
  \tab \code{affx}:\tab standard AFFX controls. \cr
  \tab \code{all}:\tab combination of above. \cr
  \tab \code{genomic}:\tab genomic background probes. \cr
  \tab \code{antigenomic}:\tab antigenomic background probes.
}
}
\value{
  A \code{\link{data.frame}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{pm}}, \code{\link{mm}}}
\keyword{methods}
