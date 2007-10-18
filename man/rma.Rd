\name{rma}
\alias{rma}
\alias{xpsRMA-methods}
\alias{xpsRMA}
\title{Robust Multi-Array Average expression measure}
\description{
 This function converts a \code{\link{DataTreeSet}} into an \code{\link{ExprTreeSet}}
 using the robust multi-array average (RMA) expression measure.
}
\usage{
rma(xps.data,
    filename   = character(0),
    filedir    = getwd(),
    tmpdir     = "",
    background = "pmonly",
    normalize  = TRUE,
    option     = "transcript",
    exonlevel  = "",
    xps.scheme = NULL,
    verbose    = TRUE)


xpsRMA(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{\link{DataTreeSet}}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{background}{probes used to compute background, one of \sQuote{pmonly}, \sQuote{mmonly}, \sQuote{both}}
  \item{normalize}{logical. If \code{TRUE} normalize data using quantile normalization.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization;
   exon/genome arrays only.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{the arguments described above.}
}
\details{
This function computes the RMA (Robust Multichip Average) expression measure described in 
Irizarry et al. for both expression arrays and exon arrays. For exon arrays it is necessary 
to supply the requested \code{option} and \code{exonlevel}.

Following \code{option}s are valid for exon arrays:
\tabular{ll}{
  \code{transcript}: \tab expression levels are computed for transcript clusters,
   i.e. probe sets containing the same \sQuote{transcript\_cluster\_id}. \cr
  \code{exon}: \tab expression levels are computed for exon clusters,
   i.e. probe sets containing the same \sQuote{exon\_id}, where each exon cluster
   consists of one or more \code{probeset}s. \cr
  \code{probeset}: \tab expression levels are computed for individual probe sets,
   i.e. for each \sQuote{probeset\_id}. \cr
}
Following \code{exonlevel} annotations are valid for exon arrays:
\tabular{lll}{
  \tab \code{core}:\tab probesets supported by RefSeq and full-length GenBank transcripts. \cr
  \tab \code{metacore}:\tab core meta-probesets. \cr
  \tab \code{extended}:\tab probesets with other cDNA support. \cr
  \tab \code{metaextended}:\tab extended meta-probesets. \cr
  \tab \code{full}:\tab probesets supported by gene predictions only. \cr
  \tab \code{metafull}:\tab full meta-probesets. \cr
  \tab \code{affx}:\tab standard AFFX controls. \cr
  \tab \code{all}:\tab combination of above.
}
Following \code{exonlevel} annotations are valid for whole genome arrays:
\tabular{lll}{
  \tab \code{core}:\tab probesets with category \sQuote{unique} and \sQuote{mixed}. \cr
  \tab \code{metacore}:\tab probesets with category \sQuote{unique} only. \cr
  \tab \code{affx}:\tab standard AFFX controls. \cr
  \tab \code{all}:\tab combination of above.
}
Exon levels can also be combined, with following combinations being most useful:
\tabular{ll}{
  \code{exonlevel="metacore+affy"}:  \tab core meta-probesets plus AFFX controls \cr
  \code{exonlevel="core+extended"}: \tab probesets with cDNA support \cr
  \code{exonlevel="core+extended+full"}:  \tab supported plus predicted probesets \cr
}

Exon level annotations are described in the Affymetrix whitepaper exon\_probeset\_trans\_clust\_whitepaper.pdf: \cr 
\dQuote{Exon Probeset Annotations and Transcript Cluster Groupings}.

In order to use an alternative \code{\link{SchemeTreeSet}} set the corresponding SchemeSet \code{xps.scheme}.

\code{xpsRMA} is the \code{DataSet} method called by function \code{rma}, 
containing the same parameters.
}
\value{
 An \code{\link{ExprTreeSet}}
}
\author{Christian Stratowa}
\note{
In contrary to other implementations of RMA the expression measure is given to you in linear scale,
analogously to the expression measures computed with \code{\link{mas5}} and \code{\link{mas4}}. 

For the analysis of many exon arrays it is strongly recommended to define a \code{tmpdir} 
unless you know that your computer has sufficient RAM.
}
\references{
 Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie M. Cope, Bridget Hobbs
 and Terence P. Speed (2003), Summaries of Affymetrix GeneChip probe level data Nucleic
 Acids Research 31(4):e15

 Bolstad, B.M., Irizarry R. A., Astrand M., and Speed, T.P. (2003), A Comparison of
 Normalization Methods for High Density Oligonucleotide Array Data Based on Bias and
 Variance. Bioinformatics 19(2):185-193

 Irizarry, RA, Hobbs, B, Collin, F, Beazer-Barclay, YD, Antonellis, KJ, Scherf, U,
 Speed, TP (2003) Exploration, Normalization, and Summaries of High Density Oligonucleotide
 Array Probe Level Data. Biostatistics .Vol. 4, Number 2: 249-264 
}
\seealso{\code{\link{express}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

data.rma <- rma(data.test3,"tmp_Test3RMA",tmpdir="",background="pmonly",normalize=TRUE,verbose=FALSE)

## get data.frame
expr.rma <- validData(data.rma)
head(expr.rma)

## plot results
if (interactive()) {
boxplot(data.rma)
boxplot(log2(expr.rma))
}

rm(scheme.test3, data.test3)
gc()
}
\keyword{manip}
