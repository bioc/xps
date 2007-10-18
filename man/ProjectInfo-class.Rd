\name{ProjectInfo-class}
\docType{class}
\alias{ProjectInfo-class}
\alias{ProjectInfo}
\alias{arrayInfo,ProjectInfo-method}
\alias{arrayInfo}
\alias{arrayInfo<-,ProjectInfo,character-method}
\alias{arrayInfo<-}
\alias{authorInfo,ProjectInfo-method}
\alias{authorInfo}
\alias{authorInfo<-,ProjectInfo,character-method}
\alias{authorInfo<-}
\alias{datasetInfo,ProjectInfo-method}
\alias{datasetInfo}
\alias{datasetInfo<-,ProjectInfo,character-method}
\alias{datasetInfo<-}
\alias{projectInfo,ProjectInfo-method}
\alias{projectInfo}
\alias{projectInfo<-,ProjectInfo,character-method}
\alias{projectInfo<-}
\alias{show,ProjectInfo-method}
\alias{sourceInfo,ProjectInfo-method}
\alias{sourceInfo}
\alias{sourceInfo<-,ProjectInfo,character-method}
\alias{sourceInfo<-}

\title{Class ProjectInfo}
\description{This class allows to save the relevant project information in the 
\code{\link{ROOT}} data file and in class \code{\link{DataTreeSet}}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \cr
\code{new("ProjectInfo", submitter=[character], laboratory=[character], contact=[character], ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{submitter}:}{Object of class \code{"character"} representing the name of the submitter.}
    \item{\code{laboratory}:}{Object of class \code{"character"} representing the laboratory of the submitter.}
    \item{\code{contact}:}{Object of class \code{"character"} representing the contact address of the submitter.}
    \item{\code{project}:}{Object of class \code{"list"} representing the project information.}
    \item{\code{author}:}{Object of class \code{"list"} representing the author information.}
    \item{\code{dataset}:}{Object of class \code{"list"} representing the dataset information.}
    \item{\code{source}:}{Object of class \code{"list"} representing the sample source information.}
    \item{\code{arraytype}:}{Object of class \code{"list"} representing the array information.}
  }
}
\section{Methods}{
  \describe{
    \item{arrayInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{arraytype}.}
    \item{arrayInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{arraytype} with character vector c(chipname,chiptype,description).}
    \item{authorInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{author}.}
    \item{authorInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{author} with character vector c(lastname,firstname,company,lab,email,phone).}
    \item{datasetInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{dataset}.}
    \item{datasetInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{dataset} with character vector c(name,type,sample,submitter,date,description).}
    \item{projectInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{project}.}
    \item{projectInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{project} with character vector c(name,date,description).}
    \item{show}{\code{signature(object = "ProjectInfo")}: shows the content of \code{ProjectInfo}.}
    \item{sourceInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{source}.}
    \item{sourceInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{source} with character vector c(name,species,subspecies,description).}
	 }
}
\author{Christian Stratowa}
\examples{
project <- new("ProjectInfo",submitter="Christian", laboratory="home",contact="email")
projectInfo(project) <- c("TestProject","20060106","use Test3 data for testing")
authorInfo(project)  <- c("Stratowa","Christian","Home","Lab","mymail","1234567")
datasetInfo(project) <- c("Test3Set","MC","Tissue","Stratowa","20060106","short description")
sourceInfo(project)  <- c("Unknown","Homo sapiens","caucasian","short description")
arrayInfo(project)   <- c("Test3","GeneChip","description of array")
show(project)
}
\keyword{classes}
