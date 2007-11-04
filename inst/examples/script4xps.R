#------------------------------------------------------------------------------#
# Script: step-by-step functions to demonstrate how to use package "xps"
#
# Copyright (c) 2007-2007 Christian Stratowa, Vienna, Austria.
# All rights reserved.
#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# 1. step: import Affymetrix chip definition and annotation files into
#          ROOT scheme files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT scheme files need only be created when new annotation files
#          are available. In order to allow all users to access the ROOT 
#          scheme files it is recommended to create all ROOT scheme files
#          in a common system directory and use function "root.scheme()"
#          to load the desired scheme into the current R session to analyse
#          a project.
#------------------------------------------------------------------------------#

### new R session: load library xps
library(xps)

### define directories:
# directory containing Affymetrix library files
libdir <- "/Volumes/GigaDrive/Affy/libraryfiles"
# directory containing Affymetrix annotation files
anndir <- "/Volumes/GigaDrive/Affy/Annotation"
# directory to store ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# note: do not separate name of ROOT files with dots, use underscores,
#       e.g. do not use "Scheme.Test3.na21" but "Scheme_Test3_na21"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na21 <- import.expr.scheme("Scheme_Test3_na21",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version06Nov/Test3.na21.annot.csv",sep="/"))
scheme.test3.na22 <- import.expr.scheme("Scheme_Test3_na22",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version07May/Test3.na22.annot.csv",sep="/"))
scheme.test3.na23 <- import.expr.scheme("Scheme_Test3_na23",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version07Jul/Test3.na23.annot.csv",sep="/"))

# Hu6800:
scheme.hu6800.na21 <- import.expr.scheme("Scheme_Hu6800_na21",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version06Nov/Hu6800.na21.annot.csv",sep="/"))
scheme.hu6800.na22 <- import.expr.scheme("Scheme_Hu6800_na22",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version07May/Hu6800.na22.annot.csv",sep="/"))
scheme.hu6800.na23 <- import.expr.scheme("Scheme_Hu6800_na23",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version07Jul/Hu6800.na23.annot.csv",sep="/"))

# HG_U95A:
scheme.hgu95a <- import.expr.scheme("Scheme_HGU95A",filedir=scmdir,paste(libdir,"HG_U95A.CDF",sep="/"),paste(libdir,"HG-U95A_probe.tab",sep="/"),annotfile="")

# HG_U95Av2:
scheme.hgu95av2.na21 <- import.expr.scheme("Scheme_HGU95Av2_na21",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG_U95Av2.na21.annot.csv",sep="/"))
scheme.hgu95av2.na22 <- import.expr.scheme("Scheme_HGU95Av2_na22",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version07May/HG_U95Av2.na22.annot.csv",sep="/"))
scheme.hgu95av2.na23 <- import.expr.scheme("Scheme_HGU95Av2_na23",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG_U95Av2.na23.annot.csv",sep="/"))

# HG-U133A:
scheme.hgu133a.na21 <- import.expr.scheme("Scheme_HGU133A_na21",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG-U133A.na21.annot.csv",sep="/"))
scheme.hgu133a.na22 <- import.expr.scheme("Scheme_HGU133A_na22",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version07May/HG-U133A.na22.annot.csv",sep="/"))
scheme.hgu133a.na23 <- import.expr.scheme("Scheme_HGU133A_na23",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG-U133A.na23.annot.csv",sep="/"))

# HG-U133B:
scheme.hgu133b.na21 <- import.expr.scheme("Scheme_HGU133B_na21",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG-U133B.na21.annot.csv",sep="/"))
scheme.hgu133b.na22 <- import.expr.scheme("Scheme_HGU133B_na22",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version07May/HG-U133B.na22.annot.csv",sep="/"))
scheme.hgu133b.na23 <- import.expr.scheme("Scheme_HGU133B_na23",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG-U133B.na23.annot.csv",sep="/"))

# HG-U133_Plus_2:
scheme.hgu133p2.na21 <- import.expr.scheme("Scheme_HGU133p2_na21",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG-U133_Plus_2.na21.annot.csv",sep="/"))
scheme.hgu133p2.na22 <- import.expr.scheme("Scheme_HGU133p2_na22",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version07May/HG-U133_Plus_2.na22.annot.csv",sep="/"))
scheme.hgu133p2.na23 <- import.expr.scheme("Scheme_HGU133p2_na23",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG-U133_Plus_2.na23.annot.csv",sep="/"))

# alternative CDF-files, e.g.:
#AffyProbeMiner
apmdir <- "/Volumes/GigaDrive/Affy/CDF_alternative/AffyProbeMiner"
scheme.hgu133p2.apm <- import.expr.scheme("Scheme_HGU133p2_apm_refseq",filedir=scmdir,paste(apmdir,"HG-U133_Plus_2_transcript_refseq/HG-U133_Plus_2_transcript_refseq.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),"","HG-U133_Plus_2_apm_refseq")
#BrainArray (UniMichigan)
umidir <- "/Volumes/GigaDrive/Affy/CDF_alternative/UniMichigan"
scheme.hgu133p2.umi <- import.expr.scheme("Scheme_HGU133p2_umi_refseq_v10",filedir=scmdir,paste(umidir,"Version10/Hs133P_Hs_REFSEQ/Hs133P_Hs_REFSEQ.cdf",sep="/"),paste(umidir,"Version10/Hs133P_Hs_REFSEQ/Hs133P_Hs_REFSEQ_probe.tab",sep="/"),"","HG-U133_Plus_2_umi_refseq")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r3:
scheme.hugene10stv1r3.na23 <- import.genome.scheme("Scheme_HuGene10stv1r3_na23",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version07Jul/HuGene-1_0-st-v1.na23.hg18.transcript.csv",sep="/"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Note: not possible on the Compaq Laptop with 512 MB RAM

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2r2.na21 <- import.exon.scheme("Scheme_HuEx10stv2r2_na21",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version06Nov/HuEx-1_0-st-v2.na21.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version06Nov/HuEx-1_0-st-v2.na21.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na22 <- import.exon.scheme("Scheme_HuEx10stv2r2_na22",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version07Mar/HuEx-1_0-st-v2.na22.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version07Mar/HuEx-1_0-st-v2.na22.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na23 <- import.exon.scheme("Scheme_HuEx10stv2r2_na23",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version07Jul/HuEx-1_0-st-v2.na23.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version07Jul/HuEx-1_0-st-v2.na23.hg18.transcript.csv",sep="/"))

# HuEx-1_0-st-v2.r2 old annotation:
scheme.huex10stv2r2.old <- import.exon.scheme("Scheme_HuEx10stv2r2_old",filedir=scmdir,
                           paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                           paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                           paste(anndir,"Version06Jun/HuEx-1_0-st-probeset-annot.csv",sep="/"),
                           paste(anndir,"Version06Jun/HuEx-1_0-st-transcript-annot.csv",sep="/"),
                           paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.controls.ps",sep="/"))

# MoEx-1_0-st-v1.r2:
scheme.moex10stv1r2.na21 <- import.exon.scheme("Scheme_MoEx10stv2r2_na21",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version06Nov/MoEx-1_0-st-v1.na21.mm8.probeset.csv",sep="/"),
                            paste(anndir,"Version06Nov/MoEx-1_0-st-v1.na21.mm8.transcript.csv",sep="/"))

scheme.moex10stv1r2.na22 <- import.exon.scheme("Scheme_MoEx10stv2r2_na22",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07May/MoEx-1_0-st-v1.na22.mm8.probeset.csv",sep="/"),
                            paste(anndir,"Version07May/MoEx-1_0-st-v1.na22.mm8.transcript.csv",sep="/"))

scheme.moex10stv1r2.na23 <- import.exon.scheme("Scheme_MoEx10stv2r2_na23",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07Jul/MoEx-1_0-st-v1.na23.mm8.probeset.csv",sep="/"),
                            paste(anndir,"Version07Jul/MoEx-1_0-st-v1.na23.mm8.transcript.csv",sep="/"))

# RaEx-1_0-st-v1.r2:
scheme.raex10stv1r2.na21 <- import.exon.scheme("Scheme_RaEx10stv2r2_na21",filedir=scmdir,
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version06Nov/RaEx-1_0-st-v1.na21.rn4.probeset.csv",sep="/"),
                            paste(anndir,"Version06Nov/RaEx-1_0-st-v1.na21.rn4.transcript.csv",sep="/"))

scheme.raex10stv1r2.na22 <- import.exon.scheme("Scheme_RaEx10stv2r2_na22",filedir=scmdir,
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07May/RaEx-1_0-st-v1.na22.rn4.probeset.csv",sep="/"),
                            paste(anndir,"Version07May/RaEx-1_0-st-v1.na22.rn4.transcript.csv",sep="/"))

scheme.raex10stv1r2.na23 <- import.exon.scheme("Scheme_RaEx10stv2r2_na23",filedir=scmdir,
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07Jul/RaEx-1_0-st-v1.na23.rn4.probeset.csv",sep="/"),
                            paste(anndir,"Version07Jul/RaEx-1_0-st-v1.na23.rn4.transcript.csv",sep="/"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# utility functions to demonstrate how to access scheme files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### export different trees from ROOT scheme file
# Test3: export as table only
export(scheme.test3.na23, treetype="idx", outfile="Test3_idx.txt")
export(scheme.test3.na23, treetype="scm", outfile="Test3_scm.txt")
export(scheme.test3.na23, treetype="prb", outfile="Test3_prb.txt")
export(scheme.test3.na23, treetype="ann", outfile="Test3_ann.txt")

# export as table and import as data.frame
idx <- export(scheme.test3.na23, treetype="idx", outfile="Test3_idx.txt",as.dataframe=T)
ann <- export(scheme.test3.na23, treetype="ann", outfile="Test3_ann.txt",as.dataframe=T)

### attach mask later: if import parameter was: as.dataframe=FALSE
scheme.test3.na23 <- attachMask(scheme.test3.na23)
str(scheme.test3.na23)
### export scheme mask
msk <- chipMask(scheme.test3.na23)
scheme.test3.na23 <- removeMask(scheme.test3.na23)
str(scheme.test3.na23)

### scheme accessors
rootFile(scheme.test3.na23)
chipName(scheme.test3.na23)
chipType(scheme.test3.na23)
probeInfo(scheme.test3.na23)

### browse ROOT scheme files
root.browser(scheme.test3.na23)


#------------------------------------------------------------------------------#
# 2. step: import CEL-files into ROOT data files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT raw data files need only be created once for each new project.
#          They can be updated later to add additional CEL-files to a project.
#          In order to allow all users to access the ROOT data files it is
#          recommended to create all ROOT data files in a common system
#          directory and use function "root.data()" to load the desired data
#          into the current R session to analyse a project.
#          To distinguish ROOT raw data files from other ROOT files, "_cel" is
#          added to each raw data file, e.g. "DataTest3_cel.root".
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 1: Test3 samples from package xps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### define directories:
# directory containing Test3 CEL files
celdir <- paste(.path.package("xps"),"raw",sep="/")

### import raw data
# first, import ROOT scheme file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
# import CEL files
data.test3 <- import.data(scheme.test3, "DataTest3", celdir=celdir)
str(data.test3)

# alternatively import CEL files and information about the project
project <- new("ProjectInfo",submitter="Christian", laboratory="home",contact="email")
projectInfo(project) <- c("TestProject","20060106","use Test3 data for testing")
authorInfo(project)  <- c("Stratowa","Christian","Home","Lab","mymail","1234567")
datasetInfo(project) <- c("Test3Set","MC","Tissue","Stratowa","20060106","description")
sourceInfo(project)  <- c("Unknown","Homo sapiens","caucasian","description")
arrayInfo(project)   <- c("Test3","GeneChip","description")
# need to delete old ROOT file "DataTest3_cel.root" first!
data.test3 <- import.data(scheme.test3,"DataTest3",filedir="./data",celdir=celdir,project=project)
str(data.test3)

### update raw data files
# import first set of CEL files
data.test3 <- import.data(scheme.test3, "NewDataTest3", celdir=celdir, celfiles=c("TestA1.CEL","TestA2.CEL"))
str(data.test3)
# import additional CEL-files and update ROOT data file
data.test3 <- addData(data.test3, celdir=celdir, celfiles=c("TestB1.CEL","TestB2.CEL"))
str(data.test3)

### to get only subset of trees from an exisiting ROOT data file
subdata.test3 <- root.data(scheme.test3, rootFile(data.test3), c("TestA1","TestB2","TestA2"))
str(subdata.test3)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration how to load raw data, access the data, and plot the data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### 1.load existing ROOT scheme file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))

### 2.load existing ROOT data file
rootfile <- paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/")
data.test3 <- root.data(scheme.test3, rootfile=rootfile)

# alternatively load only subset from a ROOT data file into a different R session
subdata.test3 <- root.data(scheme.test3, rootfile=rootfile, c("TestA1","TestB2","TestA2"))

### two possibilities to get tree data as data.frame
# 1.possibility: need to attach data first to "data.test3"
# attach intensity data
data.test3 <- attachInten(data.test3)
# get data
tmp <- intensity(data.test3)
head(tmp)
# also possible to attach only subset
data.test3 <- attachInten(data.test3, c("TestB1.cel","TestA2"))
# get data
tmp <- intensity(data.test3)
head(tmp)
# to avoid memory comsumption of R do:
data.test3 <- removeInten(data.test3)

# 2.possibility:
tmp <- getTreeData(data.test3)
head(tmp)
tmp <- getTreeData(data.test3, varlist="fStdev")
head(tmp)
tmp <- getTreeData(data.test3, varlist="fNPixels")
head(tmp)
tmp <- getTreeData(data.test3, varlist="fInten:fStdev")
head(tmp)

### export different trees from ROOT data file
# export tree data
export(data.test3, treetype="cel", outfile="TestAB_cel.txt")
export(data.test3, treetype="cel", varlist = "fInten", outfile="TestAB_int_cel.txt")
export(data.test3, treename="TestB1", treetype="cel", varlist = "fInten", outfile="TestB1_int_cel.txt")
export(data.test3, treename="TestB1", treetype="cel", varlist = "fInten:fStdev", outfile="TestB1_int_sd__cel.txt")
export(data.test3, treename="TestB1", treetype="cel", varlist = "*", outfile="TestB1_cel.txt")
export(data.test3, treename=c("TestB1","TestA2"), treetype="cel", varlist = "fInten", outfile="TestB1A2_int_cel.txt")

# export data directly from root file
schemefile <- paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/")
datafile   <- paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/")
export.root(datafile, schemefile, "DataSet", "*", "cel", "*", "DataOutFile")

# inspect ROOT file with ROOT browser (to quit ROOT, type ".q")
root.browser(data.test3)

### plot raw data
# need to attach scheme mask, since it was not attached to scheme
data.test3 <- attachMask(data.test3)
# need to attach data first
data.test3 <- attachInten(data.test3)
str(data.test3)

# plots
hist(data.test3)
image(data.test3)
boxplot(data.test3)
mboxplot(data.test3, ylim=c(-6,6))
pmplot(data.test3)

# plots to export
image.dev(data.test3)
image.dev(data.test3, dev="png", outfile="Image_DataTest3")
image.dev(data.test3, dev="png", outfile="Image_DataTest3_TestA1",names="TestA1.cel_MEAN")
boxplot.dev(data.test3)
boxplot.dev(data.test3, dev="png", outfile="Boxplot_DataTest3")
boxplot.dev(data.test3, dev="jpeg", outfile="Boxplot_DataTest3")
boxplot.dev(data.test3, dev="pdf", outfile="Boxplot_DataTest3")

# to avoid memory comsumption of R remove data:
data.test3 <- removeInten(data.test3)
data.test3 <- removeMask(data.test3)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 2: Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### define directories:
# directory of ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
# directory containing Tissues CEL files
celdir <- "/Volumes/GigaDrive/ChipData/Exon/HuMixture"
# directory to store ROOT raw data files
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"

### HG-U133_Plus_2 data: import raw data
# first, import ROOT scheme file
scheme.u133p2 <- root.scheme(paste(scmdir,"Scheme_HGU133p2_na23.root",sep="/"))

# subset of CEL files to import
celfiles <- c("u1332plus_ivt_breast_A.CEL","u1332plus_ivt_breast_B.CEL","u1332plus_ivt_breast_C.CEL",
              "u1332plus_ivt_prostate_A.CEL","u1332plus_ivt_prostate_B.CEL","u1332plus_ivt_prostate_C.CEL")
# rename CEL files
celnames <- c("BreastA","BreastB","BreastC","ProstateA","ProstateB","ProstateC")
# import CEL files
data.mix.u133p2 <- import.data(scheme.u133p2, "HuTissuesU133P2", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 3: Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# same R session for example

### HuEx-1_0-st-v2 data: import raw data
# first, import ROOT scheme file
scheme.exon <- root.scheme(paste(scmdir,"Scheme_HuEx10stv2r2_na23.root",sep="/"))

# subset of CEL files to import
celfiles <- c("huex_wta_breast_A.CEL","huex_wta_breast_B.CEL","huex_wta_breast_C.CEL",
              "huex_wta_prostate_A.CEL","huex_wta_prostate_B.CEL","huex_wta_prostate_C.CEL")
# rename CEL files
celnames <- c("BreastA","BreastB","BreastC","ProstateA","ProstateB","ProstateC")
# import CEL files
data.mix.exon <- import.data(scheme.exon, "HuTissuesExon", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 4: Tissues from Affymetrix Exon Array Dataset for HuGene-1_0-st-v1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# same R session for example

# directory containing Tissues CEL files
celdir <- "/Volumes/GigaDrive/ChipData/Exon/HuGene"

### HuGene-1_0-st-v1 data: import raw data
# first, import ROOT scheme file
scheme.genome <- root.scheme(paste(scmdir,"Scheme_HuGene10stv1r3_na23.root",sep="/"))

# subset of CEL files to import
celfiles <- c("TisMap_Breast_01_v1_WTGene1.CEL","TisMap_Breast_02_v1_WTGene1.CEL","TisMap_Breast_03_v1_WTGene1.CEL",
              "TisMap_Prostate_01_v1_WTGene1.CEL","TisMap_Prostate_02_v1_WTGene1.CEL","TisMap_Prostate_03_v1_WTGene1.CEL")
# rename CEL files
celnames <- c("Breast01","Breast02","Breast03","Prostate01","Prostate02","Prostate03")
# import CEL files
data.mix.genome <- import.data(scheme.genome, "HuTissuesGenome", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration how to access the data, and plot the data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

# import ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.u133p2 <- root.scheme(paste(scmdir,"Scheme_HGU133p2_na23.root",sep="/"))
scheme.exon   <- root.scheme(paste(scmdir,"Scheme_HuEx10stv2r2_na23.root",sep="/"))
scheme.genome <- root.scheme(paste(scmdir,"Scheme_HuGene10stv1r3_na23.root",sep="/"))

# import ROOT data files
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.u133p2 <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))
data.exon   <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))
data.genome <- root.data(scheme.genome, paste(datdir,"HuTissuesGenome_cel.root",sep="/"))

# inspect ROOT file with ROOT browser
root.browser(data.exon)


### plot raw data for HG-U133_Plus_2
# need to attach scheme mask, since it was not attached to scheme
data.u133p2 <- attachMask(data.u133p2)
# need to attach data 
data.u133p2 <- attachInten(data.u133p2)
str(data.u133p2)

# plots
hist(data.u133p2)
image(data.u133p2)
boxplot(data.u133p2)
mboxplot(data.u133p2, ylim=c(-6,6))
pmplot(data.u133p2)

# plots to export
image.dev(data.u133p2)
image.dev(data.u133p2, dev="png", col=rainbow(32), outfile="Image_DataMixU133P2")
image.dev(data.u133p2, dev="png", col=rainbow(32), outfile="Image_DataMixU133P2_BrA",names="BreastA.cel_MEAN")
boxplot.dev(data.u133p2)
boxplot.dev(data.u133p2, dev="png", w=600, h=480, outfile="Boxplot_DataMixU133P2")

# to avoid memory comsumption of R remove data:
data.u133p2 <- removeInten(data.u133p2)
data.u133p2 <- removeMask(data.u133p2)


### plot raw data for HuEx-1_0-st-v2
# On the PowerBook with 1GB RAM and the Compaq Laptop with 512 MB I did not test this
# need to attach data first
data.exon <- attachInten(data.exon)
#data.exon <- attachInten(data.exon, treenames=c("BreastA","BreastB","ProstateA","ProstateB"))
# need to attach scheme mask, since it was not attached to scheme
data.exon <- attachMask(data.exon)

# Note: On my MacBook Pro with 2 GB RAM it was necessary to use "core" values only and
#       use only a subset of size=100000 rows, otherwise R returned memory errors:
#       "Error: cannot allocate vector of size 50.0 Mb"

# plots
hist(data.exon, which="core",size=100000)
image(data.exon)
boxplot(data.exon, which="core",size=100000)
mboxplot(data.exon, which="core",size=100000, ylim=c(-6,6))
pmplot(data.exon, which="core",size=100000)

# plots to export
image.dev(data.exon)
image.dev(data.exon, dev="png",col=rainbow(32), outfile="Image_DataMixExon")
image.dev(data.exon, dev="png",col=rainbow(32), outfile="Image_DataMixExon_BreastA",names="BreastA.cel_MEAN")
boxplot.dev(data.exon, which="core",size=100000)
boxplot.dev(data.exon, which="core",size=100000, dev="png",w=600, h=480, outfile="Boxplot_DataMixExon")

# to avoid memory comsumption of R remove data:
data.exon <- removeInten(data.exon)
data.exon <- removeMask(data.exon)
gc()


### plot raw data for HuGene-1_0-st-v1
# need to attach scheme mask, since it was not attached to scheme
data.genome <- attachMask(data.genome)
# need to attach data 
data.genome <- attachInten(data.genome)

# plots
hist(data.genome)
image(data.genome)
boxplot(data.genome)
mboxplot(data.genome, ylim=c(-6,6))
pmplot(data.genome)
pmplot(data.genome, which="core")

# plots to export
image.dev(data.genome)
image.dev(data.genome, dev="png", col=rainbow(32), outfile="Image_DataMixHuGene")
image.dev(data.genome, dev="png", col=rainbow(32), outfile="Image_DataMixHuGene_Br01",names="Breast01.cel_MEAN")
boxplot.dev(data.genome)
boxplot.dev(data.genome, dev="png", w=600, h=480, outfile="Boxplot_DataMixHuGene")

# to avoid memory comsumption of R remove data:
data.genome <- removeInten(data.genome)
data.genome <- removeMask(data.genome)


#------------------------------------------------------------------------------#
# 3. step: convert raw data to expression levels
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT scheme files and ROOT raw data files are usually already stored
#          in special system directories. When a new R session is created for the
#          first time, they must fist be loaded using "root.scheme()" and "root.data()".
#          However, this is not necessary when re-opening a saved R session later.
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 1: Test3 samples from package xps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))


### preprocess raw data ###

# 1. RMA
data.rma <- rma(data.test3,"Test3RMA",tmpdir="",background="pmonly",normalize=T)

# 2. MAS5
data.mas5 <- mas5(data.test3,"Test3MAS5",,tmpdir="",normalize=T,sc=500)
# to store all trees (including e.g. background trees) in same ROOT file, use "update=T"
data.mas5 <- mas5(data.test3,"Test3MAS5All",,tmpdir="",normalize=T,sc=500, update=T)

# 3. MAS5 detection call
call.mas5 <- mas5.call(data.test3,"Test3Call",tmpdir="")

# 4. DABG detection call (yes, this is possible for expression arrays)
call.dabg <- dabg.call(data.test3,"Test3DABG")

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)


### plot results ###

# compare mas5 to rma
plot(expr.rma[,1],expr.mas5[,1])
plot(expr.rma[,1],expr.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot.dev(data.rma)
boxplot.dev(data.rma, dev="png", mar=c(4,4,1,1), w=480, h=480, outfile="BoxPlot_Test3_rma")

boxplot(data.mas5)
boxplot.dev(data.mas5)
boxplot.dev(data.mas5, dev="png", mar=c(4,4,1,1), w=480, h=480, outfile="BoxPlot_Test3_mas5")

# relative boxplots
mboxplot(data.rma)
mboxplot(data.rma, ylim=c(-2,3))
mboxplot(data.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot(data.rma, pch=20, ylim=c(-2,2))
mvaplot.dev(data.rma, pch=20, ylim=c(-2,2))
mvaplot.dev(data.rma, pch=20, ylim=c(-2,2), dev="png",outfile="MvAPlot_Test3_rma")
mvaplot.dev(data.rma, pch=20, ylim=c(-2,2), names="TestB1.mdp_LEVEL")
mvaplot.dev(data.rma, pch=20, ylim=c(-2,2), names=c("TestA2.mdp_LEVEL","TestB1.mdp_LEVEL"))
mvaplot.dev(data.rma, pch=20, ylim=c(-2,2), dev="png",outfile="MvAPlot_Test3_rma")
mvaplot.dev(data.rma, pch=20, ylim=c(-2,2), names="TestB1.mdp_LEVEL", dev="png",outfile="MvAPlot_TestB1_rma")

mvaplot(data.mas5, pch=20)
mvaplot.dev(data.mas5, pch=20)
mvaplot.dev(data.mas5, pch=20, names="TestB1.tmn_LEVEL")
mvaplot.dev(data.mas5, pch=20, names=c("TestA2.tmn_LEVEL","TestB1.tmn_LEVEL"))
mvaplot.dev(data.mas5, pch=20, dev="png",outfile="MvAPlot_Test3_mas5")
mvaplot.dev(data.mas5, pch=20, names="TestB1.tmn_LEVEL", dev="png",outfile="MvAPlot_TestB1_mas5")

# present call plots
callplot(call.mas5)
callplot(call.mas5, beside=F, ylim=c(0,125))
callplot(call.dabg)
callplot(call.dabg, beside=F, ylim=c(0,125))

# image for rma background intensities
# 1. find background treenames for data.rma
getTreeNames(rootFile(data.rma))
# 2. get background intensity for e.g. tree "TestA2.rbg"
bgrd <- export.root(rootFile(data.rma),schemeFile(data.rma),"PreprocesSet","TestA2","rbg","fBg","BgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
rbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.rma)), nrow=nrows(schemeSet(data.rma)))
# 4. create image
image(rbg)
image(log2(rbg))

# image for mas5 background intensities
# 1. find background treenames for data.mas5
getTreeNames(rootFile(data.mas5))
# 2. get background intensity for e.g. tree "TestA2.wbg"
bgrd <- export.root(rootFile(data.mas5),schemeFile(data.mas5),"PreprocesSet","TestA2","wbg","fBg","BgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
wbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.mas5)), nrow=nrows(schemeSet(data.mas5)))
# 4. create image
image(wbg)
image(log2(wbg))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 2: Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.u133p2 <- root.scheme(paste(scmdir,"Scheme_HGU133p2_na23.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.u133p2 <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))


### preprocess raw data ###

# 1. RMA
data.rma <- rma(data.u133p2,"MixU133P2RMA",tmpdir="",background="pmonly",normalize=T)

# 2. MAS5
data.mas5 <- mas5(data.u133p2,"MixU133P2MAS5",,tmpdir="",normalize=T,sc=500)
# to store all trees (including e.g. background trees) in same ROOT file, use "update=T"
data.mas5 <- mas5(data.u133p2,"MixU133P2MAS5All",,tmpdir="",normalize=T,sc=500, update=T)

# 3. MAS5 detection call
call.mas5 <- mas5.call(data.u133p2,"MixU133P2Call",tmpdir="")

# 4. DABG detection call
# Note: yes, it is possible to compute DABG call for expression arrays, but:
#       on my MacBook Pro with 2.3 GHz Intel Core 2 Duo it takes 44 min per array!!!
#call.dabg <- dabg.call(data.u133p2,"MixU133P2DABG")

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)


### plot results ###

# compare mas5 to rma
plot(expr.rma[,1],expr.mas5[,1])
plot(expr.rma[,1],expr.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot.dev(data.rma)
boxplot.dev(data.rma, dev="png", mar=c(4,4,1,1), w=480, h=480, outfile="BoxPlot_MixU133P2_rma")

boxplot(data.mas5)
boxplot.dev(data.mas5)
boxplot.dev(data.mas5, dev="png", mar=c(4,4,1,1), w=480, h=480, outfile="BoxPlot_MixU133P2_mas5")

# relative boxplots
mboxplot(data.rma)
mboxplot(data.rma, ylim=c(-2,3))
mboxplot(data.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot(data.rma, pch=20, ylim=c(-4,4))
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6))
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), dev="png",outfile="MvAPlot_MixU133P2_rma")
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), names="BreastA.mdp_LEVEL")
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), names=c("BreastA.mdp_LEVEL","BreastB.mdp_LEVEL"))
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), dev="png",outfile="MvAPlot_MixU133P2_rma")
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), names="BreastA.mdp_LEVEL", dev="png",outfile="MvAPlot_BreastA_rma")

mvaplot(data.mas5, pch=20)
mvaplot.dev(data.mas5, pch=20)
mvaplot.dev(data.mas5, pch=20, names="BreastA.tmn_LEVEL")
mvaplot.dev(data.mas5, pch=20, names=c("BreastA.tmn_LEVEL","BreastB.tmn_LEVEL"))
mvaplot.dev(data.mas5, pch=20, dev="png",outfile="MvAPlot_MixU133P2_mas5")
mvaplot.dev(data.mas5, pch=20, names="BreastA.tmn_LEVEL", dev="png",outfile="MvAPlot_BreastA_mas5")

# present call plots
callplot(call.mas5)
callplot(call.mas5, beside=F, ylim=c(0,125))
#callplot(call.dabg)

# image for rma background intensities
# 1. find background treenames for data.rma
getTreeNames(rootFile(data.rma))
# 2. get background intensity for e.g. tree "BreastA.rbg"
bgrd <- export.root(rootFile(data.rma),schemeFile(data.rma),"PreprocesSet","BreastA","rbg","fBg","BreastAU133BgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
rbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.rma)), nrow=nrows(schemeSet(data.rma)))
# 4. create image
image(rbg)
image(log2(rbg))

# image for mas5 background intensities
# 1. find background treenames for data.mas5
getTreeNames(rootFile(data.mas5))
# 2. get background intensity for e.g. tree "BreastA.wbg"
bgrd <- export.root(rootFile(data.mas5),schemeFile(data.mas5),"PreprocesSet","BreastA","wbg","fBg","BreastAU133BgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
wbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.mas5)), nrow=nrows(schemeSet(data.mas5)))
# 4. create image
image(wbg)
image(log2(wbg))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 3: Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.exon <- root.scheme(paste(scmdir,"Scheme_HuEx10stv2r2_na23.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.exon <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))


### preprocess raw data ###
datdir <- getwd()

# 1. RMA
# transcript: metacore
# ok for 6 exon arrays in RAM
data.rma <- rma(data.exon,"HuExonMixRMAMetacore",filedir=datdir,tmpdir="",background="antigenomic",
                normalize=T,option="transcript",exonlevel="metacore+affx")
# for more exon arrays better use tmpdir
tmpdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/temp"
data.rma.tmp <- rma(data.exon,"HuExonMixRMAtmpMetacore",filedir=datdir,tmpdir=tmpdir,background="antigenomic",
                    normalize=T,option="transcript",exonlevel="metacore+affx")
# probeset: metacore
data.rma.ps <- rma(data.exon,"HuExonMixRMAMetacorePS",filedir=datdir,tmpdir="",background="antigenomic",
                   normalize=T,option="probeset",exonlevel="metacore+affx")

# 2. MAS5
data.mas5 <- mas5(data.exon,"HuExonMixMAS5Metacore",filedir=datdir,tmpdir="",
                  normalize=T,sc=500,option="transcript",exonlevel="metacore+affx")
# to store all trees (including e.g. background trees) in same ROOT file, use "update=T"
data.mas5 <- mas5(data.exon,"HuExonMixExonMAS5MetacoreAll",filedir=datdir,tmpdir="",
                  normalize=T,sc=500,option="transcript",exonlevel="metacore+affx", update=T)

# 3. MAS5 detection call (yes, this is possible for exon arrays)
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls
call.mas5 <- mas5.call(data.exon,"HuExonMixCallMetacore",filedir=datdir,tmpdir="",
                       option="transcript",exonlevel="metacore+affx")

# 4. DABG detection call
# transcript: metacore
call.dabg <- dabg.call(data.exon,"HuExonMixDABGMetacore",filedir=datdir,
                       option="transcript",exonlevel="metacore+affx")
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls for transcripts
#call.dabg <- dabg.call(data.exon,"HuExonMixDABGMetacore", alpha1=???,alpha2=???)
# probeset: metacore
call.dabg.ps <- dabg.call(data.exon,"HuExonMixDABGMetacorePS",filedir=datdir,
                          option="probeset",exonlevel="metacore", alpha1=0.01,alpha2=0.015)

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)
pval.dabg <- pvalData(call.dabg)
pres.dabg <- presCall(call.dabg)


### plot results ###

# compare mas5 to rma
plot(expr.rma[,1],expr.mas5[,1])
plot(expr.rma[,1],expr.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot.dev(data.rma)
boxplot.dev(data.rma, dev="png", mar=c(4,4,1,1), w=480, h=480, outfile="BoxPlot_MixExon_rma")

boxplot(data.mas5)
boxplot.dev(data.mas5)
boxplot.dev(data.mas5, dev="png", mar=c(4,4,1,1), w=480, h=480, outfile="BoxPlot_MixExon_mas5")

# relative boxplots
mboxplot(data.rma)
mboxplot(data.rma, ylim=c(-3,3))
mboxplot(data.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot(data.rma, pch=20, ylim=c(-4,4))
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6))
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), dev="png",outfile="MvAPlot_MixExon_rma")
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), names="BreastA.mdp_LEVEL")
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), names=c("BreastA.mdp_LEVEL","BreastB.mdp_LEVEL"))
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), names="BreastA.mdp_LEVEL", dev="png",outfile="MvAPlot_MixExon_BreastA_rma")

mvaplot(data.mas5, pch=20)
mvaplot.dev(data.mas5, pch=20)
mvaplot.dev(data.mas5, pch=20, names="BreastA.tmn_LEVEL")
mvaplot.dev(data.mas5, pch=20, names=c("BreastA.tmn_LEVEL","BreastB.tmn_LEVEL"))
mvaplot.dev(data.mas5, pch=20, dev="png",outfile="MvAPlot_MixExon_mas5")
mvaplot.dev(data.mas5, pch=20, names="BreastA.tmn_LEVEL", dev="png",outfile="MvAPlot_MixExon_BreastA_mas5")

# present call plots
callplot(call.mas5)
callplot(call.mas5, beside=F, ylim=c(0,125))
callplot(call.dabg)
callplot(call.dabg, beside=F, ylim=c(0,125))
callplot(call.dabg.ps)
callplot(call.dabg.ps, beside=F, ylim=c(0,125))

# image for rma background intensities
# 1. find background treenames for data.rma
getTreeNames(rootFile(data.rma))
# 2. get background intensity for e.g. tree "BreastA.rbg"
bgrd <- export.root(rootFile(data.rma),schemeFile(data.rma),"PreprocesSet","BreastA","rbg","fBg","BreastAExonBgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
rbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.rma)), nrow=nrows(schemeSet(data.rma)))
# 4. create image
image(rbg)
image(log2(rbg))

# image for mas5 background intensities
# 1. find background treenames for data.mas5
getTreeNames(rootFile(data.mas5))
# 2. get background intensity for e.g. tree "BreastA.wbg"
bgrd <- export.root(rootFile(data.mas5),schemeFile(data.mas5),"PreprocesSet","BreastA","wbg","fBg","BreastAExonBgrdROOTOut.txt",as.dataframe=T)
rootfile <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/HuMixture/MixExon/MixMAS5Metacore_tbw.root"
bgrd <- export.root(rootfile,schemeFile(data.mas5),"PreprocesSet","BreastA","wbg","fBg","BreastAExonBgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
wbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.mas5)), nrow=nrows(schemeSet(data.mas5)))
# 4. create image
image(wbg)
image(log2(wbg))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 4: Tissues from Affymetrix Exon Array Dataset for HuGene-1_0-st-v1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.genome <- root.scheme(paste(scmdir,"Scheme_HuGene10stv1r3_na23.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.genome <- root.data(scheme.genome, paste(datdir,"HuTissuesGenome_cel.root",sep="/"))


### preprocess raw data ###
datdir <- getwd()

# 1. RMA
data.rma <- rma(data.genome,"HuGeneMixRMAMetacore",filedir=datdir,tmpdir="",
                background="antigenomic",normalize=T,exonlevel="metacore+affx")

# 2. MAS5
data.mas5 <- mas5(data.genome,"HuGeneMixMAS5Metacore",filedir=datdir,tmpdir="",
                  normalize=T,sc=500,exonlevel="metacore+affx")
# to store all trees (including e.g. background trees) in same ROOT file, use "update=T"
data.mas5 <- mas5(data.genome,"HuGeneMixExonMAS5MetacoreAll",filedir=datdir,tmpdir="",
                  normalize=T,sc=500,exonlevel="metacore+affx", update=T)

# 3. MAS5 detection call (yes, this is possible for exon arrays)
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls
call.mas5 <- mas5.call(data.genome,"HuGeneMixCallMetacore",filedir=datdir,tmpdir="",
                       exonlevel="metacore+affx")

# 4. DABG detection call
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls for transcripts
call.dabg <- dabg.call(data.genome,"HuGeneMixDABGMetacore",filedir=datdir,
                       exonlevel="metacore+affx")
#call.dabg <- dabg.call(data.genome,"HuGeneMixDABGMetacore", alpha1=???,alpha2=???)

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)
pval.dabg <- pvalData(call.dabg)
pres.dabg <- presCall(call.dabg)


### plot results ###

# compare mas5 to rma
plot(expr.rma[,1],expr.mas5[,1])
plot(expr.rma[,1],expr.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot.dev(data.rma)
boxplot.dev(data.rma, dev="png", mar=c(4,4,1,1), w=480, h=480, outfile="BoxPlot_MixHuGene_rma")

boxplot(data.mas5)
boxplot.dev(data.mas5)
boxplot.dev(data.mas5, dev="png", mar=c(4,4,1,1), w=480, h=480, outfile="BoxPlot_MixHuGene_mas5")

# relative boxplots
mboxplot(data.rma)
mboxplot(data.rma, ylim=c(-3,3))
mboxplot(data.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot(data.rma, pch=20, ylim=c(-4,4))
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6))
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), dev="png",outfile="MvAPlot_MixHuGene_rma")
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), names="Breast01.mdp_LEVEL")
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), names=c("Breast01.mdp_LEVEL","Breast02.mdp_LEVEL"))
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), names="Breast01.mdp_LEVEL", dev="png",outfile="MvAPlot_MixHuGene_Breast01_rma")

mvaplot(data.mas5, pch=20)
mvaplot.dev(data.mas5, pch=20)
mvaplot.dev(data.mas5, pch=20, names="Breast01.tmn_LEVEL")
mvaplot.dev(data.mas5, pch=20, names=c("Breast01.tmn_LEVEL","Breast02.tmn_LEVEL"))
mvaplot.dev(data.mas5, pch=20, dev="png",outfile="MvAPlot_MixHuGene_mas5")
mvaplot.dev(data.mas5, pch=20, names="Breast01.tmn_LEVEL", dev="png",outfile="MvAPlot_MixHuGene_Breast01_mas5")

# present call plots
callplot(call.mas5)
callplot(call.mas5, beside=F, ylim=c(0,125))
callplot(call.dabg)
callplot(call.dabg, beside=F, ylim=c(0,125))

# image for rma background intensities
# 1. find background treenames for data.rma
getTreeNames(rootFile(data.rma))
# 2. get background intensity for e.g. tree "Breast1.rbg"
bgrd <- export.root(rootFile(data.rma),schemeFile(data.rma),"PreprocesSet","Breast01","rbg","fBg","Breast01HuGeneBgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
rbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.rma)), nrow=nrows(schemeSet(data.rma)))
# 4. create image
image(rbg)
image(log2(rbg))

# image for mas5 background intensities
# 1. find background treenames for data.mas5
getTreeNames(rootFile(data.mas5))
# 2. get background intensity for e.g. tree "Breast01.wbg"
bgrd <- export.root(rootFile(data.mas5),schemeFile(data.mas5),"PreprocesSet","Breast01","wbg","fBg","Breast01HuGeneBgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
wbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.mas5)), nrow=nrows(schemeSet(data.mas5)))
# 4. create image
image(wbg)
image(log2(wbg))



#------------------------------------------------------------------------------#
# Demonstrations of advanced methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT scheme files and ROOT raw data files are usually already stored
#          in special system directories. When a new R session is created for the
#          first time, they must fist be loaded using "root.scheme()" and "root.data()".
#          However, this is not necessary when re-opening a saved R session later.
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration 1: compute RMA stepwise for Test3 samples
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

# directory
datdir <- getwd()

### first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(.path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(.path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

### 1.step: background - rma
data.bg.rma <- bgcorrect.rma(data.test3,"Test3RMABgrd",filedir=datdir,tmpdir="")

# attach data
data.bg.rma <- attachMask(data.bg.rma)
data.bg.rma <- attachInten(data.bg.rma)
data.bg.rma <- attachBgrd(data.bg.rma)

# plot intensities
hist(data.bg.rma)
mboxplot(data.bg.rma, ylim=c(-6,6))
pmplot(data.bg.rma)
image(data.bg.rma,col=rainbow(32))
image.dev(data.bg.rma,col=rainbow(32))
boxplot.dev(data.bg.rma)

# plot background
image.dev(data.bg.rma,bg=T,col=rainbow(32))
image.dev(data.bg.rma,bg=T,transfo=0,col=rainbow(32))

# remove data
data.bg.rma <- removeInten(data.bg.rma)
data.bg.rma <- removeBgrd(data.bg.rma)

### 2step: normalization - quantile
data.qu.rma <- normalize.quantiles(data.bg.rma,"Test3RMANorm",filedir=datdir,tmpdir="")

# plot intensities
data.qu.rma <- attachInten(data.qu.rma)

hist(data.qu.rma)
mboxplot(data.qu.rma, ylim=c(-6,6))
image.dev(data.qu.rma,col=rainbow(32))
boxplot.dev(data.qu.rma,transfo=0)

data.qu.rma <- removeInten(data.qu.rma)

### 3.step: summarization - medpol
data.mp.rma <- summarize.rma(data.qu.rma,"Test3RMAExpr",filedir=datdir,tmpdir="")

# plot expression levels
hist(data.mp.rma)
boxplot(data.mp.rma)
mboxplot(data.mp.rma)
mvaplot(data.mp.rma, pch=20, ylim=c(-4,4))
mvaplot.dev(data.mp.rma, pch=20, ylim=c(-4,4))

### alternatively save all data in same ROOT file using "update=T"
data.bg.rmall <- bgcorrect.rma(data.test3,"Test3RMAall",filedir=datdir,tmpdir="")
data.qu.rmall <- normalize.quantiles(data.bg.rmall,"Test3RMAall",filedir=datdir,tmpdir="",update=T)
data.mp.rmall <- summarize.rma(data.qu.rmall,"Test3RMAall",filedir=datdir,tmpdir="",update=T)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration 2: compute RMA for Test3 samples using function "express()"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### compute rma stepwise
expr.bg.rma <- express(data.test3,"Test3ExprsBgrd",filedir=datdir,tmpdir="",update=F,
               bgcorrect.method="rma",bgcorrect.select="none",bgcorrect.option="pmonly:epanechnikov",bgcorrect.params=c(16384))

expr.qu.rma <- express(expr.bg.rma,"Test3ExprsNorm",filedir=datdir,tmpdir="",update=F,
               normalize.method="quantile",normalize.select="pmonly",normalize.option="together:none",normalize.logbase="0",normalize.params=c(0.0))

expr.mp.rma <- express(expr.qu.rma,"Test3ExprsSum",filedir=datdir,tmpdir="",update=F,
               summarize.method="medianpolish",summarize.select="pmonly",summarize.option="transcript",summarize.logbase="log2",summarize.params=c(10, 0.01, 1.0))

# rma
data.rma <- rma(data.test3,"tmp_Test3RMA",tmpdir="",background="pmonly",normalize=T)

# compare results
expr <- exprs(data.rma)
expr.mp <- exprs(expr.mp.rma)
# plot differences
plot((expr.mp[,1] - expr[,1])/expr[,1], ylim=c(-0.0001,0.0001))

### compute rma with a single call to express()
expr.rma <- express(data.test3,"Test3Exprs",filedir=datdir,tmpdir="",update=F,
            bgcorrect.method="rma",bgcorrect.select="none",bgcorrect.option="pmonly:epanechnikov",bgcorrect.params=c(16384),
            normalize.method="quantile",normalize.select="pmonly",normalize.option="together:none",normalize.logbase="0",normalize.params=c(0.0),
            summarize.method="medianpolish",summarize.select="pmonly",summarize.option="transcript",summarize.logbase="log2",summarize.params=c(10, 0.01, 1.0))

# compare results
expr <- exprs(data.rma)
expr.mp <- exprs(expr.rma)
# plot differences
plot((expr.mp[,1] - expr[,1])/expr[,1], ylim=c(-0.0001,0.0001))











