#------------------------------------------------------------------------------#
# Script: step-by-step functions to create ROOT scheme files for package "xps"
#
# Note: please feel free to copy-paste the examples of interest and adapt the
#       examples to your own needs
#
# Copyright (c) 2010-2010 Christian Stratowa, Vienna, Austria.
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
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: do not separate name of ROOT files with dots, use underscores,
#          e.g. do not use "test3.na30" but "test3_na30"
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


#------------------------------------------------------------------------------#
# Nov 2009: Affymetrix annotation "na30"
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na30 <- import.expr.scheme("test3_na30", filedir = scmdir,
                     schemefile = file.path(libdir, "Test3.CDF"), 
                     probefile  = file.path(libdir, "Test3_probe.tab"), 
                     annotfile  = file.path(anndir, "Version09Nov", "Test3.na30.annot.csv"))

# Hu6800:
scheme.hu6800.na30 <- import.expr.scheme("hu6800_na30", filedir = scmdir,
                      schemefile = file.path(libdir, "Hu6800.CDF"), 
                      probefile  = file.path(libdir, "HuGeneFL_probe.tab"), 
                      annotfile  = file.path(anndir, "Version09Nov", "Hu6800.na30.annot.csv"))

# HG_U95A:
scheme.hgu95a.na30 <- import.expr.scheme("hgu95a_na30", filedir = scmdir,
                      schemefile = file.path(libdir, "HG_U95A.CDF"), 
                      probefile  = file.path(libdir, "HG-U95A_probe.tab"), 
                      annotfile  = file.path(anndir, "Version09Nov", "HG_U95A.na30.annot.csv"))

# HG_U95Av2:
scheme.hgu95av2.na30 <- import.expr.scheme("hgu95av2_na30", filedir = scmdir,
                        schemefile = file.path(libdir, "HG_U95Av2.CDF"), 
                        probefile  = file.path(libdir, "HG-U95Av2_probe.tab"), 
                        annotfile  = file.path(anndir, "Version09Nov", "HG_U95Av2.na30.annot.csv"))

# HG-U133A:
scheme.hgu133a.na30 <- import.expr.scheme("hgu133a_na30", filedir = scmdir,
                       schemefile = file.path(libdir, "HG-U133A.CDF"), 
                       probefile  = file.path(libdir, "HG-U133A_probe.tab"), 
                       annotfile  = file.path(anndir, "Version09Nov", "HG-U133A.na30.annot.csv"))

# HG-U133B:
scheme.hgu133b.na30 <- import.expr.scheme("hgu133b_na30", filedir = scmdir,
                       schemefile = file.path(libdir, "HG-U133B.CDF"), 
                       probefile  = file.path(libdir, "HG-U133B_probe.tab"), 
                       annotfile  = file.path(anndir, "Version09Nov", "HG-U133B.na30.annot.csv"))

# HG-U133_Plus_2:
scheme.hgu133plus2.na30 <- import.expr.scheme("hgu133plus2_na30", filedir = scmdir,
                           schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"), 
                           probefile  = file.path(libdir, "HG-U133-PLUS_probe.tab"), 
                           annotfile  = file.path(anndir, "Version09Nov", "HG-U133_Plus_2.na30.annot.csv"))

# MOE430A:
scheme.moe430a.na30 <- import.expr.scheme("moe430a_na30",filedir=scmdir,
                       schemefile = file.path(libdir,"MOE430A.CDF"),
                       probefile  = file.path(libdir,"MOE430A.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/MOE430A.na30.annot.csv"))

# MOE430B:
scheme.moe430b.na30 <- import.expr.scheme("moe430b_na30",filedir=scmdir,
                       schemefile = file.path(libdir,"MOE430B.CDF"),
                       probefile  = file.path(libdir,"MOE430B.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/MOE430B.na30.annot.csv"))

# Mouse430_2:
scheme.mouse4302.na30 <- import.expr.scheme("mouse4302_na30",filedir=scmdir,
                         schemefile = file.path(libdir,"Mouse430_2.cdf"),
                         probefile  = file.path(libdir,"Mouse430_2.probe.tab"),
                         annotfile  = file.path(anndir,"Version09Nov/Mouse430_2.na30.annot.csv"))

# RG_U34A:
scheme.rgu34a.na30 <- import.expr.scheme("rgu34a_na30",filedir=scmdir,
                       schemefile = file.path(libdir,"RG_U34A.cdf"),
                       probefile  = file.path(libdir,"RG_U34A.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RG_U34A.na30.annot.csv"))

# RG_U34B:
scheme.rgu34b.na30 <- import.expr.scheme("rgu34b_na30",filedir=scmdir,
                       schemefile = file.path(libdir,"RG_U34B.cdf"),
                       probefile  = file.path(libdir,"RG_U34B.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RG_U34B.na30.annot.csv"))

# RG_U34C:
scheme.rgu34c.na30 <- import.expr.scheme("rgu34c_na30",filedir=scmdir,
                       schemefile = file.path(libdir,"RG_U34C.cdf"),
                       probefile  = file.path(libdir,"RG_U34C.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RG_U34C.na30.annot.csv"))

# RAE230A:
scheme.rae230a.na30 <- import.expr.scheme("rae230a_na30",filedir=scmdir,
                       schemefile = file.path(libdir,"RAE230A.CDF"),
                       probefile  = file.path(libdir,"RAE230A.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RAE230A.na30.annot.csv"))

# RAE230B:
scheme.rae230b.na30 <- import.expr.scheme("rae230b_na30",filedir=scmdir,
                       schemefile = file.path(libdir,"RAE230B.CDF"),
                       probefile  = file.path(libdir,"RAE230B.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RAE230B.na30.annot.csv"))

# Rat230_2:
scheme.rat2302.na30 <- import.expr.scheme("rat2302_na30",filedir=scmdir,
                       schemefile = file.path(libdir,"Rat230_2.cdf"),
                       probefile  = file.path(libdir,"Rat230_2.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/Rat230_2.na30.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HT_HG-U133A
scheme.hthgu133a.na30 <- import.expr.scheme("hthgu133a_na30", filedir = scmdir,
                         schemefile = file.path(libdir, "HT_HG-U133A.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version09Nov", "HT_HG-U133A.na30.annot.csv"))

# HT_HG-U133_Plus_PM
scheme.hthgu133pluspm.na30 <- import.expr.scheme("hthgu133pluspm_na30", filedir = scmdir,
                              schemefile = file.path(libdir, "HT_HG-U133_Plus_PM.CDF"), 
                              probefile  = file.path(libdir, "HT_HG-U133_Plus_PM.probe.tab"), 
                              annotfile  = file.path(anndir, "Version09Nov", "HT_HG-U133_Plus_PM.na30.annot.csv"))

# HT_MG-430A
scheme.htmg430a.na30 <- import.expr.scheme("htmg430a_na30", filedir = scmdir,
                         schemefile = file.path(libdir, "HT_MG-430A.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version09Nov", "HT_MG-430A.na30.annot.csv"))

# HT_MG-430_PM
scheme.htmg430pm.na30 <- import.expr.scheme("htmg430pm_na30", filedir = scmdir,
                         schemefile = file.path(libdir, "HT_MG-430_PM.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430_PM.probe.tab"), 
                         annotfile  = file.path(anndir, "Version09Nov", "HT_MG-430_PM.na30.annot.csv"))

# HT_Rat230_PM
scheme.htrat230pm.na30 <- import.expr.scheme("htrat230pm_na30", filedir = scmdir,
                          schemefile = file.path(libdir, "HT_Rat230_PM.cdf"), 
                          probefile  = file.path(libdir, "HT_Rat230_PM.probe.tab"), 
                          annotfile  = file.path(anndir, "Version09Nov", "HT_Rat230_PM.na30.annot.csv"))

# HG-U219
scheme.hgu219.na30 <- import.expr.scheme("hgu219_na30", filedir = scmdir,
                       schemefile = file.path(libdir, "HG-U219.cdf"), 
                       probefile  = file.path(libdir, "HG-U219.probe.tab"), 
                       annotfile  = file.path(anndir, "Version09Nov", "HG-U219.na30.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1.na30 <- import.exon.scheme("hugene10stv1_na30", filedir = scmdir,
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.1.hg19.probeset.csv"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.hg19.transcript.csv"))

# MoGene-1_0-st-v1.r4: used as exon array
scheme.mogene10stv1.na30 <- import.exon.scheme("mogene10stv1_na30", filedir = scmdir,
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.probeset.csv"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.transcript.csv"))

# RaGene-1_0-st-v1.r4: used as exon array
scheme.ragene10stv1.na30 <- import.exon.scheme("ragene10stv1_na30", filedir = scmdir,
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.probeset.csv"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.transcript.csv"))

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2.na30 <- import.exon.scheme("huex10stv2_na30", filedir = scmdir,
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                          file.path(anndir, "Version09Nov", "HuEx-1_0-st-v2.na30.hg19.probeset.csv"),
                          file.path(anndir, "Version09Nov", "HuEx-1_0-st-v2.na30.hg19.transcript.csv"))

# MoEx-1_0-st-v1.r2:
scheme.moex10stv1.na30 <- import.exon.scheme("moex10stv1_na30",filedir = scmdir,
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version09Nov", "MoEx-1_0-st-v1.na30.mm9.probeset.csv"),
                          file.path(anndir, "Version09Nov", "MoEx-1_0-st-v1.na30.mm9.transcript.csv"))

# RaEx-1_0-st-v1.r2:
scheme.raex10stv1.na30 <- import.exon.scheme("raex10stv1_na30",filedir=scmdir,
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version09Nov", "RaEx-1_0-st-v1.na30.rn4.probeset.csv"),
                          file.path(anndir, "Version09Nov", "RaEx-1_0-st-v1.na30.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_1-st-v1:
scheme.hugene11stv1.na30 <- import.exon.scheme("hugene11stv1_na30", filedir = scmdir,
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.1.hg19.probeset.csv"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.hg19.transcript.csv"))

# MoGene-1_1-st-v1.r4
scheme.mogene11stv1.na30 <- import.exon.scheme("mogene11stv1_na30", filedir = scmdir,
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.probeset.csv"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.transcript.csv"))

# RaGene-1_1-st-v1.r4
scheme.ragene11stv1.na30 <- import.exon.scheme("ragene11stv1_na30", filedir = scmdir,
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.probeset.csv"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for other Affymetrix arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# E_coli_2:
scheme.ecoli2.na30 <- import.expr.scheme("ecoli2_na30", filedir = scmdir,
                      schemefile = file.path(libdir, "E_coli_2.cdf"),
                      probefile  = file.path(libdir, "E_coli_2.probe.tab"),
                      annotfile  = file.path(anndir, "Version09Nov", "E_coli_2.na30.annot.csv"))

# Ecoli_ASv2:
scheme.ecoliasv2.na30 <- import.expr.scheme("ecoliasv2_na30", filedir = scmdir,
                         schemefile = file.path(libdir, "Ecoli_ASv2.CDF"),
                         probefile  = file.path(libdir, "Ecoli_ASv2.probe.tab"),
                         annotfile  = file.path(anndir, "Version09Nov", "Ecoli_ASv2.na30.annot.csv"))

# Pae_G1a:
scheme.paeg1a.na30 <- import.expr.scheme("paeg1a_na30", filedir = scmdir,
                      schemefile = file.path(libdir, "Pae_G1a.CDF"),
                      probefile  = file.path(libdir, "Pae_G1a.probe.tab"), 
                      annotfile  = file.path(anndir, "Version09Nov", "Pae_G1a.na30.annot.csv"))

# Yeast_2:
scheme.yeast2.na30 <- import.expr.scheme("yeast2_na30", filedir = scmdir,
                      schemefile = file.path(libdir, "Yeast_2.cdf"),
                      probefile  = file.path(libdir, "Yeast_2.probe.tab"),
                      annotfile  = file.path(anndir, "Version09Nov", "Yeast_2.na30.annot.csv"))

# miRNA-1_0:
# note: you need to rename "miRNA-1_0.probe_list.20081203.txt" to "miRNA-1_0.probe.tab"
scheme.mirna10 <- import.expr.scheme("mirna10", filedir = scmdir,
                  schemefile = file.path(libdir, "miRNA-1_0.CDF"),
                  probefile  = file.path(libdir, "miRNA-1_0.probe.tab"))




#------------------------------------------------------------------------------#
