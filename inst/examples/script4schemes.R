#------------------------------------------------------------------------------#
# Script: step-by-step functions to create ROOT scheme files for package "xps"
#
# Note: please feel free to copy-paste the examples of interest and adapt the
#       examples to your own needs
#
# Copyright (c) 2010-2011 Christian Stratowa, Vienna, Austria.
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
#    Note: do not separate name of ROOT files with dots, use underscores instead.
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
# Jul 2011: Affymetrix annotation "na32"
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na32 <- import.expr.scheme("test3", filedir = file.path(scmdir, "na32"),
                     schemefile = file.path(libdir, "Test3.CDF"), 
                     probefile  = file.path(libdir, "Test3_probe.tab"), 
                     annotfile  = file.path(anndir, "Version11Jul", "Test3.na32.annot.csv"))

# Hu6800: 
# annotation file for na32 no longer provided by Affymetrix (see na30 annotations)

# HG_U95A:
scheme.hgu95a.na32 <- import.expr.scheme("hgu95a", filedir = file.path(scmdir, "na32"),
                      schemefile = file.path(libdir, "HG_U95A.CDF"), 
                      probefile  = file.path(libdir, "HG-U95A_probe.tab"), 
                      annotfile  = file.path(anndir, "Version11Jul", "HG_U95A.na32.annot.csv"))

# HG_U95Av2:
scheme.hgu95av2.na32 <- import.expr.scheme("hgu95av2", filedir = file.path(scmdir, "na32"),
                        schemefile = file.path(libdir, "HG_U95Av2.CDF"), 
                        probefile  = file.path(libdir, "HG-U95Av2_probe.tab"), 
                        annotfile  = file.path(anndir, "Version11Jul", "HG_U95Av2.na32.annot.csv"))

# HG-U133A:
scheme.hgu133a.na32 <- import.expr.scheme("hgu133a", filedir = file.path(scmdir, "na32"),
                       schemefile = file.path(libdir, "HG-U133A.CDF"), 
                       probefile  = file.path(libdir, "HG-U133A_probe.tab"), 
                       annotfile  = file.path(anndir, "Version11Jul", "HG-U133A.na32.annot.csv"))

# HG-U133B:
scheme.hgu133b.na32 <- import.expr.scheme("hgu133b", filedir = file.path(scmdir, "na32"),
                       schemefile = file.path(libdir, "HG-U133B.CDF"), 
                       probefile  = file.path(libdir, "HG-U133B_probe.tab"), 
                       annotfile  = file.path(anndir, "Version11Jul", "HG-U133B.na32.annot.csv"))

# HG-U133_Plus_2:
scheme.hgu133plus2.na32 <- import.expr.scheme("hgu133plus2", filedir = file.path(scmdir, "na32"),
                           schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"), 
                           probefile  = file.path(libdir, "HG-U133-PLUS_probe.tab"), 
                           annotfile  = file.path(anndir, "Version11Jul", "HG-U133_Plus_2.na32.annot.csv"))

# MOE430A:
scheme.moe430a.na32 <- import.expr.scheme("moe430a",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"MOE430A.CDF"),
                       probefile  = file.path(libdir,"MOE430A.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/MOE430A.na32.annot.csv"))

# MOE430B:
scheme.moe430b.na32 <- import.expr.scheme("moe430b",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"MOE430B.CDF"),
                       probefile  = file.path(libdir,"MOE430B.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/MOE430B.na32.annot.csv"))

# Mouse430_2:
scheme.mouse4302.na32 <- import.expr.scheme("mouse4302",filedir=file.path(scmdir, "na32"),
                         schemefile = file.path(libdir,"Mouse430_2.cdf"),
                         probefile  = file.path(libdir,"Mouse430_2.probe.tab"),
                         annotfile  = file.path(anndir,"Version11Jul/Mouse430_2.na32.annot.csv"))

# RG_U34A:
scheme.rgu34a.na32 <- import.expr.scheme("rgu34a",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"RG_U34A.cdf"),
                       probefile  = file.path(libdir,"RG_U34A.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/RG_U34A.na32.annot.csv"))

# RG_U34B:
scheme.rgu34b.na32 <- import.expr.scheme("rgu34b",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"RG_U34B.cdf"),
                       probefile  = file.path(libdir,"RG_U34B.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/RG_U34B.na32.annot.csv"))

# RAE230A:
scheme.rae230a.na32 <- import.expr.scheme("rae230a",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"RAE230A.CDF"),
                       probefile  = file.path(libdir,"RAE230A.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/RAE230A.na32.annot.csv"))

# RAE230B:
scheme.rae230b.na32 <- import.expr.scheme("rae230b",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"RAE230B.CDF"),
                       probefile  = file.path(libdir,"RAE230B.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/RAE230B.na32.annot.csv"))

# Rat230_2:
scheme.rat2302.na32 <- import.expr.scheme("rat2302",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"Rat230_2.cdf"),
                       probefile  = file.path(libdir,"Rat230_2.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/Rat230_2.na32.annot.csv"))

# Bovine:
scheme.bovine.na32 <- import.expr.scheme("bovine",filedir=file.path(scmdir, "na32"),
                      schemefile = file.path(libdir,"Bovine.cdf"),
                      probefile  = file.path(libdir,"Bovine.probe.tab"),
                      annotfile  = file.path(anndir,"Version11Jul/Bovine.na32.annot.csv"))

# Porcine:
scheme.porcine.na32 <- import.expr.scheme("porcine",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"Porcine.cdf"),
                       probefile  = file.path(libdir,"Porcine.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/Porcine.na32.annot.csv"))

# Rhesus:
scheme.rhesus.na32 <- import.expr.scheme("rhesus",filedir=file.path(scmdir, "na32"),
                      schemefile = file.path(libdir,"Rhesus.cdf"),
                      probefile  = file.path(libdir,"Rhesus.probe.tab"),
                      annotfile  = file.path(anndir,"Version11Jul/Rhesus.na32.annot.csv"))

# Zebrafish:
scheme.zebrafish.na32 <- import.expr.scheme("zebrafish",filedir=file.path(scmdir, "na32"),
                         schemefile = file.path(libdir,"Zebrafish.cdf"),
                         probefile  = file.path(libdir,"Zebrafish.probe.tab"),
                         annotfile  = file.path(anndir,"Version11Jul/Zebrafish.na32.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HT_HG-U133A
scheme.hthgu133a.na32 <- import.expr.scheme("hthgu133a", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "HT_HG-U133A.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "HT_HG-U133A.na32.annot.csv"))

# HT_HG-U133B
scheme.hthgu133b.na32 <- import.expr.scheme("hthgu133b", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "HT_HG-U133B.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133B.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "HT_HG-U133B.na32.annot.csv"))

# HT_HG-U133_Plus_PM
scheme.hthgu133pluspm.na32 <- import.expr.scheme("hthgu133pluspm", filedir = file.path(scmdir, "na32"),
                              schemefile = file.path(libdir, "HT_HG-U133_Plus_PM.CDF"), 
                              probefile  = file.path(libdir, "HT_HG-U133_Plus_PM.probe.tab"), 
                              annotfile  = file.path(anndir, "Version11Jul", "HT_HG-U133_Plus_PM.na32.annot.csv"))

# HT_MG-430A
scheme.htmg430a.na32 <- import.expr.scheme("htmg430a", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "HT_MG-430A.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "HT_MG-430A.na32.annot.csv"))

# HT_MG-430B
scheme.htmg430b.na32 <- import.expr.scheme("htmg430b", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "HT_MG-430B.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430B.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "HT_MG-430B.na32.annot.csv"))

# HT_MG-430_PM
scheme.htmg430pm.na32 <- import.expr.scheme("htmg430pm", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "HT_MG-430_PM.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430_PM.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "HT_MG-430_PM.na32.annot.csv"))

# HT_Rat230_PM
scheme.htrat230pm.na32 <- import.expr.scheme("htrat230pm", filedir = file.path(scmdir, "na32"),
                          schemefile = file.path(libdir, "HT_Rat230_PM.cdf"), 
                          probefile  = file.path(libdir, "HT_Rat230_PM.probe.tab"), 
                          annotfile  = file.path(anndir, "Version11Jul", "HT_Rat230_PM.na32.annot.csv"))

# HG-U219
scheme.hgu219.na32 <- import.expr.scheme("hgu219", filedir = file.path(scmdir, "na32"),
                       schemefile = file.path(libdir, "HG-U219.cdf"), 
                       probefile  = file.path(libdir, "HG-U219.probe.tab"), 
                       annotfile  = file.path(anndir, "Version11Jul", "HG-U219.na32.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1.na32 <- import.exon.scheme("hugene10stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "HuGene-1_0-st-v1.na32.hg19.probeset.csv"),
                            file.path(anndir, "Version11Jul", "HuGene-1_0-st-v1.na32.hg19.transcript.csv"))

   <257385> probeset tree entries read...Finished
Warning: number of  probes <33> for probeset <8180192> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180193> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180194> is not equal to annotated probe_count <0>
Warning: number of  probes <37> for probeset <8180195> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180196> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180197> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <8180198> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180199> is not equal to annotated probe_count <0>
Warning: number of  probes <9> for probeset <8180200> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180201> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180202> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180203> is not equal to annotated probe_count <0>
Warning: number of  probes <37> for probeset <8180204> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180205> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180206> is not equal to annotated probe_count <0>
Warning: number of  probes <17> for probeset <8180207> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180208> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180209> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180210> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180211> is not equal to annotated probe_count <0>
Warning: number of  probes <44> for probeset <8180212> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <8180213> is not equal to annotated probe_count <0>
Warning: number of  probes <17> for probeset <8180214> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <8180215> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180216> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180217> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180218> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180219> is not equal to annotated probe_count <0>
Warning: number of  probes <18> for probeset <8180220> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180221> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180222> is not equal to annotated probe_count <0>
Warning: number of  probes <18> for probeset <8180223> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180224> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <8180225> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180226> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180227> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180228> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180229> is not equal to annotated probe_count <0>
Warning: number of  probes <17> for probeset <8180230> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180231> is not equal to annotated probe_count <0>
Warning: number of  probes <18> for probeset <8180232> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180233> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180234> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180235> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180236> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180237> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180238> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180239> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180240> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180241> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180242> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180243> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180244> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180245> is not equal to annotated probe_count <0>
Warning: number of  probes <12> for probeset <8180246> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180247> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180248> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180249> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180250> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <8180251> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180252> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180253> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180254> is not equal to annotated probe_count <0>
Warning: number of  probes <6> for probeset <8180255> is not equal to annotated probe_count <0>
Warning: number of  probes <42> for probeset <8180256> is not equal to annotated probe_count <0>
Warning: number of  probes <43> for probeset <8180257> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180258> is not equal to annotated probe_count <0>
Warning: number of  probes <8> for probeset <8180259> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180260> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180261> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180262> is not equal to annotated probe_count <0>
Warning: number of  probes <51> for probeset <8180263> is not equal to annotated probe_count <0>
Warning: number of  probes <6> for probeset <8180264> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180265> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180266> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180267> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180268> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180269> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180270> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180271> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180272> is not equal to annotated probe_count <0>
Warning: number of  probes <15> for probeset <8180273> is not equal to annotated probe_count <0>
Warning: number of  probes <9> for probeset <8180274> is not equal to annotated probe_count <0>
Warning: number of  probes <10> for probeset <8180275> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180276> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180277> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180278> is not equal to annotated probe_count <0>
Warning: number of  probes <10> for probeset <8180279> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180280> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180281> is not equal to annotated probe_count <0>
Warning: number of  probes <8> for probeset <8180282> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180283> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180284> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180285> is not equal to annotated probe_count <0>
Warning: number of  probes <18> for probeset <8180286> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180287> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180288> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180289> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <8180290> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <8180291> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <8180292> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180293> is not equal to annotated probe_count <0>
Warning: number of  probes <7> for probeset <8180294> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180295> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <8180296> is not equal to annotated probe_count <0>
Warning: number of  probes <109> for probeset <8180297> is not equal to annotated probe_count <0>
Warning: number of  probes <37> for probeset <8180298> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180299> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180300> is not equal to annotated probe_count <0>
Warning: number of  probes <48> for probeset <8180301> is not equal to annotated probe_count <0>
Warning: number of  probes <47> for probeset <8180302> is not equal to annotated probe_count <0>
Warning: number of  probes <15> for probeset <8180303> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180304> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180305> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180306> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180307> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180308> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180309> is not equal to annotated probe_count <0>
Warning: number of  probes <49> for probeset <8180310> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180311> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180312> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180313> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180314> is not equal to annotated probe_count <0>
Warning: number of  probes <17> for probeset <8180315> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180316> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180317> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180318> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180319> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180320> is not equal to annotated probe_count <0>
Warning: number of  probes <6> for probeset <8180321> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180322> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180323> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180324> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180325> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180326> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180327> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180328> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180329> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180330> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180331> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180332> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180333> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180334> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180335> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <8180336> is not equal to annotated probe_count <0>
Warning: number of  probes <40> for probeset <8180337> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180338> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180339> is not equal to annotated probe_count <0>
Warning: number of  probes <242> for probeset <8180340> is not equal to annotated probe_count <0>
Warning: number of  probes <42> for probeset <8180341> is not equal to annotated probe_count <0>
Warning: number of  probes <40> for probeset <8180342> is not equal to annotated probe_count <0>
Warning: number of  probes <42> for probeset <8180343> is not equal to annotated probe_count <0>
Warning: number of  probes <39> for probeset <8180344> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180345> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180346> is not equal to annotated probe_count <0>
Warning: number of  probes <80> for probeset <8180347> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180348> is not equal to annotated probe_count <0>
Warning: number of  probes <37> for probeset <8180349> is not equal to annotated probe_count <0>
Warning: number of  probes <8> for probeset <8180350> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180351> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180352> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180353> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180354> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180355> is not equal to annotated probe_count <0>
Warning: number of  probes <49> for probeset <8180356> is not equal to annotated probe_count <0>
Warning: number of  probes <47> for probeset <8180357> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180358> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180359> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180360> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180361> is not equal to annotated probe_count <0>
Warning: number of  probes <35> for probeset <8180362> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180363> is not equal to annotated probe_count <0>
Warning: number of  probes <65> for probeset <8180364> is not equal to annotated probe_count <0>
Warning: number of  probes <65> for probeset <8180365> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180366> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180367> is not equal to annotated probe_count <0>
Warning: number of  probes <48> for probeset <8180368> is not equal to annotated probe_count <0>
Warning: number of  probes <42> for probeset <8180369> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <8180370> is not equal to annotated probe_count <0>
Warning: number of  probes <40> for probeset <8180371> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180372> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <8180373> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180374> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180375> is not equal to annotated probe_count <0>
Warning: number of  probes <34> for probeset <8180376> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180377> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180378> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180379> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180380> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180381> is not equal to annotated probe_count <0>
Warning: number of  probes <36> for probeset <8180382> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180383> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180384> is not equal to annotated probe_count <0>
Warning: number of  probes <7> for probeset <8180385> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180386> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180387> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180388> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180389> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180390> is not equal to annotated probe_count <0>
Warning: number of  probes <56> for probeset <8180391> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180392> is not equal to annotated probe_count <0>
Warning: number of  probes <37> for probeset <8180393> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <8180394> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180395> is not equal to annotated probe_count <0>
Warning: number of  probes <35> for probeset <8180396> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180397> is not equal to annotated probe_count <0>
Warning: number of  probes <48> for probeset <8180398> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180399> is not equal to annotated probe_count <0>
Warning: number of  probes <6> for probeset <8180400> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180401> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180402> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180403> is not equal to annotated probe_count <0>
Warning: number of  probes <34> for probeset <8180404> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180405> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180406> is not equal to annotated probe_count <0>
Warning: number of  probes <111> for probeset <8180407> is not equal to annotated probe_count <0>
Warning: number of  probes <121> for probeset <8180408> is not equal to annotated probe_count <0>
Warning: number of  probes <98> for probeset <8180409> is not equal to annotated probe_count <0>
Warning: number of  probes <101> for probeset <8180410> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180411> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180412> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180413> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180414> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180415> is not equal to annotated probe_count <0>
Warning: number of  probes <39> for probeset <8180416> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180417> is not equal to annotated probe_count <0>
Warning: number of  probes <6> for probeset <8180418> is not equal to annotated probe_count <0>
Note: Number of exons imported <210414> is not equal to number of annotated exons <210415>.


# MoGene-1_0-st-v1.r4: used as exon array
scheme.mogene10stv1.na32 <- import.exon.scheme("mogene10stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "MoGene-1_0-st-v1.na32.mm9.probeset.csv"),
                            file.path(anndir, "Version11Jul", "MoGene-1_0-st-v1.na32.mm9.transcript.csv"))

   <241531> probeset tree entries read...Finished
Warning: number of  probes <34> for probeset <10608634> is not equal to annotated probe_count <0>
Warning: number of  probes <79> for probeset <10608635> is not equal to annotated probe_count <0>
Warning: number of  probes <57> for probeset <10608636> is not equal to annotated probe_count <0>
Warning: number of  probes <12> for probeset <10608637> is not equal to annotated probe_count <0>
Warning: number of  probes <15> for probeset <10608638> is not equal to annotated probe_count <0>
Warning: number of  probes <49> for probeset <10608639> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608640> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <10608641> is not equal to annotated probe_count <0>
Warning: number of  probes <61> for probeset <10608642> is not equal to annotated probe_count <0>
Warning: number of  probes <75> for probeset <10608643> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <10608644> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <10608645> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <10608646> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <10608647> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <10608648> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <10608649> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <10608650> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <10608651> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <10608652> is not equal to annotated probe_count <0>
Warning: number of  probes <100> for probeset <10608653> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <10608654> is not equal to annotated probe_count <0>
Warning: number of  probes <51> for probeset <10608655> is not equal to annotated probe_count <0>
Warning: number of  probes <49> for probeset <10608656> is not equal to annotated probe_count <0>
Warning: number of  probes <36> for probeset <10608657> is not equal to annotated probe_count <0>
Warning: number of  probes <17> for probeset <10608658> is not equal to annotated probe_count <0>
Warning: number of  probes <15> for probeset <10608659> is not equal to annotated probe_count <0>
Warning: number of  probes <40> for probeset <10608660> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <10608661> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <10608662> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <10608663> is not equal to annotated probe_count <0>
Warning: number of  probes <35> for probeset <10608664> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <10608665> is not equal to annotated probe_count <0>
Warning: number of  probes <9> for probeset <10608666> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <10608667> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <10608668> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608669> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <10608670> is not equal to annotated probe_count <0>
Warning: number of  probes <35> for probeset <10608671> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <10608672> is not equal to annotated probe_count <0>
Warning: number of  probes <85> for probeset <10608674> is not equal to annotated probe_count <0>
Warning: number of  probes <15> for probeset <10608675> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <10608676> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <10608677> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <10608678> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <10608679> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <10608680> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <10608681> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608682> is not equal to annotated probe_count <0>
Warning: number of  probes <34> for probeset <10608683> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608684> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <10608685> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <10608686> is not equal to annotated probe_count <0>
Warning: number of  probes <95> for probeset <10608687> is not equal to annotated probe_count <0>
Warning: number of  probes <105> for probeset <10608688> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <10608689> is not equal to annotated probe_count <0>
Warning: number of  probes <35> for probeset <10608690> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <10608691> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <10608692> is not equal to annotated probe_count <0>
Warning: number of  probes <43> for probeset <10608695> is not equal to annotated probe_count <0>
Warning: number of  probes <88> for probeset <10608696> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608697> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <10608698> is not equal to annotated probe_count <0>
Warning: number of  probes <74> for probeset <10608699> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608700> is not equal to annotated probe_count <0>
Warning: number of  probes <74> for probeset <10608701> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608702> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608703> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608704> is not equal to annotated probe_count <0>
Warning: number of  probes <133> for probeset <10608705> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608706> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608707> is not equal to annotated probe_count <0>
Warning: number of  probes <77> for probeset <10608708> is not equal to annotated probe_count <0>
Warning: number of  probes <56> for probeset <10608709> is not equal to annotated probe_count <0>
Warning: number of  probes <100> for probeset <10608710> is not equal to annotated probe_count <0>
Warning: number of  probes <94> for probeset <10608711> is not equal to annotated probe_count <0>
Warning: number of  probes <87> for probeset <10608712> is not equal to annotated probe_count <0>
Warning: number of  probes <49> for probeset <10608713> is not equal to annotated probe_count <0>
Warning: number of  probes <86> for probeset <10608714> is not equal to annotated probe_count <0>
Warning: number of  probes <70> for probeset <10608715> is not equal to annotated probe_count <0>
Warning: number of  probes <98> for probeset <10608716> is not equal to annotated probe_count <0>
Warning: number of  probes <80> for probeset <10608717> is not equal to annotated probe_count <0>
Warning: number of  probes <165> for probeset <10608718> is not equal to annotated probe_count <0>
Warning: number of  probes <97> for probeset <10608719> is not equal to annotated probe_count <0>
Warning: number of  probes <60> for probeset <10608720> is not equal to annotated probe_count <0>
Warning: number of  probes <65> for probeset <10608721> is not equal to annotated probe_count <0>
Warning: number of  probes <41> for probeset <10608722> is not equal to annotated probe_count <0>
Warning: number of  probes <83> for probeset <10608723> is not equal to annotated probe_count <0>
Warning: number of  probes <50> for probeset <10608724> is not equal to annotated probe_count <0>
Note: Number of exons imported <207392> is not equal to number of annotated exons <207393>.


# RaGene-1_0-st-v1.r4: used as exon array
scheme.ragene10stv1.na32 <- import.exon.scheme("ragene10stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "RaGene-1_0-st-v1.na32.rn4.probeset.csv"),
                            file.path(anndir, "Version11Jul", "RaGene-1_0-st-v1.na32.rn4.transcript.csv"))

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2.na32 <- import.exon.scheme("huex10stv2", filedir = file.path(scmdir, "na32"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                          file.path(anndir, "Version11Jul", "HuEx-1_0-st-v2.na32.hg19.probeset.csv"),
                          file.path(anndir, "Version11Jul", "HuEx-1_0-st-v2.na32.hg19.transcript.csv"))

# MoEx-1_0-st-v1.r2:
scheme.moex10stv1.na32 <- import.exon.scheme("moex10stv1",filedir = file.path(scmdir, "na32"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version11Jul", "MoEx-1_0-st-v1.na32.mm9.probeset.csv"),
                          file.path(anndir, "Version11Jul", "MoEx-1_0-st-v1.na32.mm9.transcript.csv"))

# RaEx-1_0-st-v1.r2:
scheme.raex10stv1.na32 <- import.exon.scheme("raex10stv1",filedir=file.path(scmdir, "na32"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version11Jul", "RaEx-1_0-st-v1.na32.rn4.probeset.csv"),
                          file.path(anndir, "Version11Jul", "RaEx-1_0-st-v1.na32.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_1-st-v1:
scheme.hugene11stv1.na32 <- import.exon.scheme("hugene11stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "HuGene-1_1-st-v1.na32.hg19.probeset.csv"),
                            file.path(anndir, "Version11Jul", "HuGene-1_1-st-v1.na32.hg19.transcript.csv"))

   <257385> probeset tree entries read...Finished
Warning: number of  probes <33> for probeset <8180192> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180193> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180194> is not equal to annotated probe_count <0>
Warning: number of  probes <37> for probeset <8180195> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180196> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180197> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <8180198> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180199> is not equal to annotated probe_count <0>
Warning: number of  probes <9> for probeset <8180200> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180201> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180202> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180203> is not equal to annotated probe_count <0>
Warning: number of  probes <37> for probeset <8180204> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180205> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180206> is not equal to annotated probe_count <0>
Warning: number of  probes <17> for probeset <8180207> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180208> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180209> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180210> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180211> is not equal to annotated probe_count <0>
Warning: number of  probes <44> for probeset <8180212> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <8180213> is not equal to annotated probe_count <0>
Warning: number of  probes <17> for probeset <8180214> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <8180215> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180216> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180217> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180218> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180219> is not equal to annotated probe_count <0>
Warning: number of  probes <18> for probeset <8180220> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180221> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180222> is not equal to annotated probe_count <0>
Warning: number of  probes <18> for probeset <8180223> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180224> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <8180225> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180226> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180227> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180228> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180229> is not equal to annotated probe_count <0>
Warning: number of  probes <17> for probeset <8180230> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180231> is not equal to annotated probe_count <0>
Warning: number of  probes <18> for probeset <8180232> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180233> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180234> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180235> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180236> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180237> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180238> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180239> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180240> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180241> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180242> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180243> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180244> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180245> is not equal to annotated probe_count <0>
Warning: number of  probes <12> for probeset <8180246> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180247> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180248> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180249> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180250> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <8180251> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180252> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180253> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180254> is not equal to annotated probe_count <0>
Warning: number of  probes <6> for probeset <8180255> is not equal to annotated probe_count <0>
Warning: number of  probes <42> for probeset <8180256> is not equal to annotated probe_count <0>
Warning: number of  probes <43> for probeset <8180257> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180258> is not equal to annotated probe_count <0>
Warning: number of  probes <8> for probeset <8180259> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180260> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180261> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180262> is not equal to annotated probe_count <0>
Warning: number of  probes <51> for probeset <8180263> is not equal to annotated probe_count <0>
Warning: number of  probes <6> for probeset <8180264> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180265> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180266> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180267> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180268> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180269> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180270> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180271> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180272> is not equal to annotated probe_count <0>
Warning: number of  probes <15> for probeset <8180273> is not equal to annotated probe_count <0>
Warning: number of  probes <9> for probeset <8180274> is not equal to annotated probe_count <0>
Warning: number of  probes <10> for probeset <8180275> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180276> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180277> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180278> is not equal to annotated probe_count <0>
Warning: number of  probes <10> for probeset <8180279> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180280> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180281> is not equal to annotated probe_count <0>
Warning: number of  probes <8> for probeset <8180282> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180283> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180284> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180285> is not equal to annotated probe_count <0>
Warning: number of  probes <18> for probeset <8180286> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180287> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180288> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180289> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <8180290> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <8180291> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <8180292> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180293> is not equal to annotated probe_count <0>
Warning: number of  probes <7> for probeset <8180294> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180295> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <8180296> is not equal to annotated probe_count <0>
Warning: number of  probes <109> for probeset <8180297> is not equal to annotated probe_count <0>
Warning: number of  probes <37> for probeset <8180298> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180299> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180300> is not equal to annotated probe_count <0>
Warning: number of  probes <48> for probeset <8180301> is not equal to annotated probe_count <0>
Warning: number of  probes <47> for probeset <8180302> is not equal to annotated probe_count <0>
Warning: number of  probes <15> for probeset <8180303> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180304> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180305> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180306> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180307> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180308> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180309> is not equal to annotated probe_count <0>
Warning: number of  probes <49> for probeset <8180310> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180311> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180312> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180313> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180314> is not equal to annotated probe_count <0>
Warning: number of  probes <17> for probeset <8180315> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180316> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180317> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <8180318> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180319> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180320> is not equal to annotated probe_count <0>
Warning: number of  probes <6> for probeset <8180321> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180322> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180323> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180324> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180325> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180326> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180327> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180328> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180329> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180330> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180331> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180332> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180333> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180334> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180335> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <8180336> is not equal to annotated probe_count <0>
Warning: number of  probes <40> for probeset <8180337> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180338> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <8180339> is not equal to annotated probe_count <0>
Warning: number of  probes <242> for probeset <8180340> is not equal to annotated probe_count <0>
Warning: number of  probes <42> for probeset <8180341> is not equal to annotated probe_count <0>
Warning: number of  probes <40> for probeset <8180342> is not equal to annotated probe_count <0>
Warning: number of  probes <42> for probeset <8180343> is not equal to annotated probe_count <0>
Warning: number of  probes <39> for probeset <8180344> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180345> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180346> is not equal to annotated probe_count <0>
Warning: number of  probes <80> for probeset <8180347> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180348> is not equal to annotated probe_count <0>
Warning: number of  probes <37> for probeset <8180349> is not equal to annotated probe_count <0>
Warning: number of  probes <8> for probeset <8180350> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180351> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180352> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180353> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180354> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180355> is not equal to annotated probe_count <0>
Warning: number of  probes <49> for probeset <8180356> is not equal to annotated probe_count <0>
Warning: number of  probes <47> for probeset <8180357> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <8180358> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180359> is not equal to annotated probe_count <0>
Warning: number of  probes <14> for probeset <8180360> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180361> is not equal to annotated probe_count <0>
Warning: number of  probes <35> for probeset <8180362> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180363> is not equal to annotated probe_count <0>
Warning: number of  probes <65> for probeset <8180364> is not equal to annotated probe_count <0>
Warning: number of  probes <65> for probeset <8180365> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180366> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180367> is not equal to annotated probe_count <0>
Warning: number of  probes <48> for probeset <8180368> is not equal to annotated probe_count <0>
Warning: number of  probes <42> for probeset <8180369> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <8180370> is not equal to annotated probe_count <0>
Warning: number of  probes <40> for probeset <8180371> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180372> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <8180373> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180374> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <8180375> is not equal to annotated probe_count <0>
Warning: number of  probes <34> for probeset <8180376> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180377> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <8180378> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180379> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180380> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180381> is not equal to annotated probe_count <0>
Warning: number of  probes <36> for probeset <8180382> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180383> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180384> is not equal to annotated probe_count <0>
Warning: number of  probes <7> for probeset <8180385> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180386> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <8180387> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180388> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180389> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180390> is not equal to annotated probe_count <0>
Warning: number of  probes <56> for probeset <8180391> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180392> is not equal to annotated probe_count <0>
Warning: number of  probes <37> for probeset <8180393> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <8180394> is not equal to annotated probe_count <0>
Warning: number of  probes <21> for probeset <8180395> is not equal to annotated probe_count <0>
Warning: number of  probes <35> for probeset <8180396> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180397> is not equal to annotated probe_count <0>
Warning: number of  probes <48> for probeset <8180398> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180399> is not equal to annotated probe_count <0>
Warning: number of  probes <6> for probeset <8180400> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180401> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <8180402> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <8180403> is not equal to annotated probe_count <0>
Warning: number of  probes <34> for probeset <8180404> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <8180405> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <8180406> is not equal to annotated probe_count <0>
Warning: number of  probes <111> for probeset <8180407> is not equal to annotated probe_count <0>
Warning: number of  probes <121> for probeset <8180408> is not equal to annotated probe_count <0>
Warning: number of  probes <98> for probeset <8180409> is not equal to annotated probe_count <0>
Warning: number of  probes <101> for probeset <8180410> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180411> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <8180412> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <8180413> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <8180414> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180415> is not equal to annotated probe_count <0>
Warning: number of  probes <39> for probeset <8180416> is not equal to annotated probe_count <0>
Warning: number of  probes <33> for probeset <8180417> is not equal to annotated probe_count <0>
Warning: number of  probes <6> for probeset <8180418> is not equal to annotated probe_count <0>
Note: Number of exons imported <210414> is not equal to number of annotated exons <210415>.



# MoGene-1_1-st-v1.r4
scheme.mogene11stv1.na32 <- import.exon.scheme("mogene11stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "MoGene-1_1-st-v1.na32.mm9.probeset.csv"),
                            file.path(anndir, "Version11Jul", "MoGene-1_1-st-v1.na32.mm9.transcript.csv"))

   <241531> probeset tree entries read...Finished
Warning: number of  probes <34> for probeset <10608634> is not equal to annotated probe_count <0>
Warning: number of  probes <79> for probeset <10608635> is not equal to annotated probe_count <0>
Warning: number of  probes <57> for probeset <10608636> is not equal to annotated probe_count <0>
Warning: number of  probes <12> for probeset <10608637> is not equal to annotated probe_count <0>
Warning: number of  probes <15> for probeset <10608638> is not equal to annotated probe_count <0>
Warning: number of  probes <49> for probeset <10608639> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608640> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <10608641> is not equal to annotated probe_count <0>
Warning: number of  probes <61> for probeset <10608642> is not equal to annotated probe_count <0>
Warning: number of  probes <75> for probeset <10608643> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <10608644> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <10608645> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <10608646> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <10608647> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <10608648> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <10608649> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <10608650> is not equal to annotated probe_count <0>
Warning: number of  probes <13> for probeset <10608651> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <10608652> is not equal to annotated probe_count <0>
Warning: number of  probes <100> for probeset <10608653> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <10608654> is not equal to annotated probe_count <0>
Warning: number of  probes <51> for probeset <10608655> is not equal to annotated probe_count <0>
Warning: number of  probes <49> for probeset <10608656> is not equal to annotated probe_count <0>
Warning: number of  probes <36> for probeset <10608657> is not equal to annotated probe_count <0>
Warning: number of  probes <17> for probeset <10608658> is not equal to annotated probe_count <0>
Warning: number of  probes <15> for probeset <10608659> is not equal to annotated probe_count <0>
Warning: number of  probes <40> for probeset <10608660> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <10608661> is not equal to annotated probe_count <0>
Warning: number of  probes <27> for probeset <10608662> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <10608663> is not equal to annotated probe_count <0>
Warning: number of  probes <35> for probeset <10608664> is not equal to annotated probe_count <0>
Warning: number of  probes <29> for probeset <10608665> is not equal to annotated probe_count <0>
Warning: number of  probes <9> for probeset <10608666> is not equal to annotated probe_count <0>
Warning: number of  probes <32> for probeset <10608667> is not equal to annotated probe_count <0>
Warning: number of  probes <11> for probeset <10608668> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608669> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <10608670> is not equal to annotated probe_count <0>
Warning: number of  probes <35> for probeset <10608671> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <10608672> is not equal to annotated probe_count <0>
Warning: number of  probes <85> for probeset <10608674> is not equal to annotated probe_count <0>
Warning: number of  probes <15> for probeset <10608675> is not equal to annotated probe_count <0>
Warning: number of  probes <28> for probeset <10608676> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <10608677> is not equal to annotated probe_count <0>
Warning: number of  probes <19> for probeset <10608678> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <10608679> is not equal to annotated probe_count <0>
Warning: number of  probes <31> for probeset <10608680> is not equal to annotated probe_count <0>
Warning: number of  probes <22> for probeset <10608681> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608682> is not equal to annotated probe_count <0>
Warning: number of  probes <34> for probeset <10608683> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608684> is not equal to annotated probe_count <0>
Warning: number of  probes <30> for probeset <10608685> is not equal to annotated probe_count <0>
Warning: number of  probes <23> for probeset <10608686> is not equal to annotated probe_count <0>
Warning: number of  probes <95> for probeset <10608687> is not equal to annotated probe_count <0>
Warning: number of  probes <105> for probeset <10608688> is not equal to annotated probe_count <0>
Warning: number of  probes <24> for probeset <10608689> is not equal to annotated probe_count <0>
Warning: number of  probes <35> for probeset <10608690> is not equal to annotated probe_count <0>
Warning: number of  probes <26> for probeset <10608691> is not equal to annotated probe_count <0>
Warning: number of  probes <16> for probeset <10608692> is not equal to annotated probe_count <0>
Warning: number of  probes <43> for probeset <10608695> is not equal to annotated probe_count <0>
Warning: number of  probes <88> for probeset <10608696> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608697> is not equal to annotated probe_count <0>
Warning: number of  probes <20> for probeset <10608698> is not equal to annotated probe_count <0>
Warning: number of  probes <74> for probeset <10608699> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608700> is not equal to annotated probe_count <0>
Warning: number of  probes <74> for probeset <10608701> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608702> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608703> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608704> is not equal to annotated probe_count <0>
Warning: number of  probes <133> for probeset <10608705> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608706> is not equal to annotated probe_count <0>
Warning: number of  probes <25> for probeset <10608707> is not equal to annotated probe_count <0>
Warning: number of  probes <77> for probeset <10608708> is not equal to annotated probe_count <0>
Warning: number of  probes <56> for probeset <10608709> is not equal to annotated probe_count <0>
Warning: number of  probes <100> for probeset <10608710> is not equal to annotated probe_count <0>
Warning: number of  probes <94> for probeset <10608711> is not equal to annotated probe_count <0>
Warning: number of  probes <87> for probeset <10608712> is not equal to annotated probe_count <0>
Warning: number of  probes <49> for probeset <10608713> is not equal to annotated probe_count <0>
Warning: number of  probes <86> for probeset <10608714> is not equal to annotated probe_count <0>
Warning: number of  probes <70> for probeset <10608715> is not equal to annotated probe_count <0>
Warning: number of  probes <98> for probeset <10608716> is not equal to annotated probe_count <0>
Warning: number of  probes <80> for probeset <10608717> is not equal to annotated probe_count <0>
Warning: number of  probes <165> for probeset <10608718> is not equal to annotated probe_count <0>
Warning: number of  probes <97> for probeset <10608719> is not equal to annotated probe_count <0>
Warning: number of  probes <60> for probeset <10608720> is not equal to annotated probe_count <0>
Warning: number of  probes <65> for probeset <10608721> is not equal to annotated probe_count <0>
Warning: number of  probes <41> for probeset <10608722> is not equal to annotated probe_count <0>
Warning: number of  probes <83> for probeset <10608723> is not equal to annotated probe_count <0>
Warning: number of  probes <50> for probeset <10608724> is not equal to annotated probe_count <0>
Note: Number of exons imported <207392> is not equal to number of annotated exons <207393>.



# RaGene-1_1-st-v1.r4
# use updated Affymetrix annotation files from 09/17/10
scheme.ragene11stv1.na32 <- import.exon.scheme("ragene11stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "RaGene-1_1-st-v1.na32.rn4.probeset.csv"),
                            file.path(anndir, "Version11Jul", "RaGene-1_1-st-v1.na32.rn4.transcript.csv"))



#------------------------------------------------------------------------------#
# Sep 2010: Affymetrix annotation "na31"
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na31 <- import.expr.scheme("test3", filedir = file.path(scmdir, "na31"),
                     schemefile = file.path(libdir, "Test3.CDF"), 
                     probefile  = file.path(libdir, "Test3_probe.tab"), 
                     annotfile  = file.path(anndir, "Version10Sep", "Test3.na31.annot.csv"))

# Hu6800: 
# annotation file for na31 no longer provided by Affymetrix (see na30 annotations)

# HG_U95A:
# annotation file for na31 no longer provided by Affymetrix (see na30 annotations)

# HG_U95Av2:
# annotation file for na31 no longer provided by Affymetrix (see na30 annotations)

# HG-U133A:
scheme.hgu133a.na31 <- import.expr.scheme("hgu133a", filedir = file.path(scmdir, "na31"),
                       schemefile = file.path(libdir, "HG-U133A.CDF"), 
                       probefile  = file.path(libdir, "HG-U133A_probe.tab"), 
                       annotfile  = file.path(anndir, "Version10Sep", "HG-U133A.na31.annot.csv"))

# HG-U133B:
scheme.hgu133b.na31 <- import.expr.scheme("hgu133b", filedir = file.path(scmdir, "na31"),
                       schemefile = file.path(libdir, "HG-U133B.CDF"), 
                       probefile  = file.path(libdir, "HG-U133B_probe.tab"), 
                       annotfile  = file.path(anndir, "Version10Sep", "HG-U133B.na31.annot.csv"))

# HG-U133_Plus_2:
scheme.hgu133plus2.na31 <- import.expr.scheme("hgu133plus2", filedir = file.path(scmdir, "na31"),
                           schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"), 
                           probefile  = file.path(libdir, "HG-U133-PLUS_probe.tab"), 
                           annotfile  = file.path(anndir, "Version10Sep", "HG-U133_Plus_2.na31.annot.csv"))

# MOE430A:
scheme.moe430a.na31 <- import.expr.scheme("moe430a",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"MOE430A.CDF"),
                       probefile  = file.path(libdir,"MOE430A.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/MOE430A.na31.annot.csv"))

# MOE430B:
scheme.moe430b.na31 <- import.expr.scheme("moe430b",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"MOE430B.CDF"),
                       probefile  = file.path(libdir,"MOE430B.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/MOE430B.na31.annot.csv"))

# Mouse430_2:
scheme.mouse4302.na31 <- import.expr.scheme("mouse4302",filedir=file.path(scmdir, "na31"),
                         schemefile = file.path(libdir,"Mouse430_2.cdf"),
                         probefile  = file.path(libdir,"Mouse430_2.probe.tab"),
                         annotfile  = file.path(anndir,"Version10Sep/Mouse430_2.na31.annot.csv"))

# RG_U34A:
scheme.rgu34a.na31 <- import.expr.scheme("rgu34a",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"RG_U34A.cdf"),
                       probefile  = file.path(libdir,"RG_U34A.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/RG_U34A.na31.annot.csv"))

# RG_U34B:
scheme.rgu34b.na31 <- import.expr.scheme("rgu34b",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"RG_U34B.cdf"),
                       probefile  = file.path(libdir,"RG_U34B.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/RG_U34B.na31.annot.csv"))

# RG_U34C:
scheme.rgu34c.na31 <- import.expr.scheme("rgu34c",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"RG_U34C.cdf"),
                       probefile  = file.path(libdir,"RG_U34C.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/RG_U34C.na31.annot.csv"))

# RAE230A:
scheme.rae230a.na31 <- import.expr.scheme("rae230a",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"RAE230A.CDF"),
                       probefile  = file.path(libdir,"RAE230A.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/RAE230A.na31.annot.csv"))

# RAE230B:
scheme.rae230b.na31 <- import.expr.scheme("rae230b",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"RAE230B.CDF"),
                       probefile  = file.path(libdir,"RAE230B.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/RAE230B.na31.annot.csv"))

# Rat230_2:
scheme.rat2302.na31 <- import.expr.scheme("rat2302",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"Rat230_2.cdf"),
                       probefile  = file.path(libdir,"Rat230_2.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/Rat230_2.na31.annot.csv"))

# Bovine:
scheme.bovine.na31 <- import.expr.scheme("bovine",filedir=file.path(scmdir, "na31"),
                      schemefile = file.path(libdir,"Bovine.cdf"),
                      probefile  = file.path(libdir,"Bovine.probe.tab"),
                      annotfile  = file.path(anndir,"Version10Sep/Bovine.na31.annot.csv"))

# Porcine:
scheme.porcine.na31 <- import.expr.scheme("porcine",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"Porcine.cdf"),
                       probefile  = file.path(libdir,"Porcine.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/Porcine.na31.annot.csv"))

# Rhesus:
scheme.rhesus.na31 <- import.expr.scheme("rhesus",filedir=file.path(scmdir, "na31"),
                      schemefile = file.path(libdir,"Rhesus.cdf"),
                      probefile  = file.path(libdir,"Rhesus.probe.tab"),
                      annotfile  = file.path(anndir,"Version10Sep/Rhesus.na31.annot.csv"))

# Zebrafish:
scheme.zebrafish.na31 <- import.expr.scheme("zebrafish",filedir=file.path(scmdir, "na31"),
                         schemefile = file.path(libdir,"Zebrafish.cdf"),
                         probefile  = file.path(libdir,"Zebrafish.probe.tab"),
                         annotfile  = file.path(anndir,"Version10Sep/Zebrafish.na31.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HT_HG-U133A
scheme.hthgu133a.na31 <- import.expr.scheme("hthgu133a", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "HT_HG-U133A.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version10Sep", "HT_HG-U133A.na31.annot.csv"))

# HT_HG-U133B
scheme.hthgu133b.na31 <- import.expr.scheme("hthgu133b", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "HT_HG-U133B.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133B.probe.tab"), 
                         annotfile  = file.path(anndir, "Version10Sep", "HT_HG-U133B.na31.annot.csv"))

# HT_HG-U133_Plus_PM
scheme.hthgu133pluspm.na31 <- import.expr.scheme("hthgu133pluspm", filedir = file.path(scmdir, "na31"),
                              schemefile = file.path(libdir, "HT_HG-U133_Plus_PM.CDF"), 
                              probefile  = file.path(libdir, "HT_HG-U133_Plus_PM.probe.tab"), 
                              annotfile  = file.path(anndir, "Version10Sep", "HT_HG-U133_Plus_PM.na31.annot.csv"))

# HT_MG-430A
scheme.htmg430a.na31 <- import.expr.scheme("htmg430a", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "HT_MG-430A.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version10Sep", "HT_MG-430A.na31.annot.csv"))

# HT_MG-430B
scheme.htmg430b.na31 <- import.expr.scheme("htmg430b", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "HT_MG-430B.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430B.probe.tab"), 
                         annotfile  = file.path(anndir, "Version10Sep", "HT_MG-430B.na31.annot.csv"))

# HT_MG-430_PM
scheme.htmg430pm.na31 <- import.expr.scheme("htmg430pm", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "HT_MG-430_PM.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430_PM.probe.tab"), 
                         annotfile  = file.path(anndir, "Version10Sep", "HT_MG-430_PM.na31.annot.csv"))

# HT_Rat230_PM
scheme.htrat230pm.na31 <- import.expr.scheme("htrat230pm", filedir = file.path(scmdir, "na31"),
                          schemefile = file.path(libdir, "HT_Rat230_PM.cdf"), 
                          probefile  = file.path(libdir, "HT_Rat230_PM.probe.tab"), 
                          annotfile  = file.path(anndir, "Version10Sep", "HT_Rat230_PM.na31.annot.csv"))

# HG-U219
scheme.hgu219.na31 <- import.expr.scheme("hgu219", filedir = file.path(scmdir, "na31"),
                       schemefile = file.path(libdir, "HG-U219.cdf"), 
                       probefile  = file.path(libdir, "HG-U219.probe.tab"), 
                       annotfile  = file.path(anndir, "Version10Sep", "HG-U219.na31.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1.na31 <- import.exon.scheme("hugene10stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "HuGene-1_0-st-v1.na31.hg19.probeset.csv"),
                            file.path(anndir, "Version10Sep", "HuGene-1_0-st-v1.na31.hg19.transcript.csv"))

# MoGene-1_0-st-v1.r4: used as exon array
scheme.mogene10stv1.na31 <- import.exon.scheme("mogene10stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "MoGene-1_0-st-v1.na31.mm9.probeset.csv"),
                            file.path(anndir, "Version10Sep", "MoGene-1_0-st-v1.na31.mm9.transcript.csv"))

# RaGene-1_0-st-v1.r4: used as exon array
scheme.ragene10stv1.na31 <- import.exon.scheme("ragene10stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "RaGene-1_0-st-v1.na31.rn4.probeset.csv"),
                            file.path(anndir, "Version10Sep", "RaGene-1_0-st-v1.na31.rn4.transcript.csv"))

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2.na31 <- import.exon.scheme("huex10stv2", filedir = file.path(scmdir, "na31"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                          file.path(anndir, "Version10Sep", "HuEx-1_0-st-v2.na31.hg19.probeset.csv"),
                          file.path(anndir, "Version10Sep", "HuEx-1_0-st-v2.na31.hg19.transcript.csv"))

# MoEx-1_0-st-v1.r2:
# use updated Affymetrix annotation files from 09/08/10
scheme.moex10stv1.na31 <- import.exon.scheme("moex10stv1",filedir = file.path(scmdir, "na31"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version10Sep", "MoEx-1_0-st-v1.na31.mm9.probeset.csv"),
                          file.path(anndir, "Version10Sep", "MoEx-1_0-st-v1.na31.mm9.transcript.csv"))

# original annotation files from 08/30/10: need to delete control->affx which are neg_control
source(paste(.path.package("xps"),"examples/updateAnnotation.R",sep="/"))
deleteNegControlFromAffxControl("MoEx-1_0-st-v1.na31.mm9.probeset.csv", "MoEx-1_0-st-v1.na31.mm9.probeset.fixed.csv", eol="\n")
deleteNegControlFromAffxControl("MoEx-1_0-st-v1.na31.mm9.transcript.csv", "MoEx-1_0-st-v1.na31.mm9.transcript.fixed.csv", eol="\n")
# use fixed annotation files
scheme.moex10stv1.na31 <- import.exon.scheme("moex10stv1",filedir = file.path(scmdir, "na31"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version10Sep", "MoEx-1_0-st-v1.na31.mm9.probeset.fixed.csv"),
                          file.path(anndir, "Version10Sep", "MoEx-1_0-st-v1.na31.mm9.transcript.fixed.csv"))

# RaEx-1_0-st-v1.r2:
scheme.raex10stv1.na31 <- import.exon.scheme("raex10stv1",filedir=file.path(scmdir, "na31"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version10Sep", "RaEx-1_0-st-v1.na31.rn4.probeset.csv"),
                          file.path(anndir, "Version10Sep", "RaEx-1_0-st-v1.na31.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_1-st-v1:
scheme.hugene11stv1.na31 <- import.exon.scheme("hugene11stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "HuGene-1_1-st-v1.na31.hg19.probeset.csv"),
                            file.path(anndir, "Version10Sep", "HuGene-1_1-st-v1.na31.hg19.transcript.csv"))

# MoGene-1_1-st-v1.r4
scheme.mogene11stv1.na31 <- import.exon.scheme("mogene11stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "MoGene-1_1-st-v1.na31.mm9.probeset.csv"),
                            file.path(anndir, "Version10Sep", "MoGene-1_1-st-v1.na31.mm9.transcript.csv"))

# RaGene-1_1-st-v1.r4
# use updated Affymetrix annotation files from 09/17/10
scheme.ragene11stv1.na31 <- import.exon.scheme("ragene11stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "RaGene-1_1-st-v1.na31.rn4.probeset.csv"),
                            file.path(anndir, "Version10Sep", "RaGene-1_1-st-v1.na31.rn4.transcript.csv"))

# original annotation files from 08/27/10: need to add missing header lines
scheme.ragene11stv1.na31 <- import.exon.scheme("ragene11stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "RaGene-1_1-st-v1.na31.rn4.probeset.fixed.csv"),
                            file.path(anndir, "Version10Sep", "RaGene-1_1-st-v1.na31.rn4.transcript.fixed.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for other Affymetrix arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# E_coli_2:
scheme.ecoli2.na31 <- import.expr.scheme("ecoli2", filedir = file.path(scmdir, "na31"),
                      schemefile = file.path(libdir, "E_coli_2.cdf"),
                      probefile  = file.path(libdir, "E_coli_2.probe.tab"),
                      annotfile  = file.path(anndir, "Version10Sep", "E_coli_2.na31.annot.csv"))

# Ecoli_ASv2:
scheme.ecoliasv2.na31 <- import.expr.scheme("ecoliasv2", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "Ecoli_ASv2.CDF"),
                         probefile  = file.path(libdir, "Ecoli_ASv2.probe.tab"),
                         annotfile  = file.path(anndir, "Version10Sep", "Ecoli_ASv2.na31.annot.csv"))

# Pae_G1a:
scheme.paeg1a.na31 <- import.expr.scheme("paeg1a", filedir = file.path(scmdir, "na31"),
                      schemefile = file.path(libdir, "Pae_G1a.CDF"),
                      probefile  = file.path(libdir, "Pae_G1a.probe.tab"), 
                      annotfile  = file.path(anndir, "Version10Sep", "Pae_G1a.na31.annot.csv"))

# Yeast_2:
scheme.yeast2.na31 <- import.expr.scheme("yeast2", filedir = file.path(scmdir, "na31"),
                      schemefile = file.path(libdir, "Yeast_2.cdf"),
                      probefile  = file.path(libdir, "Yeast_2.probe.tab"),
                      annotfile  = file.path(anndir, "Version10Sep", "Yeast_2.na31.annot.csv"))

# miRNA-1_0:
# note: you need to rename "miRNA-1_0.probe_list.20081203.txt" to "miRNA-1_0.probe.tab"
scheme.mirna10 <- import.expr.scheme("mirna10", filedir = file.path(scmdir, "na31"),
                  schemefile = file.path(libdir, "miRNA-1_0.CDF"),
                  probefile  = file.path(libdir, "miRNA-1_0.probe.tab"))



#------------------------------------------------------------------------------#
# Nov 2009: Affymetrix annotation "na30"
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na30 <- import.expr.scheme("test3", filedir = file.path(scmdir, "na30"),
                     schemefile = file.path(libdir, "Test3.CDF"), 
                     probefile  = file.path(libdir, "Test3_probe.tab"), 
                     annotfile  = file.path(anndir, "Version09Nov", "Test3.na30.annot.csv"))

# Hu6800:
scheme.hu6800.na30 <- import.expr.scheme("hu6800", filedir = file.path(scmdir, "na30"),
                      schemefile = file.path(libdir, "Hu6800.CDF"), 
                      probefile  = file.path(libdir, "HuGeneFL_probe.tab"), 
                      annotfile  = file.path(anndir, "Version09Nov", "Hu6800.na30.annot.csv"))

# HG_U95A:
scheme.hgu95a.na30 <- import.expr.scheme("hgu95a", filedir = file.path(scmdir, "na30"),
                      schemefile = file.path(libdir, "HG_U95A.CDF"), 
                      probefile  = file.path(libdir, "HG-U95A_probe.tab"), 
                      annotfile  = file.path(anndir, "Version09Nov", "HG_U95A.na30.annot.csv"))

# HG_U95Av2:
scheme.hgu95av2.na30 <- import.expr.scheme("hgu95av2", filedir = file.path(scmdir, "na30"),
                        schemefile = file.path(libdir, "HG_U95Av2.CDF"), 
                        probefile  = file.path(libdir, "HG-U95Av2_probe.tab"), 
                        annotfile  = file.path(anndir, "Version09Nov", "HG_U95Av2.na30.annot.csv"))

# HG-U133A:
scheme.hgu133a.na30 <- import.expr.scheme("hgu133a", filedir = file.path(scmdir, "na30"),
                       schemefile = file.path(libdir, "HG-U133A.CDF"), 
                       probefile  = file.path(libdir, "HG-U133A_probe.tab"), 
                       annotfile  = file.path(anndir, "Version09Nov", "HG-U133A.na30.annot.csv"))

# HG-U133B:
scheme.hgu133b.na30 <- import.expr.scheme("hgu133b", filedir = file.path(scmdir, "na30"),
                       schemefile = file.path(libdir, "HG-U133B.CDF"), 
                       probefile  = file.path(libdir, "HG-U133B_probe.tab"), 
                       annotfile  = file.path(anndir, "Version09Nov", "HG-U133B.na30.annot.csv"))

# HG-U133_Plus_2:
scheme.hgu133plus2.na30 <- import.expr.scheme("hgu133plus2", filedir = file.path(scmdir, "na30"),
                           schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"), 
                           probefile  = file.path(libdir, "HG-U133-PLUS_probe.tab"), 
                           annotfile  = file.path(anndir, "Version09Nov", "HG-U133_Plus_2.na30.annot.csv"))

# MOE430A:
scheme.moe430a.na30 <- import.expr.scheme("moe430a",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"MOE430A.CDF"),
                       probefile  = file.path(libdir,"MOE430A.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/MOE430A.na30.annot.csv"))

# MOE430B:
scheme.moe430b.na30 <- import.expr.scheme("moe430b",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"MOE430B.CDF"),
                       probefile  = file.path(libdir,"MOE430B.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/MOE430B.na30.annot.csv"))

# Mouse430_2:
scheme.mouse4302.na30 <- import.expr.scheme("mouse4302",filedir=file.path(scmdir, "na30"),
                         schemefile = file.path(libdir,"Mouse430_2.cdf"),
                         probefile  = file.path(libdir,"Mouse430_2.probe.tab"),
                         annotfile  = file.path(anndir,"Version09Nov/Mouse430_2.na30.annot.csv"))

# RG_U34A:
scheme.rgu34a.na30 <- import.expr.scheme("rgu34a",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"RG_U34A.cdf"),
                       probefile  = file.path(libdir,"RG_U34A.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RG_U34A.na30.annot.csv"))

# RG_U34B:
scheme.rgu34b.na30 <- import.expr.scheme("rgu34b",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"RG_U34B.cdf"),
                       probefile  = file.path(libdir,"RG_U34B.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RG_U34B.na30.annot.csv"))

# RG_U34C:
scheme.rgu34c.na30 <- import.expr.scheme("rgu34c",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"RG_U34C.cdf"),
                       probefile  = file.path(libdir,"RG_U34C.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RG_U34C.na30.annot.csv"))

# RAE230A:
scheme.rae230a.na30 <- import.expr.scheme("rae230a",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"RAE230A.CDF"),
                       probefile  = file.path(libdir,"RAE230A.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RAE230A.na30.annot.csv"))

# RAE230B:
scheme.rae230b.na30 <- import.expr.scheme("rae230b",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"RAE230B.CDF"),
                       probefile  = file.path(libdir,"RAE230B.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RAE230B.na30.annot.csv"))

# Rat230_2:
scheme.rat2302.na30 <- import.expr.scheme("rat2302",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"Rat230_2.cdf"),
                       probefile  = file.path(libdir,"Rat230_2.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/Rat230_2.na30.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HT_HG-U133A
scheme.hthgu133a.na30 <- import.expr.scheme("hthgu133a", filedir = file.path(scmdir, "na30"),
                         schemefile = file.path(libdir, "HT_HG-U133A.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version09Nov", "HT_HG-U133A.na30.annot.csv"))

# HT_HG-U133_Plus_PM
scheme.hthgu133pluspm.na30 <- import.expr.scheme("hthgu133pluspm", filedir = file.path(scmdir, "na30"),
                              schemefile = file.path(libdir, "HT_HG-U133_Plus_PM.CDF"), 
                              probefile  = file.path(libdir, "HT_HG-U133_Plus_PM.probe.tab"), 
                              annotfile  = file.path(anndir, "Version09Nov", "HT_HG-U133_Plus_PM.na30.annot.csv"))

# HT_MG-430A
scheme.htmg430a.na30 <- import.expr.scheme("htmg430a", filedir = file.path(scmdir, "na30"),
                         schemefile = file.path(libdir, "HT_MG-430A.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version09Nov", "HT_MG-430A.na30.annot.csv"))

# HT_MG-430_PM
scheme.htmg430pm.na30 <- import.expr.scheme("htmg430pm", filedir = file.path(scmdir, "na30"),
                         schemefile = file.path(libdir, "HT_MG-430_PM.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430_PM.probe.tab"), 
                         annotfile  = file.path(anndir, "Version09Nov", "HT_MG-430_PM.na30.annot.csv"))

# HT_Rat230_PM
scheme.htrat230pm.na30 <- import.expr.scheme("htrat230pm", filedir = file.path(scmdir, "na30"),
                          schemefile = file.path(libdir, "HT_Rat230_PM.cdf"), 
                          probefile  = file.path(libdir, "HT_Rat230_PM.probe.tab"), 
                          annotfile  = file.path(anndir, "Version09Nov", "HT_Rat230_PM.na30.annot.csv"))

# HG-U219
scheme.hgu219.na30 <- import.expr.scheme("hgu219", filedir = file.path(scmdir, "na30"),
                       schemefile = file.path(libdir, "HG-U219.cdf"), 
                       probefile  = file.path(libdir, "HG-U219.probe.tab"), 
                       annotfile  = file.path(anndir, "Version09Nov", "HG-U219.na30.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1.na30 <- import.exon.scheme("hugene10stv1", filedir = file.path(scmdir, "na30"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.1.hg19.probeset.csv"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.hg19.transcript.csv"))

# MoGene-1_0-st-v1.r4: used as exon array
scheme.mogene10stv1.na30 <- import.exon.scheme("mogene10stv1", filedir = file.path(scmdir, "na30"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.probeset.csv"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.transcript.csv"))

# RaGene-1_0-st-v1.r4: used as exon array
scheme.ragene10stv1.na30 <- import.exon.scheme("ragene10stv1", filedir = file.path(scmdir, "na30"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.probeset.csv"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.transcript.csv"))

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2.na30 <- import.exon.scheme("huex10stv2", filedir = file.path(scmdir, "na30"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                          file.path(anndir, "Version09Nov", "HuEx-1_0-st-v2.na30.hg19.probeset.csv"),
                          file.path(anndir, "Version09Nov", "HuEx-1_0-st-v2.na30.hg19.transcript.csv"))

# MoEx-1_0-st-v1.r2:
scheme.moex10stv1.na30 <- import.exon.scheme("moex10stv1",filedir = file.path(scmdir, "na30"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version09Nov", "MoEx-1_0-st-v1.na30.mm9.probeset.csv"),
                          file.path(anndir, "Version09Nov", "MoEx-1_0-st-v1.na30.mm9.transcript.csv"))

# RaEx-1_0-st-v1.r2:
scheme.raex10stv1.na30 <- import.exon.scheme("raex10stv1",filedir=file.path(scmdir, "na30"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version09Nov", "RaEx-1_0-st-v1.na30.rn4.probeset.csv"),
                          file.path(anndir, "Version09Nov", "RaEx-1_0-st-v1.na30.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_1-st-v1:
scheme.hugene11stv1.na30 <- import.exon.scheme("hugene11stv1", filedir = file.path(scmdir, "na30"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.1.hg19.probeset.csv"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.hg19.transcript.csv"))

# MoGene-1_1-st-v1.r4
scheme.mogene11stv1.na30 <- import.exon.scheme("mogene11stv1", filedir = file.path(scmdir, "na30"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.probeset.csv"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.transcript.csv"))

# RaGene-1_1-st-v1.r4
scheme.ragene11stv1.na30 <- import.exon.scheme("ragene11stv1", filedir = scmdir,
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.probeset.csv"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for other Affymetrix arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# E_coli_2:
scheme.ecoli2.na30 <- import.expr.scheme("ecoli2", filedir = file.path(scmdir, "na30"),
                      schemefile = file.path(libdir, "E_coli_2.cdf"),
                      probefile  = file.path(libdir, "E_coli_2.probe.tab"),
                      annotfile  = file.path(anndir, "Version09Nov", "E_coli_2.na30.annot.csv"))

# Ecoli_ASv2:
scheme.ecoliasv2.na30 <- import.expr.scheme("ecoliasv2", filedir = file.path(scmdir, "na30"),
                         schemefile = file.path(libdir, "Ecoli_ASv2.CDF"),
                         probefile  = file.path(libdir, "Ecoli_ASv2.probe.tab"),
                         annotfile  = file.path(anndir, "Version09Nov", "Ecoli_ASv2.na30.annot.csv"))

# Pae_G1a:
scheme.paeg1a.na30 <- import.expr.scheme("paeg1a", filedir = file.path(scmdir, "na30"),
                      schemefile = file.path(libdir, "Pae_G1a.CDF"),
                      probefile  = file.path(libdir, "Pae_G1a.probe.tab"), 
                      annotfile  = file.path(anndir, "Version09Nov", "Pae_G1a.na30.annot.csv"))

# Yeast_2:
scheme.yeast2.na30 <- import.expr.scheme("yeast2", filedir = file.path(scmdir, "na30"),
                      schemefile = file.path(libdir, "Yeast_2.cdf"),
                      probefile  = file.path(libdir, "Yeast_2.probe.tab"),
                      annotfile  = file.path(anndir, "Version09Nov", "Yeast_2.na30.annot.csv"))

# miRNA-1_0:
# note: you need to rename "miRNA-1_0.probe_list.20081203.txt" to "miRNA-1_0.probe.tab"
scheme.mirna10 <- import.expr.scheme("mirna10", filedir = file.path(scmdir, "na30"),
                  schemefile = file.path(libdir, "miRNA-1_0.CDF"),
                  probefile  = file.path(libdir, "miRNA-1_0.probe.tab"))



#------------------------------------------------------------------------------#
# schemes for older Affymetrix annotation files
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na21 <- import.expr.scheme("Scheme_Test3_na21",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version06Nov/Test3.na21.annot.csv",sep="/"))
scheme.test3.na22 <- import.expr.scheme("Scheme_Test3_na22",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version07May/Test3.na22.annot.csv",sep="/"))
scheme.test3.na23 <- import.expr.scheme("Scheme_Test3_na23",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version07Jul/Test3.na23.annot.csv",sep="/"))
scheme.test3.na24 <- import.expr.scheme("Scheme_Test3_na24",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version07Nov/Test3.na24.annot.csv",sep="/"))
scheme.test3.na25 <- import.expr.scheme("Scheme_Test3_na25",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version08Mar/Test3.na25.annot.csv",sep="/"))
scheme.test3.na26 <- import.expr.scheme("Scheme_Test3_na26",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version08Jul/Test3.na26.annot.csv",sep="/"))
scheme.test3.na27 <- import.expr.scheme("Scheme_Test3_na27",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version08Nov/Test3.na27.annot.csv",sep="/"))
scheme.test3.na29 <- import.expr.scheme("Scheme_Test3_na29",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version09Jul/Test3.na29.annot.csv",sep="/"))

# Hu6800:
scheme.hu6800.na21 <- import.expr.scheme("Scheme_Hu6800_na21",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version06Nov/Hu6800.na21.annot.csv",sep="/"))
scheme.hu6800.na22 <- import.expr.scheme("Scheme_Hu6800_na22",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version07May/Hu6800.na22.annot.csv",sep="/"))
scheme.hu6800.na23 <- import.expr.scheme("Scheme_Hu6800_na23",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version07Jul/Hu6800.na23.annot.csv",sep="/"))
scheme.hu6800.na24 <- import.expr.scheme("Scheme_Hu6800_na24",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version07Nov/Hu6800.na24.annot.csv",sep="/"))
scheme.hu6800.na25 <- import.expr.scheme("Scheme_Hu6800_na25",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version08Mar/Hu6800.na25.annot.csv",sep="/"))

# HG_U95A:
scheme.hgu95a <- import.expr.scheme("Scheme_HGU95A",filedir=scmdir,paste(libdir,"HG_U95A.CDF",sep="/"),paste(libdir,"HG-U95A_probe.tab",sep="/"),annotfile="")
scheme.hgu95a.na23 <- import.expr.scheme("Scheme_HGU95A_na23",filedir=scmdir,paste(libdir,"HG_U95A.CDF",sep="/"),paste(libdir,"HG-U95A_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG_U95A.na23.annot.csv",sep="/"))
scheme.hgu95a.na24 <- import.expr.scheme("Scheme_HGU95A_na24",filedir=scmdir,paste(libdir,"HG_U95A.CDF",sep="/"),paste(libdir,"HG-U95A_probe.tab",sep="/"),paste(anndir,"Version07Nov/HG_U95A.na24.annot.csv",sep="/"))
scheme.hgu95a.na25 <- import.expr.scheme("Scheme_HGU95A_na25",filedir=scmdir,paste(libdir,"HG_U95A.CDF",sep="/"),paste(libdir,"HG-U95A_probe.tab",sep="/"),paste(anndir,"Version08Mar/HG_U95A.na25.annot.csv",sep="/"))

# HG_U95Av2:
scheme.hgu95av2.na21 <- import.expr.scheme("Scheme_HGU95Av2_na21",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG_U95Av2.na21.annot.csv",sep="/"))
scheme.hgu95av2.na22 <- import.expr.scheme("Scheme_HGU95Av2_na22",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version07May/HG_U95Av2.na22.annot.csv",sep="/"))
scheme.hgu95av2.na23 <- import.expr.scheme("Scheme_HGU95Av2_na23",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG_U95Av2.na23.annot.csv",sep="/"))
scheme.hgu95av2.na24 <- import.expr.scheme("Scheme_HGU95Av2_na24",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version07Nov/HG_U95Av2.na24.annot.csv",sep="/"))
scheme.hgu95av2.na25 <- import.expr.scheme("Scheme_HGU95Av2_na25",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version08Mar/HG_U95Av2.na25.annot.csv",sep="/"))

# HG-U133A:
scheme.hgu133a.na21 <- import.expr.scheme("Scheme_HGU133A_na21",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG-U133A.na21.annot.csv",sep="/"))
scheme.hgu133a.na22 <- import.expr.scheme("Scheme_HGU133A_na22",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version07May/HG-U133A.na22.annot.csv",sep="/"))
scheme.hgu133a.na23 <- import.expr.scheme("Scheme_HGU133A_na23",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG-U133A.na23.annot.csv",sep="/"))
scheme.hgu133a.na24 <- import.expr.scheme("Scheme_HGU133A_na24",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version07Nov/HG-U133A.na24.annot.csv",sep="/"))
scheme.hgu133a.na25 <- import.expr.scheme("Scheme_HGU133A_na25",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version08Mar/HG-U133A.na25.annot.csv",sep="/"))

# HG-U133B:
scheme.hgu133b.na21 <- import.expr.scheme("Scheme_HGU133B_na21",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG-U133B.na21.annot.csv",sep="/"))
scheme.hgu133b.na22 <- import.expr.scheme("Scheme_HGU133B_na22",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version07May/HG-U133B.na22.annot.csv",sep="/"))
scheme.hgu133b.na23 <- import.expr.scheme("Scheme_HGU133B_na23",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG-U133B.na23.annot.csv",sep="/"))
scheme.hgu133b.na24 <- import.expr.scheme("Scheme_HGU133B_na24",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version07Nov/HG-U133B.na24.annot.csv",sep="/"))
scheme.hgu133b.na25 <- import.expr.scheme("Scheme_HGU133B_na25",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version08Mar/HG-U133B.na25.annot.csv",sep="/"))

# HG-U133_Plus_2:
scheme.hgu133p2.na21 <- import.expr.scheme("Scheme_HGU133p2_na21",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG-U133_Plus_2.na21.annot.csv",sep="/"))
scheme.hgu133p2.na22 <- import.expr.scheme("Scheme_HGU133p2_na22",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version07May/HG-U133_Plus_2.na22.annot.csv",sep="/"))
scheme.hgu133p2.na23 <- import.expr.scheme("Scheme_HGU133p2_na23",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG-U133_Plus_2.na23.annot.csv",sep="/"))
scheme.hgu133p2.na24 <- import.expr.scheme("Scheme_HGU133p2_na24",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version07Nov/HG-U133_Plus_2.na24.annot.csv",sep="/"))
scheme.hgu133p2.na25 <- import.expr.scheme("Scheme_HGU133p2_na25",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version08Mar/HG-U133_Plus_2.na25.annot.csv",sep="/"))
scheme.hgu133p2.na26 <- import.expr.scheme("Scheme_HGU133p2_na26",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version08Jul/HG-U133_Plus_2.na26.annot.csv",sep="/"))
scheme.hgu133p2.na27 <- import.expr.scheme("Scheme_HGU133p2_na27",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version08Nov/HG-U133_Plus_2.na27.annot.csv",sep="/"))
scheme.hgu133p2.na28 <- import.expr.scheme("Scheme_HGU133p2_na28",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version09Mar/HG-U133_Plus_2.na28.annot.csv",sep="/"))
scheme.hgu133p2.na29 <- import.expr.scheme("Scheme_HGU133p2_na29",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version09Jul/HG-U133_Plus_2.na29.annot.csv",sep="/"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HT_HG-U133A
scheme.hthgu133A.na28 <- import.expr.scheme("Scheme_HTHGU133A_na28",filedir=scmdir,paste(libdir,"HT_HG-U133A.cdf",sep="/"),paste(libdir,"HT_HG-U133A.probe.tab",sep="/"),paste(anndir,"Version09Mar/HT_HG-U133A.na28.annot.csv",sep="/"))

# HT_HG-U133_Plus_PM
scheme.hthgu133ppm.na27 <- import.expr.scheme("Scheme_HTHGU133pPM_na27",filedir=scmdir,paste(libdir,"HT_HG-U133_Plus_PM.CDF",sep="/"),paste(libdir,"HT_HG-U133_Plus_PM.probe.tab",sep="/"),paste(anndir,"Version09Feb/HT_HG-U133_Plus_PM.na27.1.annot.csv",sep="/"))
scheme.hthgu133ppm.na28 <- import.expr.scheme("Scheme_HTHGU133pPM_na28",filedir=scmdir,paste(libdir,"HT_HG-U133_Plus_PM.CDF",sep="/"),paste(libdir,"HT_HG-U133_Plus_PM.probe.tab",sep="/"),paste(anndir,"Version09Mar/HT_HG-U133_Plus_PM.na28.annot.csv",sep="/"))

# HT_MG-430_PM
scheme.htmg430pm.na28 <- import.expr.scheme("Scheme_HTMG430PM_na28",filedir=scmdir,paste(libdir,"HT_MG-430_PM.cdf",sep="/"),paste(libdir,"HT_MG-430_PM.probe.tab",sep="/"),paste(anndir,"Version09Mar/HT_MG-430_PM.na28.annot.csv",sep="/"))

# HT_Rat230_PM
scheme.htrat230.na28 <- import.expr.scheme("Scheme_HTRat230PM_na28",filedir=scmdir,paste(libdir,"HT_Rat230_PM.cdf",sep="/"),paste(libdir,"HT_Rat230_PM.probe.tab",sep="/"),paste(anndir,"Version09Mar/HT_Rat230_PM.na28.annot.csv",sep="/"))

# HG-U219
scheme.hgu219.na29 <- import.expr.scheme("Scheme_HGU219_na29",filedir=scmdir,paste(libdir,"HG-U219.cdf",sep="/"),paste(libdir,"HG-U219.probe.tab",sep="/"),paste(anndir,"Version09Oct/HG-U219.na29.annot.csv",sep="/"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for alternative CDF-files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

# HuGene-1_0-st-v1.r3: used as whole genome array
scheme.hugene10stv1r3.na23 <- import.genome.scheme("Scheme_HuGene10stv1r3_na23",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version07Jul/HuGene-1_0-st-v1.na23.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r3.na24 <- import.genome.scheme("Scheme_HuGene10stv1r3_na24",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version07Nov/HuGene-1_0-st-v1.na24.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r3.na25 <- import.genome.scheme("Scheme_HuGene10stv1r3_na25",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Mar/HuGene-1_0-st-v1.na25.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r3.na26 <- import.genome.scheme("Scheme_HuGene10stv1r3_na26",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Jul/HuGene-1_0-st-v1.na26.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r3.na27 <- import.genome.scheme("Scheme_HuGene10stv1r3_na27",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Nov/HuGene-1_0-st-v1.na27.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r3.na28 <- import.genome.scheme("Scheme_HuGene10stv1r3_na28",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version09Mar/HuGene-1_0-st-v1.na28.hg18.transcript.csv",sep="/"))

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1r4.na27 <- import.exon.scheme("Scheme_HuGene10stv1r4_na27",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Feb/HuGene-1_0-st-v1.na27.2.hg18.probeset.csv",sep="/"),
                              paste(anndir,"Version09Feb/HuGene-1_0-st-v1.na27.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r4.na28 <- import.exon.scheme("Scheme_HuGene10stv1r4_na28",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Mar/HuGene-1_0-st-v1.na28.hg18.probeset.csv",sep="/"),
                              paste(anndir,"Version09Mar/HuGene-1_0-st-v1.na28.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r4.na29 <- import.exon.scheme("Scheme_HuGene10stv1r4_na29",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Jul/HuGene-1_0-st-v1.na29.hg18.probeset.csv",sep="/"),
                              paste(anndir,"Version09Jul/HuGene-1_0-st-v1.na29.hg18.transcript.csv",sep="/"))

# note: need to use corrected the probeset annotation file "na30.1"
scheme.hugene10stv1r4.na30 <- import.exon.scheme("Scheme_HuGene10stv1r4_na30_hg19",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Nov/HuGene-1_1-st-v1.na30.1.hg19.probeset.csv",sep="/"),
                              paste(anndir,"Version09Nov/HuGene-1_1-st-v1.na30.hg19.transcript.csv",sep="/"))

# MoGene-1_0-st-v1.r3:
scheme.mogene10stv1r3.na24 <- import.genome.scheme("Scheme_MoGene10stv1r3_na24",filedir=scmdir,
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version07Nov/MoGene-1_0-st-v1.na24.mm8.transcript.csv",sep="/"))

scheme.mogene10stv1r3.na25 <- import.genome.scheme("Scheme_MoGene10stv1r3_na25",filedir=scmdir,
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Mar/MoGene-1_0-st-v1.na25.mm9.transcript.csv",sep="/"))

scheme.mogene10stv1r3.na27 <- import.genome.scheme("Scheme_MoGene10stv1r3_na27",filedir=scmdir,
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Nov/MoGene-1_0-st-v1.na27.mm9.transcript.csv",sep="/"))

# MoGene-1_0-st-v1.r4: used as exon array
# need to fix problem with annotation files first, see:
# https://www.stat.math.ethz.ch/pipermail/bioconductor/2009-August/029049.html
source(paste(.path.package("xps"),"examples/updateAnnotation.R",sep="/"))
updateAnnotation("MoGene-1_0-st-v1.na29.mm9.probeset.csv", "MoGene-1_0-st-v1.na29.mm9.probeset.fixed.csv", probeset="10338063", skip=18, eol="\n")
updateAnnotation("MoGene-1_0-st-v1.na29.mm9.transcript.csv", "MoGene-1_0-st-v1.na29.mm9.transcript.fixed.csv", probeset="10338063", skip=19, eol="\n")

scheme.mogene10stv1r4.na29 <- import.exon.scheme("Scheme_MoGene10stv1r4_na29",filedir=scmdir,
                              paste(libdir,"MoGene-1_0-st-v1.r4.analysis-lib-files/MoGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"MoGene-1_0-st-v1.r4.analysis-lib-files/MoGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Jul/MoGene-1_0-st-v1.na29.mm9.probeset.fixed.csv",sep="/"),
                              paste(anndir,"Version09Jul/MoGene-1_0-st-v1.na29.mm9.transcript.fixed.csv",sep="/"))

# in Dec 2009 Affymetrix has corrected the annotation files for release "na30.1"
scheme.mogene10stv1r4.na30 <- import.exon.scheme("Scheme_MoGene10stv1r4_na30_1",filedir=scmdir,
                              paste(libdir,"MoGene-1_0-st-v1.r4.analysis-lib-files/MoGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"MoGene-1_0-st-v1.r4.analysis-lib-files/MoGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Dec/MoGene-1_1-st-v1.na30.1.mm9.probeset.csv",sep="/"),
                              paste(anndir,"Version09Dec/MoGene-1_1-st-v1.na30.1.mm9.transcript.csv",sep="/"))


# RaGene-1_0-st-v1.r3:
scheme.ragene10stv1r3.na24 <- import.genome.scheme("Scheme_RaGene10stv1r3_na24",filedir=scmdir,
                              paste(libdir,"RaGene-1_0-st-v1.r3.analysis-lib-files/RaGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"RaGene-1_0-st-v1.r3.analysis-lib-files/RaGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version07Nov/RaGene-1_0-st-v1.na24.rn4.transcript.csv",sep="/"))

scheme.ragene10stv1r3.na25 <- import.genome.scheme("Scheme_RaGene10stv1r3_na25",filedir=scmdir,
                              paste(libdir,"RaGene-1_0-st-v1.r3.analysis-lib-files/RaGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"RaGene-1_0-st-v1.r3.analysis-lib-files/RaGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Mar/RaGene-1_0-st-v1.na25.rn4.transcript.csv",sep="/"))

# RaGene-1_0-st-v1.r4: used as exon array
# need to fix problem with annotation files first, see:
# https://www.stat.math.ethz.ch/pipermail/bioconductor/2009-August/029049.html
source(paste(.path.package("xps"),"examples/updateAnnotation.R",sep="/"))
updateAnnotation("RaGene-1_0-st-v1.na29.rn4.probeset.csv", "RaGene-1_0-st-v1.na29.rn4.probeset.fixed.csv", probeset="10700063", skip=18, eol="\n")
updateAnnotation("RaGene-1_0-st-v1.na29.rn4.transcript.csv", "RaGene-1_0-st-v1.na29.rn4.transcript.fixed.csv", probeset="10700063", skip=19, eol="\n")

scheme.ragene10stv1r4.na29 <- import.exon.scheme("Scheme_RaGene10stv1r4_na29",filedir=scmdir,
                              paste(libdir,"RaGene-1_0-st-v1.r4.analysis-lib-files/RaGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"RaGene-1_0-st-v1.r4.analysis-lib-files/RaGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Jul/RaGene-1_0-st-v1.na29.rn4.probeset.fixed.csv",sep="/"),
                              paste(anndir,"Version09Jul/RaGene-1_0-st-v1.na29.rn4.transcript.fixed.csv",sep="/"))

# in Dec 2009 Affymetrix has corrected the annotation files for release "na30.1"
scheme.ragene10stv1r4.na30 <- import.exon.scheme("Scheme_RaGene10stv1r4_na30_1",filedir=scmdir,
                              paste(libdir,"RaGene-1_0-st-v1.r4.analysis-lib-files/RaGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"RaGene-1_0-st-v1.r4.analysis-lib-files/RaGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Dec/RaGene-1_1-st-v1.na30.1.rn4.probeset.csv",sep="/"),
                              paste(anndir,"Version09Dec/RaGene-1_1-st-v1.na30.1.rn4.transcript.csv",sep="/"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

scheme.huex10stv2r2.na24 <- import.exon.scheme("Scheme_HuEx10stv2r2_na24",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version07Nov/HuEx-1_0-st-v2.na24.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version07Nov/HuEx-1_0-st-v2.na24.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na25 <- import.exon.scheme("Scheme_HuEx10stv2r2_na25",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version08Mar/HuEx-1_0-st-v2.na25.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version08Mar/HuEx-1_0-st-v2.na25.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na26 <- import.exon.scheme("Scheme_HuEx10stv2r2_na26",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version08Jul/HuEx-1_0-st-v2.na26.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version08Jul/HuEx-1_0-st-v2.na26.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na27 <- import.exon.scheme("Scheme_HuEx10stv2r2_na27",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version08Nov/HuEx-1_0-st-v2.na27.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version08Nov/HuEx-1_0-st-v2.na27.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na28 <- import.exon.scheme("Scheme_HuEx10stv2r2_na28",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version09Mar/HuEx-1_0-st-v2.na28.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version09Mar/HuEx-1_0-st-v2.na28.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na29 <- import.exon.scheme("Scheme_HuEx10stv2r2_na29",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version09Jul/HuEx-1_0-st-v2.na29.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version09Jul/HuEx-1_0-st-v2.na29.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na30 <- import.exon.scheme("Scheme_HuEx10stv2r2_na30",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version09Jul/HuEx-1_0-st-v2.na30.hg19.probeset.csv",sep="/"),
                            paste(anndir,"Version09Jul/HuEx-1_0-st-v2.na30.hg19.transcript.csv",sep="/"))

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

scheme.moex10stv1r2.na24 <- import.exon.scheme("Scheme_MoEx10stv2r2_na24",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07Nov/MoEx-1_0-st-v1.na24.mm8.probeset.csv",sep="/"),
                            paste(anndir,"Version07Nov/MoEx-1_0-st-v1.na24.mm8.transcript.csv",sep="/"))

scheme.moex10stv1r2.na25 <- import.exon.scheme("Scheme_MoEx10stv2r2_na25",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version08Mar/MoEx-1_0-st-v1.na25.mm9.probeset.csv",sep="/"),
                            paste(anndir,"Version08Mar/MoEx-1_0-st-v1.na25.mm9.transcript.csv",sep="/"))

scheme.moex10stv1r2.na30 <- import.exon.scheme("Scheme_MoEx10stv2r2_na30",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version09Nov/MoEx-1_0-st-v1.na30.mm9.probeset.csv",sep="/"),
                            paste(anndir,"Version09Nov/MoEx-1_0-st-v1.na30.mm9.transcript.csv",sep="/"))

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

scheme.raex10stv1r2.na24 <- import.exon.scheme("Scheme_RaEx10stv2r2_na24",filedir=scmdir,
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07Nov/RaEx-1_0-st-v1.na24.rn4.probeset.csv",sep="/"),
                            paste(anndir,"Version07Nov/RaEx-1_0-st-v1.na24.rn4.transcript.csv",sep="/"))

scheme.raex10stv1r2.na25 <- import.exon.scheme("Scheme_RaEx10stv2r2_na25",filedir=scmdir,
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version08Mar/RaEx-1_0-st-v1.na25.rn4.probeset.csv",sep="/"),
                            paste(anndir,"Version08Mar/RaEx-1_0-st-v1.na25.rn4.transcript.csv",sep="/"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for other Affymetrix arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Yeast_2:
scheme.yeast2.na28 <- import.expr.scheme("Scheme_Yeast2_na27",filedir=scmdir,
                      paste(libdir,"Yeast_2.cdf",sep="/"),
                      paste(libdir,"Yeast_2.probe.tab",sep="/"),
                      paste(anndir,"Version09Mar/Yeast_2.na28.annot.csv",sep="/"))

# Pae_G1a:
scheme.paeg1a.na28 <- import.expr.scheme("Scheme_PaeG1a_na28",filedir=scmdir,
                      paste(libdir,"Pae_G1a.CDF",sep="/"),
                      paste(libdir,"Pae_G1a.probe.tab",sep="/"), 
                      paste(anndir,"Version09Mar/Pae_G1a.na28.annot.csv",sep="/"))

# Mouse430_2:
scheme.mouse430.2.na28 <- import.expr.scheme("Scheme_Mouse430_2_na28",filedir=scmdir,
                          paste(libdir,"Mouse430_2.cdf",sep="/"),
                          paste(libdir,"Mouse430_2.probe.tab",sep="/"),
                          paste(anndir,"Version09Mar/Mouse430_2.na28.annot.csv",sep="/"))

# miRNA-1_0:
# note: you need to rename "miRNA-1_0.probe_list.20081203.txt" to "miRNA-1_0.probe.tab"
scheme.mirna <- import.expr.scheme("Scheme_miRNA_1_0",filedir=scmdir,
                paste(libdir,"miRNA-1_0.CDF",sep="/"),
                paste(libdir,"miRNA-1_0.probe.tab",sep="/"))



#------------------------------------------------------------------------------#
