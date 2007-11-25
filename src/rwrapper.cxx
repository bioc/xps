/*
 *******************************************************************************
 *
 * File: rwrapper.cpp
 *
 * Implementation by: Christian Stratowa
 *
 * Copyright (C) Christian Stratowa 2002-2007
 *
 * A wrapper for the XPS libraries
 *
 *******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TKey.h"
#include "TString.h"
#include "TSystem.h"

#include "rwrapper.h"

#include "XPSBase.h"
#include "XPSSchemes.h"
#include "XPSDataTypes.h"
#include "XPSData.h"
#include "XPSPreProcessing.h"
#include "XPSNormation.h"
//#include "XPSAnalysis.h"


/*//////////////////////////////////////////////////////////////////////////////
//                                                                            //
// XSchemeManager: import/export chip definitions and chip annotations        //
//                                                                            //
//////////////////////////////////////////////////////////////////////////////*/

/*____________________________________________________________________________*/
void ImportExprSchemes(char **filename, char **dirname, char **chipname,
                       char **scheme, char **probe, char **annot,
                       int *verbose, int *err)
{
// Import Affymetrix chip definition file *.CDF files into XPS

// create new scheme manager
   XSchemeManager *manager = new XSchemeManager("SchemeManager","", *verbose);

// create new root schemes file
   int r = 0;
   r += manager->New(filename[0], dirname[0], "GeneChip", "Schemes");

// store chip definitions,  probe sequences and annotations
   r += manager->NewScheme(chipname[0], scheme[0], "scheme");
   r += manager->NewScheme(chipname[0], probe[0], "probe");
   if (strcmp(annot[0], "") != 0) {
      r += manager->NewAnnotation(chipname[0], annot[0]);
   }//if
   *err = (int)r;

// cleanup
   manager->Close();
   delete manager;
}//ImportExprSchemes

/*____________________________________________________________________________*/
void ImportExonSchemes(char **filename, char **dirname, char **chipname,
                       char **layout, char **scheme, char **probeset,
                       char **transcript, char **control,
                       int *verbose, int *err)
{
// Import Affymetrix chip definition file *.CDF files into XPS

// create new scheme manager
   XSchemeManager *manager = new XSchemeManager("SchemeManager","", *verbose);

// create new root schemes file
   int r = 0;
   r += manager->New(filename[0], dirname[0], "ExonChip", "Schemes");

// store chip definitions,  probe sequences and annotations
   r += manager->NewScheme(chipname[0], layout[0], "layout");
   if (strcmp(control[0], "") != 0) {
      r += manager->NewAnnotation(chipname[0], control[0],"control");
   }//if
   r += manager->NewAnnotation(chipname[0], probeset[0],"probeset");
   r += manager->NewScheme(chipname[0], scheme[0], "scheme");
   r += manager->NewAnnotation(chipname[0], transcript[0],"transcript");
   *err = (int)r;

// cleanup
   manager->Close();
   delete manager;
}//ImportExonSchemes

/*____________________________________________________________________________*/
void ImportGenomeSchemes(char **filename, char **dirname, char **chipname,
                         char **layout, char **scheme, char **transcript,
                         int *verbose, int *err)
{
// Import Affymetrix chip definition file *.CDF files into XPS

// create new scheme manager
   XSchemeManager *manager = new XSchemeManager("SchemeManager","", *verbose);

// create new root schemes file
   int r = 0;
   r += manager->New(filename[0], dirname[0], "GenomeChip", "Schemes");

// store chip definitions,  probe sequences and annotations
   r += manager->NewScheme(chipname[0], layout[0], "layout");
   r += manager->NewAnnotation(chipname[0], transcript[0],"transcript");
   r += manager->NewScheme(chipname[0], scheme[0], "scheme");
   *err = (int)r;

// cleanup
   manager->Close();
   delete manager;
}//ImportGenomeSchemes

/*____________________________________________________________________________*/
void ExportScheme(char **schemefile, char **treename, char **varlist,
                  char **outfile, char **sep, int *verbose)
{
// Export data from scheme tree

// create new scheme manager
   XSchemeManager *manager = new XSchemeManager("SchemeManager","", *verbose);

// open root schemes file
   manager->Open(schemefile[0]);

// export tree
   manager->Export(treename[0], varlist[0], outfile[0], sep[0]);

// cleanup
   manager->Close();
   delete manager;
}//ExportScheme


/*//////////////////////////////////////////////////////////////////////////////
//                                                                            //
// XDataManager: import/export raw chip data, i.e. *.CEL files                //
//                                                                            //
//////////////////////////////////////////////////////////////////////////////*/

/*____________________________________________________________________________*/
void ImportData(char **filename, char **dirname, char **chiptype,
                char **chipname, char **schemefile, char **treeset,
                char **celfiles, char **celnames, int *numdata,
                char **project, int *nproject, char **author, int *nauthor,
                char **dataset, int *ndataset, char **source, int *nsource,
                char **sample, int *nsample, char **cell, int *ncell,
                char **pcell, int *npcell, char **tissue, int *ntissue,
                char **biopsy, int *nbiopsy, char **array, int *narray,
                char **hyb, int *nhyb, char **treat, int *ntreat, 
                int *replace, int *update, int *verbose, int *err)
{
// Import Affymetrix *.CEL files into XPS

// create new data manager
   XDataManager *manager = new XDataManager("DataManager","", *verbose);

// initialize chip type and variable list
   int r = 0;
   r += manager->Initialize(chiptype[0]);
   r += manager->InitInput(chipname[0], "cel", "MEAN/D:STDV/D:NPIXELS/I", "RawData");

// open root scheme file
   r += manager->OpenSchemes(schemefile[0]);

// create/update root data file 
   if (*update == 1) {
      r += manager->Update(filename[0], "RawData", "R");
      // need to force to delete root file from memory to prevent R using former root file!
      manager->SetFileOwner(kTRUE);
   } else {
      r += manager->New(filename[0], dirname[0], chiptype[0], "RawData");
   }//if

// optional: add database info
   int nh = *nhyb/9;    //number of hybridization infos
   int nt = *ntreat/8;  //number of treatment infos
   int db = *nproject + *nauthor + *ndataset + *nsource + *nsample 
          + *ncell + *npcell + *ntissue + *nbiopsy + *narray + nh + nt;
   if (db > 0) {
      manager->BeginTransaction();
   }//if
   if (*nproject == 5) {
      manager->ProjectInfo(project[0], atol(project[1]), project[2],
                           project[3], project[4], *replace);
   }//if
   if (*nauthor == 8) {
      manager->AuthorInfo(author[0], author[1], author[2], author[3], author[4],
                          author[5], author[6], author[7], *replace);
   }//if
   if (*ndataset == 7) {
      manager->DatasetInfo(dataset[0], dataset[1], dataset[2], dataset[3],
                           atol(dataset[4]), dataset[5], dataset[6], *replace);
   }//if
   if (*nsource == 6) {
      manager->SourceInfo(source[0], source[1], source[2], source[3],
                          source[4], source[5], *replace);
   }//if
   if (*nsample == 12) {
      manager->SampleInfo(sample[0], sample[1], sample[2], sample[3],
                          sample[4], sample[5], atoi(sample[6]), sample[7],
                          sample[8], atof(sample[9]), sample[10], sample[11],
                          *replace);
   }//if
   if (*ncell == 15) {
      manager->CellLineInfo(cell[0], cell[1], cell[2], cell[3], cell[4],
                            cell[5], cell[6], cell[7], cell[8], atoi(cell[9]),
                            cell[10], cell[11], atof(cell[12]), cell[13],
                            cell[14], *replace);
   }//if
   if (*npcell == 14) {
      manager->PrimaryCellInfo(pcell[0], pcell[1], atol(pcell[2]), pcell[3], 
                               pcell[4], pcell[5], pcell[6], pcell[7],
                               atoi(pcell[8]), pcell[9], pcell[10],
                               atof(pcell[11]), pcell[12], pcell[13], *replace);
   }//if
   if (*ntissue == 19) {
      manager->TissueInfo(tissue[0], tissue[1], tissue[2], tissue[3], tissue[4],
                          tissue[5], atof(tissue[6]), tissue[7], tissue[8],
                          tissue[9], tissue[10], tissue[11], tissue[12],
                          atoi(tissue[13]), tissue[14], tissue[15],
                          atof(tissue[16]), tissue[17], tissue[18], *replace);
   }//if
   if (*nbiopsy == 19) {
      manager->BiopsyInfo(biopsy[0], biopsy[1], biopsy[2], biopsy[3], biopsy[4],
                          atof(biopsy[5]), biopsy[6], biopsy[7], biopsy[8],
                          biopsy[9], biopsy[10], biopsy[11], atoi(biopsy[12]),
                          biopsy[13], biopsy[14], atof(biopsy[15]),
                          biopsy[16], biopsy[17], *replace);
   }//if
   if (*narray == 4) {
      manager->ArrayInfo(array[0], array[1], array[2], array[3], *replace);
   }//if
   for (int i=0; i<nh; i++) {
      manager->HybridizationInfo(hyb[0+9*i], hyb[1+9*i], hyb[2+9*i],
                                 atol(hyb[3+9*i]), hyb[4+9*i], hyb[5+9*i],
                                 hyb[6+9*i], atoi(hyb[7+9*i]), hyb[8+9*i],
                                 *replace);
   }//for_i
   for (int i=0; i<nt; i++) {
      manager->TreatmentInfo(treat[0+8*i], treat[1+8*i], atof(treat[2+8*i]),
                             treat[3+8*i], atof(treat[4+8*i]), treat[5+8*i],
                             treat[6+8*i], treat[7+8*i], *replace);
   }//for_i
   if (db > 0) {
      manager->CommitTransaction();
   }//if

// store *.CEL data as trees in data file
   for (int i=0; i<*numdata; i++) {
      if (strcmp(celnames[0], "") == 0) {
         r += manager->Import(treeset[0], celfiles[i], "");
      } else {
         r += manager->Import(treeset[0], celfiles[i], celnames[i]);
      }//if
   }//for_i
   *err = (int)r;

// cleanup
   manager->Close();
   delete manager;
}//ImportData

/*____________________________________________________________________________*/
void ExportData(char **filename, char **schemefile, char **chiptype,
                char **datatype, char **treenames, int *ntrees, char **exten,
                char **varlist, char **outfile, char **sep,
                int *verbose, int *err)
{
// Export data from scheme tree

// create new manager
   int r = 0;
   XManager *manager = 0;
   if (strcmp(datatype[0], "rawdata") == 0) {
      manager = new XDataManager("DataManager","", *verbose);

      r += manager->Initialize(chiptype[0], "rawdata", "R");
      r += ((XDataManager*)manager)->OpenSchemes(schemefile[0]);
   } else if (strcmp(datatype[0], "preprocess") == 0) {
      manager = new XPreProcessManager("PreProcessManager","", *verbose);

      r += manager->Initialize(chiptype[0], "preprocess", "R");
      r += ((XPreProcessManager*)manager)->OpenSchemes(schemefile[0]);
   } else if (strcmp(datatype[0], "normation") == 0) {
      manager = new XNormationManager("NormationManager","", *verbose);

      r += manager->Initialize(chiptype[0], "normation", "R");
      r += ((XNormationManager*)manager)->OpenSchemes(schemefile[0]);
   }//if

// open root data file
   r += manager->Open(filename[0]);

// export tree
   if (*ntrees == 1) {
      r += manager->Export(treenames[0], varlist[0], outfile[0], sep[0]);
   } else if (*ntrees > 1) {
      for (int i=0; i<*ntrees; i++) {
         r += manager->AddTree("ExportSet", treenames[i]);
      }//for_i

      r += manager->ExportSet("ExportSet", exten[0], varlist[0], outfile[0], sep[0]);
   }//if
   *err = (int)r;

// cleanup
   manager->Close();
   delete manager;
}//ExportData


/*////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// XPreProcessManager: preprocess data, e.g. RMA, MAS5, DABG                    //
//                                                                              //
////////////////////////////////////////////////////////////////////////////////*/

/*____________________________________________________________________________*/
void PreprocessRMA(char **filename, char **dirname, char **chipname,
                   char **chiptype, char **schemefile, char **tmpdir,
                   char **bgrdoption, char **exproption, char **treeset,
                   char **treenames, int *ntrees, int *normalize,
                   int *level, int *verbose, int *err)
{
// Preprocess trees using RMA

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager","", *verbose);

// initialize chip type
   int r = 0;
   r += manager->Initialize(chiptype[0]);

// temporary files for background, normalization, expression algorithms
   char *bgrdfile = new char[strlen(tmpdir[0]) + 22];
   if (strcmp(tmpdir[0], "") != 0) {
      bgrdfile = strcpy(bgrdfile, tmpdir[0]);
      bgrdfile = strcat(bgrdfile, "/tmp_bgrd_310151.root");
   } else {
      bgrdfile = "";
   }//if

   char *normfile = new char[strlen(tmpdir[0]) + 14];
   if (strcmp(tmpdir[0], "") != 0) {
      normfile = strcpy(normfile, tmpdir[0]);
      normfile = strcat(normfile, "/tmp_rkq.root");
   } else {
      normfile = "";
   }//if

   char *exprfile = new char[strlen(tmpdir[0]) + 22];
   if (strcmp(tmpdir[0], "") != 0) {
      exprfile = strcpy(exprfile, tmpdir[0]);
      exprfile = strcat(exprfile, "/tmp_expr_310151.root");
   } else {
      exprfile = "";
   }//if

// initialize backgrounder
   char *bgrdopt = new char[strlen(bgrdoption[0]) + 14];
   if (strcmp(bgrdoption[0], "pmonly") == 0 || strcmp(bgrdoption[0], "mmonly") == 0) {
      bgrdopt = strcpy(bgrdopt, bgrdoption[0]);
      bgrdopt = strcat(bgrdopt, ":epanechnikov");

      r += manager->InitAlgorithm("selector", "probe", bgrdoption[0], 0);
      r += manager->InitAlgorithm("backgrounder", "rma", bgrdopt, bgrdfile, 1, 16384);
   } else if (strcmp(bgrdoption[0], "genomic") == 0 || strcmp(bgrdoption[0], "antigenomic") == 0) {
//??      int bgrdlevel = 16383;
      int bgrdlevel = *level; //use only those probes which are also used for expression
      int bgrdtype  = (strcmp(bgrdoption[0], "genomic") == 0) ? -1 : -2;

      r += manager->InitAlgorithm("selector", "probe", "exon", 0, 2, bgrdlevel, bgrdtype);
      r += manager->InitAlgorithm("backgrounder", "rma", "pmonly:epanechnikov", bgrdfile, 1, 16384);
   }//if

// initialize normalizer
   if (normalize) {
      if (strcmp(chiptype[0], "GeneChip") == 0) {
         r += manager->InitAlgorithm("selector", "probe", "pmonly", 0);
      } else if (strcmp(chiptype[0], "GenomeChip") == 0) {
//??         int normlevel = 16383; 
         int normlevel = *level; //use only those probes which are also used for expression
         r += manager->InitAlgorithm("selector", "probe", "genome", 0, 1, normlevel);
      } else if (strcmp(chiptype[0], "ExonChip") == 0) {
//??         int normlevel = 16383; 
         int normlevel = *level; //use only those probes which are also used for expression
         r += manager->InitAlgorithm("selector", "probe", "exon", 0, 1, normlevel);
      }//if

      r += manager->InitAlgorithm("normalizer", "quantile", "together:none:0", normfile, 1, 0.0);
   }//if

// initialize expressor
   if (strcmp(chiptype[0], "GeneChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", "pmonly", 0);
   } else if (strcmp(chiptype[0], "GenomeChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", "genome", 0, 1, *level);
   } else if (strcmp(chiptype[0], "ExonChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", "exon", 0, 1, *level);
   }//if

   char *expropt = new char[strlen(exproption[0]) + 6];
   expropt = strcpy(expropt, exproption[0]);
   expropt = strcat(expropt, ":log2");
   manager->InitAlgorithm("expressor", "medianpolish", expropt, exprfile, 3, 10, 0.01, 1.0);

// create new root data file 
   r += manager->New(filename[0], dirname[0], chiptype[0], "preprocess");

// open root scheme file
   r += manager->OpenSchemes(schemefile[0], chipname[0]);

// add trees for rma
   for (int i=0; i<*ntrees; i++) {
      r += manager->AddTree(treeset[0], treenames[i]);
   }//for_i

// preprocess expression values and store as trees in new file
   r += manager->Preprocess(treeset[0], "preprocess");
   *err = (int)r;

// cleanup
   delete [] expropt;
   delete [] bgrdopt;
   if (strcmp(tmpdir[0], "") != 0) {
      delete [] exprfile;
      delete [] normfile;
      delete [] bgrdfile;
   }//if

   manager->Close();
   delete manager;
}//PreprocessRMA

/*____________________________________________________________________________*/
void PreprocessMAS4(char **filename, char **dirname, char **chipname,
                    char **chiptype, char **schemefile, char **tmpdir,
                    char **exproption, char **treeset, char **treenames,
                    int *ntrees, int *level, int *verbose, int *err)
{
// Preprocess trees using MAS4

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager","", *verbose);

// initialize chip type
   int r = 0;
   r += manager->Initialize(chiptype[0]);

// temporary files for background and expression algorithms
   char *tmpfile = "";
   if (strcmp(tmpdir[0], "") != 0) {
      tmpfile = new char[strlen(tmpdir[0]) + 22];
      tmpfile = strcpy(tmpfile, tmpdir[0]);
      tmpfile = strcat(tmpfile, "/tmp_mas4_310151.root");
   }//if

// initialize algorithms
   char *expropt = 0;
   if (strcmp(chiptype[0], "GeneChip") == 0) {
   // initialize backgrounder
      r += manager->InitAlgorithm("selector","probe","all", 0, 0);
      r += manager->InitAlgorithm("backgrounder","sector","subtractbg",0, 4,0.02,4,4,0);

   // initialize expressor
      r += manager->InitAlgorithm("selector","default","none", 0, 0);
      r += manager->InitAlgorithm("expressor","avgdiff","0",tmpfile,1,3.0);
   } else if ((strcmp(chiptype[0], "GenomeChip") == 0) ||
              (strcmp(chiptype[0], "ExonChip")   == 0)) {
   // initialize backgrounder
      r += manager->InitAlgorithm("selector","probe","exon",0, 1, *level);
      r += manager->InitAlgorithm("backgrounder","sector","subtractbg",0, 4,0.02,4,4,0);

   // initialize expressor
      expropt = new char[strlen(exproption[0]) + 6];
      expropt = strcpy(expropt, exproption[0]);
      expropt = strcat(expropt, ":0");

      r += manager->InitAlgorithm("selector","probe","exon",0, 2, *level, 1); //pm=level; mm=affx
      r += manager->InitAlgorithm("expressor","avgdiff",expropt,tmpfile,1,3.0);
   }//if

// create new root data file 
   r += manager->New(filename[0], dirname[0], chiptype[0], "preprocess");

// open root scheme file
   r += manager->OpenSchemes(schemefile[0], chipname[0]);

// add trees for rma
   for (int i=0; i<*ntrees; i++) {
      r += manager->AddTree(treeset[0], treenames[i]);
   }//for_i

// preprocess expression values and store as trees in new file
   r += manager->Preprocess(treeset[0], "preprocess");
   *err = (int)r;

// cleanup
   if (expropt) delete [] expropt;
   if (tmpfile && strcmp(tmpfile, "") != 0) delete [] tmpfile;

   manager->Close();
   delete manager;
}//PreprocessMAS4

/*____________________________________________________________________________*/
void PreprocessMAS5(char **filename, char **dirname, char **chipname,
                    char **chiptype, char **schemefile, char **tmpdir,
                    char **exproption, char **treeset, char **treenames,
                    int *ntrees, int *level, int *verbose, int *err)
{
// Preprocess trees using MAS5

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager","", *verbose);

// initialize chip type
   int r = 0;
   r += manager->Initialize(chiptype[0]);

// temporary files for background and expression algorithms
   char *tmpfile = "";
   if (strcmp(tmpdir[0], "") != 0) {
      tmpfile = new char[strlen(tmpdir[0]) + 22];
      tmpfile = strcpy(tmpfile, tmpdir[0]);
      tmpfile = strcat(tmpfile, "/tmp_mas5_310151.root");
   }//if

// initialize algorithms
   char *expropt = 0;
   if (strcmp(chiptype[0], "GeneChip") == 0) {
   // initialize backgrounder
      r += manager->InitAlgorithm("selector","probe","both", 0, 0);
      r += manager->InitAlgorithm("backgrounder","weightedsector","correctbg", tmpfile, 6,0.02,4,4,0,100,0.5);

   // initialize expressor
      r += manager->InitAlgorithm("selector","probe","none", 0, 0);
      r += manager->InitAlgorithm("expressor","TukeyBiweight","log2",tmpfile,7,0.03,10.0,2.0e-20,5.0,0.0001,1.0,0.5);
   } else if ((strcmp(chiptype[0], "GenomeChip") == 0) ||
              (strcmp(chiptype[0], "ExonChip")   == 0)) {
   // initialize backgrounder
      r += manager->InitAlgorithm("selector","probe","exon",0, 1, *level);
      r += manager->InitAlgorithm("backgrounder","weightedsector","correctbg",tmpfile, 6,0.02,4,4,0,100,0.5);

   // initialize expressor
      expropt = new char[strlen(exproption[0]) + 6];
      expropt = strcpy(expropt, exproption[0]);
      expropt = strcat(expropt, ":log2");

      r += manager->InitAlgorithm("selector","probe","exon",0, 1, *level);
      r += manager->InitAlgorithm("expressor","TukeyBiweight",expropt,tmpfile,7,0.03,10.0,2.0e-20,5.0,0.0001,1.0,0.5);
   }//if

// create new root data file 
   r += manager->New(filename[0], dirname[0], chiptype[0], "preprocess");

// open root scheme file
   r += manager->OpenSchemes(schemefile[0], chipname[0]);

// add trees for rma
   for (int i=0; i<*ntrees; i++) {
      r += manager->AddTree(treeset[0], treenames[i]);
   }//for_i

// preprocess expression values and store as trees in new file
   r += manager->Preprocess(treeset[0], "preprocess");
   *err = (int)r;

// cleanup
   if (expropt) delete [] expropt;
   if (tmpfile && strcmp(tmpfile, "") != 0) delete [] tmpfile;

   manager->Close();
   delete manager;
}//PreprocessMAS5

/*____________________________________________________________________________*/
void PreprocessMAS5Call(char **filename, char **dirname, char **chipname,
                        char **chiptype, char **schemefile, char **tmpdir,
                        char **calloption, char **treeset, char **treenames,
                        int *ntrees, double *tau, double *alpha1, double *alpha2,
                        int *ignore, int *level, int *verbose, int *err)
{
// Preprocess trees using MAS5 call

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager","", *verbose);

// initialize chip type
   int r = 0;
   r += manager->Initialize(chiptype[0]);

// temporary file for background
   char *bgrdfile = "";
   if (strcmp(tmpdir[0], "") != 0) {
      bgrdfile = new char[strlen(tmpdir[0]) + 22];
      bgrdfile = strcpy(bgrdfile, tmpdir[0]);
      bgrdfile = strcat(bgrdfile, "/tmp_call_310151.root");
   }//if

// initialize algorithms
   char *callopt = 0;
   if (strcmp(chiptype[0], "GeneChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", "none", 0, 0);
      r += manager->InitAlgorithm("calldetector", "dc5", "raw", 0, 6, *tau, *alpha1, *alpha2, *ignore, 0, 0); //tau,alpha1,alpha2,ignore,exact,correct
   } else if (strcmp(chiptype[0], "GenomeChip") == 0) {
      callopt = new char[strlen(calloption[0]) + 9];
      callopt = strcpy(callopt, calloption[0]);
      callopt = strcat(callopt, ":adjusted");

      // initialize backgrounder
      r += manager->InitAlgorithm("selector", "probe", "genome", 0, 1, *level);
      // the following setting calculates bg but does not subtract bg from intensity:
      r += manager->InitAlgorithm("backgrounder", "weightedsector", "correctbg", bgrdfile, 6, 0.02, 4, 4, 0, 100, 0.5);
//?      r += manager->InitAlgorithm("backgrounder", "weightedsector", "none", bgrdfile, 6, 0.02, 4, 4, 0, 100, -1.0);

      // initialize call detector
      r += manager->InitAlgorithm("selector", "probe", "genome", 0, 2, *level, -2); //??mm=antigenommic??
      r += manager->InitAlgorithm("calldetector", "dc5", callopt, 0, 6, *tau, *alpha1, *alpha2, *ignore, 0, 0);
   } else if (strcmp(chiptype[0], "ExonChip") == 0) {
      callopt = new char[strlen(calloption[0]) + 9];
      callopt = strcpy(callopt, calloption[0]);
      callopt = strcat(callopt, ":adjusted");

      // initialize backgrounder
      r += manager->InitAlgorithm("selector", "probe", "exon", 0, 1, *level);
      // the following setting calculates bg but does not subtract bg from intensity:
      r += manager->InitAlgorithm("backgrounder", "weightedsector", "none", bgrdfile, 6, 0.02, 4, 4, 0, 100, -1.0);

      // initialize call detector
      r += manager->InitAlgorithm("selector", "probe", "exon", 0, 2, *level, -2); //??mm=antigenommic??
      r += manager->InitAlgorithm("calldetector", "dc5", callopt, 0, 6, *tau, *alpha1, *alpha2, *ignore, 0, 0);
   }//if

// create new root data file 
   r += manager->New(filename[0], dirname[0], chiptype[0], "preprocess");

// open root scheme file
   r += manager->OpenSchemes(schemefile[0], chipname[0]);

// add trees for rma
   for (int i=0; i<*ntrees; i++) {
      r += manager->AddTree(treeset[0], treenames[i]);
   }//for_i

// preprocess expression values and store as trees in new file
   r += manager->Preprocess(treeset[0], "preprocess");
   *err = (int)r;

// cleanup
   if (callopt  && strcmp(callopt,  "") != 0) delete [] callopt;
   if (bgrdfile && strcmp(bgrdfile, "") != 0) delete [] bgrdfile;

   manager->Close();
   delete manager;
}//PreprocessMAS5Call

/*____________________________________________________________________________*/
void PreprocessDABGCall(char **filename, char **dirname, char **chipname,
                        char **chiptype, char **schemefile, char **calloption,
                        char **treeset, char **treenames, int *ntrees,
                        double *alpha1, double *alpha2, int *level,
                        int *verbose, int *err)
{
// Preprocess trees using DABG call

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager","", *verbose);

// initialize chip type
   int r = 0;
   r += manager->Initialize(chiptype[0]);

// initialize algorithms
   char *callopt = 0;
   if (strcmp(chiptype[0], "GeneChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", "none", 0, 0);
      r += manager->InitAlgorithm("calldetector", "dab", "raw", 0, 3, 2.0, *alpha1, *alpha2); //cut, alpha1, alpha2
   } else if ((strcmp(chiptype[0], "GenomeChip") == 0) ||
              (strcmp(chiptype[0], "ExonChip")   == 0)) {
      callopt = new char[strlen(calloption[0]) + 5];
      callopt = strcpy(callopt, calloption[0]);
      callopt = strcat(callopt, ":raw");

      r += manager->InitAlgorithm("selector","probe","exon",0, 2, *level, -2); //??mm=antigenommic??
      r += manager->InitAlgorithm("calldetector", "dab", callopt, 0, 3, 2.0, *alpha1, *alpha2); //cut, alpha1, alpha2
   }//if

// create new root data file 
   r += manager->New(filename[0], dirname[0], chiptype[0], "preprocess");

// open root scheme file
   r += manager->OpenSchemes(schemefile[0], chipname[0]);

// add trees for rma
   for (int i=0; i<*ntrees; i++) {
      r += manager->AddTree(treeset[0], treenames[i]);
   }//for_i

// preprocess expression values and store as trees in new file
   r += manager->Preprocess(treeset[0], "DetectCall");
   *err = (int)r;

// cleanup
   if (callopt) delete [] callopt;

   manager->Close();
   delete manager;
}//PreprocessDABGCall

/*____________________________________________________________________________*/
void Preprocess(char **filename, char **dirname, char **chipname, char **chiptype,
                char **schemefile, char **tmpdir, int *update,
                char **bgrdtype, char **bgrdselection, char **bgrdoption,
                int *nbgrdpar, double *bgrdpars,
                char **normtype, char **normselection, char **normoption,
                int *nnormpar, double *normpars,
                char **exprtype, char **exprselection, char **exproption,
                int *nexprpar, double *exprpars,
                char **reftree, char **refmethod, double *refparam,
                char **treeset, char **treenames, int *ntrees,
                int *level, int *verbose, int *err)
{
// Preprocess trees
   int    r = 0;
   int    n;
   double p0, p1, p2, p3, p4, p5, p6, p7;

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager","", *verbose);

// initialize chip type
   r += manager->Initialize(chiptype[0]);

// temporary files for background, normalization, expression algorithms
   char *bgrdfile = new char[strlen(tmpdir[0]) + 22];
   if (strcmp(tmpdir[0], "") != 0) {
      bgrdfile = strcpy(bgrdfile, tmpdir[0]);
      bgrdfile = strcat(bgrdfile, "/tmp_bgrd_310151.root");
   } else {
      bgrdfile = "";
   }//if

   char *normfile = new char[strlen(tmpdir[0]) + 14];
   if (strcmp(tmpdir[0], "") != 0) {
      normfile = strcpy(normfile, tmpdir[0]);
      normfile = strcat(normfile, "/tmp_rkq.root");
   } else {
      normfile = "";
   }//if

   char *exprfile = new char[strlen(tmpdir[0]) + 22];
   if (strcmp(tmpdir[0], "") != 0) {
      exprfile = strcpy(exprfile, tmpdir[0]);
      exprfile = strcat(exprfile, "/tmp_expr_310151.root");
   } else {
      exprfile = "";
   }//if

// initialize backgrounder
   if (strcmp(bgrdtype[0], "none") != 0) {
      n  = *nbgrdpar;
      p0 = p1 = p2 = p3 = p4 = p5 = p6 = p7 = 0.0;
      if (n > 0) p0 = bgrdpars[0];
      if (n > 1) p1 = bgrdpars[1];
      if (n > 2) p2 = bgrdpars[2];
      if (n > 3) p3 = bgrdpars[3];
      if (n > 4) p4 = bgrdpars[4];
      if (n > 5) p5 = bgrdpars[5];
      if (n > 6) p6 = bgrdpars[6];
      if (n > 7) p7 = bgrdpars[7];

      if (strcmp(bgrdselection[0], "pmonly") == 0 ||
          strcmp(bgrdselection[0], "mmonly") == 0 ||
          strcmp(bgrdselection[0], "both")   == 0 ||
          strcmp(bgrdselection[0], "all")    == 0 ||
          strcmp(bgrdselection[0], "none")   == 0) {
         r += manager->InitAlgorithm("selector", "probe", bgrdselection[0], 0, 0);
         r += manager->InitAlgorithm("backgrounder", bgrdtype[0], bgrdoption[0],
                                     bgrdfile, n, p0, p1, p2, p3, p4, p5, p6, p7);
      } else if (strcmp(bgrdselection[0], "genomic")     == 0 ||
                 strcmp(bgrdselection[0], "antigenomic") == 0) {
//??         int bgrdlevel = 16383;
         int bgrdlevel = *level; //use only those probes which are also used for expression
         int bgrdopt   = (strcmp(bgrdselection[0], "genomic") == 0) ? -1 : -2;

         r += manager->InitAlgorithm("selector", "probe", "exon", 0, 2,
                                     bgrdlevel, bgrdopt);
         r += manager->InitAlgorithm("backgrounder", bgrdtype[0], bgrdoption[0],
                                     bgrdfile, n, p0, p1, p2, p3, p4, p5, p6, p7);
      }//if
   }//if

// initialize normalizer
   if (strcmp(normtype[0], "none") != 0) {
      n  = *nnormpar;
      p0 = p1 = p2 = 0.0;
      if (n > 0) p0 = normpars[0];
      if (n > 1) p1 = normpars[1];
      if (n > 2) p2 = normpars[2];

      if (strcmp(chiptype[0], "GeneChip") == 0) {
         r += manager->InitAlgorithm("selector", "probe", normselection[0], 0);
      } else if (strcmp(chiptype[0], "GenomeChip") == 0) {
         r += manager->InitAlgorithm("selector", "probe", "genome", 0, 1, *level);
      } else if (strcmp(chiptype[0], "ExonChip") == 0) {
         r += manager->InitAlgorithm("selector", "probe", "exon", 0, 1, *level);
      }//if

      r += manager->InitAlgorithm("normalizer", normtype[0], normoption[0],
                                  normfile, n, p0, p1, p2);
   }//if

// initialize expressor
   if (strcmp(exprtype[0], "none") != 0) {
      n  = *nexprpar;
      p0 = p1 = p2 = p3 = p4 = p5 = p6 = p7 = 0.0;
      if (n > 0) p0 = exprpars[0];
      if (n > 1) p1 = exprpars[1];
      if (n > 2) p2 = exprpars[2];
      if (n > 3) p3 = exprpars[3];
      if (n > 4) p4 = exprpars[4];
      if (n > 5) p5 = exprpars[5];
      if (n > 6) p6 = exprpars[6];
      if (n > 7) p7 = exprpars[7];

      if (strcmp(chiptype[0], "GeneChip") == 0) {
         r += manager->InitAlgorithm("selector", "probe", exprselection[0], 0, 0);
      } else if (strcmp(chiptype[0], "GenomeChip") == 0) {
         r += manager->InitAlgorithm("selector", "probe", "genome", 0, 1, *level);
      } else if (strcmp(chiptype[0], "ExonChip") == 0) {
         r += manager->InitAlgorithm("selector", "probe", "exon", 0, 1, *level);
      }//if
      r += manager->InitAlgorithm("expressor", exprtype[0], exproption[0],
                                  exprfile, n, p0, p1, p2, p3, p4, p5, p6, p7);
   }//if

// open root scheme file
   r += manager->OpenSchemes(schemefile[0], chipname[0]);

// create/update root data file 
   if (*update == 1) {
      r += manager->Update(filename[0], "preprocess", "R");
      // need to force to delete root file from memory to prevent R using former root file!
      manager->SetFileOwner(kTRUE);
   } else {
      r += manager->New(filename[0], dirname[0], chiptype[0], "preprocess");
   }//if

// add trees
   for (int i=0; i<*ntrees; i++) {
      r += manager->AddTree(treeset[0], treenames[i]);
   }//for_i

// set reference tree
   if (strcmp(normtype[0], "none") != 0) {
      r += manager->SetReference(reftree[0], refmethod[0], *refparam);
   }//if

// preprocess expression values and store as trees in new file
   r += manager->Preprocess(treeset[0], "preprocess");
   *err = (int)r;

// cleanup
   if (strcmp(tmpdir[0], "") != 0) {
      delete [] exprfile;
      delete [] normfile;
      delete [] bgrdfile;
   }//if

   manager->Close();
   delete manager;
}//Preprocess

/*____________________________________________________________________________*/
void BgCorrect(char **filename, char **dirname, char **chiptype,
               char **schemefile, char **tmpdir, char **seloption,
               char **bgrdtype, char **bgrdoption, int *npar, double *pars,
               char **treeset, char **treenames, int *ntrees,
               int *update, int *level, int *verbose, int *err)
{
// Correct for background

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager","", *verbose);

// initialize chip type
   int r = 0;
   r += manager->Initialize(chiptype[0]);

// temporary file for background algorithms
   char *tmpfile = new char[strlen(tmpdir[0]) + 22];
   if (strcmp(tmpdir[0], "") != 0) {
      tmpfile = strcpy(tmpfile, tmpdir[0]);
      tmpfile = strcat(tmpfile, "/tmp_bgrd_310151.root");
   } else {
      tmpfile = "";
   }//if

// copy parameters
   double p0, p1, p2, p3, p4, p5, p6, p7;
   p0 = p1 = p2 = p3 = p4 = p5 = p6 = p7 = 0.0;
   if (*npar > 0) p0 = pars[0];
   if (*npar > 1) p1 = pars[1];
   if (*npar > 2) p2 = pars[2];
   if (*npar > 3) p3 = pars[3];
   if (*npar > 4) p4 = pars[4];
   if (*npar > 5) p5 = pars[5];
   if (*npar > 6) p6 = pars[6];
   if (*npar > 7) p7 = pars[7];

// initialize backgrounder
   if (strcmp(seloption[0], "pmonly") == 0 ||
       strcmp(seloption[0], "mmonly") == 0 ||
       strcmp(seloption[0], "both")   == 0 ||
       strcmp(seloption[0], "all")    == 0 ||
       strcmp(seloption[0], "none")   == 0) {
      r += manager->InitAlgorithm("selector", "probe", seloption[0], 0, 0);
      r += manager->InitAlgorithm("backgrounder", bgrdtype[0], bgrdoption[0], tmpfile,
                                  *npar, p0, p1, p2, p3, p4, p5, p6, p7);
//mas4
//   manager->InitAlgorithm("selector","probe","all", 0, 0);
//   manager->InitAlgorithm("backgrounder","sector","subtractbg",0, 4,0.02,4,4,0);
//mas5
//   manager->InitAlgorithm("selector","probe","both", 0, 0);
//   manager->InitAlgorithm("backgrounder","weightedsector","correctbg", "", 6,0.02,4,4,0,100,0.5);
//rma
//   manager->InitAlgorithm("selector","probe","none", 0, 0);
//   manager->InitAlgorithm("backgrounder","rma","pmonly:epanechnikov","tmp", 1,16384);
//gc
//   manager->InitAlgorithm("selector","probe","none", 0, 0);
//   manager->InitAlgorithm("backgrounder","gccontent","attenuatebg", "", 3,0.4,0.005,-1.0);
   } else if (strcmp(seloption[0], "genomic")     == 0 ||
              strcmp(seloption[0], "antigenomic") == 0) {
//      int bgrdlevel = 16383;
      int bgrdlevel = *level; //use only those probes which are also used for expression
      int bgrdopt   = (strcmp(seloption[0], "genomic") == 0) ? -1 : -2;

      r += manager->InitAlgorithm("selector", "probe", "exon", 0, 2, bgrdlevel, bgrdopt);
      r += manager->InitAlgorithm("backgrounder", bgrdtype[0], bgrdoption[0], tmpfile,
                                  *npar, p0, p1, p2, p3, p4, p5, p6, p7);

//mas4
//   manager->InitAlgorithm("selector","probe","exon",0, 2, 1016, -2);
//   manager->InitAlgorithm("backgrounder","sector","subtractbg",0, 4,0.02,4,4,0);
//mas5
//   manager->InitAlgorithm("selector","probe","exon",0, 2, 1016, -2);
//   manager->InitAlgorithm("backgrounder","weightedsector","correctbg",0, 6,0.02,4,4,0,100,0.5);
//rma
//   manager->InitAlgorithm("selector","probe","exon",0, 2, 1016, -2);
//   manager->InitAlgorithm("backgrounder","rma","pmonly:epanechnikov",0, 1,16384);
//gc
//   manager->InitAlgorithm("selector","probe","exon",0, 2, 1016, -2);
//   manager->InitAlgorithm("backgrounder","gccontent","attenuatebg", "", 3,0.4,0.005,-1.0);
   }//if

// open root scheme file
   r += manager->OpenSchemes(schemefile[0]);

// create/update root data file 
   if (*update == 1) {
      r += manager->Update(filename[0], "preprocess", "R");
      // need to force to delete root file from memory to prevent R using former root file!
      manager->SetFileOwner(kTRUE);
   } else {
      r += manager->New(filename[0], dirname[0], chiptype[0], "preprocess");
   }//if

// add trees for rma
   for (int i=0; i<*ntrees; i++) {
      r += manager->AddTree(treeset[0], treenames[i]);
   }//for_i

// preprocess expression values and store as trees in new file
   r += manager->Preprocess(treeset[0], "AdjustBgrd");
   *err = (int)r;

// cleanup
   if (strcmp(tmpdir[0], "") != 0) {
      delete [] tmpfile;
   }//if

   manager->Close();
   delete manager;
}//BgCorrect

/*____________________________________________________________________________*/
void Normalize(char **filename, char **dirname, char **chiptype,
               char **schemefile, char **tmpdir, char **seloption,
               char **type, char **option, int *npar, double *pars,
               int *level, char **treeset, char **treenames, int *ntrees,
               char **reftree, char **refmethod, int *update,
               int *verbose, int *err)
{
// Normalize trees

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager","", *verbose);

// initialize chip type
   int r = 0;
   r += manager->Initialize(chiptype[0]);

// temporary file
   char *tmpfile = new char[strlen(tmpdir[0]) + 14];
   if (strcmp(tmpdir[0], "") != 0) {
      tmpfile = strcpy(tmpfile, tmpdir[0]);
      tmpfile = strcat(tmpfile, "/tmp_rkq.root");
   } else {
      tmpfile = "";
   }//if

// copy parameters
   double p0, p1;
   p0 = p1 = 0.0;
   if (*npar > 0) p0 = pars[0];
   if (*npar > 1) p1 = pars[1];

// initialize normalizer
   if (strcmp(chiptype[0], "GeneChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", seloption[0], 0);
   } else if (strcmp(chiptype[0], "GenomeChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", "genome", 0, 1, *level);
   } else if (strcmp(chiptype[0], "ExonChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", "exon", 0, 1, *level);
   }//if
   r += manager->InitAlgorithm("normalizer",type[0],option[0],tmpfile, *npar,p0,p1);
//NEED TO TEST:
   if (strcmp(type[0], "lowess") == 0 || strcmp(type[0], "supsmu") == 0) {
      r += manager->InitAlgorithm("normalizer", "approx", "linear:mean", "", 2, 0.0, 0.0);
   }//if


// open root scheme file
   r += manager->OpenSchemes(schemefile[0]);

// create/update root data file 
   if (*update == 1) {
      r += manager->Update(filename[0], "preprocess", "R");
      // need to force to delete root file from memory to prevent R using former root file!
      manager->SetFileOwner(kTRUE);
   } else {
      r += manager->New(filename[0], dirname[0], chiptype[0], "preprocess");
   }//if

// add trees for rma
   for (int i=0; i<*ntrees; i++) {
      r += manager->AddTree(treeset[0], treenames[i]);
   }//for_i

// set reference tree
   r += manager->SetReference(reftree[0], refmethod[0], 0.0);

// preprocess expression values and store as trees in new file
   r += manager->Preprocess(treeset[0], "preprocess");
   *err = (int)r;

// cleanup
   if (strcmp(tmpdir[0], "") != 0) {
      delete [] tmpfile;
   }//if

   manager->Close();
   delete manager;
}//Normalize

/*____________________________________________________________________________*/
void Summarize(char **filename, char **dirname, char **chipname, char **chiptype,
               char **schemefile, char **tmpdir, char **seloption, char **type,
               char **option, int *npar, double *pars, int *level,
               char **treeset, char **treenames, int *ntrees,
               int *update, int *verbose, int *err)
{
// Summarize trees

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager","", *verbose);

// initialize chip type
   int r = 0;
   r += manager->Initialize(chiptype[0]);

// temporary file for expression algorithms
   char *tmpfile = new char[strlen(tmpdir[0]) + 22];
   if (strcmp(tmpdir[0], "") != 0) {
      tmpfile = strcpy(tmpfile, tmpdir[0]);
      tmpfile = strcat(tmpfile, "/tmp_expr_310151.root");
   } else {
      tmpfile = "";
   }//if

// copy parameters
   double p0, p1, p2, p3, p4, p5, p6, p7;
   p0 = p1 = p2 = p3 = p4 = p5 = p6 = p7 = 0.0;
   if (*npar > 0) p0 = pars[0];
   if (*npar > 1) p1 = pars[1];
   if (*npar > 2) p2 = pars[2];
   if (*npar > 3) p3 = pars[3];
   if (*npar > 4) p4 = pars[4];
   if (*npar > 5) p5 = pars[5];
   if (*npar > 6) p6 = pars[6];
   if (*npar > 7) p7 = pars[7];

// initialize expressor
   if (strcmp(chiptype[0], "GeneChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", seloption[0], 0, 0);
   } else if (strcmp(chiptype[0], "GenomeChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", "genome", 0, 1, *level);
   } else if (strcmp(chiptype[0], "ExonChip") == 0) {
      r += manager->InitAlgorithm("selector", "probe", "exon", 0, 1, *level);
   }//if
   r += manager->InitAlgorithm("expressor", type[0], option[0], tmpfile,
                               *npar, p0, p1, p2, p3, p4, p5, p6, p7);

// open root scheme file
   r += manager->OpenSchemes(schemefile[0], chipname[0]);

// create/update root data file 
   if (*update == 1) {
      r += manager->Update(filename[0], "preprocess", "R");
      // need to force to delete root file from memory to prevent R using former root file!
      manager->SetFileOwner(kTRUE);
   } else {
      r += manager->New(filename[0], dirname[0], chiptype[0], "preprocess");
   }//if

// add trees for rma
   for (int i=0; i<*ntrees; i++) {
      r += manager->AddTree(treeset[0], treenames[i]);
   }//for_i

// preprocess expression values and store as trees in new file
   r += manager->Preprocess(treeset[0], "preprocess");
   *err = (int)r;

// cleanup
   if (strcmp(tmpdir[0], "") != 0) {
      delete [] tmpfile;
   }//if

   manager->Close();
   delete manager;
}//Summarize


/*////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// XNormationManager: normalize data, e.g. mean, median, lowess, supsmu         //
//                                                                              //
////////////////////////////////////////////////////////////////////////////////*/

/*____________________________________________________________________________*/
void Normxpress(char **filename, char **dirname, char **chiptype,
                char **schemefile, char **tmpdir, char **seloption, double *pc,
                char **type, char **option, int *npar, double *pars,
                int *level, char **treeset, char **treenames, int *ntrees,
                char **reftree, char **refmethod, int *update,
                int *verbose, int *err)
{
// Normalize trees using method: mean, median, lowess, supsmu

// create new normation manager
   XNormationManager *manager = new XNormationManager("NormationManager","", *verbose);

// initialize chip type
   int r = 0;
   r += manager->Initialize(chiptype[0]);

// temporary file
   char *tmpfile = new char[strlen(tmpdir[0]) + 14];
   if (strcmp(tmpdir[0], "") != 0) {
      tmpfile = strcpy(tmpfile, tmpdir[0]);
      tmpfile = strcat(tmpfile, "/tmp_rkq.root");
   } else {
      tmpfile = "";
   }//if

// copy parameters
   double p0, p1;
   p0 = p1 = 0.0;
   if (*npar > 0) p0 = pars[0];
   if (*npar > 1) p1 = pars[1];

// initialize algorithms
   r += manager->InitAlgorithm("selector","rank", seloption[0],"", 4, 0, 0.3, *pc, 0);
   r += manager->InitAlgorithm("normalizer", type[0], option[0], tmpfile, 2, p0, p1);
   if (strcmp(type[0], "lowess") == 0 || strcmp(type[0], "supsmu") == 0) {
      r += manager->InitAlgorithm("normalizer", "approx", "linear:mean", "", 2, 0.0, 0.0);
   }//if

// open root scheme file
   r += manager->OpenSchemes(schemefile[0]);

// create/update root data file 
   if (*update == 1) {
//x      r += manager->Update(filename[0]);
      r += manager->Update(filename[0], "normalize", "R");
      // need to force to delete root file from memory to prevent R using former root file!
      manager->SetFileOwner(kTRUE);
   } else {
      r += manager->New(filename[0], dirname[0], chiptype[0]);
   }//if

// add trees for rma
   for (int i=0; i<*ntrees; i++) {
      r += manager->AddTree(treeset[0], treenames[i]);
   }//for_i

// set reference tree
   r += manager->SetReference(reftree[0], refmethod[0], 0.0);

// preprocess expression values and store as trees in new file
   r += manager->Normalize(treeset[0]);
   *err = (int)r;

// cleanup
   if (strcmp(tmpdir[0], "") != 0) {
      delete [] tmpfile;
   }//if

// cleanup
   manager->Close();
   delete manager;
}//Normxpress


/*//////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Utility functions for utils.R                                              //
//                                                                            //
//////////////////////////////////////////////////////////////////////////////*/

/*____________________________________________________________________________*/
void ChipNameType(char **filename, char **nametype)
{
// open chip type from scheme file

// Open file
   TFile *file = TFile::Open(filename[0], "READ");
   if (!file || file->IsZombie()) {
      printf("Could not open file <%s>", filename[0]);
      return;
   }//if

// Get content from file
   XFolder *content = (XFolder*)(file->Get("Content"));
   if (!content) {
      printf("Content for file <%s> not found.", filename[0]);
      nametype[0] = nametype[1] = "";
      return;
   }//if

// Check if file is schemefile
   TString datatype = content->GetTitle();
   if (strcmp(content->GetTitle(), "Schemes") != 0) {
      printf("File <%s> is not a scheme file.", filename[0]);
      nametype[0] = nametype[1] = "";
      return;
   }//if

// Get chip type
   int i = 0;
   TIter next(content->GetListOfFolders());
   XDNAChip *chip = 0;
   TString chip_type;
   while ((chip = (XDNAChip*)next())) {
//      nametype[0] = (char*)(chip->GetName());
//      nametype[1] = (char*)(chip->GetTitle());

      chip_type = chip->GetName();
      nametype[0] = new char[chip_type.Length() + 1];
      strcpy(nametype[0], (char*)chip_type.Data());

      chip_type = chip->GetTitle();
      nametype[1] = new char[chip_type.Length() + 1];
      strcpy(nametype[1], (char*)chip_type.Data());
      i++;
   }//while
}//ChipNameType

/*____________________________________________________________________________*/
void GeneChipProbeInfo(char **filename, int *value)
{
// Return parameter of type "integer" as value

// Open file
   TFile *file = TFile::Open(filename[0], "READ");
   if (!file || file->IsZombie()) {
      printf("Could not open file <%s>", filename[0]);
      return;
   }//if

// Get content from file
   XFolder *content = (XFolder*)(file->Get("Content"));
   if (!content) {
      printf("Content for file <%s> not found.", filename[0]);
      for (int i=0; i<8; i++) value[i] = 0;
      return;
   }//if

// Check if file is schemefile
   TString datatype = content->GetTitle();
   if (strcmp(content->GetTitle(), "Schemes") != 0) {
      printf("File <%s> is not a scheme file.", filename[0]);
      for (int i=0; i<8; i++) value[i] = 0;
      return;
   }//if

// Get chip integer values
   int i = 0;
   TIter next(content->GetListOfFolders());
   XGeneChip *chip = 0;
   while ((chip = (XGeneChip*)next())) {
      value[0] = chip->GetNumRows();
      value[1] = chip->GetNumColumns();
      value[2] = chip->GetNumProbes();
      value[3] = chip->GetNumControls();
      value[4] = chip->GetNumGenes();
      value[5] = chip->GetNumUnits();
      value[6] = chip->GetNumProbesets();
      value[7] = chip->GetNumAffx();
      i++;
   }//while
}//GeneChipProbeInfo

/*____________________________________________________________________________*/
void GetNumberOfTrees(char **filename, char **setname, char **exten, int *numtrees)
{
// get number of trees stored in root file

   TString savedir = gSystem->WorkingDirectory();

// Open file
   TFile *file = TFile::Open(filename[0], "READ");
   if (!file || file->IsZombie()) {
      printf("Could not open file <%s>", filename[0]);
      return;
   }//if

// Change directory
   TDirectory *dir = file->GetDirectory(setname[0]);
   if (dir) {
      file->cd(setname[0]);
   } else {
      printf("Could not open file directory <%s>", setname[0]);
//?      *numtrees = 0;
      return;
   }//if

// Get number of trees
   int   ntrees = 0;
   TKey *key    = 0;
   TIter next(gDirectory->GetListOfKeys());
   while ((key = (TKey*)next())) {
      if (strcmp(key->GetClassName(), "TTree") != 0) continue;
      if (!(strcmp((Path2Name(key->GetName(),".",";")).Data(), exten[0]) == 0 ||
            strcmp(exten[0], "*") == 0)) continue;
      ntrees++;
   }//while
   *numtrees = ntrees;

// Close file
   file->Close();

   gSystem->ChangeDirectory(savedir);
}//GetNumberOfTrees

/*____________________________________________________________________________*/
void GetNumberOfTrees4Exten(char **filename, char **exten, int *numtrees)
{
// get number of trees stored in root file

   TString savedir = gSystem->WorkingDirectory();

// Open file
   TFile *file = TFile::Open(filename[0], "READ");
   if (!file || file->IsZombie()) {
      printf("Could not open file <%s>", filename[0]);
      return;
   }//if

// Get content from file
   XFolder *content = (XFolder*)(file->Get("Content"));
   if (!content) {
      printf("Content for file <%s> not found.", filename[0]);
      numtrees = 0;
      return;
   }//if

// Get number of trees
   int   ntrees = 0;
   TKey *key    = 0;
   TIter next1(content->GetListOfFolders());
   XTreeSet *set = 0;
   TString setname = "";
   while ((set = (XTreeSet*)next1())) {
      setname = set->GetName();

      // change directory
      TDirectory *dir = file->GetDirectory(setname);
      if (dir) {
         file->cd(setname);
      } else {
         printf("Could not open file directory <%s>", setname.Data());
         return;
      }//if

      // number of trees
      TIter *next2 = new TIter(gDirectory->GetListOfKeys());
      while ((key = (TKey*)next2->Next())) {
         if (strcmp(key->GetClassName(), "TTree") != 0) continue;
         if (!(strcmp((Path2Name(key->GetName(),".",";")).Data(), exten[0]) == 0 ||
               strcmp(exten[0], "*") == 0)) continue;
         ntrees++;
      }//while
      delete next2;
   }//while
   *numtrees = ntrees;

// Close file
   file->Close();

   gSystem->ChangeDirectory(savedir);
}//GetNumberOfTrees4Exten

/*____________________________________________________________________________*/
void GetTreeNames(char **filename, char **setname, char **exten, int *gettitle,
                  char **treenames)
{
// get tree names stored in root file

   TString savedir = gSystem->WorkingDirectory();

// Open file
   TFile *file = TFile::Open(filename[0], "READ");
   if (!file || file->IsZombie()) {
      printf("Could not open file <%s>", filename[0]);
      return;
   }//if

// Change directory
   TDirectory *dir = file->GetDirectory(setname[0]);
   if (dir) {
      file->cd(setname[0]);
   } else {
      printf("Could not open file directory <%s>", setname[0]);
      return;
   }//if

// Get number of trees
   int   ntrees = 0;
   TKey *key    = 0;
   TIter next(gDirectory->GetListOfKeys());
   while ((key = (TKey*)next())) {
      if (strcmp(key->GetClassName(), "TTree") != 0) continue;
      if (!(strcmp((Path2Name(key->GetName(),".",";")).Data(), exten[0]) == 0 ||
            strcmp(exten[0], "*") == 0)) continue;
      ntrees++;
   }//while

// Get tree names
   TString names[ntrees];
//   TString *names = new TString[ntrees];
   ntrees = 0;
   key    = 0;
   next.Reset();
   while ((key = (TKey*)next())) {
      if (strcmp(key->GetClassName(), "TTree") != 0) continue;
      if (!(strcmp((Path2Name(key->GetName(),".",";")).Data(), exten[0]) == 0 ||
            strcmp(exten[0], "*") == 0)) continue;
      names[ntrees] = (*gettitle == 0) ? key->GetName() : key->GetTitle();
      ntrees++;
   }//while

   for (int i=0; i<ntrees; i++) {
//      treenames[i] = (char*)(names[i].Data());
      treenames[i] = new char[names[i].Length() + 1];
      strcpy(treenames[i], (char*)names[i].Data());
   }//for_i

// Close file
   file->Close();

   gSystem->ChangeDirectory(savedir);
}//GetTreeNames

/*____________________________________________________________________________*/
void GetTreeNames4Exten(char **filename, char **exten, int *gettitle,
                        char **treenames)
{
// get tree names stored in root file

   TString savedir = gSystem->WorkingDirectory();

// Open file
   TFile *file = TFile::Open(filename[0], "READ");
   if (!file || file->IsZombie()) {
      printf("Could not open file <%s>", filename[0]);
      return;
   }//if

// Get content from file
   XFolder *content = (XFolder*)(file->Get("Content"));
   if (!content) {
      printf("Content for file <%s> not found.", filename[0]);
      treenames[0] = "";
      return;
   }//if

// Get number of trees
   int   ntrees = 0;
   TKey *key    = 0;
   TIter next1(content->GetListOfFolders());
   XTreeSet *set = 0;
   TString setname = "";
   while ((set = (XTreeSet*)next1())) {
      setname = set->GetName();

      // change directory
      TDirectory *dir = file->GetDirectory(setname);
      if (dir) {
         file->cd(setname);
      } else {
         printf("Could not open file directory <%s>", setname.Data());
         return;
      }//if

      // number of trees
      TIter *next2 = new TIter(gDirectory->GetListOfKeys());
      while ((key = (TKey*)next2->Next())) {
         if (strcmp(key->GetClassName(), "TTree") != 0) continue;
         if (!(strcmp((Path2Name(key->GetName(),".",";")).Data(), exten[0]) == 0 ||
               strcmp(exten[0], "*") == 0)) continue;
         ntrees++;
      }//while
      delete next2;
   }//while

// Get tree names
   TString names[ntrees];
//   TString *names = new TString[ntrees];
   ntrees = 0;
   key    = 0;
   next1.Reset();
   while ((set = (XTreeSet*)next1())) {
      setname = set->GetName();

      // change directory
      TDirectory *dir = file->GetDirectory(setname);
      if (dir) {
         file->cd(setname);
      } else {
         printf("Could not open file directory <%s>", setname.Data());
         return;
      }//if

      // number of trees
      TIter *next2 = new TIter(gDirectory->GetListOfKeys());
      while ((key = (TKey*)next2->Next())) {
         if (strcmp(key->GetClassName(), "TTree") != 0) continue;
         if (!(strcmp((Path2Name(key->GetName(),".",";")).Data(), exten[0]) == 0 ||
               strcmp(exten[0], "*") == 0)) continue;
         names[ntrees] = (*gettitle == 0) ? key->GetName() : key->GetTitle();
         ntrees++;
      }//while
      delete next2;
   }//while

   for (int i=0; i<ntrees; i++) {
//      treenames[i] = (char*)(names[i].Data());
      treenames[i] = new char[names[i].Length() + 1];
      strcpy(treenames[i], (char*)names[i].Data());
   }//for_i

// Close file
   file->Close();

   gSystem->ChangeDirectory(savedir);
}//GetTreeNames4Exten


/*////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// XPlot: functions for plotting data using ROOT Canvas                         //
//                                                                              //
////////////////////////////////////////////////////////////////////////////////*/

// NOTE: It is not possible to call this function from within R to display a
//       ROOT canvas with an image file since the canvas will not be opened.
//       We need to call macro "macroOpenImage.C" instead!

//______________________________________________________________________________
void PlotImage(char **filename, char **canvasname, char **treename,
               char **varlist, char **logbase, char **option, int *err)
{
// macro to draw image

// create new manager
   XPlot *plotter = new XPlot("Plot");

// open root file
   int r = 0;
   r += plotter->Open(filename[0]);

   plotter->NewCanvas(canvasname[0],"title",10,10,500,500);
//?   plotter->NewCanvas(canvasname[0],"title",10,10,*w+10,*h+10);

// draw in separate canvases
   r += plotter->DrawImage(canvasname[0],treename[0],varlist[0],logbase[0],option[0],0,0,"U");

   *err = (int)r;

// cleanup
   delete plotter;
}//PlotImage


/*____________________________________________________________________________*/

