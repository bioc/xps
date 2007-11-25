/*
 ***************************************************************
 *
 * File: rwrapper.h
 *
 * Implementation by: Christian Stratowa
 *
 * Copyright (C) Christian Stratowa 2002-2007
 *
 * A wrapper for the XPS libraries
 *
 ******************************************************************
*/

#ifndef R_WRAPPER
#define R_WRAPPER 1

extern "C" {
void ImportExprSchemes(char **filename, char **dirname, char **chipname,
                       char **scheme, char **probe, char **annot,
                       int *verbose, int *err);
void ImportExonSchemes(char **filename, char **dirname, char **chipname,
                       char **layout, char **scheme, char **probeset,
                       char **transcript, char **control,
                       int *verbose, int *err);
void ImportGenomeSchemes(char **filename, char **dirname, char **chipname,
                         char **layout, char **scheme, char **transcript,
                         int *verbose, int *err);
void ExportScheme(char **filename, char **treename, char **varlist,
                  char **outfile, char **sep, int *verbose);

void ImportData(char **filename, char **dirname, char **chiptype,
                char **chipname, char **schemefile, char **treeset,
                char **celfiles, char **celnames, int *numdata,
                char **project, int *nproject, char **author, int *nauthor,
                char **dataset, int *ndataset, char **source, int *nsource,
                char **sample, int *nsample, char **cell, int *ncell,
                char **pcell, int *npcell, char **tissue, int *ntissue,
                char **biopsy, int *nbiopsy, char **array, int *narray,
                char **hyb, int *nhyb, char **treat, int *ntreat, 
                int *replace, int *update, int *verbose, int *err);
void ExportData(char **filename, char **schemefile, char **chiptype,
                char **datatype, char **treenames, int *ntrees, char **exten,
                char **varlist, char **outfile, char **sep,
                int *verbose, int *err);

void PreprocessRMA(char **filename, char **dirname, char **chipname,
                   char **chiptype, char **schemefile, char **tmpdir,
                   char **bgrdoption, char **exproption, char **treeset,
                   char **treenames, int *ntrees, int *normalize,
                   int *level, int *verbose, int *err);
void PreprocessMAS4(char **filename, char **dirname, char **chipname,
                    char **chiptype, char **schemefile, char **tmpdir,
                    char **exproption, char **treeset, char **treenames,
                    int *ntrees, int *level, int *verbose, int *err);
void PreprocessMAS5(char **filename, char **dirname, char **chipname,
                    char **chiptype, char **schemefile, char **tmpdir,
                    char **exproption, char **treeset, char **treenames,
                    int *ntrees, int *level, int *verbose, int *err);

void PreprocessMAS5Call(char **filename, char **dirname, char **chipname,
                        char **chiptype, char **schemefile, char **tmpdir,
                        char **calloption, char **treeset, char **treenames,
                        int *ntrees, double *tau, double *alpha1, double *alpha2,
                        int *ignore, int *level, int *verbose, int *err);
void PreprocessDABGCall(char **filename, char **dirname, char **chipname,
                        char **chiptype, char **schemefile, char **calloption,
                        char **treeset, char **treenames, int *ntrees,
                        double *alpha1, double *alpha2, int *level,
                        int *verbose, int *err);

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
                int *level, int *verbose, int *err);
void BgCorrect(char **filename, char **dirname, char **chiptype,
               char **schemefile, char **tmpdir, char **seloption,
               char **bgrdtype, char **bgrdoption, int *npars, double *pars,
               char **treeset, char **treenames, int *ntrees,
               int *update, int *level, int *verbose, int *err);
void Normalize(char **filename, char **dirname, char **chiptype,
               char **schemefile, char **tmpdir, char **seloption,
               char **type, char **option, int *npar, double *pars,
               int *level, char **treeset, char **treenames, int *ntrees,
               char **reftree, char **refmethod, int *update,
               int *verbose, int *err);
void Summarize(char **filename, char **dirname, char **chipname, char **chiptype,
               char **schemefile, char **tmpdir, char **seloption, char **type,
               char **option, int *npar, double *pars, int *level,
               char **treeset, char **treenames, int *ntrees,
               int *update, int *verbose, int *err);

void Normxpress(char **filename, char **dirname, char **chiptype,
                char **schemefile, char **tmpdir, char **seloption, double *pc,
                char **type, char **option, int *npar, double *pars,
                int *level, char **treeset, char **treenames, int *ntrees,
                char **reftree, char **refmethod, int *update,
                int *verbose, int *err);

void ChipNameType(char **filename, char **nametype);
void GeneChipProbeInfo(char **filename, int *value);
void GetNumberOfTrees(char **filename, char **setname, char **exten, int *numtrees);
void GetNumberOfTrees4Exten(char **filename, char **exten, int *numtrees);
void GetTreeNames(char **filename, char **setname, char **exten, int *gettitle,
                  char **treenames);
void GetTreeNames4Exten(char **filename, char **exten, int *gettitle,
                        char **treenames);

void PlotImage(char **filename, char **canvasname, char **treename,
               char **varlist, char **logbase, char **option, int *err);
}

#endif
