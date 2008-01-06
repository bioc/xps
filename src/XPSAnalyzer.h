// File created: 12/16/2002                          last modified: 11/11/2007
// Author: Christian Stratowa 06/18/2000

/******************************************************************************
* Copyright(c) 2000-2006, Dr. Christian Stratowa, Vienna, Austria.            *
* All rights reserved.                                                        *
* Author: Christian Stratowa.                                                 *
*                                                                             *
*******************************************************************************
*********************  XPS - eXpression Profiling System  *********************
*******************************************************************************
*                                                                             *
* Based on: "The ROOT System", All rights reserved.                           *
* Authors: Rene Brun and Fons Rademakers.                                     *
* For the licensing terms of "The ROOT System" see $ROOTSYS/AA_LICENSE.       *
* For the list of contributors to "The ROOT System" see $ROOTSYS/AA_CREDITS.  *
******************************************************************************/

#ifndef __XPSAnalyzer__
#define __XPSAnalyzer__

#include <Riostream.h>
#include <cmath>

#ifndef __XPSProcessing__
#include "XPSProcessing.h"
#endif

class TUnivariateTest;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAnalyser                                                            //
//                                                                      //
// Base class for analysis of data                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XAnalyser: public XAlgorithm {

   protected:

   public:
      XAnalyser()                {}
      XAnalyser(const char *name, const char *type)
         :XAlgorithm(name, type) {}
      virtual ~XAnalyser()       {}

      ClassDef(XAnalyser,1) //Analyser
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUniTester                                                           //
//                                                                      //
// Class for univariate analysis of normalized data                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XUniTester: public XAnalyser {

   protected:
      TUnivariateTest  *fUniTest;     //! univariate test

   public:
      XUniTester();
      XUniTester(const char *name, const char *type);
      virtual ~XUniTester();

      virtual Int_t InitType(const char *type, Option_t *options,
                       Int_t npars, Double_t *pars);
      virtual Int_t Analyse(const char *infile, const char *outfile,
                       const char *varlist, Int_t nrows,
                       Int_t nperm = -1, Double_t mu = 0,
                       const char *sepi = "\t", const char *sepo = "\t",
                       char delim = '\n', Int_t linebuf = 16635);
      virtual Int_t Analyse(Int_t n, Int_t *gid, TTree **intree, const char *leafname,
                       TTree *outtree, const char *varlist = "*",
                       Int_t nperm = -1, Double_t mu = 0);

      TUnivariateTest *GetUniTest() const {return fUniTest;}

      ClassDef(XUniTester,1) //UniTester
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultiTester                                                         //
//                                                                      //
// Class for multivariate analysis of normalized data                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMultiTester: public XAnalyser {

   protected:
//      TMultivariateTest  *fMultiTest;     //! multivariate test

   public:
      XMultiTester();
      XMultiTester(const char *name, const char *type);
      virtual ~XMultiTester();

      virtual Int_t InitType(const char *type, Option_t *options,
                       Int_t npars, Double_t *pars);
      virtual Int_t Analyse(const char *infile, const char *outfile,
                       const char *varlist, Int_t nrows,
                       Int_t nperm = -1, Double_t mu = 0,
                       const char *sepi = "\t", const char *sepo = "\t",
                       char delim = '\n', Int_t linebuf = 16635);
      virtual Int_t Analyse(Int_t n, Int_t *gid, TTree **intree, const char *leafname,
                       TTree *outtree, const char *varlist = "*",
                       Int_t nperm = -1, Double_t mu = 0);

//      TMultivariateTest *GetMultiTest() const {return fMultiTest;}

      ClassDef(XMultiTester,1) //XultiTester
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XClusterizer                                                         //
//                                                                      //
// Class for multivariate analysis of normalized data                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XClusterizer: public XAnalyser {

   protected:
//      TClusterAnalysis  *fCluster;     //! cluster analysis

   public:
      XClusterizer();
      XClusterizer(const char *name, const char *type);
      virtual ~XClusterizer();

      virtual Int_t InitType(const char *type, Option_t *options,
                       Int_t npars, Double_t *pars);
      virtual Int_t Analyse(const char *infile, const char *outfile,
                       const char *varlist, Int_t nrows,
                       Int_t nperm = -1, Double_t mu = 0,
                       const char *sepi = "\t", const char *sepo = "\t",
                       char delim = '\n', Int_t linebuf = 16635);
      virtual Int_t Analyse(Int_t n, Int_t *gid, TTree **intree, const char *leafname,
                       TTree *outtree, const char *varlist = "*",
                       Int_t nperm = -1, Double_t mu = 0);

//      TClusterAnalysis *GetClusterAnalysis() const {return fCluster;}

      ClassDef(XClusterizer,1) //Clusterizer
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRegressor                                                           //
//                                                                      //
// Class for regression analysis of normalized data                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XRegressor: public XAnalyser {

   protected:
//      TRegression  *fRegression;     //! regression analysis

   public:
      XRegressor();
      XRegressor(const char *name, const char *type);
      virtual ~XRegressor();

      virtual Int_t InitType(const char *type, Option_t *options,
                       Int_t npars, Double_t *pars);
      virtual Int_t Analyse(const char *infile, const char *outfile,
                       const char *varlist, Int_t nrows,
                       Int_t nperm = -1, Double_t mu = 0,
                       const char *sepi = "\t", const char *sepo = "\t",
                       char delim = '\n', Int_t linebuf = 16635);
      virtual Int_t Analyse(Int_t n, Int_t *gid, TTree **intree, const char *leafname,
                       TTree *outtree, const char *varlist = "*",
                       Int_t nperm = -1, Double_t mu = 0);

//      TRegression *GetRegression() const {return fRegression;}

      ClassDef(XRegressor,1) //Regressor
};

#endif

