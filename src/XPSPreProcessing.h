// File created: 08/05/2002                          last modified: 10/03/2008
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2008 Dr. Christian Stratowa
 *
 *  Written by: Christian Stratowa, Vienna, Austria <cstrato@aon.at>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy of the GNU General Public
 *  License is available at http://www.gnu.org/copyleft/gpl.html. You
 *  can also obtain it by writing to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA.
 *
 *******************************************************************************
 * Based on: "The ROOT System", http://root.cern.ch/                           *
 * ROOT:     An Object-Oriented Data Analysis Framework                        *
 * Authors:  Rene Brun and Fons Rademakers.                                    *
 * For the licensing terms of "The ROOT System" see $ROOTSYS/LICENSE.          *
 * For the list of contributors to "The ROOT System" see http://root.cern.ch/  *
 *******************************************************************************
 */

#ifndef __XPSPreProcessing__
#define __XPSPreProcessing__

#ifndef __XPSProcessing__
#include "XPSProcessing.h"
#endif

class XBackgrounder;
class XSelector;
class XNormalizer;
class XExpressor;
class XCallDetector;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreProcessManager                                                   //
//                                                                      //
// Class for managing pre-processing of microarray data                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XPreProcessManager: public XProcessManager {

   private:

   protected:
      virtual Int_t     ImportDefaults(const char *infile);
      virtual Int_t     InitDefaults();
      virtual XSetting *NewSetting(const char *type, const char *infile);
      virtual XTreeSet *NewTreeSet(const char *type);

   public:
      XPreProcessManager();
      XPreProcessManager(const char *name, const char *arraytype = "",
                         Int_t verbose = kTRUE);
      virtual ~XPreProcessManager();

      Int_t Preprocess(const char *setname, const char *method = "preprocess");
      Int_t Ratio(const char *intree, const char *reftree,
               const char *outtree = "");
//test (see answer Rene)
      void  MacroTest(const char *name);
      
      ClassDef(XPreProcessManager,1) //PreProcessManager
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCallTreeInfo                                                        //
//                                                                      //
// Class containing info about present call tree stored in fUserInfo    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XCallTreeInfo: public XTreeInfo {

   protected:
      Int_t      fNUnits;       //number of units
      Int_t      fNAbsent;      //number of absent call units
      Int_t      fNMarginal;    //number of marginal call units
      Int_t      fNPresent;     //number of present call units
      Double_t   fPcAbsent;     //percent absent call units
      Double_t   fPcMarginal;   //percent marginal call units
      Double_t   fPcPresent;    //percent present call units
      Double_t   fMinPValue;    //minimal detection call p-value
      Double_t   fMaxPValue;    //maximal detection call p-value

   public :
      XCallTreeInfo();
      XCallTreeInfo(const char *name, const char *title);
      virtual ~XCallTreeInfo();

      using XTreeInfo::AddUserInfo;
      virtual void     AddUserInfo(Int_t nunits, Int_t nabsent, Int_t nmarginal,
                          Int_t npresent, Double_t minpval, Double_t maxpval);
      virtual Double_t GetValue(const char *name);

      ClassDef(XCallTreeInfo,1) //CallTreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreProcesSetting                                                    //
//                                                                      //
// Class for initialization of pre-processing parameter settings        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XPreProcesSetting: public XProcesSetting {

   protected:
      XSelector     *fSelector;      //temporary/current selector
      XSelector     *fBgrdSelector;  //selector for background algorithm
      XBackgrounder *fBackgrounder;  //background algorithm
      XSelector     *fNormSelector;  //selector for normalization algorithm
      XNormalizer   *fNormalizer;    //normalization algorithm
      XSelector     *fExprSelector;  //selector for expression algorithm
      XExpressor    *fExpressor;     //expression algorithm
      XSelector     *fCallSelector;  //selector for call algorithm
      XCallDetector *fCaller;        //present call algorithm
//      XRatioAlgorithm *fRatio;        //ratio algorithm

   protected:
      Int_t InitBackgrounder(const char *type, Option_t *options,
               const char *filename, Int_t npars, Double_t *pars);
      Int_t InitSelector(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);
      Int_t InitNormalizer(const char *type, Option_t *options,
               const char *filename, Int_t npars, Double_t *pars);
      Int_t InitApprox(Option_t *options, Int_t npars, Double_t *pars);
      Int_t InitExpressor(const char *type, Option_t *options,
               const char *filename, Int_t npars, Double_t *pars);
      Int_t InitCallDetector(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);

   public:
      XPreProcesSetting();
      XPreProcesSetting(const char *arraytype, const char *infile);
      virtual ~XPreProcesSetting();

      virtual Int_t InitAlgorithm(const char *name, const char *type,
                       Option_t *options, const char *filename,
                       Int_t npars, Double_t *pars);

      XSelector     *GetSelector()     const {return fSelector;}
      XSelector     *GetBgrdSelector() const {return fBgrdSelector;}
      XBackgrounder *GetBackgrounder() const {return fBackgrounder;}
      XSelector     *GetNormSelector() const {return fNormSelector;}
      XNormalizer   *GetNormalizer()   const {return fNormalizer;}
      XSelector     *GetExprSelector() const {return fExprSelector;}
      XExpressor    *GetExpressor()    const {return fExpressor;}
      XSelector     *GetCallSelector() const {return fCallSelector;}
      XCallDetector *GetCallDetector() const {return fCaller;}

      ClassDef(XPreProcesSetting,1) //PreProcesSetting
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreProcesSet                                                        //
//                                                                      //
// Base class for microarray pre-processing type                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XPreProcesSet: public XProcesSet {

   protected:
      XSelector     *fBgrdSelector;  //! selector for background algorithm
      XBackgrounder *fBackgrounder;  //! background algorithm
      XSelector     *fNormSelector;  //! selector for normalization algorithm
      XNormalizer   *fNormalizer;    //! normalization algorithm
      XSelector     *fExprSelector;  //! selector for expression algorithm
      XExpressor    *fExpressor;     //! expression algorithm
      XSelector     *fCallSelector;  //! selector for call algorithm
      XCallDetector *fCaller;        //! present call algorithm
//      XRatioAlgorithm *fRatio;        //! ratio algorithm

   protected:
      virtual Int_t AdjustBackground(Int_t /*numdata*/, TTree ** /*datatree*/,
                       Int_t &/*numbgrd*/, TTree ** /*bgrdtree*/) {return 0;}
      virtual Int_t Normalize(Int_t /*numdata*/, TTree ** /*datatree*/,
                       Int_t &/*numbgrd*/, TTree ** /*bgrdtree*/) {return 0;}
      virtual Int_t Express(Int_t /*numdata*/, TTree ** /*datatree*/,
                       Int_t &/*numbgrd*/, TTree ** /*bgrdtree*/) {return 0;}
      virtual Int_t DetectCall(Int_t /*numdata*/, TTree ** /*datatree*/,
                       Int_t &/*numbgrd*/, TTree ** /*bgrdtree*/) {return 0;}

      virtual void  AddDataTreeInfo(TTree *tree, const char *name, Option_t *option,
                       Int_t nrows, Int_t ncols, Int_t nmin, Double_t min, Int_t nmax,
                       Double_t max, Int_t maxnpix);
      virtual void  AddMaskTreeInfo(TTree *tree, const char *name, Option_t *option,
                       Int_t nrows, Int_t ncols, Int_t nflags);
      virtual void  AddCallTreeInfo(TTree *tree, const char *name, Option_t *option,
                       Int_t nunits, Int_t nabsent, Int_t nmarginal, Int_t npresent,
                       Double_t minpval, Double_t maxpval);

      Int_t XY2Index(Int_t x, Int_t y, Int_t ncol) {return (x + y*ncol);}
      Int_t Index2X(Int_t index, Int_t ncol)       {return (index % ncol);}
      Int_t Index2Y(Int_t index, Int_t ncol)       {return (index / ncol);}

   public:
      XPreProcesSet();
      XPreProcesSet(const char *name, const char *type);
      virtual ~XPreProcesSet();

      virtual Int_t Initialize(TFile *file, XSetting *setting,
                       const char *infile = "", const char *treename = "");
      virtual Int_t Preprocess(const char * /*method*/)          {return 0;}

      ClassDef(XPreProcesSet,2) //PreProcesSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCProcesSet                                                         //
//                                                                      //
// Class for GeneChip oligonucleotide array pre-processing              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGCProcesSet: public XPreProcesSet {

   protected:
      Int_t     fNBgPar;        //number of parameters
      Double_t *fBgPars;        //[fNBgPar] Array of fNBgPar parameters

   protected:
      virtual Int_t AdjustBackground(Int_t numdata, TTree **datatree,
                       Int_t &numbgrd, TTree **bgrdtree);
      virtual Int_t Normalize(Int_t numdata, TTree **datatree,
                       Int_t &numbgrd, TTree **bgrdtree);
      virtual Int_t DetectCall(Int_t numdata, TTree **datatree,
                       Int_t &numbgrd, TTree **bgrdtree);
      virtual Int_t DoCall(Int_t numdata, TTree **datatree,
                       Int_t &numbgrd, TTree **bgrdtree);
      virtual Int_t DoMultichipCall(Int_t numdata, TTree **datatree,
                       Int_t &numbgrd, TTree **bgrdtree, TFile *file);
      virtual Int_t Express(Int_t numdata, TTree **datatree,
                       Int_t &numbgrd, TTree **bgrdtree);
      virtual Int_t DoExpress(Int_t numdata, TTree **datatree,
                       Int_t numbgrd, TTree **bgrdtree);
      virtual Int_t DoMultichipExpress(Int_t numdata, TTree **datatree,
                       Int_t numbgrd, TTree **bgrdtree, TFile *file);
      virtual Int_t ExportBgrdTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportIntnTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportNormTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportExprTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportCallTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t SchemeMask(XDNAChip *chip, Int_t level, Int_t n, Int_t *msk);

      virtual Int_t *FillMaskArray(XDNAChip *chip, TTree *scmtree, XScheme *scheme,
                       Int_t level, Int_t n, Int_t *msk);
      virtual Int_t *FillUnitArray(TTree *idxtree, XGCUnit *unit,
                         Int_t n, Int_t *arr, Int_t *msk);
      virtual void   FillProbeSets(Int_t &p, Int_t &m, Int_t &idx, Double_t *pm, 
                        Double_t *mm, Double_t *sp, Double_t *sm, Int_t *xp, 
                        Int_t *xm, Int_t msk, Double_t inten, Double_t stdev,
                        Double_t bgrd, Double_t bgdev, Int_t npix);
      virtual Int_t  MaxNumberCells(TTree *idxtree);
      virtual TTree *SchemeTree(XAlgorithm *algorithm, void *scheme, TLeaf **scmleaf);
      virtual TTree *UnitTree(XAlgorithm *algorithm, void *unit, Int_t &numunits);

      Double_t AdjustIntensity(Double_t inten, Double_t bgrd, Double_t stdv);
      Bool_t   BackgroundParameters(XAlgorithm *algorithm, const char *option);
      TString  ChipType(const char *type);
      Int_t    FillDataArrays(TTree *datatree, TTree *bgrdtree, Bool_t doBg,
                  Int_t nrow, Int_t ncol, Double_t *inten, Double_t *stdev,
                  Int_t *npix, Double_t *bgrd, Double_t *bgdev);
      Int_t    FillDataArrays(TTree *datatree, TTree *bgrdtree, Bool_t doBg,
                  Int_t nrow, Int_t ncol, Double_t *inten, Double_t *stdev, Int_t *npix);
      Int_t    FillDataArrays(TTree *datatree, Int_t nrow, Int_t ncol,
                  Double_t *inten, Double_t *stdev, Int_t *npix);
      Int_t    FillBgrdArrays(TTree *bgrdtree, Int_t nrow, Int_t ncol,
                  Double_t *inten, Double_t *stdev);
      TTree   *FillDataTree(TTree *oldtree, const char *exten, XAlgorithm *algorithm,
                  Int_t nrow, Int_t ncol, Double_t *inten, Double_t *stdev);
      TTree   *FillDataTree(const char *name, XAlgorithm *algorithm, Int_t nrow,
                  Int_t ncol, Double_t *arr);
      Int_t    FillMaskTree(const char *name, XAlgorithm *algorithm, Int_t nrow,
                  Int_t ncol, Int_t *arr);
      void     FillProbeSets(Int_t &p, Int_t &idx, Double_t *pm, Double_t *sp, 
                  Int_t *xp, Int_t msk, Double_t inten, Double_t stdev);
      Int_t    MaskArray2GC(XDNAChip *chip, Int_t *msk);
      Int_t    MeanReference(Int_t numdata, TTree **datatree,
                  Int_t numbgrd, TTree **bgrdtree,
                  Int_t nrow, Int_t ncol, Double_t *arr, Bool_t doBg);
      Int_t    MedianReference(Int_t numdata, TTree **datatree,
                  Int_t numbgrd, TTree **bgrdtree,
                  Int_t nrow, Int_t ncol, Double_t *arr, Bool_t doBg);

   public:
      XGCProcesSet();
      XGCProcesSet(const char *name, const char *type);
      virtual ~XGCProcesSet();

      virtual Int_t Preprocess(const char *method);
      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);

      ClassDef(XGCProcesSet,2) //GCProcesSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeProcesSet                                                     //
//                                                                      //
// Class for GenomeChip oligonucleotide array pre-processing            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGenomeProcesSet: public XGCProcesSet {

   protected:

   protected:
      virtual Int_t *FillMaskArray(XDNAChip *chip, TTree *scmtree, XScheme *scheme,
                        Int_t level, Int_t n, Int_t *msk);
      virtual Int_t *FillUnitArray(TTree *idxtree, XGCUnit *unit,
                         Int_t n, Int_t *arr, Int_t *msk);
      virtual void   FillProbeSets(Int_t &p, Int_t &m, Int_t &idx, Double_t *pm, 
                        Double_t *mm, Double_t *sp, Double_t *sm, Int_t *xp, 
                        Int_t *xm, Int_t msk, Double_t inten, Double_t stdev,
                        Double_t bgrd, Double_t bgdev, Int_t npix);
      virtual Int_t  MaxNumberCells(TTree *idxtree);

   public:
      XGenomeProcesSet();
      XGenomeProcesSet(const char *name, const char *type);
      virtual ~XGenomeProcesSet();

      ClassDef(XGenomeProcesSet,2) //GenomeProcesSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonProcesSet                                                       //
//                                                                      //
// Class for ExonChip oligonucleotide array pre-processing              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExonProcesSet: public XGenomeProcesSet {

   protected:

   protected:
      virtual Int_t *FillMaskArray(XDNAChip *chip, TTree *scmtree, XScheme *scheme,
                        Int_t level, Int_t n, Int_t *msk);
      virtual TTree *SchemeTree(XAlgorithm *algorithm, void *scheme, TLeaf **scmleaf);
      virtual TTree *UnitTree(XAlgorithm *algorithm, void *unit, Int_t &numunits);

      virtual const char *GetTranscriptID(XUnit *unit, XTransAnnotation *annot, Int_t type);

   public:
      XExonProcesSet();
      XExonProcesSet(const char *name, const char *type);
      virtual ~XExonProcesSet();

      ClassDef(XExonProcesSet,2) //ExonProcesSet
};

#endif

