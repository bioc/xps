// File created: 08/05/2002                          last modified: 10/12/2007
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2007 Dr. Christian Stratowa
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


/******************************************************************************
* Major Revision History:
* May 2002 - Initial versions finished
* Jun 2002 - Remove InitXXX() (replaced with InitAlgorithm() in XManager)
* Dec 2002 - Completely reorganize source file structure
* Aug 2004 - Rename XxxxAlgorithm to XBackgrounder, XExpressor, XCallDetector
*          - Use XNormalizer, add additional method Preprocess()
* Sep 2004 - Use XSelector and derived classes
* Oct 2004 - Add classes XExpressionTreeInfo, XCallTreeInfo
*          - Implement algorithms for RMA and MAS5
* May 2005 - Replace XContent with XFolder.
* Jul 2005 - Replace TTree::AddFriend() with TList fTrees.
*          - Store trees for every XTreeSet in separate TDirectory.
* Jan 2006 - To decrease size of saved trees replace Int_t with Short_t for 
*            mask/flag arrays
*          - Add parameter filename to InitAlgorithm() and field fFile and
*            method NewFile() to XAlgorithm to allow saving certain data such
*            as e.g. temporary trees in a (temporary) file
* Mar 2006 - Add support for alternative schemes, i.e. CDF files.
* Apr 2006 - Optionally reduce memory consumption of MedianPolish().
* Oct 2006 - Add support for exon arrays, class XExonProcesSet
* Dec 2006 - Implement background class XGCBackground and call class XDABGCall
*          - Add method AdjustIntensity() to allow different background corrections
* Aug 2007 - Add support for whole genome arrays, class XGenomeProcesSet
*
******************************************************************************/

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include <cfloat>
#include <new>  //needed for new (nothrow)


#include "TBranch.h"
#include "TFriendElement.h"
#include "THashTable.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TSystem.h"

//for Benchmark test only:
#include <TBenchmark.h>

#include "XPSPreProcessing.h"
#include "XPSHybridizer.h"
#include "XPSSelector.h"
#include "XPSNormalizer.h"
#include "TStat.h"

//debug: print function names
const Bool_t  kCS  = 0; 
const Bool_t  kCSa = 0; //debug: print function names in loops

ClassImp(XPreProcessManager);
ClassImp(XCallTreeInfo);
ClassImp(XPreProcesSetting);
ClassImp(XPreProcesSet);
ClassImp(XGCProcesSet);
ClassImp(XGenomeProcesSet);
ClassImp(XExonProcesSet);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreProcessManager                                                   //
//                                                                      //
// Class for managing pre-processing of microarray data                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XPreProcessManager::XPreProcessManager()
                   :XProcessManager()
{
   // Default PreProcessManager constructor
   if(kCS) cout << "---XPreProcessManager::XProcessingManager(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XPreProcessManager::XPreProcessManager(const char *name, const char *arraytype,
                    Int_t verbose)
                   :XProcessManager(name, arraytype, verbose)
{
   // Normal PreProcessManager constructor
   // If initfile is given, import default values
   if(kCS) cout << "---XPreProcessManager::XPreProcessManager------" << endl;

//TEST Benchmark
//   if (!gBenchmark) gBenchmark = new TBenchmark();
}//Constructor

//______________________________________________________________________________
XPreProcessManager::~XPreProcessManager()
{
   // PreProcessManager destructor
   if(kCS) cout << "---XPreProcessManager::~XPreProcessManager------" << endl;

//TEST Benchmark
//   SafeDelete(gBenchmark); 
}//Destructor

//______________________________________________________________________________
Int_t XPreProcessManager::Preprocess(const char *setname, const char *method)
{
   // Preprocess all trees added to tree set "setname" with AddTree().
   // Calculate background, intensity, expression level and/or present call.
   // The following methods are available:
   // method = "adjustbgrd": calculate background trees used for bgrd adjustment
   // method = "normalize":  normalize raw or background corrected intensities
   // method = "express":    convert intensities to expression levels
   // method = "detectcall": calculate present calls
   // The following methods can be used to preprocess all initialized algorithms:
   // method = "preprocess": preprocess all initialized algorithms (default)
   // method = "rma":        use the rma method of Irizarry et al. for preprocessing
   // method = "mas5":       preprocess data according to Affymetrix MAS5
   // method = "mas4":       preprocess data according to Affymetrix MAS4
   // Note: Methods "rma", "mas5", "mas4" will automatically initialize the correct
   //       algorithms by calling the corresponding InitAlgorithm() methods
   if(kCS) cout << "------XPreProcessManager::Preprocess------" << endl;

   if (fAbort) return errAbort;

// Create directory with name of treeset
   TDirectory *dir = fFile->GetDirectory(setname);
   if (!dir) fFile->mkdir(setname, fDataType);
   fFile->cd();

// Get treeset "setname" (set created in AddTree())
   XPreProcesSet *set = 0;
//??   set = (XPreProcesSet*)fContent->FindObject(setname, "XPreProcesSet");
   set = (XPreProcesSet*)fContent->FindObject(setname);
   if (!set) {
      return HandleError(errGetTreeSet, setname);
   }//if

// Set must be inherited from XPreProcesSet
   if (!set->InheritsFrom("XPreProcesSet")) {
      return HandleError(errClassTreeSet, setname, set->ClassName());
   }//if

// to lower
   TString smethod = method; smethod.ToLower();

// Initialize algorithms for methods "rma" or "mas"
   Int_t err = errNoErr;
   if (strcmp(smethod.Data(), "rma") == 0) {
      if (!err) err = InitAlgorithm("selector","probe","none", 0, 0);
      if (!err) err = InitAlgorithm("backgrounder","rma","pmonly:epanechnikov", 0, 1,16384);
      if (!err) err = InitAlgorithm("selector","probe","pmonly", 0);
      if (!err) err = InitAlgorithm("normalizer","quantile","together:none:0", 0, 1,0.0);
      if (!err) err = InitAlgorithm("selector","probe","pmonly", 0);
      if (!err) err = InitAlgorithm("expressor","medianpolish","log2",0,3,10,0.01,1.0);
   } else if (strcmp(smethod.Data(), "mas5") == 0) {
      if (!err) err = InitAlgorithm("selector","probe","both",0, 0);
      if (!err) err = InitAlgorithm("backgrounder","weightedsector","correctbg",0, 6,0.02,4,4,0,100,0.5);
      if (!err) err = InitAlgorithm("selector","probe","none",0, 0);
      if (!err) err = InitAlgorithm("expressor","TukeyBiweight","log2",0, 7,0.03,10.0,2.0e-20,5.0,0.0001,1.0,0.5);
      if (!err) err = InitAlgorithm("selector","default","none",0, 0);
      if (!err) err = InitAlgorithm("calldetector","dc5","raw",0, 6,0.015,0.04,0.06,1,0,0);
   } else if (strcmp(smethod.Data(), "mas4") == 0) {
      if (!err) err = InitAlgorithm("selector","probe","all",0, 0);
      if (!err) err = InitAlgorithm("backgrounder","sector","subtractbg",0, 4,0.02,4,4,0);
      if (!err) err = InitAlgorithm("selector","default","none",0, 0);
      if (!err) err = InitAlgorithm("expressor","avgdiff","0",0, 1,3.0);
      if (!err) err = InitAlgorithm("selector","default","none",0, 0);
      if (!err) err = InitAlgorithm("calldetector","dc5","raw",0, 6,0.015,0.04,0.06,1,0,0);
   } else if ((strcmp(smethod.Data(), "preprocess") == 0) ||
              (strcmp(smethod.Data(), "adjustbgrd") == 0) ||
              (strcmp(smethod.Data(), "normalize")  == 0) ||
              (strcmp(smethod.Data(), "express")    == 0) ||
              (strcmp(smethod.Data(), "detectcall") == 0)) {
      /*ok*/;
   } else {
      return HandleError(errAlgorithm, "Preprocessing", method);
   }//if

// Preprocess trees using method smethod
//ev?   if ((set->GetNumTrees() > 1) && (set->GetNumSelections() > 1)) {
   if (set->GetNumSelections() > 1) {
      if (!err) err = set->Initialize(fFile, fSetting);
      if (!err) err = set->Preprocess(smethod.Data());

      HandleError(err, "in XPreProcessManager::Preprocess", "???");
   } else {
//better in XGCProcesSet::Preprocess()!!
      cerr << "Error: At least two trees need to be selected." << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

//ev??   SafeDelete(set);  //no???

   return err;
}//Preprocess

//______________________________________________________________________________
Int_t XPreProcessManager::Ratio(const char *intree, const char *reftree,
                          const char *outtree)
{
   // Calculate ratio intree versus reference tree. Reference can be another
   // data tree (e.g. at time 0) or "mean" (mean of all data trees).
   // Store result in outtree. If outtree = "" use "name" part from intree.
   if(kCS) cout << "------XPreProcessManager::Ratio------" << endl;

// NOT YET USED - to prevent compiler warnings:
   intree = 0; reftree = 0; outtree = 0;

   if (fAbort) return errAbort;

   Int_t err = errNoErr;

   cerr << "Error: XPreProcessManager::Ratio() is not yet implemented." << endl;

   if (err) fAbort = kTRUE;
   return err;
}//Ratio

//______________________________________________________________________________
void XPreProcessManager::MacroTest(const char *name)
{
   // Testing only
   if(kCS) cout << "------XPreProcessManager::MacroTest------" << endl;

// NOT YET USED - to prevent compiler warnings:
   name = 0;

}//MacroTest

//______________________________________________________________________________
Int_t XPreProcessManager::ImportDefaults(const char *infile)
{
   // Import default values from infile
   if(kCS) cout << "------XPreProcessManager::ImportDefaults------" << endl;

// NOT YET USED - to prevent compiler warnings:
   infile = 0;

   return errNoErr;
}//ImportDefaults

//______________________________________________________________________________
Int_t XPreProcessManager::InitDefaults()
{
   // Initialize default values
   if(kCS) cout << "------XPreProcessManager::InitDefaults------" << endl;

   //IMPORTANT: do not initialize default algorithms since method = "preprocess"
   //           will preprocess all initialized algorithms! In this way the
   //           user can choose the combination of algorithms to be processed.
   return errNoErr;
}//InitDefaults

//______________________________________________________________________________
XSetting *XPreProcessManager::NewSetting(const char *type, const char *infile)
{
   // Create new setting for type
   if(kCS) cout << "------XPreProcessManager::NewSetting------" << endl;

   XPreProcesSetting *setting = new XPreProcesSetting(type, infile);
   return setting;
}//NewSetting

//______________________________________________________________________________
XTreeSet *XPreProcessManager::NewTreeSet(const char *type)
{
   // Create new tree set of type
   if(kCS) cout << "------XPreProcessManager::NewTreeSet------" << endl;

   XTreeSet *set = 0;
   if (strcmp(type,"GeneChip") == 0) {
      set = new XGCProcesSet("GeneChipProcesSet", type);
   } else if (strcmp(type,"SNPChip") == 0) {
//      set = new XSNPProcesSet("SNPChipProcesSet", type);
cout << "Note: to be done in the future: SNPChip analysis" << endl;
   } else if (strcmp(type,"GenomeChip") == 0) {
      set = new XGenomeProcesSet("GenomeChipProcesSet", type);
   } else if (strcmp(type,"ExonChip") == 0) {
      set = new XExonProcesSet("ExonChipProcesSet", type);
   } else if (strcmp(type,"GenePix") == 0) {
//      set = new XGPProcesSet("GenePixProcesSet", type);
cout << "Note: to be done in the future: GenePix analysis" << endl;
   } else {
      HandleError(errUnknownType, "Array", type);
   }//if

   return set;
}//NewTreeSet


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCallTreeInfo                                                        //
//                                                                      //
// Class containing info about present call tree stored in fUserInfo    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XCallTreeInfo::XCallTreeInfo() 
              :XTreeInfo()
{
   // Default CallTreeInfo constructor
   if(kCS) cout << "---XCallTreeInfo::XCallTreeInfo(default)------" << endl;

   fNUnits     = 0;
   fNAbsent    = 0;
   fNMarginal  = 0;
   fNPresent   = 0;
   fPcAbsent   = 0;
   fPcMarginal = 0;
   fPcPresent  = 0;
   fMinPValue  = -1;
   fMaxPValue  = -1;
}//Constructor

//______________________________________________________________________________
XCallTreeInfo::XCallTreeInfo(const char *name, const char *title) 
              :XTreeInfo(name, title)
{
   // Normal CallTreeInfo constructor
   if(kCS) cout << "---XCallTreeInfo::XCallTreeInfo------" << endl;

   fNUnits     = 0;
   fNAbsent    = 0;
   fNMarginal  = 0;
   fNPresent   = 0;
   fPcAbsent   = 0;
   fPcMarginal = 0;
   fPcPresent  = 0;
   fMinPValue  = -1;
   fMaxPValue  = -1;
}//Constructor

//______________________________________________________________________________
XCallTreeInfo::~XCallTreeInfo()
{
   // CallTreeInfo destructor
   if(kCS) cout << "---XCallTreeInfo::~XCallTreeInfo------" << endl;

}//Destructor

//______________________________________________________________________________
void XCallTreeInfo::AddUserInfo(Int_t nunits, Int_t nabsent, Int_t nmarginal,
                    Int_t npresent, Double_t minpval, Double_t maxpval)
{
   // Add user info from tree set
   if(kCS) cout << "------XCallTreeInfo::AddUserInfo------" << endl;

   fNUnits     = nunits;
   fNAbsent    = nabsent;
   fNMarginal  = nmarginal;
   fNPresent   = npresent;
   fPcAbsent   = 100.0*(Double_t)nabsent/(Double_t)nunits;
   fPcMarginal = 100.0*(Double_t)nmarginal/(Double_t)nunits;
   fPcPresent  = 100.0*(Double_t)npresent/(Double_t)nunits;
   fMinPValue  = minpval;
   fMaxPValue  = maxpval;
}//AddUserInfo

//______________________________________________________________________________
Double_t XCallTreeInfo::GetValue(const char *name)
{
   // Return value for class member field name
   if(kCS) cout << "------XCallTreeInfo::GetValue------" << endl;

   if (strcmp(name, "fNUnits") == 0) {
      return fNUnits;
   } else if (strcmp(name, "fNAbsent")    == 0) {
      return fNAbsent;
   } else if (strcmp(name, "fNMarginal")  == 0) {
      return fNMarginal;
   } else if (strcmp(name, "fNPresent")   == 0) {
      return fNPresent;
   } else if (strcmp(name, "fPcAbsent")   == 0) {
      return fPcAbsent;
   } else if (strcmp(name, "fPcMarginal") == 0) {
      return fPcMarginal;
   } else if (strcmp(name, "fPcPresent")  == 0) {
      return fPcPresent;
   } else if (strcmp(name, "fMinPValue")  == 0) {
      return fMinPValue;
   } else if (strcmp(name, "fMaxPValue")  == 0) {
      return fMaxPValue;
   }//if
   return 0;
}//GetValue


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreProcesSetting                                                    //
//                                                                      //
// Class for initialization of pre-processing parameter settings        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XPreProcesSetting::XPreProcesSetting()
                  :XProcesSetting()
{
   // Default PreProcesSetting constructor
   if(kCS) cout << "---XPreProcesSetting::XPreProcesSetting(default)------" << endl;

   fSelector     = 0;
   fBgrdSelector = 0;
   fBackgrounder = 0;
   fNormSelector = 0;
   fNormalizer   = 0;
   fExprSelector = 0;
   fExpressor    = 0;
   fCallSelector = 0;
   fCaller       = 0;
}//Constructor

//______________________________________________________________________________
XPreProcesSetting::XPreProcesSetting(const char *arraytype, const char *infile)
                  :XProcesSetting(arraytype, infile)
{
   // Normal PreProcesSetting constructor
   if(kCS) cout << "---XPreProcesSetting::XPreProcesSetting------" << endl;

   fSelector     = 0;
   fBgrdSelector = 0;
   fBackgrounder = 0;
   fNormSelector = 0;
   fNormalizer   = 0;
   fExprSelector = 0;
   fExpressor    = 0;
   fCallSelector = 0;
   fCaller       = 0;
}//Constructor

//______________________________________________________________________________
XPreProcesSetting::~XPreProcesSetting()
{
   // PreProcesSetting destructor
   if(kCS) cout << "---XPreProcesSetting::~XPreProcesSetting------" << endl;

   SafeDelete(fCaller);
   SafeDelete(fCallSelector);
   SafeDelete(fExpressor);
   SafeDelete(fExprSelector);
   SafeDelete(fNormalizer);
   SafeDelete(fNormSelector);
   SafeDelete(fBackgrounder);
   SafeDelete(fBgrdSelector);
   SafeDelete(fSelector);
}//Destructor

//______________________________________________________________________________
Int_t XPreProcesSetting::InitAlgorithm(const char *name, const char *type,
                         Option_t *options, const char *filename,
                         Int_t npars, Double_t *pars)
{
   // Initialize algorithm "name" with "type"
   if(kCS) cout << "------XPreProcesSetting::InitAlgorithm------" << endl;

   if (strcmp(name, "selector") == 0) {
      return this->InitSelector(type, options, npars, pars);
   } else if (strcmp(name, "backgrounder") == 0) {
      return this->InitBackgrounder(type, options, filename, npars, pars);
   } else if ((strcmp(name, "normalizer") == 0) &&
              (strcmp(type, "approx")     != 0)) {
      return this->InitNormalizer(type, options, filename, npars, pars);
   } else if ((strcmp(name, "normalizer") == 0) &&
              (strcmp(type, "approx")     == 0)) {
      return this->InitApprox(options, npars, pars);
   } else if (strcmp(name, "expressor") == 0) {
      return this->InitExpressor(type, options, filename, npars, pars);
   } else if (strcmp(name, "calldetector") == 0) {
      return this->InitCallDetector(type, options, npars, pars);
   } else {
      cerr << "Error: Algorithm <" << name << "> is not known." << endl;
   }//if

   return errInitSetting;
}//InitAlgorithm

//______________________________________________________________________________
Int_t XPreProcesSetting::InitSelector(const char *type, Option_t *options,
                         Int_t npars, Double_t *pars)
{
   // Initialize selector algorithm of type with options:
   // type = "default": default selector 
   //      - "none":    mask is not changed by selector
   //      - "all":     all probes are selected
   //    parameters are: numpars
   //    - numpars:  number of other parameters as integer, i.e. numpars = 0:
   //
   // type = "probe": probes are selected using following options:
   //      - "all":    all probes are selected
   //      - "pmonly": only PM probes for genes are selected
   //      - "mmonly": only MM probes for genes are selected
   //      - "both":   both PM and MM probes for genes are selected
   //      - "exon":   probes are selected according to following parameters
   //    parameters are: numpars
   //    - numpars:  number of other parameters as integer, 
   //                numpars = 0: all options except "exon"
   //                numpars = 1 or 2: for option is "exon" (pm or pm,mm type)
   //
   // type = "rank": selector based on difference of ranks, with options
   //      - "separate": use separate mask for each normalization
   //      - "together": use numflags masks for normalization
   //    parameters are: numpars, cutoff, trim, minunits, minflags
   //    - numpars:  number of other parameters as integer, i.e. numpars = 4:
   //    - cutoff:   cutoff value for selection of units, if = 0 then calculated
   //    - trim:     trim value for the fitting line to calculate cutoff
   //    - minunits: minimum number of units to be selected (as percentage)
   //    - numflags: number of flags, e.g. 9999 (=all) for intersect, 1 for union
   //
   // type = "user": user can either select "all" entries or supply a mask file
   //      for the selection of tree entries, with following options:
   //      - "all:separate" or "all:together", all units are used
   //      - "fullname.txt:separate" or "fullname.txt:together", where
   //        "fullname.txt" is the filename of the mask to import, and mask
   //         need only contain unitIDs and flag = 1 of selected units
   //    parameters are: numpars
   //    - numpars:  number of other parameters as integer, i.e. numpars = 0:
   if(kCS) cout << "------XPreProcesSetting::InitSelector------" << endl;

// Delete default selector setting
   SafeDelete(fSelector);

   TString exten = Type2Extension(type, kTypeSlct, kExtenSlct);
   TString stype = Extension2Type(type, kTypeSlct, kExtenSlct);

   if (strcmp(exten.Data(), kExtenSlct[0]) == 0) {
      fSelector = new XSelector(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenSlct[1]) == 0) {
      fSelector = new XRankSelector(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenSlct[2]) == 0) {
      fSelector = new XProbeSelector(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenSlct[3]) == 0) {
      fSelector = new XUnitSelector(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenSlct[4]) == 0) {
      fSelector = new XUserSelector(stype.Data(), exten.Data());
   } else {
      cerr << "Error: Selector <" << type << "> is not known." << endl;
      return errInitSetting;
   }//if
   if (fSelector == 0) return errInitMemory;

   fSelector->SetOptions(options);

   return fSelector->InitParameters(npars, pars);
}//InitSelector

//______________________________________________________________________________
Int_t XPreProcesSetting::InitBackgrounder(const char *type, Option_t *options,
                         const char *filename, Int_t npars, Double_t *pars)
{
   // Initialize background algorithm and parameters
   //
   // type = "sector": constant bgrd is calculated for each sector separately,
   //      and subtracted from intensities according to options:
   //      - "subtractbg":  subtract bgrd from intensity - result can be negative
   //      - "correctbg":   correct bgrd with noise fraction to avoid negative results
   //      - "attenuatebg": use generalized log-transform to avoid negative results
   //    parameters are: numpars, pcntcells, secrows, seccols, smoothiter, (nfrac, l, h)
   //    - numpars:    number of other parameters as integer, i.e. numpars = 4-6:
   //    - pcntcells:  precent of cells with lowest intensities in each sector
   //    - secrows:    number of sector rows
   //    - seccols:    number of sector columns
   //    - smoothiter: number of iterations used for smoothing of sector bgrd values
   //    - nfrac:      noise fraction for bgrd option "correctbg", or
   //    - l:          tunable parameter, 0<=l<=1 (default is 0.005) for "attenuatebg", and     
   //    - h:          parameter (default is -1) for "attenuatebg"
   //
   // type = "weightedsector": first, a constant bgrd is calculated for each sector
   //      separately, then a weight is calculated and each bgrd value multiplied
   //      with the weight to smooth transition between sector, and finally the
   //      smoothed bgrd is subtracted from intensities according to options:
   //      - "subtractbg":  subtract bgrd from intensity - result can be negative
   //      - "correctbg":   correct bgrd with noise fraction to avoid negative results
   //      - "attenuatebg": use generalized log-transform to avoid negative results
   //    parameters are: numpars, pcntcells, secrows, seccols, smoothiter, smooth, (nfrac, l, h)
   //    - numpars:    number of other parameters as integer, i.e. numpars = 5-7:
   //    - pcntcells:  precent of cells with lowest intensities in each sector
   //    - secrows:    number of sector rows
   //    - seccols:    number of sector columns
   //    - smoothiter: number of iterations used for smoothing of sector bgrd values
   //    - smooth:     smoothing parameter added to avoid infinite weights
   //    - nfrac:      noise fraction for bgrd option "correctbg", or
   //    - l:          optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:          optional parameter (default is -1)
   //
   // type = "rma":    bgrd is calculated using a global model for the distribution
   //      of probe intensities using the following options:
   //      - "pmonly": only PM probes for genes are used for bgrd calculation ("rma")
   //      - "mmonly": only MM probes for genes are used for bgrd calculation
   //      - "both":   both PM and MM probes are used for bgrd calculation ("rma2")
   //      Note: The default kernel is "epanechnikov"; if a different kernel is
   //            used than options must be e.g. "pmonly:gaussian"
   //    parameters are: numpars, numpoints
   //    - numpars:   number of other parameters as integer, i.e. numpars = 1:
   //    - numpoints: size of density array, preferrable power of 2, e.g. 16384
   //
   // type = "gccontent": bgrd is calculated based on the median intensity of probes
   //      with similar GC content.
   //    parameters are: numpars, trim, (nfrac, l, h)
   //    - trim:  trim value for trimmed mean (default is median, i.e. trim=0.5)
   //    - nfrac: noise fraction for bgrd option "correctbg", or
   //    - l:     optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:     optional parameter (default is -1)
   //
   // filename = "":    all trees will be written to main file 
   //          - "tmp": optional filename to create temporary file "tmp_exten"
   if(kCS) cout << "------XPreProcesSetting::InitBackgrounder------" << endl;

// Init (default) background selector
   Int_t err = errNoErr;
   if (!fSelector) err = this->InitSelector("probe", "none", 0, 0);
   if (err != errNoErr) return err;
   fBgrdSelector = fSelector;
   fSelector = 0; //free temporary selector
//   fBgrdSelector = new XProbeSelector(*(XProbeSelector*)fSelector); //need to if(exten==)
//   SafeDelete(fSelector); //not necessary ? done in InitSelector()

// Delete default background setting
   SafeDelete(fBackgrounder);

   TString exten = Type2Extension(type, kTypeBgrd, kExtenBgrd);
   TString stype = Extension2Type(type, kTypeBgrd, kExtenBgrd);

   if (strcmp(exten.Data(), kExtenBgrd[0]) == 0) {
      fBackgrounder = new XSectorBackground(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenBgrd[1]) == 0) {
      fBackgrounder = new XWeightedBackground(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenBgrd[2]) == 0) {
      fBackgrounder = new XRMABackground(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenBgrd[3]) == 0) {
      fBackgrounder = new XGCBackground(stype.Data(), exten.Data());
   } else {
      cerr << "Error: Backgrounder <" << type << "> is not known." << endl;
      return errInitSetting;
   }//if
   if (fBackgrounder == 0) return errInitMemory;

   fBackgrounder->SetOptions(options);
   fBackgrounder->NewFile(filename, exten.Data());

   return fBackgrounder->InitParameters(npars, pars);
}//InitBackgrounder

//______________________________________________________________________________
Int_t XPreProcesSetting::InitNormalizer(const char *type, Option_t *options,
                         const char *filename, Int_t npars, Double_t *pars)
{
   // Initialize normalizer algorithm and parameters
   // All normalization algorithms have the following option settings:
   //    "logbase" or "option:logbase" or "option:bgrdoption:logbase" with:
   //    option;
   //    - "transcript": use unit tree for transcripts (default)
   //    - "exon":       use unit tree for exons (exon array only)
   //    - "probeset":   use unit tree for probesets
   //    bgrdoption;
   //    - "none":        no background subtraction
   //    - "subtractbg":  subtract bgrd from intensity - result can be negative
   //    - "correctbg":   correct bgrd with noise fraction to avoid negative results
   //    - "attenuatebg": use generalized log-transform to avoid negative results
   //    logbase:
   //    - "0" - linear,
   //    - "log", "log2", "log10" - log with base e, 2, 10
   //
   // type = "mean": scale trimmed mean of data with options "option:logbase"
   //      - "all":     use all data to calculate trimmed mean
   //      - "sel":     use data selected with selector to calculate trimmed mean
   //      - logbase:    "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //    parameters are: numpars, trim, targetinten
   //    - numpars:     number of other parameters as integer, i.e. numpars = 2:
   //    - trim:        trim value for mean, in range [0, 0.5]
   //    - targetinten: if > 0 then scale trimmed mean to value of target intensity
   //
   // type = "median": scale median of data with options "option:logbase"
   //      - "all":     use all data to calculate median
   //      - "sel":     use data selected with selector to calculate median
   //      - logbase:    "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //    parameters are: numpars,  targetinten
   //    - numpars:     number of other parameters as integer, i.e. numpars = 1:
   //    - targetinten: if > 0 then scale median to value of target intensity
   //
   // type = "lowess": apply Lowess smoother with options:
   //      - logbase:    "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //    parameters are: numpars,  span, iter
   //    - numpars: number of other parameters as integer, i.e. numpars = 2:
   //    - span:    smoother span
   //    - iter:    number of robustifying iterations
   //
   // type = "supsmu": apply super smoother with options:
   //      - logbase:    "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //    parameters are: numpars,  bass, span
   //    - numpars: number of other parameters as integer, i.e. numpars = 2:
   //    - bass:    controls smoothness of fitted curve
   //    - span:    fraction of observations in span
   //
   // type = "quantile": apply quantile normalization with options:
   //      - "together": quantile will be applied to PM and MM together
   //      - "separate": quantile will be applied to PM and MM separately
   //      - logbase:    "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //      Note: if background should be subtracted first, then options must be
   //            "option:bgrdoption:logbase" where bgrdoption is:
   //      - "subtractbg": subtract bgrd from intensity - result can be negative
   //      - "correctbg":  correct bgrd with noise fraction to avoid negative results
   //    parameters are: numpars, trim, (nfrac, l, h)
   //    - numpars: number of other parameters as integer, i.e. numpars = 1-3:
   //    - trim:    trim value for mean, in range [0, 0.5]
   //    - nfrac:   noise fraction for bgrd option "correctbg", or
   //    - l:       optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:       optional parameter (default is -1)
   //
   // filename = "":    all trees will be written to main file; exception:
   //                   for "quantile" tmp_rkq.root will be created in current dir
   //          - "tmp": optional filename to create temporary file "tmp_exten"
   if(kCS) cout << "------XPreProcesSetting::InitNormalizer------" << endl;

// Init (default) normalization selector
   Int_t err = errNoErr;
   if (!fSelector) err = this->InitSelector("probe","both", 0, 0);
   if (err != errNoErr) return err;
   fNormSelector = fSelector;
   fSelector = 0; //free temporary selector
//   fNormSelector = new XProbeSelector(*(XProbeSelector*)fSelector); //need to if(exten==)
//   SafeDelete(fSelector);

// Delete default normalizer setting
   SafeDelete(fNormalizer);

   TString exten = Type2Extension(type, kTypeCNrm, kExtenCNrm);
   TString stype = Extension2Type(type, kTypeCNrm, kExtenCNrm);

   if (strcmp(exten.Data(), kExtenCNrm[4]) == 0) {
      fNormalizer = new XQuantileNormalizer(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenCNrm[0]) == 0) {
      fNormalizer = new XMeanNormalizer(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenCNrm[1]) == 0) {
      fNormalizer = new XMedianNormalizer(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenCNrm[2]) == 0) {
      fNormalizer = new XLowessNormalizer(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenCNrm[3]) == 0) {
      fNormalizer = new XSuperNormalizer(stype.Data(), exten.Data());
   } else {
      cerr << "Error: Normalizer <" << type << "> is not known/applicable." << endl;
      return errInitSetting;
   }//if
   if (fNormalizer == 0) return errInitMemory;

   fNormalizer->SetOptions(options);
   fNormalizer->NewFile(filename, exten.Data());

   return fNormalizer->InitParameters(npars, pars);
}//InitNormalizer

//______________________________________________________________________________
Int_t XPreProcesSetting::InitApprox(Option_t *options, Int_t npar, Double_t *pars)
{
   // Initialize normalizer algorithm with type = "approx" and parameters
   // and with options of the form "method:ties"
   //    - method: "linear" or "constant" (interpolation method used)
   //    - ties:   "mean", "min", "max", "ordered" (handling of tied values)
   //    parameters are: numpars, rule, f
   //    - numpars: number of other parameters as integer, i.e. numpars = 2:
   //    - rule:    how interpolation should be done at ends
   //    - f:       for method "constant" a value between 0 and 1
   // Note: approx must only be initialized for normalizers "lowess" and "supsmu".
   if(kCS) cout << "------XPreProcesSetting::InitApprox------" << endl;

   if (fNormalizer == 0) {
      cerr << "Error: Need to initialize Normalizer first" << endl;
      return errGeneral;
   }//if

   TString optcpy = options;
   TString method = "linear";
   TString ties   = "mean";

	char *opt = (char*)optcpy.Data();
   if (NumSeparators(options, ":") == 1) {
      method = strtok(opt, ":");
      ties   = strtok(NULL,":");
   } else {
      cout << "Warning: InitAlgorithm() must have two options for approx." << endl;
      cout << "         Using default options." << endl;
   }//if

   return fNormalizer->InitApprox(method.Data(), ties.Data(), npar, pars);
}//InitApprox

//______________________________________________________________________________
Int_t XPreProcesSetting::InitExpressor(const char *type, Option_t *options,
                         const char *filename, Int_t npars, Double_t *pars)
{
   // Initialize expression algorithm and parameters
   // All expression algorithms have the following option settings:
   //    "logbase" or "option:logbase" or "option:bgrdoption:logbase" with:
   //    option;
   //    - "transcript": use unit tree for transcripts (default)
   //    - "exon":       use unit tree for exons (exon array only)
   //    - "probeset":   use unit tree for probesets
   //    bgrdoption;
   //    - "none":        no background subtraction
   //    - "subtractbg":  subtract bgrd from intensity - result can be negative
   //    - "correctbg":   correct bgrd with noise fraction to avoid negative results
   //    - "attenuatebg": use generalized log-transform to avoid negative results
   //    logbase:
   //    - "0" - linear,
   //    - "log", "log2", "log10" - log with base e, 2, 10
   //
   // type = "ArithmeticMean": arithmetic mean,
   // type = "GeometricMean": geometric mean, each with parameters:
   //    parameters are: numpars, trim, (nfrac, l, h)
   //    - numpars: number of other parameters as integer, i.e. numpars = 1-3:
   //    - trim:    trim value for mean, in range [0, 0.5]
   //    - nfrac:   noise fraction for bgrd option "correctbg", or
   //    - l:       optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:       optional parameter (default is -1) for "attenuatebg"
   //
   // type = "WeightedMean": weighted mean, 
   // type = "WeightedDiff": weighted difference, each with parameters:
   //    parameters are: numpars, wmaxinten, (nfrac, l, h)
   //    - numpars:   number of other parameters as integer, i.e. numpars = 1-3:
   //    - wmaxinten: weight for maximal (=saturation) intensity
   //    - nfrac:   noise fraction for bgrd option "correctbg", or
   //    - l:       optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:       optional parameter (default is -1) for "attenuatebg"
   //
   // type = "GCcorrectedMean": GC-corrected mean, with parameters:
   //    Note: not yet implemented
   //
   // type = "AvgDiff": average difference of MAS4, with parameters:
   //    parameters are: numpars, STP, (nfrac, l, h)
   //    - numpars: number of other parameters as integer, i.e. numpars = 1-3:
   //    - STP:     default of MAS4 is STP = 3
   //    - nfrac:   noise fraction for bgrd option "correctbg", or
   //    - l:       optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:       optional parameter (default is -1) for "attenuatebg"
   //
   // type = "TukeyBiweight": Tukey biweight of MAS5, as described in the
   //        sadd_whitepaper.pdf from Affymetrix, with parameters:
   //    parameters are: numpars, tau, scaletau, delta, c, eps, neglog, (nfrac, l, h)
   //    - numpars:  number of other parameters as integer, i.e. numpars = 6-8:
   //    - tau:      contrast tau, default is 0.03
   //    - scaletau: scale tau, default is 10
   //    - delta:    delta for probe value calculation, default is 2.0e-20
   //    - c:        a tuning constant, default is 5
   //    - eps:      small value to avoid zeros in division, default is 0.0001
   //    - neglog:   substitution for logarithm of negative values
   //    - nfrac:   noise fraction for bgrd option "correctbg", or
   //    - l:       optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:       optional parameter (default is -1) for "attenuatebg"
   //
   // type = "MedianPolish": median polish for RMA, with parameters:
   //    parameters are: numpars, maxiter, eps, neglog, (nfrac, l, h)
   //    - numpars: number of other parameters as integer, i.e. numpars = 3-5:
   //    - maxiter: maximal number of iterations, default is 10
   //    - eps:     epsilon of test for convergence, default is 0.01
   //    - neglog:  substitution for logarithm of negative values
   //    - nfrac:   noise fraction for bgrd option "correctbg", or
   //    - l:       optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:       optional parameter (default is -1) for "attenuatebg"
   //
   //    filename = "": data for "medianpolish" will be stored as table in RAM  
   //             = "tmp": optional filename to create temporary file "tmp_exten",
   //               where data will be stored temporarily for "medianpolish"
   if(kCS) cout << "------XPreProcesSetting::InitExpressor------" << endl;

// Init (default) expression selector
   Int_t err = errNoErr;
   if (!fSelector) err = this->InitSelector("probe","none", 0, 0);
   if (err != errNoErr) return err;
   fExprSelector = fSelector;
   fSelector = 0; //free temporary selector
//   fExprSelector = new XProbeSelector(*(XProbeSelector*)fSelector); //need to if(exten==)
//   SafeDelete(fSelector);

// Delete default expression setting
   SafeDelete(fExpressor);

   TString exten = Type2Extension(type, kTypeExpr, kExtenExpr);
   TString stype = Extension2Type(type, kTypeExpr, kExtenExpr);

   if (strcmp(exten.Data(), kExtenExpr[0]) == 0) {
      fExpressor = new XArithmeticMean(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenExpr[1]) == 0) {
      fExpressor = new XGeometricMean(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenExpr[2]) == 0) {
      fExpressor = new XWeightedMean(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenExpr[3]) == 0) {
      fExpressor = new XGCCorrectedMean(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenExpr[4]) == 0) {
      fExpressor = new XWeightedDiff(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenExpr[5]) == 0) {
      fExpressor = new XAvgDif(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenExpr[6]) == 0) {
      fExpressor = new XTukeyBiweight(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenExpr[7]) == 0) {
      fExpressor = new XMedianPolish(stype.Data(), exten.Data());
      if (fExpressor) fExpressor->NewFile(filename, exten.Data());
   } else {
      cerr << "Error: Expressor <" << type << "> is not known." << endl;
      return errInitSetting;
   }//if
   if (fExpressor == 0) return errInitMemory;

   fExpressor->SetOptions(options);

   return fExpressor->InitParameters(npars, pars);
}//InitExpressor

//______________________________________________________________________________
Int_t XPreProcesSetting::InitCallDetector(const char *type, Option_t *options,
                         Int_t npars, Double_t *pars)
{
   // Initialize present call algorithm and parameters
   // All present call algorithms have the following option settings:
   //    "dataoption", "option:dataoption" or "option:dataoption:bgrdoption" with:
   //    option;
   //    - "transcript": use unit tree for transcripts (default)
   //    - "exon":       use unit tree for exons (exon array only)
   //    - "probeset":   use unit tree for probesets
   //   dataoption determines to which data algorithm should be applied:
   //    - "raw":        apply to raw data
   //    - "adjusted":   apply to background-adjusted data
   //    - "normalized": apply to normalized data
   //   bgrdoption determines if background should be subtracted:
   //    - "none":        no background subtraction (default)
   //    - "subtractbg":  subtract bgrd from intensity - result can be negative
   //    - "correctbg":   correct bgrd with noise fraction to avoid negative results
   //    - "attenuatebg": use generalized log-transform to avoid negative results
   //
   // type = "MeanDifferenceCall": mean difference call, with parameters:
   //    parameters are: numpars, trim, trim1, trim2, percent, (nfrac, l, h)
   //    - numpars: number of other parameters as integer, i.e. numpars = 4-6:
   //    - trim:    trim value for mean of PMs
   //    - trim1:   trim value for mean1 of PMs and mean1 of MMs
   //    - trim2:   trim value for mean2 of MMs and mean2 of MMs
   //    - percent: percent of means for PM1 or PM2 as test condition
   //    - nfrac:   noise fraction for bgrd option "correctbg", or
   //    - l:       optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:       optional parameter (default is -1)
   //
   // type = "DetectionCall": detection call of MAS5, with parameters:
   //    parameters are: numpars, tau, alpha1, alpha2, ignore, exact, correct, (nfrac, l, h)
   //    - numpars: number of other parameters as integer, i.e. numpars = 6-8:
   //    - tau:     threshold for discrimination score, default is 0.015
   //    - alpha1:  significance level between P and M call, default is 0.04
   //    - alpha2:  significance level between M and A call, default is 0.06
   //    - ignore:  if > 0 then ignore probes with saturated MMs
   //    - exact:   set to one to compute exact p-value for Wilcox test
   //    - correct: set to one to apply continuity correction for Wilcox test
   //    - nfrac:   noise fraction for bgrd option "correctbg", or
   //    - l:       optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:       optional parameter (default is -1)
   //
   // type = "Mas4Call": detection call of MAS4, with parameters:
   //    parameters are: numpars, ????
   //    NEED TO BE IMPLEMENTED!
   //
   // type = "DABGCall": detection above background call of APT, with parameters:
   //    parameters are: numpars, cut, alpha1, alpha2
   //    - numpars: number of other parameters as integer, i.e. numpars = 3:
   //    - cut:     cut>1 for Fisher's chi-squared method, 0<=cut<=1 for percentile
   //    - alpha1:  significance level between P and M call, default is 0.01 (t.b.d.)
   //    - alpha2:  significance level between M and A call, default is 0.015 (t.b.d.)
   if(kCS) cout << "------XPreProcesSetting::InitCallDetector------" << endl;

// Init (default) call selector
   Int_t err = errNoErr;
   if (!fSelector) err = this->InitSelector("probe","none", 0, 0);
   if (err != errNoErr) return err;
   fCallSelector = fSelector;
   fSelector = 0; //free temporary selector
//   fCallSelector = new XProbeSelector(*(XProbeSelector*)fSelector); //need to if(exten==)
//   SafeDelete(fSelector);

// Delete default present call setting
   SafeDelete(fCaller);

   TString exten = Type2Extension(type, kTypeCall, kExtenCall);
   TString stype = Extension2Type(type, kTypeCall, kExtenCall);

   if (strcmp(exten.Data(), kExtenCall[0]) == 0) {
      fCaller = new XMeanDifferenceCall(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenCall[1]) == 0) {
      fCaller = new XDetectionCall(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenCall[2]) == 0) {
      fCaller = new XMAS4Call(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenCall[3]) == 0) {
      fCaller = new XDABGCall(stype.Data(), exten.Data());
   } else {
      cerr << "Error: Call detector <" << type << "> is not known." << endl;
      return errInitSetting;
   }//if
   if (fCaller == 0) return errInitMemory;

   fCaller->SetOptions(options);

   return fCaller->InitParameters(npars, pars);
}//InitCallDetector


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreProcesSet                                                        //
//                                                                      //
// Base class for microarray analysis type                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XPreProcesSet::XPreProcesSet()
              :XProcesSet()
{
   // Default PreProcesSet constructor
   if(kCS) cout << "---XPreProcesSet::XPreProcesSet(default)------" << endl;

   fBgrdSelector = 0;
   fBackgrounder = 0;
   fNormSelector = 0;
   fNormalizer   = 0;
   fExprSelector = 0;
   fExpressor    = 0;
   fCallSelector = 0;
   fCaller       = 0;
}//Constructor

//______________________________________________________________________________
XPreProcesSet::XPreProcesSet(const char *name, const char *type)
              :XProcesSet(name, type)
{
   // Normal PreProcesSet constructor
   if(kCS) cout << "---XPreProcesSet::XPreProcesSet------" << endl;

   fBgrdSelector = 0;
   fBackgrounder = 0;
   fNormSelector = 0;
   fNormalizer   = 0;
   fExprSelector = 0;
   fExpressor    = 0;
   fCallSelector = 0;
   fCaller       = 0;
}//Constructor

//______________________________________________________________________________
XPreProcesSet::~XPreProcesSet()
{
   // PreProcesSet destructor
   if(kCS) cout << "---XPreProcesSet::~XPreProcesSet------" << endl;

   // will be deleted in destructor of XPreProcesSetting
   fBgrdSelector = 0;
   fBackgrounder = 0;
   fNormSelector = 0;
   fNormalizer   = 0;
   fExprSelector = 0;
   fExpressor    = 0;
   fCallSelector = 0;
   fCaller       = 0;
}//Destructor

//______________________________________________________________________________
Int_t XPreProcesSet::Initialize(TFile *file, XSetting *setting,
                     const char *infile, const char *treename)
{
   // Initialize algorithms
   if(kCS) cout << "------XPreProcesSet::Initialize------" << endl;

   if ((file == 0) || (setting == 0)) return errAbort;

   fFile     = file;
   fSetting  = (XPreProcesSetting*)setting;
   fInfile   = infile;
   fTreeName = treename;

// Get files from settings
   fDataFile   = ((XPreProcesSetting*)fSetting)->GetDataFile();
   fSchemeFile = ((XPreProcesSetting*)fSetting)->GetSchemeFile();
   fSchemeName = ((XPreProcesSetting*)fSetting)->GetSchemeName();

// Get algorithms from settings
   fBgrdSelector = ((XPreProcesSetting*)fSetting)->GetBgrdSelector();
   fBackgrounder = ((XPreProcesSetting*)fSetting)->GetBackgrounder();
   fNormSelector = ((XPreProcesSetting*)fSetting)->GetNormSelector();
   fNormalizer   = ((XPreProcesSetting*)fSetting)->GetNormalizer();
   fExprSelector = ((XPreProcesSetting*)fSetting)->GetExprSelector();
   fExpressor    = ((XPreProcesSetting*)fSetting)->GetExpressor();
   fCallSelector = ((XPreProcesSetting*)fSetting)->GetCallSelector();
   fCaller       = ((XPreProcesSetting*)fSetting)->GetCallDetector();

   return errNoErr;
}//Initialize

//______________________________________________________________________________
void XPreProcesSet::AddDataTreeInfo(TTree *tree, const char *name, Option_t *option,
                    Int_t nrows, Int_t ncols, Int_t nmin, Double_t min, Int_t nmax,
                    Double_t max, Int_t maxnpix)
{
   // Add background/intensity tree info to list fUserInfo of tree
   if(kCS) cout << "------XPreProcesSet::AddDataTreeInfo------" << endl;

   XDataTreeInfo *info = new XDataTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info
   info->AddUserInfo(nrows, ncols, nmin, min, nmax, max, maxnpix);

   tree->GetUserInfo()->Add(info);
}//AddDataTreeInfo

//______________________________________________________________________________
void XPreProcesSet::AddMaskTreeInfo(TTree *tree, const char *name, Option_t *option,
                    Int_t nrows, Int_t ncols, Int_t nflags)
{
   // Add mask tree info to list fUserInfo of tree
   if(kCS) cout << "------XPreProcesSet::AddMaskTreeInfo------" << endl;

// store name of tree set as title
   XMaskTreeInfo *info = new XMaskTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info
   info->AddUserInfo(nrows, ncols, nflags);

   tree->GetUserInfo()->Add(info);
}//AddMaskTreeInfo

//______________________________________________________________________________
void XPreProcesSet::AddCallTreeInfo(TTree *tree, const char *name, Option_t *option,
                    Int_t nunits, Int_t nabsent, Int_t nmarginal, Int_t npresent,
                    Double_t minpval, Double_t maxpval)
{
   // Add detection call tree info to list fUserInfo of tree
   if(kCS) cout << "------XPreProcesSet::AddCallTreeInfo------" << endl;

// store name of tree set as title
   XCallTreeInfo *info = new XCallTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info
   info->AddUserInfo(nunits, nabsent, nmarginal, npresent, minpval, maxpval);

   tree->GetUserInfo()->Add(info);
}//AddCallTreeInfo


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCProcesSet                                                         //
//                                                                      //
// Class for GeneChip oligonucleotide array pre-processing              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGCProcesSet::XGCProcesSet()
             :XPreProcesSet()
{
   // Default GeneChip Processing constructor
   if(kCS) cout << "---XGCProcesSet::XGCProcesSet(default)------" << endl;

   fNBgPar = 0;
   fBgPars = 0;
}//Constructor

//______________________________________________________________________________
XGCProcesSet::XGCProcesSet(const char *name, const char *type)
             :XPreProcesSet(name, type)
{
   // Normal GeneChip Processing constructor
   if(kCS) cout << "---XGCProcesSet::XGCProcesSet------" << endl;

   fNBgPar = 0;
   fBgPars = 0;
}//Constructor

//______________________________________________________________________________
XGCProcesSet::~XGCProcesSet()
{
   // GeneChip Processing destructor
   if(kCS) cout << "---XGCProcesSet::~XGCProcesSet------" << endl;

   if (fBgPars) {delete [] fBgPars; fBgPars = 0;}
   fNBgPar = 0;
}//Destructor

//______________________________________________________________________________
Int_t XGCProcesSet::Preprocess(const char *method)
{
   // Preprocess microarray data
   if(kCS) cout << "------XGCProcesSet::Preprocess------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "Preprocessing data using method <" << method << ">..." << endl;
   }//if

   Int_t err = errNoErr;

// Test if algorithms are initialized
   if (!(fBackgrounder || fNormalizer || fExpressor || fCaller)) {
      cerr << "Error: At least one algorithm must be initialized!" << endl;
      return errAbort;
   }//if

// Flags for algorithms to be processed
   Bool_t doAll  = (strcmp(method, kProcessMethod[0]) == 0);
   Bool_t doRMA  = (strcmp(method, kProcessMethod[1]) == 0);
   Bool_t doMAS4 = (strcmp(method, kProcessMethod[2]) == 0);
   Bool_t doMAS5 = (strcmp(method, kProcessMethod[3]) == 0);
   Bool_t doBgrd = (strcmp(method, kProcessMethod[4]) == 0) || doRMA || doMAS4 || doMAS5;
   Bool_t doNorm = (strcmp(method, kProcessMethod[5]) == 0) || doRMA;
   Bool_t doExpr = (strcmp(method, kProcessMethod[6]) == 0) || doRMA || doMAS4 || doMAS5;
   Bool_t doCall = (strcmp(method, kProcessMethod[7]) == 0) ||          doMAS4 || doMAS5;
//??   Bool_t doCall = (strcmp(method, kProcessMethod[7]) == 0) || doRMA || doMAS4 || doMAS5;

// Flags to determine kind of datatree to be used for detection call
   Bool_t atRaw  = kTRUE;
   Bool_t atBgrd = kFALSE;
   Bool_t atNorm = kFALSE;
   if (fCaller) {
      atRaw  = (strcmp(fCaller->GetDataOption(), kCallOption[0]) == 0);
      atBgrd = (strcmp(fCaller->GetDataOption(), kCallOption[1]) == 0);
      atNorm = (strcmp(fCaller->GetDataOption(), kCallOption[2]) == 0);

      // set default to raw data
      if (!(atRaw || atBgrd || atNorm)) atRaw = kTRUE;
   }//if

// Get schemes from scheme file
   fSchemeFile->cd();
   fSchemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!fSchemes) {
      return fManager->HandleError(errMissingContent, "Scheme", kContent);
   }//if

// Get data from data file
   fDataFile->cd();
   fData = (XFolder*)(fDataFile->Get(kContent));
   if (!fData) {
      return fManager->HandleError(errMissingContent, "Data", kContent);
   }//if

// Get number of data trees and background trees
   Int_t numdata  = 0;
   Int_t numbgrd  = 0;
   Int_t numtrees = fTrees->GetSize();
//or?   Int_t numtrees = fSelections->GetSize();

   for (Int_t k=0; k<numtrees; k++) {
      TTree* tree = (TTree*)fTrees->At(k);
//or?      TTree* tree = (TTree*)fSelections->At(k);
      if      ((tree->GetBranch("DataBranch")) != 0) numdata++;
      else if ((tree->GetBranch("BgrdBranch")) != 0) numbgrd++;
   }//for_k

// Need to add equal number of bgrd and data trees
   if ((numbgrd > 0) && (numbgrd != numdata)) {
      cerr << "Error: Number of background trees <" << numbgrd
           << "> is not equal to number of data trees <" << numdata
           << ">!" << endl;
      return errAbort;
   }//if

// Create hash table to store datatree names
   TString     tmpstr;
   XIdxString *idxstr = 0;
   THashTable *htable = 0;
   if (!(htable = new THashTable(2*numdata))) return errInitMemory;

// Initialize data trees and background trees
   TTree *datatree[numdata+1]; //for possible Reference tree for normalization
   TTree *bgrdtree[numdata];
   TTree *temptree[numdata];
   for (Int_t k=0; k<numdata; k++) datatree[k] = bgrdtree[k] = temptree[k] = 0;

   numdata = 0;
   numbgrd = 0;
   for (Int_t k=0; k<numtrees; k++) {
      TTree* tree = (TTree*)fTrees->At(k);
//or?      TTree* tree = (TTree*)fSelections->At(k);
      if ((tree->GetBranch("DataBranch")) != 0) {
         datatree[numdata] = tree;
         tmpstr = Path2Name(datatree[numdata]->GetName(), "", ".");
         idxstr = new XIdxString(numdata, tmpstr.Data());
         htable->Add(idxstr);
         numdata++;
      } else if ((tree->GetBranch("BgrdBranch")) != 0) {
         temptree[numbgrd++] = tree;
      }//if
   }//for_k

// Sort bgrdtrees in order of datatrees
   if (numbgrd > 0) {
      numbgrd = 0;
      for (Int_t k=0; k<numdata; k++) {
         tmpstr = Path2Name(temptree[k]->GetName(), "", ".");
         idxstr = (XIdxString*)(htable->FindObject(tmpstr.Data()));
         if (idxstr) {
            bgrdtree[idxstr->GetIndex()] = temptree[k];
            numbgrd++;
         }//if
      }//for_k
   }//if

   // delete idxstr before deleting htable
   if (htable) {htable->Delete(); delete htable; htable = 0;}

// Check again for equal number of bgrd trees and corresponding data trees
   if ((numbgrd > 0) && (numbgrd != numdata)) {
      cerr << "Error: <" << numdata-numbgrd 
           << "> data trees have no corresponding background tree!" << endl;
      err = errAbort; goto cleanup;
   }//if

// Detect present call based on raw datatree
   if (fCaller && fCallSelector && atRaw && (doCall || doAll)) {
      err = this->DetectCall(numdata, datatree, numbgrd, bgrdtree);
      if (err != errNoErr) goto cleanup;
   } else if (doCall && !fCaller) {
      cerr << "Error: CallDetector algorithm is not initialized!" << endl;
      err = errAbort; goto cleanup;
   }//if

// Adjust background
   if (fBackgrounder && fBgrdSelector && (doBgrd || doAll)) {
      err = this->AdjustBackground(numdata, datatree, numbgrd, bgrdtree);
      if (err != errNoErr) goto cleanup;
   } else if (doBgrd && !fBackgrounder) {
      cerr << "Error: Backgrounder algorithm is not initialized!" << endl;
      err = errAbort; goto cleanup;
   }//if

// Detect present call based on background-adjusted datatree
   if (fCaller && fCallSelector && atBgrd && (doCall || doAll)) {
      err = this->DetectCall(numdata, datatree, numbgrd, bgrdtree);
      if (err != errNoErr) goto cleanup;
   } else if (doCall && !fCaller) {
      cerr << "Error: CallDetector algorithm is not initialized!" << endl;
      err = errAbort; goto cleanup;
   }//if

// Normalize
   if (fNormalizer && fNormSelector && (doNorm || doAll)) {
//? if () cerr << "Error: At least two trees need to be selected." << endl;
      err = this->Normalize(numdata, datatree, numbgrd, bgrdtree);
      if (err != errNoErr) goto cleanup;
   } else if (doNorm && !fNormalizer) {
      cerr << "Error: Normalizer or selector algorithm is not initialized!" << endl;
      err = errAbort; goto cleanup;
   }//if

// Detect present call based on normalized datatree
   if (fCaller && fCallSelector && atNorm && (doCall || doAll)) {
      err = this->DetectCall(numdata, datatree, numbgrd, bgrdtree);
      if (err != errNoErr) goto cleanup;
   } else if (doCall && !fCaller) {
      cerr << "Error: CallDetector algorithm is not initialized!" << endl;
      err = errAbort; goto cleanup;
   }//if

// Condense
   if (fExpressor && fExprSelector && (doExpr || doAll)) {
//? if () cerr << "Error: At least two trees need to be selected." << endl;
      err = this->Express(numdata, datatree, numbgrd, bgrdtree);
      if (err != errNoErr) goto cleanup;
   } else if (doExpr && !fExpressor) {
      cerr << "Error: Expressor algorithm is not initialized!" << endl;
      err = errAbort; goto cleanup;
   }//if

// Cleanup
cleanup:
   SafeDelete(fSchemes);
   SafeDelete(fData);

   return err;
}//Preprocess

//______________________________________________________________________________
Int_t XGCProcesSet::AdjustBackground(Int_t numdata, TTree **datatree,
                    Int_t &numbgrd, TTree **bgrdtree)
{
   // Adjust background, i.e. calculate background and subtact from intensity.
   // All datatrees will be replaced with bg-corrected datatrees
   // Note: each datatree is processed individually and thus can have different
   //       number of entries (i.e. belong to different chip types)
   // Note: Parameter numbgrd will be set to numbgrd=0 to prevent further use
   //       of externally added bgrdtrees, since datatrees are already corrected!
   if(kCS) cout << "------XGCProcesSet::AdjustBackground------" << endl;

   Int_t err   = errNoErr;
   Int_t split = 99;

   TDirectory *savedir = gDirectory;
   TFile      *tmpfile = fBackgrounder->GetFile();

// Check for background trees
   if (numbgrd > 0) {
      cout << "Warning: Externally added background trees will be ignored."
           << endl;
      // to prevent bg-subtraction from bg-corrected data trees!
      numbgrd = 0;
   }//if

// Extension for background-corrected datatrees
   TString exten = kExtenIntn[0];

// Init local arrays to store data from trees
   Int_t    *arrMask  = 0;  //mask for selection of probes to be used
   Double_t *arrInten = 0;  //initial/bg-corrected intensity of probes
   Double_t *arrStdev = 0;  //standard deviation of probe intensities
   Double_t *arrBgrd  = 0;  //background signal of probes
   Double_t *arrNoise = 0;  //standard deviation of probe background signal

// Calculate background
   Int_t i, j, ij, x, y;
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) {err = errGetTree; break;}

   // Informing user
      if (XManager::fgVerbose) {
         cout << "   Calculating background for <" << datatree[k]->GetName() << ">..."
              << endl;
      }//if

   // Init data tree
      XGCCell *gccell = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell);

   // Get chip parameters from scheme file (also for alternative CDFs)
      if (strcmp(fSchemeName.Data(), "") == 0) {
         fSchemeName = datatree[k]->GetTitle();
      } else if (!fSchemeName.Contains(datatree[k]->GetTitle())) {
         return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[k]->GetTitle());
      }//if
      if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

      XDNAChip *chip = (XDNAChip*)fSchemes->FindObject(fSchemeName, kTRUE);
      if (!chip) return fManager->HandleError(errGetScheme, fSchemeName);
      Int_t numrows = chip->GetNumRows();
      Int_t numcols = chip->GetNumColumns();

   // Init counters for min/max values
      Double_t min    = DBL_MAX;  //defined in float.h
      Double_t max    = 0;
      Int_t    nummin = 0;
      Int_t    nummax = 0;

   // Get exon level of annotation
      Int_t level = 0;
      if (fBgrdSelector->GetNumParameters() > 0) {
         level = (Int_t)(fBgrdSelector->GetParameters())[0];
      }//if

   // Initialize memory for mask array
      Int_t size = numrows*numcols;
      if (!(arrMask  = new (nothrow) Int_t[size])) return errInitMemory;
      //not =0,since msk=1 and msk=0 must be determined!
      for (i=0; i<size; i++) arrMask[i] = eINITMASK;

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

   // Get mask for PM/MM from scheme tree and store in array 
      err = this->SchemeMask(chip, level, size, arrMask);
      if (err != errNoErr) return err;

   // Calculate mask for background
      err = fBgrdSelector->Calculate(size, 0, 0, arrMask);
      if (err != errNoErr) return err;

   // For GCBackground only, fill arrMask with GC content:
   // for PM set GC >= 0 , i.e. GC = 0...kProbeLength,
   // for MM set GC < eINITMASK, i.e. GC = eINITMASK - (1...kProbeLength+1)
      if (strcmp(fBackgrounder->GetName(), "gccontent") == 0) {
      // Get probe tree for scheme
         XGCProbe *probe = 0;
         TTree *probetree = (TTree*)gDirectory->Get(chip->GetProbeTree()); 
         if (probetree == 0) return errGetTree;
         probetree->SetBranchAddress("PrbBranch", &probe);

      // Get GC content from probe tree and store in array
         for (i=0; i<size; i++) {
            probetree->GetEntry(i);

            x  = probe->GetX();
            y  = probe->GetY();
            ij = XY2Index(x, y, numcols);

            if (arrMask[ij] == 1) {
               arrMask[ij] = probe->GetNumberGC();
            } else if (arrMask[ij] == 0) {
               //need to use (numberGC + 1) to avoid setting arrMask=eINITMASK for GC=0!!
               arrMask[ij] = eINITMASK - (probe->GetNumberGC() + 1);
//no!!               arrMask[ij] = eINITMASK - probe->GetNumberGC();
            }//if
         }//for_i
      }//if

   // Change directory
      if (tmpfile != 0) tmpfile->cd();
      else if (!fFile->cd(fName)) return errGetDir;

   // Create new tree bgrdtree
      TString dataname = Path2Name(datatree[k]->GetName(),"/",".");
      TString bgrdname = dataname + "." + fBackgrounder->GetTitle();
      bgrdtree[k] = new TTree(bgrdname, fSchemeName);
      if (bgrdtree[k] == 0) return errCreateTree;
      XBgCell *bgcell = new XBgCell();
      bgrdtree[k]->Branch("BgrdBranch", "XBgCell", &bgcell, 64000, split);

   // Initialize memory for data arrays
      if (!(arrInten = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrStdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrBgrd  = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrNoise = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      for (i=0; i<size; i++) {
         arrInten[i] = arrStdev[i] = arrBgrd[i] = arrNoise[i] = 0.0;
      }//for_i

   // Get data from data tree and store in array
      for (i=0; i<size; i++) {
         datatree[k]->GetEntry(i);

         x  = gccell->GetX();
         y  = gccell->GetY();
         ij = XY2Index(x, y, numcols);

         arrInten[ij] = gccell->GetIntensity();
         arrStdev[ij] = gccell->GetStdev();
      }//for_i

   // Calculate background and stdev
      fBackgrounder->SetNumRows(numrows);
      fBackgrounder->SetNumCols(numcols);
      err = fBackgrounder->Calculate(size, arrInten, arrBgrd, arrNoise, arrMask);
      if (err != errNoErr) goto cleanup;

   // Fill background tree 
      for (j=0; j<numrows; j++) {
         for (i=0; i<numcols; i++) {
            ij = XY2Index(i, j, numcols);

            // number of cells with minimal intensity
            if      (arrBgrd[ij] <  min) {min = arrBgrd[ij]; nummin = 1;}
            else if (arrBgrd[ij] == min) {nummin++;}

            // number of cells with maximal intensity
            if      (arrBgrd[ij] >  max) {max = arrBgrd[ij]; nummax = 1;}
            else if (arrBgrd[ij] == max) {nummax++;}

            bgcell->SetX(i);
            bgcell->SetY(j);
            bgcell->SetBackground(arrBgrd[ij]);
            bgcell->SetStdev(arrNoise[ij]);
            bgrdtree[k]->Fill(); 
         }//for_j
      }//for_i

      if (XManager::fgVerbose) {
         cout << "      background statistics: " << endl;
         cout << "         " << nummin << " cells with minimal intensity " << min << endl;
         cout << "         " << nummax << " cells with maximal intensity " << max << endl;
      }//if

   // Add tree info to tree
      AddDataTreeInfo(bgrdtree[k], bgrdtree[k]->GetName(), fBackgrounder->GetOption(),
                      numrows, numcols, nummin, min, nummax, max, -1);

   // Write background tree to file 
      if ((err = WriteTree(bgrdtree[k], TObject::kOverwrite)) != errNoErr) goto cleanup;
//??      if ((err = WriteTree(bgrdtree[k], 0)) != errNoErr) goto cleanup;

   // Add background tree header to list
//tmp      Int_t useTmpFile = fBackgrounder->UseTemporaryFile(); //??
//tmp      if useTmpFile then fSelections->Add()??? else fHeaders->Add() ??
//tmp      AddTreeHeader(bgrdtree[k]->GetName(), "Bgrd", useTmpFile, fBackgrounder->GetNumParameters(),
//////////////////////
//tmp or?:  if (tmpfile == 0) AddTreeHeader(...); //do NOT add to fSelections->Add()??
//or?:  if (tmpfile == 0) treeid = 0 else treeid = 1; //do NOt add to fHeaders->Add()??
//PROBLEM: for tmpfile!=0  fSize of fHeaders too large!!!
//////////////////////
      AddTreeHeader(bgrdtree[k]->GetName(), "Bgrd", 0, fBackgrounder->GetNumParameters(),
                    fBackgrounder->GetParameters());

   // Write background-corrected data tree to file 
      arrStdev = fBackgrounder->AdjustError(size, arrStdev, arrNoise);
      datatree[k] = this->FillDataTree(datatree[k], exten, fBackgrounder,
                          numrows, numcols, arrInten, arrStdev);
      if (datatree[k] == 0) {err = errCreateTree; goto cleanup;}

   cleanup:
      // delete arrays
      if (arrNoise) {delete [] arrNoise; arrNoise = 0;}
      if (arrBgrd)  {delete [] arrBgrd;  arrBgrd  = 0;}
      if (arrStdev) {delete [] arrStdev; arrStdev = 0;}
      if (arrInten) {delete [] arrInten; arrInten = 0;}
      if (arrMask)  {delete [] arrMask;  arrMask  = 0;}
//?      // delete scheme tree from RAM
//?      if (scmtree)  {scmtree->Delete(""); scmtree = 0;}
      // note: keep datatree[k] and bgrdtree[k] for later use!

      if (err != errNoErr) break;
   }//for_k

   savedir->cd();

   return err;
}//AdjustBackground

//______________________________________________________________________________
Int_t XGCProcesSet::Normalize(Int_t numdata, TTree **datatree,
                    Int_t &numbgrd, TTree **bgrdtree)
{
   // Normalize data trees after subtraction of background from corresponding
   // background trees (if present).
   // Note: all trees must have same number of entries.
   // Note: pointers **datatree will return normalized data trees which are 
   //       already background subtracted, thus &numbgrd returns numbgrd=0,
   //       to avoid background subtraction in method Express()
   if(kCS) cout << "------XGCProcesSet::Normalize------" << endl;

   Int_t err = errNoErr;

   TString treename;
   TString refstr;
   TObjString *selstr;

   TDirectory *savedir = gDirectory;

// Get chip parameters from scheme file (also alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XDNAChip *chip = (XDNAChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) {
      return fManager->HandleError(errGetScheme, fSchemeName);
   }//if
   Int_t numrows = chip->GetNumRows();
   Int_t numcols = chip->GetNumColumns();
   Int_t size    = numrows*numcols;

// Change directory
   TFile *tmpfile = fNormalizer->GetFile();
   if (tmpfile != 0) tmpfile->cd();
   else if (!fFile->cd(fName)) return errGetDir;

// Check for equal number of data tree entries
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;

      if ((Int_t)(datatree[k]->GetEntries()) != size) {
         TString str = ""; str += size;
         return fManager->HandleError(errNumTreeEntries, datatree[k]->GetName(), str);
      }//if
   }//for_k

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fNormalizer, fNormalizer->GetBgrdOption());

// Check for equal number of background tree entries
   if (numbgrd > 0) {
      for (Int_t k=0; k<numdata; k++) {
         if (bgrdtree[k] == 0) return errGetTree;

         if ((Int_t)(bgrdtree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, bgrdtree[k]->GetName(), str);
         }//if
      }//for_k
   } else if (doBg == kTRUE) {
      cout << "Warning: No background trees available for background subtraction."
           << endl;
      doBg = kFALSE;
   } else {
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Initialize reference trees
   Int_t  numsels = fSelections ? fSelections->GetSize() : 0;
   Int_t  numrefs = fReferences ? fReferences->GetSize() : 0;
   TTree *reftree[numrefs];
   for (Int_t j=0; j<numrefs; j++) reftree[j] = 0;

///////////
//?? Get reference trees (need not be contained in list of trees)
//?? Problem: for trees from other root files no bgrd subtraction allowed
///////////
// Get reference trees (must be contained in list of trees)
   Int_t refid = 0;
   if (numrefs == 1) {
      // get reference tree
      reftree[0] = (TTree*)fReferences->At(0); 
      if (reftree[0] == 0) return errGetTree;

      // get id of reference datatree
      refstr = reftree[0]->GetName();
//??      for (Int_t k=0; k<numsels; k++) {  //??since numsels > numdata ??
      for (Int_t k=0; k<numdata; k++) {
         selstr = (TObjString*)(fSelections->At(k));
         if (strcmp(selstr->GetName(),reftree[0]->GetName()) == 0) {
            refid = k;
            break;
         }//if
      }//for_k
   } else if (numrefs > 1) {
      for (Int_t j=0; j<numrefs; j++) {
         reftree[j] = (TTree*)fReferences->At(j); 
         if (reftree[j] == 0) return errGetTree;
      }//for_j
   }//if

// Get exon level of annotation
   Int_t level = 0;
   if (fNormSelector->GetNumParameters() > 0) {
      level = (Int_t)(fNormSelector->GetParameters())[0];
   }//if

// Init local arrays to store data from trees
   Int_t    *arrMask = 0;  //mask data
   Double_t *arrIntx = 0;  //intensities for tree x (e.g. raw intensities)
   Double_t *arrInty = 0;  //intensities for tree y (e.g. output intensities)

// Initialize memory for local arrays
   if (!(arrMask = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
   if (!(arrIntx = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
   if (!(arrInty = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<size; i++) {
      arrMask[i] = eINITMASK;
      arrIntx[i] = arrInty[i] = 0.0;
   }//for_i

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

// Get mask for PM/MM from scheme tree and store in array 
   err = this->SchemeMask(chip, level, size, arrMask);
   if (err != errNoErr) goto cleanup;

   // informing user
   if (XManager::fgVerbose) {
      cout << "   Normalizing data using method <" << fNormalizer->GetName() << ">..." << endl;
   }//if

// Compute reference arrays
   if (strcmp(fNormalizer->GetName(), "quantile") == 0) {

      // get mask for data to be used for normalization
      err = fNormSelector->Calculate(size, 0, 0, arrMask);
      if (err != errNoErr) goto cleanup;
/*
//?? save mask tree???
      TString selname = Path2Name(datatree[k]->GetName(),"",".");
      treename = selname + "." + kExtenMask[0];
      err = FillMaskTree(treename, fNormSelector, numrows, numcols, arrMask);
      if (err != errNoErr) break;
*/

      // sort arrays for quantile normalization
      for (Int_t k=0; k<numdata; k++) {
         treename = Path2Name(datatree[k]->GetName(),"",".");

         // informing user
         if (XManager::fgVerbose) {
            cout << "      filling array <" << treename << ">..." << endl;
         }//if

         err = this->FillDataArrays(datatree[k], bgrdtree[k], doBg,
                     numrows, numcols, arrIntx, 0, 0);
         if (err != errNoErr) goto cleanup;

         err = fNormalizer->AddArray(size, arrIntx, arrMask, treename);
         if (err != errNoErr) goto cleanup;
      }//for_k
   } else if (numrefs == 1) {
//////////
//TO DO: if reftree is not one of datatrees, ev. from different root file
//       datatree[refid] not possible, need to use reftree[0] (w/o bgrd??)
//////////
      // informing user
      if (XManager::fgVerbose) {
         cout << "      filling array <" << datatree[refid]->GetName() << ">..." << endl;
      }//if

      err = this->FillDataArrays(datatree[refid], bgrdtree[refid], doBg,
                  numrows, numcols, arrIntx, 0, 0);
      if (err != errNoErr) goto cleanup;

      TString exten = fNormalizer->GetTitle();
      datatree[refid] = this->FillDataTree(datatree[refid], exten,
                              fNormalizer, numrows, numcols, arrIntx, 0);
      if (datatree[refid] == 0) {err = errCreateTree; goto cleanup;}
   } else if (numrefs > 1) {
      // informing user
      if (XManager::fgVerbose) {
         cout << "      filling array <" << kReference << ">..." << endl;
      }//if

      TTree *rbgtree[numrefs];
      for (Int_t j=0; j<numrefs; j++) rbgtree[j] = 0;

      // find correct bgrdtree for each reftree
      Int_t numrbgs = 0;
      if (numbgrd > 0) {
         // create hash table to store datatree names
         TString     tmpstr;
         XIdxString *idxstr = 0;
         THashTable *htable = 0;
         if (!(htable = new THashTable(2*numdata))) {err = errInitMemory; goto cleanup;}

         // fill htable with index and name (w/o extension) of reftrees
         for (Int_t j=0; j<numrefs; j++) {
            tmpstr = Path2Name(reftree[j]->GetName(), "", ".");
            idxstr = new XIdxString(j, tmpstr.Data());
            htable->Add(idxstr);
         }//for_j

         // initialize rbgtrees in order of reftrees
         for (Int_t j=0; j<numrefs; j++) {
            tmpstr = Path2Name(bgrdtree[j]->GetName(), "", ".");
            idxstr = (XIdxString*)(htable->FindObject(tmpstr.Data()));
            if (idxstr) {
               rbgtree[idxstr->GetIndex()] = bgrdtree[j];
               numrbgs++;
            }//if
         }//for_j

         // check for equal number of bgrd trees and corresponding reference trees
         if ((numrbgs > 0) && (numrbgs != numrefs)) {
            cerr << "Error: <" << numrefs - numrbgs 
                 << "> reference trees have no corresponding background tree!"
                 << endl;
            err = errAbort;
            goto cleanup;
         }//if

         // delete idxstr before deleting htable
         if (htable) {htable->Delete(); delete htable; htable = 0;}
      }//if

      // get array of mean/median intensities
      if (strcmp(fRefOpt.Data(), "mean") == 0) {
         err = MeanReference(numrefs, reftree, numrbgs, rbgtree, numrows, numcols, arrIntx, doBg);
         if (err != errNoErr) goto cleanup;
      } else if (strcmp(fRefOpt.Data(), "median") == 0) {
         err = MedianReference(numrefs, reftree, numrbgs, rbgtree, numrows, numcols, arrIntx, doBg);
         if (err != errNoErr) goto cleanup;
      }//if

      treename = TString(kReference) + "." + TString(fNormalizer->GetTitle());
      datatree[numsels] = FillDataTree(treename, fNormalizer, numrows, numcols, arrIntx);
      if (datatree[numsels] == 0) {err = errCreateTree; goto cleanup;}

      // add reference tree to selections list
      Select(kReference, 1);
      refstr = ((TObjString*)fSelections->Last())->GetString();
//PROBLEM?? see XPSNormation.cxx ca line 1030!! need to create NEW refstrg??:
//??refstr = new TObjString(kReference);
      numsels++;
   } else {
         cerr << "Error: No reference tree is selected." << endl;
         err = errGeneral;
         goto cleanup;
   }//if

// Normalization
   if (strcmp(fNormalizer->GetName(), "quantile") == 0) {
      // quantile normalization
      if ((err = fNormalizer->Calculate(size, arrIntx, arrInty, arrMask))) goto cleanup;

      // get normalized arrays and fill data trees
      for (Int_t k=0; k<numdata; k++) {
         treename = Path2Name(datatree[k]->GetName(),"",".");
         TString exten = fNormalizer->GetTitle();

         // informing user
         if (XManager::fgVerbose) {
            cout << "      filling tree <" << (treename + "." + exten) << ">..." << endl;
         }//if

         arrInty = fNormalizer->GetArray(size, arrInty, arrMask, treename);
         if (arrInty == 0) {err = errAbort; goto cleanup;}

         datatree[k] = this->FillDataTree(datatree[k], exten, fNormalizer,
                             numrows, numcols, arrInty, 0);
         if (datatree[k] == 0) {err = errCreateTree; goto cleanup;}
      }//for_k
   } else if ((strcmp(fNormSelector->GetName(),   "probe")    == 0) ||
              ((strcmp(fNormSelector->GetName(),  "rank")     == 0) &&
               (strcmp(fNormSelector->GetOption(),"separate") == 0))) {
      for (Int_t k=0; k<numsels; k++) {
         // informing user
         if (XManager::fgVerbose) {
            cout << "      normalizing <" << datatree[k]->GetName() << ">..." << endl;
         }//if

         selstr = (TObjString*)(fSelections->At(k));

//not necessary?          // if not mean/median fill reference tree only
         // for mean/median the tree used as reference must be normalized, too!
         // e.g. to scale tree to targetinten
         if (strcmp(selstr->GetName(),refstr.Data()) == 0) {
            // skip if reference tree is kReference (numrefs > 1)
            if (strcmp(selstr->GetName(), kReference) == 0) {
               continue;
            }//if

            // skip reference tree if not mean/median 
            if (!((strcmp(fNormalizer->GetName(),"mean")   == 0) ||
                  (strcmp(fNormalizer->GetName(),"median") == 0))) {
               continue;
            }//if
         }//if

         err = FillDataArrays(datatree[k], numrows, numcols, arrInty, 0, 0);
         if (err != errNoErr) break;

//ev??// reset mask to initial value for each k
//for (Int_t i=0; i<size; i++) arrMask[i] = arrTmp[i];
         // select non-variant units, calculate mask and fill mask tree
         if ((err = fNormSelector->Calculate(size, arrIntx, arrInty, arrMask))) break;

         TString selname = Path2Name(datatree[k]->GetName(),"",".");
         treename = selname + "." + kExtenMask[0];
         err = FillMaskTree(treename, fNormSelector, numrows, numcols, arrMask);
         if (err != errNoErr) break;

         // normalize expression levels arrY and fill expression tree
         if ((err = fNormalizer->Calculate(size, arrIntx, arrInty, arrMask))) break;

         TString exten = fNormalizer->GetTitle();
         datatree[k] = this->FillDataTree(datatree[k], exten, fNormalizer,
                             numrows, numcols, arrInty, 0);
         if (datatree[k] == 0) {err = errCreateTree; goto cleanup;}
      }//for_k
   } else if ((strcmp(fNormSelector->GetName(),   "rank")     == 0) &&
              (strcmp(fNormSelector->GetOption(), "together") == 0)) {
      cerr << "Note: Rank selector option <together> is not yet implemented."
           << endl;
      err = errGeneral;
   } else if (strcmp(fNormSelector->GetName(), "user") == 0) {
      cerr << "Note: User selector is not yet implemented." << endl;
      err = errGeneral;
   } else {
      cerr << "Error: Selector option <" << fNormSelector->GetOption()
           << "> is not known/applicable." << endl;
      err = errGeneral;
   }//if

// IMPORTANT: set numbgrd=0 since background already subtracted!!!
   numbgrd = 0;

// Cleanup
cleanup:
   // delete arrays
   if (arrInty) {delete [] arrInty; arrInty = 0;}
   if (arrIntx) {delete [] arrIntx; arrIntx = 0;}
   if (arrMask) {delete [] arrMask; arrMask = 0;}
   // delete scheme tree from RAM
//   if (scmtree) {scmtree->Delete(""); scmtree = 0;}

   savedir->cd();

   return err;
}//Normalize

//______________________________________________________________________________
Int_t XGCProcesSet::Express(Int_t numdata, TTree **datatree,
                    Int_t &numbgrd, TTree **bgrdtree)
{
   // Convert intensity data to expression levels
   if(kCS) cout << "------XGCProcesSet::Express------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "   Converting raw data to expression levels..." << endl;
   }//if

   Int_t err = errNoErr;

   if (strcmp(fExpressor->GetName(), "medianpolish") == 0) {
      if (fExpressor->GetFile()) {
         err = this->DoMedianPolish(numdata, datatree, numbgrd, bgrdtree,
                                    fExpressor->GetFile());
      } else {
         err = this->DoMedianPolish(numdata, datatree, numbgrd, bgrdtree);
      }//if
   } else {
      err = this->DoExpress(numdata, datatree, numbgrd, bgrdtree);
   }//if

   return err;
}//Express

//______________________________________________________________________________
Int_t XGCProcesSet::DetectCall(Int_t numdata, TTree **datatree,
                    Int_t &numbgrd, TTree **bgrdtree)
{
   // Detect presence call
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   // Note: each datatree is processed individually and thus can have different
   //       number of entries (i.e. belong to different chip types)
   if(kCS) cout << "------XGCProcesSet::DetectCall------" << endl;

   Int_t err= errNoErr;

   fFile->cd();

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fCaller, fCaller->GetBgrdOption());

// Check for presence of background trees
   if (numbgrd == 0) {
      if (doBg == kTRUE) {
         cout << "Warning: No background trees available for background subtraction."
              << endl;
      }//if
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Init local arrays to store data from trees
   Int_t    *arrMask  = 0;
   Double_t *arrInten = 0;
   Double_t *arrStdev = 0;
   Int_t    *arrNPix  = 0;
   Double_t *arrPM    = 0; 
   Double_t *arrMM    = 0;
   Double_t *arrSP    = 0; 
   Double_t *arrSM    = 0;
   Int_t    *arrXP    = 0; 
   Int_t    *arrXM    = 0;

   TTree  *calltree = 0;
   XPCall *call     = 0;
   Int_t   split    = 99;

// Calculate detection call
   Bool_t doGC = (Bool_t)(strcmp(fCaller->GetName(), "dabgcall") == 0);
   Int_t  i, j, ij, idx, x, y, start, end;
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;

   // Informing user
      TString name = datatree[k]->GetName();
      if (XManager::fgVerbose) {
         cout << "   Calculating present call for <" << name << ">..." << endl;
      }//if

   // Get tree info for datatree name.exten
      XDataTreeInfo *info = 0;
      info = (XDataTreeInfo*)datatree[k]->GetUserInfo()->FindObject(name);
      if (!info) {
         cerr << "Error: Could not get tree info for <" << name << ">." << endl;
         return errGeneral;
      }//if

   // Get chip parameters from scheme file (also alternative CDFs)
      if (strcmp(fSchemeName.Data(), "") == 0) {
         fSchemeName = datatree[k]->GetTitle();
      } else if (!fSchemeName.Contains(datatree[k]->GetTitle())) {
         return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[k]->GetTitle());
      }//if
      if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

      XDNAChip *chip = (XDNAChip*)fSchemes->FindObject(fSchemeName, kTRUE);
      if (!chip) {
         return fManager->HandleError(errGetScheme, fSchemeName);
      }//if
      Int_t numrows  = chip->GetNumRows();
      Int_t numcols  = chip->GetNumColumns();
      Int_t numctrls = chip->GetNumControls();
      Int_t numunits = chip->GetNumUnits();

      Int_t numabsent  = 0;
      Int_t numarginal = 0;
      Int_t numpresent = 0;
      Double_t minpval = 1.0;
      Double_t maxpval = 0.0;

   // Get scheme tree for scheme
      XScheme *scheme = 0;
      TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
      if (scmtree == 0) return errGetTree;
      scmtree->SetBranchAddress("ScmBranch", &scheme);

   // Get unit tree for scheme
      XGCUnit *unit = 0;
      TTree *idxtree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
      if (idxtree == 0) return errGetTree;
      idxtree->SetBranchAddress("IdxBranch", &unit);

   // Get maximum number of pairs from tree info for unit tree
      XUnitTreeInfo *idxinfo = 0;
      idxinfo = (XUnitTreeInfo*)idxtree->GetUserInfo()->FindObject(idxtree->GetName());
      if (!idxinfo) {
         cerr << "Error: Could not get tree info for <" << idxtree->GetName() << ">." << endl;
         return errGeneral;
      }//if
      Int_t maxnumpairs = (Int_t)idxinfo->GetValue("fMaxNPairs");

   // Get exon level of annotation
      Int_t level = 0;
      if (fCallSelector->GetNumParameters() > 0) {
         level = (Int_t)(fCallSelector->GetParameters())[0];
      }//if

   // Initialize memory for data arrays
      Int_t size     = numrows*numcols;
      Int_t nentries = scmtree->GetEntries();
      if (!(arrMask  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
      if (!(arrInten = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrStdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrNPix  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}

      for (i=0; i<size; i++) {
         arrMask[i]  = eINITMASK;  //not =0,since msk=1 and msk=0 must be determined!
         arrInten[i] = arrStdev[i] = 0.0;
         arrNPix[i]  = 0;
      }//for_i

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

   // Get mask for PM/MM from scheme tree and store in array 
//?0?      arrMask = this->FillMaskArray(chip, scmtree, scheme, 0, size, arrMask);
      arrMask = this->FillMaskArray(chip, scmtree, scheme, level, size, arrMask);

   // Calculate mask for detection call
      err = fCallSelector->Calculate(nentries, 0, 0, arrMask);
      if (err != errNoErr) goto cleanup;

   // Get data from datatree and store in arrays
      err = this->FillDataArrays(datatree[k], bgrdtree[k], doBg,
                       numrows, numcols, arrInten, arrStdev, arrNPix);
      if (err != errNoErr) goto cleanup;

   // For DABG only, fill arrMask with GC content (GC>=0 for PM)
   // for PM set GC >= 0 , i.e. GC = 0...kProbeLength,
   // for MM set GC < eINITMASK, i.e. GC = eINITMASK - (1...kProbeLength+1)
      if (doGC) {
         // get probe tree for scheme
         XGCProbe *probe = 0;
         TTree *probetree = (TTree*)gDirectory->Get(chip->GetProbeTree()); 
         if (probetree == 0) {err = errGetTree; goto cleanup;}
         probetree->SetBranchAddress("PrbBranch", &probe);

         // get GC content from probe tree and store in array
         for (i=0; i<size; i++) {
            probetree->GetEntry(i);

            x  = probe->GetX();
            y  = probe->GetY();
            ij = XY2Index(x, y, numcols);

            if (arrMask[ij] == 1) {
               arrMask[ij] = probe->GetNumberGC();
            } else if (arrMask[ij] == 0) {
               //need to use (numberGC + 1) to avoid setting arrMask=eINITMASK for GC=0!!
               arrMask[ij] = eINITMASK - (probe->GetNumberGC() + 1);
//no!!               arrMask[ij] = eINITMASK - probe->GetNumberGC();
            }//if
         }//for_i

         // get intensities and save in table sorted for GC-content
         if ((err = fCaller->Calculate(size, arrInten, arrMask))) goto cleanup;
      }//if

   // Create new tree calltree
      if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

      name = Path2Name(datatree[k]->GetName(),"/",".") + "." + fCaller->GetTitle();
      calltree = new TTree(name, fSchemeName);
      if (calltree == 0) {err = errCreateTree; goto cleanup;}
      call = new XPCall();
      calltree->Branch("CallBranch", "XPCall", &call, 64000, split);

   // Initialize maximum memory for PM/MM arrays (maxnumpairs+1 to avoid potential buffer overflow)
      if (!(arrPM = new (nothrow) Double_t[maxnumpairs+1])) {err = errInitMemory; goto cleanup;}
      if (!(arrMM = new (nothrow) Double_t[maxnumpairs+1])) {err = errInitMemory; goto cleanup;}
      if (!(arrSP = new (nothrow) Double_t[maxnumpairs+1])) {err = errInitMemory; goto cleanup;}
      if (!(arrSM = new (nothrow) Double_t[maxnumpairs+1])) {err = errInitMemory; goto cleanup;}
      if (!(arrXP = new (nothrow) Int_t[maxnumpairs+1]))    {err = errInitMemory; goto cleanup;}
      if (!(arrXM = new (nothrow) Int_t[maxnumpairs+1]))    {err = errInitMemory; goto cleanup;}

   // Calculate detection call values
      start = 0;
      end   = 0;
      idx   = 0;
      for (Int_t id=0; id<numunits; id++) { 
         idxtree->GetEntry(id);

         Int_t unitID   = id - numctrls;
         Int_t numcells = unit->GetNumCells();
         // skip negative unit entries (controls)
         if (unit->GetUnitID() < 0) {
            start += numcells;
            end = start;
            continue;
         }//if

         Int_t p = 0;
         Int_t m = 0;
         end += numcells;
         for (j=start; j<end; j++) {
            scmtree->GetEntry(j);

            if ((scheme->GetUnitID()) != unitID) {
               cerr << "Error: unitID is not equal to: " << unitID << endl;
               err = errAbort;
               goto cleanup;
            }//if

            x  = scheme->GetX();
            y  = scheme->GetY();
            ij = XY2Index(x, y, numcols);

            if (doGC) {
               if (arrMask[ij] >= 0) {
                  if (p == 0) idx++;  //count number of units
                  arrPM[p] = arrInten[ij];
                  arrSP[p] = arrStdev[ij];
                  arrXP[p] = arrMask[ij];
                  p++;
               }//if
            } else {
               if (arrMask[ij] == 1) {
                  if (p == 0) idx++;  //count number of units
                  arrPM[p] = arrInten[ij];
                  arrSP[p] = arrStdev[ij];
                  arrXP[p] = arrNPix[ij];
                  p++;
               } else if (arrMask[ij] == 0) {
                  arrMM[m] = arrInten[ij];
                  arrSM[m] = arrStdev[ij];
                  arrXM[p] = arrNPix[ij];
                  m++;
               }//if
            }//if

            if (p > maxnumpairs || m > maxnumpairs) {
               cerr << "Error: unitID <" << unitID << "> exceeds maximum number of pairs <"
                    << maxnumpairs << ">. Buffer overflow!" << endl;
               err = errAbort;
               goto cleanup;
            }//if
         }//for_j
         start += numcells;

         // continue if arrays arrPM etc are not filled
         if (p == 0) continue;
         if (!doGC && (p != m)) {
            cerr << "Error: UnitID <" << unitID << "> has different numbers of PM <"
                 << p << "> and MM <" << m << "> data." << endl;
//          continue;
            err = errAbort;
            goto cleanup;
         }//if

         // calculate detection call
         Int_t    dummy    = 0;
         Double_t prescall = 0.0;
         Double_t pvalue   = 1.0;
         fCaller->InitTreeInfo(info);
         fCaller->InitArrays(p, arrPM, arrSP, arrXP, arrMM, arrSM, arrXM);
         if ((err = fCaller->Calculate(prescall, pvalue, dummy))) goto cleanup;

         // number of present/absent calls
         if      (prescall == 2.0) numpresent++;
         else if (prescall == 1.0) numarginal++;
         else if (prescall == 0.0) numabsent++;

         // minimal/maximal detection call p-value
         if      (pvalue < minpval) minpval = pvalue;
         else if (pvalue > maxpval) maxpval = pvalue;

         // fill call tree
         call->SetUnitID(unitID);
         call->SetCall((Short_t)prescall);
         call->SetPValue(pvalue);
         calltree->Fill();

         if (XManager::fgVerbose && id%1000 == 0) {
            cout << "      <" << id+1 << "> of <" << numunits << "> calls processed...\r" << flush;
         }//if
      }//for_id
      if (XManager::fgVerbose) {
         cout << "      <" << numunits << "> calls processed." << endl;
      }//if

      if (XManager::fgVerbose) {
         cout << "      detection call statistics: " << endl;
         cout << "         minimum detection p-value = " << minpval << endl;
         cout << "         maximum detection p-value = " << maxpval << endl;
         cout << "         P: <" << numpresent*100.0/idx << "> percent units present."  << endl;
         cout << "         M: <" << numarginal*100.0/idx << "> percent units marginal." << endl;
         cout << "         A: <" << numabsent*100.0/idx  << "> percent units absent."   << endl;
      }//if

   // Add tree info to tree
      AddCallTreeInfo(calltree, calltree->GetName(), fCaller->GetOption(),
                      numunits-numctrls, numabsent, numarginal, numpresent,
                      minpval, maxpval);

   // Write call tree to file 
      if ((err = WriteTree(calltree, TObject::kOverwrite)) == errNoErr) {
//?      if ((err = WriteTree(calltree, 0)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(calltree->GetName(), "Call", 0, fCaller->GetNumParameters(),
                       fCaller->GetParameters());
      }//if

   cleanup:
      // delete arrays
      if (arrPM)    {delete [] arrPM;    arrPM    = 0;}
      if (arrMM)    {delete [] arrMM;    arrMM    = 0;}
      if (arrSP)    {delete [] arrSP;    arrSP    = 0;}
      if (arrSM)    {delete [] arrSM;    arrSM    = 0;}
      if (arrXP)    {delete [] arrXP;    arrXP    = 0;}
      if (arrXM)    {delete [] arrXM;    arrXM    = 0;}
      if (arrNPix)  {delete [] arrNPix;  arrNPix  = 0;}
      if (arrStdev) {delete [] arrStdev; arrStdev = 0;}
      if (arrInten) {delete [] arrInten; arrInten = 0;}
      if (arrMask)  {delete [] arrMask;  arrMask  = 0;}
//??      // delete scheme tree from RAM
//??      if (scmtree)  {scmtree->Delete(""); scmtree = 0;}

      if (err != errNoErr) break;
   }//for_k

   return err;
}//DetectCall

//______________________________________________________________________________
Int_t XGCProcesSet::ExportTreeType(const char *exten, Int_t n, TString *names, 
                    const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output
   if(kCS) cout << "------XGCProcesSet::ExportTreeType------" << endl;

// Set scheme file to be able to access scheme data for exporting
   if (fSetting) {
      fSchemeFile = ((XPreProcesSetting*)fSetting)->GetSchemeFile();
   }//if

   if (HasExtension(exten, kExtenBgrd)) {
      return this->ExportBgrdTrees(n, names, varlist, output, sep);
   } else if (HasExtension(exten, kExtenIntn)) {
      return this->ExportIntnTrees(n, names, varlist, output, sep);
   } else if (HasExtension(exten, kExtenCNrm)) {
      return this->ExportNormTrees(n, names, varlist, output, sep);
   } else if (HasExtension(exten, kExtenExpr)) {
      return this->ExportExprTrees(n, names, varlist, output, sep);
   } else if (HasExtension(exten, kExtenCall)) {
      return this->ExportCallTrees(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if

   return errNoErr;
}//ExportTreeType

//______________________________________________________________________________
Int_t XGCProcesSet::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                    const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output as XML-file
   if(kCS) cout << "------XGCProcesSet::ExportTreeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   exten = 0; n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTreeXML

//______________________________________________________________________________
Int_t XGCProcesSet::ExportBgrdTrees(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export data stored in background tree to file output
   if(kCS) cout << "------XGCProcesSet::ExportBgrdTrees------" << endl;

// Decompose varlist
   Bool_t hasBgrd = kFALSE;
   Bool_t hasStdv = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasBgrd = kTRUE;
      hasStdv = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fBg")    == 0) {hasBgrd = kTRUE;}
         if (strcmp(name,"fStdev") == 0) {hasStdv = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get trees
   TTree   *tree[n];
   XBgCell *cell[n];
   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         cell[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("BgrdBranch", &cell[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         cell[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("BgrdBranch", &cell[k]);
      }//for_k
   }//if

// Output header
   output << "X" << sep << "Y";
   if (n > 1) {
      for (Int_t i=0; i<n; i++) {
         if (hasBgrd) output << sep << (names[i] + "_BGRD");
         if (hasStdv) output << sep << (names[i] + "_STDV");
      }//for_i
   } else {
      if (hasBgrd) output << sep << "BGRD";
      if (hasStdv) output << sep << "STDV";
   }//if
   output << endl;

// Loop over tree entries and trees
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);
         if (k == 0)  output << cell[k]->GetX() << sep << cell[k]->GetY();
         if (hasBgrd) output << sep << cell[k]->GetBackground();
         if (hasStdv) output << sep << cell[k]->GetStdev();
      }//for_k
      output << endl;
   }//for_i

   return errNoErr;
}//ExportBgrdTrees

//______________________________________________________________________________
Int_t XGCProcesSet::ExportIntnTrees(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export data stored in background-corrected intensity tree to file output
   if(kCS) cout << "------XGCProcesSet::ExportIntnTrees------" << endl;

// Decompose varlist
   Bool_t hasMean = kFALSE;
   Bool_t hasStdv = kFALSE;
   Bool_t hasNPix = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasMean = kTRUE;
      hasStdv = kTRUE;
      hasNPix = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fInten")   == 0) {hasMean = kTRUE;}
         if (strcmp(name,"fStdev")   == 0) {hasStdv = kTRUE;}
         if (strcmp(name,"fNPixels") == 0) {hasNPix = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get trees
   TTree   *tree[n];
   XGCCell *cell[n];
   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         cell[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("DataBranch", &cell[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         cell[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("DataBranch", &cell[k]);
      }//for_k
   }//if

// Output header
   output << "X" << sep << "Y";
   if (n > 1) {
      for (Int_t i=0; i<n; i++) {
         if (hasMean) output << sep << (names[i] + "_MEAN");
         if (hasStdv) output << sep << (names[i] + "_STDV");
         if (hasNPix) output << sep << (names[i] + "_NPIXELS");
      }//for_i
   } else {
      if (hasMean) output << sep << "MEAN";
      if (hasStdv) output << sep << "STDV";
      if (hasNPix) output << sep << "NPIXELS";
   }//if
   output << endl;

// Loop over tree entries and tree branches
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);
         if (k == 0)  output << cell[k]->GetX() << sep << cell[k]->GetY();
         if (hasMean) output << sep << cell[k]->GetIntensity();
         if (hasStdv) output << sep << cell[k]->GetStdev();
         if (hasNPix) output << sep << cell[k]->GetNumPixels();
      }//for_k
      output << endl;
   }//for_i

   return errNoErr;
}//ExportIntnTrees

//______________________________________________________________________________
Int_t XGCProcesSet::ExportNormTrees(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export data stored in normalized tree to file output
   if(kCS) cout << "------XGCProcesSet::ExportNormTrees------" << endl;

   return this->ExportIntnTrees(n, names, varlist, output, sep);
}//ExportNormTrees

//______________________________________________________________________________
Int_t XGCProcesSet::ExportExprTrees(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export variables from varlist for expression tree(s) to file output
   if(kCS) cout << "------XGCProcesSet::ExportExprTrees------" << endl;

// Decompose varlist
   Bool_t hasUnit   = kFALSE;
   Bool_t hasName   = kFALSE;
   Bool_t hasSymbol = kFALSE;
   Bool_t hasCyto   = kFALSE;
   Bool_t hasAnnot  = kFALSE;
   Bool_t hasLevel  = kFALSE;
   Bool_t hasStdev  = kFALSE;
   Bool_t hasNPairs = kFALSE;
////////////////
//TO DO
//fTranscriptID
//fAccession
//fEntrezID
//fChromosome
//fStart
//fStop
//fStrand
////////////////

   if (strcmp(varlist,"*")  == 0) {
      hasUnit   = kTRUE;
      hasName   = kTRUE;
      hasSymbol = kTRUE;
      hasCyto   = kTRUE;
      hasLevel  = kTRUE;
      hasStdev  = kTRUE;
      hasNPairs = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName") == 0) {hasUnit   = kTRUE;}
         if (strcmp(name,"fName")     == 0) {hasName   = kTRUE;}
         if (strcmp(name,"fSymbol")   == 0) {hasSymbol = kTRUE;}
         if (strcmp(name,"fCytoBand") == 0) {hasCyto   = kTRUE;}
         if (strcmp(name,"fLevel")    == 0) {hasLevel  = kTRUE;}
         if (strcmp(name,"fStdev")    == 0) {hasStdev  = kTRUE;}
         if (strcmp(name,"fNPairs")   == 0) {hasNPairs = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if
   hasAnnot = (hasName || hasSymbol || hasCyto);

// Get trees
   TTree       *tree[n];
   XExpression *expr[n];

   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         expr[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("ExprBranch", &expr[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         expr[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("ExprBranch", &expr[k]);
      }//for_k
   }//if

// Get scheme name (also for alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = tree[0]->GetTitle();
   } else if (!fSchemeName.Contains(tree[0]->GetTitle())) {
      cerr << "Error: Scheme <" << fSchemeName << "> is not derived from <"
           << tree[0]->GetTitle() << ">." << endl;
      hasUnit  = kFALSE;
      hasAnnot = kFALSE;
   }//if

// If scheme file is present get unit tree
   Int_t    numgenes = 0;
   Int_t    numctrls = 0;
   XFolder *schemes  = 0;
   XGCUnit *unit     = 0;
   TTree   *unittree = 0; 
   TTree   *anntree  = 0; 
   XTransAnnotation *annot = 0;
   if ((hasUnit || hasAnnot) && fSchemeFile) {
   // Get chip from scheme file
      fSchemeFile->cd();
      schemes = (XFolder*)(fSchemeFile->Get(kContent));
      if (!schemes) {
         return fManager->HandleError(errMissingContent, "Scheme", kContent);
      }//if

      XDNAChip *chip = 0;
      chip = (XDNAChip*)schemes->FindObject(fSchemeName, kTRUE);
      if (chip) {
         numgenes = chip->GetNumGenes();
         numctrls = chip->GetNumControls();

         if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

      // Get unit tree for scheme
         if (hasUnit) {
            unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
            if (unittree == 0) return errGetTree;
            unittree->SetBranchAddress("IdxBranch", &unit);
         }//if

      // Get annotation tree for scheme
         if (hasAnnot) {
            anntree = (TTree*)gDirectory->Get(chip->GetAnnotTree()); 
            if (anntree) {
               anntree->SetBranchAddress("AnnBranch", &annot);
            } else {
               cout << "Warning: Missing annotation, gene info not exported."
                    << endl;
               hasAnnot = kFALSE;
            }//if
         }//if
      } else {
         cerr << "Error: Could not find scheme <" << fSchemeName << ">." << endl;
         hasUnit  = kFALSE;
         hasAnnot = kFALSE;
      }//if
   } else if (hasUnit || hasAnnot) {
      cout << "Warning: Missing scheme file, unit names not exported." << endl;
      hasUnit  = kFALSE;
      hasAnnot = kFALSE;
   }//if

////////////
//better: see XExonProcesSet!!!
////////////
// Check if tree entries is equal to number of genes in unittree
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   if ((hasUnit || hasAnnot) && (entries != numgenes)) {
      if (XManager::fgVerbose) {
         cout << "Note: Number of tree entries <" << entries  
              << "> is not equal to number of units <" << numgenes << ">." << endl;
//?         cout << "         Unit data will not be exported." << endl;
      }//if
//?      hasUnit  = kFALSE;
//?      hasAnnot = kFALSE;
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit)      output << sep << "UNIT_NAME";
   if (hasAnnot) {
      if (hasName)   output << sep << "GENE_NAME";
      if (hasSymbol) output << sep << "GENE_SYMBOL";
      if (hasCyto)   output << sep << "CYTOBAND";
   }//if
   if (n == 1) {
      if (hasLevel)  output << sep << "LEVEL";
      if (hasStdev)  output << sep << "STDEV";
      if (hasNPairs) output << sep << "NUMBER_PAIRS";
   } else {
      for (Int_t k=0; k<n; k++) {
         if (hasLevel)  output << sep << (names[k] + "_LEVEL");
         if (hasStdev)  output << sep << (names[k] + "_STDEV");
         if (hasNPairs) output << sep << (names[k] + "_NUMBER_PAIRS");
      }//for_k
   }//if
   output << endl;

// Loop over tree entries and trees
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);

         // export annotation
         if (k == 0) {
            Int_t unitID = expr[k]->GetUnitID();
            output << unitID;

            if (hasUnit) {
               unittree->GetEntry(unitID + numctrls); //unit names for genes only
               output << sep << unit->GetUnitName();
            }//if

            if (hasAnnot) {
               anntree->GetEntry(unitID + numctrls);
               if (hasName)   output << sep << annot->GetName();
               if (hasSymbol) output << sep << annot->GetSymbol();
               if (hasCyto)   output << sep << annot->GetCytoBand();
            }//if
         }//if

         // export data
         if (hasLevel)  output << sep << expr[k]->GetLevel();
         if (hasStdev)  output << sep << ((XGCExpression*)expr[k])->GetStdev();
         if (hasNPairs) output << sep << ((XGCExpression*)expr[k])->GetNumPairs();
      }//for_j
      output << endl;
   }//for_i

//Cleanup
   // remove trees from RAM
   if (anntree)  {anntree->Delete("");  anntree  = 0;}
   if (unittree) {unittree->Delete(""); unittree = 0;}
   SafeDelete(schemes);

   return errNoErr;
}//ExportExprTrees

//______________________________________________________________________________
Int_t XGCProcesSet::ExportCallTrees(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export data stored in call tree to file output
   if(kCS) cout << "------XGCProcesSet::ExportCallTrees------" << endl;

// Decompose varlist
   Bool_t hasUnit   = kFALSE;
   Bool_t hasName   = kFALSE;
   Bool_t hasSymbol = kFALSE;
   Bool_t hasCyto   = kFALSE;
   Bool_t hasAnnot  = kFALSE;
   Bool_t hasCall   = kFALSE;
   Bool_t hasPVal   = kFALSE;
////////////////
//TO DO
//fTranscriptID
//fAccession
//fEntrezID
//fChromosome
//fStart
//fStop
//fStrand
////////////////

   if (strcmp(varlist,"*")  == 0) {
      hasUnit   = kTRUE;
      hasName   = kTRUE;
      hasSymbol = kTRUE;
      hasCyto   = kTRUE;
      hasCall   = kTRUE;
      hasPVal   = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while (name) {
         if (strcmp(name,"fUnitName") == 0) {hasUnit   = kTRUE;}
         if (strcmp(name,"fName")     == 0) {hasName   = kTRUE;}
         if (strcmp(name,"fSymbol")   == 0) {hasSymbol = kTRUE;}
         if (strcmp(name,"fCytoBand") == 0) {hasCyto   = kTRUE;}
         if (strcmp(name,"fCall")     == 0) {hasCall   = kTRUE;}
         if (strcmp(name,"fPValue")   == 0) {hasPVal   = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if
   hasAnnot = (hasName || hasSymbol || hasCyto);

// Get trees
   TTree  *tree[n];
   XPCall *call[n];
   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         call[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("CallBranch", &call[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         call[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("CallBranch", &call[k]);
      }//for_k
   }//if

// Get scheme name
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = tree[0]->GetTitle();
   } else if (!fSchemeName.Contains(tree[0]->GetTitle())) {
      cerr << "Error: Scheme <" << fSchemeName << "> is not derived from <"
           << tree[0]->GetTitle() << ">." << endl;
      hasUnit  = kFALSE;
      hasAnnot = kFALSE;
   }//if

// If scheme file is present get unit tree
   Int_t        numgenes = 0;
   Int_t        numctrls = 0;
   XFolder     *schemes  = 0;
   XGCUnit     *unit     = 0;
   TTree       *unittree = 0; 
   TTree       *anntree  = 0; 
   XTransAnnotation *annot = 0;
   if ((hasUnit || hasAnnot) && fSchemeFile) {
   // Get chip from scheme file
      fSchemeFile->cd();
      schemes = (XFolder*)(fSchemeFile->Get(kContent));
      if (!schemes) {
         return fManager->HandleError(errMissingContent, "Scheme", kContent);
      }//if

      XDNAChip *chip = 0;
      chip = (XDNAChip*)schemes->FindObject(fSchemeName, kTRUE);
      if (chip) {
         numgenes = chip->GetNumGenes();
         numctrls = chip->GetNumControls();

         if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

      // Get unit tree for scheme
         if (hasUnit) {
            unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
            if (unittree == 0) return errGetTree;
            unittree->SetBranchAddress("IdxBranch", &unit);
         }//if

      // Get annotation tree for scheme
         if (hasAnnot) {
            anntree = (TTree*)gDirectory->Get(chip->GetAnnotTree()); 
            if (anntree) {
               anntree->SetBranchAddress("AnnBranch", &annot);
            } else {
               cout << "Warning: Missing annotation, gene info not exported."
                    << endl;
               hasAnnot = kFALSE;
            }//if
         }//if
      } else {
         cerr << "Error: Could not find scheme <" << fSchemeName << ">." << endl;
         hasUnit  = kFALSE;
         hasAnnot = kFALSE;
      }//if
   } else if (hasUnit || hasAnnot) {
      cout << "Warning: Missing scheme file, unit names not exported." << endl;
      hasUnit  = kFALSE;
      hasAnnot = kFALSE;
   }//if

////////////
//better: see XExonProcesSet!!!
////////////
// Check if tree entries is equal to number of genes in unittree
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   if ((hasUnit || hasAnnot) && (entries != numgenes)) {
      if (XManager::fgVerbose) {
         cout << "Note: Number of tree entries <" << entries  
              << "> is not equal to number of units <" << numgenes << ">" << endl;
//?         cout << "         Unit data will not be exported." << endl;
      }//if
//?      hasUnit  = kFALSE;
//?      hasAnnot = kFALSE;
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit)      output << sep << "UNIT_NAME";
   if (hasAnnot) {
      if (hasName)   output << sep << "GENE_NAME";
      if (hasSymbol) output << sep << "GENE_SYMBOL";
      if (hasCyto)   output << sep << "CYTOBAND";
   }//if
   if (n == 1) {
      if (hasCall)   output << sep << "CALL";
      if (hasPVal)   output << sep << "PVALUE";
   } else {
      for (Int_t i=0; i<n; i++) {
         if (hasCall)  output << sep << (names[i] + "_CALL");
         if (hasPVal)  output << sep << (names[i] + "_PVALUE");
      }//for_i
   }//if
   output << endl;

// Loop over tree entries and tree branches
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);

         // export annotation
         if (k == 0) {
            Int_t unitID = call[k]->GetUnitID();
            output << unitID;

            if (hasUnit) {
               unittree->GetEntry(unitID + numctrls); //unit names for genes only
               output << sep << unit->GetUnitName();
            }//if

            if (hasAnnot) {
               anntree->GetEntry(unitID + numctrls);
               if (hasName)   output << sep << annot->GetName();
               if (hasSymbol) output << sep << annot->GetSymbol();
               if (hasCyto)   output << sep << annot->GetCytoBand();
            }//if
         }//if

         // export data
         if (hasCall) {
            Int_t cl = call[k]->GetCall();
            char *ch = "NA";
            if      (cl == 2) ch = "P";
            else if (cl == 0) ch = "A";
            else if (cl == 1) ch = "M";
            output << sep << ch;
         }//if

         if (hasPVal) {
            output << sep << call[k]->GetPValue();
         }//if
      }//for_j

      output << endl;
   }//for_i

//Cleanup
   // remove trees from RAM
   if (anntree)  {anntree->Delete("");  anntree  = 0;}
   if (unittree) {unittree->Delete(""); unittree = 0;}
   SafeDelete(schemes);

   return errNoErr;
}//ExportCallTrees

//______________________________________________________________________________
Int_t XGCProcesSet::SchemeMask(XDNAChip *chip, Int_t level, Int_t n, Int_t *msk)
{
   // Get name of scheme tree from "chip"
   // "level=0" is cutoff to prevent PM to be redefined as MM for alternative CDFs
   // Fill array "msk" with mask from scheme tree
   if(kCS) cout << "------XGCProcesSet::SchemeMask------" << endl;

   Int_t err = errNoErr;

// Change dir to scheme file
   TDirectory *savedir = gDirectory;
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

// Get scheme tree for scheme
   XScheme *scheme = 0;
   TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
   if (scmtree == 0) return errGetTree;
   scmtree->SetBranchAddress("ScmBranch", &scheme);

// Get mask for PM/MM from scheme tree and store in array 
   msk = this->FillMaskArray(chip, scmtree, scheme, level, n, msk);

   // delete scheme tree from RAM
   if (scmtree) {scmtree->Delete(""); scmtree = 0;}

   savedir->cd();

   return err;
}//SchemeMask

//______________________________________________________________________________
Int_t *XGCProcesSet::FillMaskArray(XDNAChip *chip, TTree *scmtree, XScheme *scheme,
                     Int_t level, Int_t /*n*/, Int_t *msk)
{
   // "level=0" is cutoff to prevent PM to be redefined as MM for alternative CDFs
   // Fill array "msk" with mask from scheme tree
   if(kCS) cout << "------XGCProcesSet::FillMaskArray------" << endl;

// Get mask for PM/MM from scheme tree and store in array 
   Int_t x, y, ij;
   Int_t numcols = chip->GetNumColumns();
   for (Int_t i=0; i<scmtree->GetEntries(); i++) {
      scmtree->GetEntry(i);

      x  = scheme->GetX();
      y  = scheme->GetY();
      ij = XY2Index(x, y, numcols);

      // once an oligo is defined as PM it cannot be redefined as MM:
      msk[ij] = (msk[ij] > level) ? msk[ij] : scheme->GetMask();
   }//for_i

   return msk;
}//FillMaskArray

//______________________________________________________________________________
Bool_t XGCProcesSet::BackgroundParameters(XAlgorithm *algorithm, const char *option)
{
   // Get parameters fBgPars for background subtraction for the following options
   // - "none":        no background subtraction
   // - "subtractbg":  subtract bgrd from intensity - result can be negative
   // - "correctbg":   correct bgrd with noise fraction to avoid negative results
   //    parameters are:
   //    - nfrac:      fBgPars[1] - noise fraction nfrac > 0 then replace negative
   //                  intensities with fraction nfrac of background intensity
   // - "attenuatebg": use generalized log-transform to avoid negative results
   //    parameters are:
   //    - l:          fBgPars[1] - tunable parameter, 0<=l<=1 (default is 0.005)    
   //    - h:          fBgPars[2] - parameter (default is -1)
   //
   // Note: First parameter  fBgPars[0] gives the background option, i.e.:
   //       - "subtractbg":  fBgPars[0] = 1;
   //       - "correctbg":   fBgPars[0] = 2;
   //       - "attenuatebg": fBgPars[0] = 3;
   //
   // Note: Set/add last parameter of/to InitAlgorithm() to "-100" to prevent background
   //       subtraction even if bgrdtrees exist. This allows e.g. calculation of bgrdtrees
   //       w/o subtracting background from datatrees
   if(kCS) cout << "------XGCProcesSet::BackgroundParameters------" << endl;

// Get algorithm parameters
   Int_t    npars = algorithm->GetNumParameters();
   Double_t *pars = algorithm->GetParameters();

// Reset background parameters
   if (fBgPars) {delete [] fBgPars; fBgPars = 0;}
   fNBgPar = 0;

// Set background parameters
   if ((strcmp(option, "none") == 0) || (strcmp(option, "") == 0)) {
      fNBgPar = 1;
      fBgPars = new (nothrow) Double_t[fNBgPar];

      fBgPars[0] = 0; //"none"
      return kFALSE;
   } else if (strcmp(option, "subtractbg") == 0) {
      fNBgPar = 1;
      fBgPars = new (nothrow) Double_t[fNBgPar];

      fBgPars[0] = 1; //"subtractbg"
   } else if (strcmp(option, "correctbg") == 0) {
      fNBgPar = 2;
      fBgPars = new (nothrow) Double_t[fNBgPar];

      fBgPars[0] = 2; //"correctbg"
      fBgPars[1] = (npars > 0) ? pars[npars-1] : 0.5;  //nfrac
   } else if (strcmp(option, "attenuatebg") == 0) {
      fNBgPar = 3;
      fBgPars = new (nothrow) Double_t[fNBgPar];

      fBgPars[0] = 3; //"attenuatebg"
      fBgPars[1] = (npars > 1) ? pars[npars-2] : 0.005;  //l
      fBgPars[2] = (npars > 1) ? pars[npars-1] : -1.0;   //h
   }//if

// Do bgrd subtraction?
   Bool_t doBg = kTRUE;
   if (npars > 0 && pars[npars-1] == -100) doBg = kFALSE;

   return doBg;
}//BackgroundParameters

//______________________________________________________________________________
Double_t XGCProcesSet::AdjustIntensity(Double_t inten, Double_t bgrd, Double_t stdv)
{
   // Adjust intensity based on fOption
   if(kCSa) cout << "------XGCProcesSet::AdjustIntensity------" << endl;

   if ((Int_t)fBgPars[0] == 0) { //"none"
      return inten;
   } else if ((Int_t)fBgPars[0] == 1) { //"subtractbg"
      inten = inten - bgrd;
   } else if ((Int_t)fBgPars[0] == 2) { //"correctbg"
      inten = inten - bgrd;
      inten = TMath::Max(inten, fBgPars[1]*stdv);
   } else if ((Int_t)fBgPars[0] == 3) { //"attenuatebg"
      Double_t hh = (fBgPars[2] < 0) ? 4*inten*bgrd*fBgPars[1] : fBgPars[2];
      Double_t xx = inten - bgrd;
      inten = (xx + TMath::Sqrt(TMath::Power(xx, 2) + hh))/2.0;
   }//if

   return inten;
}//AdjustIntensity

//______________________________________________________________________________
Int_t XGCProcesSet::FillDataArrays(TTree *datatree, TTree *bgrdtree, Bool_t doBg,
                    Int_t nrow, Int_t ncol, Double_t *inten, Double_t *stdev,
                    Int_t *npix)
{
   // Fill arrays inten, stdev, npix with data from datatree.
   // If bgrdtree != 0 and doBg == kTRUE then subtract background from intensity
   if(kCS) cout << "------XGCProcesSet::FillDataArrays------" << endl;

// Init datatree
   XGCCell *gccell = 0;
   datatree->SetBranchAddress("DataBranch", &gccell);

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

// Get data from datatree and store in arrays
//   Int_t size = nrow*ncol;
   Int_t size = datatree->GetEntries();
   Int_t x, y, ij;
   for (Int_t i=0; i<size; i++) {
      datatree->GetEntry(i);

      x  = gccell->GetX();
      y  = gccell->GetY();
      ij = XY2Index(x, y, ncol);

      if (inten) inten[ij] = gccell->GetIntensity();
      if (stdev) stdev[ij] = gccell->GetStdev();
      if (npix)  npix[ij]  = gccell->GetNumPixels();
   }//for_i

// Get background from bgrdtree and subtract from intensity
   if (bgrdtree && inten && doBg) {
      // init background tree
      XBgCell *bgcell = 0;
      bgrdtree->SetBranchAddress("BgrdBranch", &bgcell);

      // compare entries with number of entries in data tree
      Int_t numentries = (Int_t)(bgrdtree->GetEntries());
      if (numentries != size) {
         TString str = ""; str += size;
         return fManager->HandleError(errNumTreeEntries, bgrdtree->GetName(), str);
      }//if

      for (Int_t i=0; i<size; i++) {
         bgrdtree->GetEntry(i);

         x  = bgcell->GetX();
         y  = bgcell->GetY();
         ij = XY2Index(x, y, ncol);

         // subtract background from intensity
         inten[ij] = this->AdjustIntensity(inten[ij], bgcell->GetBackground(), bgcell->GetStdev());
      }//for_i
   }//if

   return errNoErr;
}//FillDataArrays

//______________________________________________________________________________
Int_t XGCProcesSet::FillDataArrays(TTree *datatree, Int_t nrow, Int_t ncol,
                    Double_t *inten, Double_t *stdev, Int_t *npix)
{
   // Fill arrays inten, stdev, npix with data from datatree.
   if(kCS) cout << "------XGCProcesSet::FillDataArrays------" << endl;

// Init datatree
   XGCCell *gccell = 0;
   datatree->SetBranchAddress("DataBranch", &gccell);

// Get data from datatree and store in arrays
//   Int_t size = nrow*ncol;
   Int_t size = datatree->GetEntries();
   Int_t x, y, ij;
   for (Int_t i=0; i<size; i++) {
      datatree->GetEntry(i);

      x  = gccell->GetX();
      y  = gccell->GetY();
      ij = XY2Index(x, y, ncol);

      if (inten) inten[ij] = gccell->GetIntensity();
      if (stdev) stdev[ij] = gccell->GetStdev();
      if (npix)  npix[ij]  = gccell->GetNumPixels();
   }//for_i

   return errNoErr;
}//FillDataArrays

//______________________________________________________________________________
Int_t XGCProcesSet::FillBgrdArrays(TTree *bgrdtree, Int_t nrow, Int_t ncol,
                    Double_t *inten, Double_t *stdev)
{
   // Fill arrays inten and stdev with data from bgrdtree.
   if(kCS) cout << "------XGCProcesSet::FillBgrdArrays------" << endl;

// Init bgrdtree
   XBgCell *bgcell = 0;
   bgrdtree->SetBranchAddress("BgrdBranch", &bgcell);

// Get data from bgrdtree and store in arrays
//   Int_t size = nrow*ncol;
   Int_t size = bgrdtree->GetEntries();
   Int_t x, y, ij;
   for (Int_t i=0; i<size; i++) {
      bgrdtree->GetEntry(i);

      x  = bgcell->GetX();
      y  = bgcell->GetY();
      ij = XY2Index(x, y, ncol);

      if (inten) inten[ij] = bgcell->GetBackground();
      if (stdev) stdev[ij] = bgcell->GetStdev();
   }//for_i

   return errNoErr;
}//FillBgrdArrays

//______________________________________________________________________________
TTree *XGCProcesSet::FillDataTree(TTree *oldtree, const char *exten,
                     XAlgorithm *algorithm, Int_t nrow, Int_t ncol,
                     Double_t *inten, Double_t *stdev)
{
   // Fill data tree with array arr, write to file and return new datatree
   // Note: Get (X,Y), StdDev and NumPix from oldtree
   if(kCS) cout << "------XGCProcesSet::FillDataTree(newtree)------" << endl;

   if (oldtree == 0) return 0;

   Int_t size = nrow*ncol;

// Init branch address for oldtree
   XGCCell *oldcell = 0;
   oldtree->SetBranchAddress("DataBranch", &oldcell);

// Create tree newtree
   TString name = Path2Name(oldtree->GetName(),"",".");
   name = name + "." + exten;
   TTree *newtree = new TTree(name, fSchemeName.Data());
   if (newtree == 0) return 0;

   Int_t    split   = 99;
   XGCCell *newcell = 0;
   newcell = new XGCCell();
   newtree->Branch("DataBranch", "XGCCell", &newcell, 64000, split);

   Double_t min    = DBL_MAX;  //defined in float.h
   Double_t max    = 0;
   Int_t    nummin = 0;
   Int_t    nummax = 0;
   Int_t    numpix = 0;
   Int_t    maxpix = 0;

   Int_t x, y, ij;
   for (Int_t i=0; i<size; i++) {
      oldtree->GetEntry(i);

      x  = oldcell->GetX();
      y  = oldcell->GetY();
      ij = XY2Index(x, y, ncol);

      // number of cells with minimal intensity
      if      (inten[ij] <  min) {min = inten[ij]; nummin = 1;}
      else if (inten[ij] == min) {nummin++;}

      // number of cells with maximal intensity
      if      (inten[ij] >  max) {max = inten[ij]; nummax = 1;}
      else if (inten[ij] == max) {nummax++;}

      // maximal pixel number
      numpix = oldcell->GetNumPixels();
      if (numpix > maxpix) {maxpix = numpix;}

      newcell->SetX(x);
      newcell->SetY(y);
      newcell->SetIntensity(inten[ij]);
      newcell->SetStdev((stdev == 0) ? oldcell->GetStdev() : stdev[ij]);
      newcell->SetNumPixels((Short_t)numpix);

      newtree->Fill();
   }//for_i

// Add tree info to tree
   AddDataTreeInfo(newtree, newtree->GetName(), algorithm->GetOption(),
                   nrow, ncol, nummin, min, nummax, max, maxpix);

//////////////////
// to do: check if exten is kExtenBgrd or kExtenIntn or kExtenCNrm and save 
// tree only for e.g.: save=fSaveBgrdTree(==TRUE) or safe=fSaveIntenTree(==TRUE)
// or save=fSaveCNrmTree(==TRUE)
// option if (save==kTRUE) WriteTree()
// ev if (save2tmpfile) WriteTree(to file tmp_data.root)
//////////////////

// Write expression tree to file 
//   if (newtree->Write("", TObject::kOverwrite) > 0) {
   if (WriteTree(newtree, TObject::kOverwrite) == errNoErr) {
   // Delete tree header in case of overwrite
      XTreeHeader *header = GetTreeHeader(name);
      if (header) {fHeaders->Remove(header); delete header;}  //????

   // Add tree header to list
//tmp      if tmpfile then fSelections->Add()??? else fHeaders->Add() ??
//tmp      AddTreeHeader(newtree->GetName(), algorithm->GetName(),, useTmpFile,
///////////////
//tmp: TEST:
      if (algorithm->GetFile() == 0)
///////////////
      AddTreeHeader(newtree->GetName(), algorithm->GetName(), 0,
                    algorithm->GetNumParameters(), algorithm->GetParameters());
   }//if

   return newtree;
}//FillDataTree

//______________________________________________________________________________
TTree *XGCProcesSet::FillDataTree(const char *name, XAlgorithm *algorithm, 
                     Int_t nrow, Int_t ncol, Double_t *arr)
{
   // Fill data tree with array arr and write to file
   if(kCS) cout << "------XGCProcesSet::FillDataTree------" << endl;

   TTree *tree = new TTree(name, fSchemeName.Data());
   if (tree == 0) return 0;

   Int_t    split = 99;
   XGCCell *cell  = 0;
   cell = new XGCCell();
   tree->Branch("DataBranch", "XGCCell", &cell, 64000, split);

   Double_t min    = DBL_MAX;  //defined in float.h
   Double_t max    = 0;
   Int_t    nummin = 0;
   Int_t    nummax = 0;

// Fill data tree
   Int_t i, j, ij;
   for (j=0; j<nrow; j++) {
      for (i=0; i<ncol; i++) {
         ij = XY2Index(i, j, ncol);

         // number of cells with minimal intensity
         if      (arr[ij] <  min) {min = arr[ij]; nummin = 1;}
         else if (arr[ij] == min) {nummin++;}

         // number of cells with maximal intensity
         if      (arr[ij] >  max) {max = arr[ij]; nummax = 1;}
         else if (arr[ij] == max) {nummax++;}

         cell->SetX(i);
         cell->SetY(j);
         cell->SetIntensity(arr[ij]);
         cell->SetStdev(0.0);
         cell->SetNumPixels(0);

         tree->Fill();
      }//for_i
   }//for_j

// Add tree info to tree
   AddDataTreeInfo(tree, name, algorithm->GetOption(), nrow, ncol,
                   nummin, min, nummax, max, 0);

//////////////////
// to do: check if exten(from name) is kExtenBgrd or kExtenIntn or kExtenCNrm and save 
// tree only for e.g.: save=fSaveBgrdTree(==TRUE) or safe=fSaveIntenTree(==TRUE)
// or save=fSaveCNrmTree(==TRUE)
// option if (save==kTRUE) WriteTree()
// ev if (save2tmpfile) WriteTree(to file tmp_data.root)
//////////////////

// Write data tree to file 
   if (WriteTree(tree, TObject::kOverwrite) == errNoErr) {
   // Delete tree header in case of overwrite
      XTreeHeader *header = GetTreeHeader(name);
      if (header) {fHeaders->Remove(header); delete header;}  //????

   // Write tree header to list
//tmp      if tmpfile then fSelections->Add()??? else fHeaders->Add() ??
//tmp      AddTreeHeader(tree->GetName(), algorithm->GetName(),, useTmpFile,
///////////////
//tmp: test:
      if (algorithm->GetFile() == 0)
///////////////
      AddTreeHeader(tree->GetName(), algorithm->GetName(), 0,
                    algorithm->GetNumParameters(), algorithm->GetParameters());
   }//if

   return tree;
}//FillDataTree

//______________________________________________________________________________
Int_t XGCProcesSet::FillMaskTree(const char *name, XAlgorithm *algorithm, 
                    Int_t nrow, Int_t ncol, Int_t *arr)
{
   // Fill mask tree with array arr and write to file
   if(kCS) cout << "------XGCProcesSet::FillMaskTree------" << endl;

   Int_t err = errNoErr;

   TTree *tree = new TTree(name, fSchemeName.Data());
   if (tree == 0) return errCreateTree;

   Int_t     split = 99;
   XCellMask *mask = 0;
   mask = new XCellMask();
   tree->Branch("MaskBranch", "XCellMask", &mask, 64000, split);

// Fill mask tree
   Int_t i, j, ij;
   Int_t nflags = 0;
   for (j=0; j<nrow; j++) {
      for (i=0; i<ncol; i++) {
         ij = XY2Index(i, j, ncol);

         // get number of positive flags
         if (arr[ij] > 0) nflags++;

         mask->SetX(i);
         mask->SetY(j);
         mask->SetFlag((Short_t)arr[ij]); //convert to Short_t for fFlag
         tree->Fill();
      }//for_i
   }//for_j

// Add tree info to tree
   AddMaskTreeInfo(tree, name, algorithm->GetOption(), nrow, ncol, nflags);

// Write mask tree to file 
   if ((err = WriteTree(tree, TObject::kOverwrite)) == errNoErr) {
   // Delete tree header in case of overwrite
      XTreeHeader *header = GetTreeHeader(name);
      if (header) {fHeaders->Remove(header); delete header;}  //????

   // Write tree header to list
///////////////
//tmp: test:
      if (algorithm->GetFile() == 0)
///////////////
      AddTreeHeader(tree->GetName(), algorithm->GetName(), 0,
                    algorithm->GetNumParameters(), algorithm->GetParameters());
   }//if

// Cleanup
   tree->Delete(""); tree = 0; //delete tree from heap
   delete mask;

   return err;
}//FillMaskTree

//______________________________________________________________________________
Int_t XGCProcesSet::MeanReference(Int_t numdata, TTree **datatree, Int_t numbgrd,
                    TTree **bgrdtree, Int_t nrow, Int_t ncol, Double_t *arr, Bool_t doBg)
{
   // Fill array with trimmed mean values from numdata datatrees
   if(kCS) cout << "------XGCProcesSet::MeanReference------" << endl;

// Init branch addresses
   XBgCell *bgcell[numdata];
   XGCCell *gccell[numdata];
   for (Int_t k=0; k<numdata; k++) {
      bgcell[k] = 0;
      gccell[k] = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell[k]);
      if (numbgrd > 0) bgrdtree[k]->SetBranchAddress("BgrdBranch", &bgcell[k]);
   }//for_k

   Double_t *arrRef = 0;  //to store entries of tree branches
   if (!(arrRef = new Double_t[numdata])) return errInitMemory; 
   for (Int_t i=0; i<numdata;  i++) arrRef[i] = 0.0;

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

// Get data from datatrees (and bgrdtrees) and fill array with mean
   Int_t size = nrow*ncol;
   Int_t x, y, ij;
   if ((numbgrd > 0) && (doBg == kTRUE)) {
      for (Int_t i=0; i<size; i++) {
         datatree[0]->GetEntry(i);

         x  = gccell[0]->GetX();
         y  = gccell[0]->GetY();
         ij = XY2Index(x, y, ncol);

         for (Int_t k=0; k<numdata; k++) {
            datatree[k]->GetEntry(i);
            bgrdtree[k]->GetEntry(i);

            // subtract background from intensity
            arrRef[k] = this->AdjustIntensity(gccell[k]->GetIntensity(),
                                              bgcell[k]->GetBackground(),
                                              bgcell[k]->GetStdev());
         }//for_k

         arr[ij] = TStat::Mean(numdata, arrRef, fRefTrim);
      }//for_i
   } else {
      for (Int_t i=0; i<size; i++) {
         datatree[0]->GetEntry(i);

         x  = gccell[0]->GetX();
         y  = gccell[0]->GetY();
         ij = XY2Index(x, y, ncol);

         arrRef[0] = gccell[0]->GetIntensity();
         for (Int_t k=1; k<numdata; k++) {
            datatree[k]->GetEntry(i);
            arrRef[k] = gccell[k]->GetIntensity();
         }//for_k

         arr[ij] = TStat::Mean(numdata, arrRef, fRefTrim);
      }//for_i
   }//if

   delete [] arrRef;
   return errNoErr;
}//MeanReference

//______________________________________________________________________________
Int_t XGCProcesSet::MedianReference(Int_t numdata, TTree **datatree, Int_t numbgrd,
                    TTree **bgrdtree, Int_t nrow, Int_t ncol, Double_t *arr, Bool_t doBg)
{
   // Fill array with median values from numdata datatrees
   if(kCS) cout << "------XGCProcesSet::MedianReference------" << endl;

   fRefTrim = 0.5;
   return this->MeanReference(numdata, datatree, numbgrd, bgrdtree,  nrow, ncol, arr, doBg);
}//MedianReference

//______________________________________________________________________________
Int_t XGCProcesSet::DoExpress(Int_t numdata, TTree **datatree,
                    Int_t numbgrd, TTree **bgrdtree)
{
   // Condense data from  datatrees to expression values
   // Note: each datatree is processed individually and thus can have different
   //       number of entries (i.e. belong to different chip types)
   if(kCS) cout << "------XGCProcesSet::DoExpress------" << endl;

   Int_t err = errNoErr;
   Int_t split = 99;

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fExpressor, fExpressor->GetBgrdOption());

// Check for presence of background trees
   if (numbgrd == 0) {
      if (doBg == kTRUE) {
         cout << "Warning: No background trees available for background subtraction."
              << endl;
      }//if
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Init local arrays to store data from trees
   Int_t    *arrMask  = 0;
   Double_t *arrInten = 0;
   Double_t *arrStdev = 0;
   Int_t    *arrNPix  = 0;
   Double_t *arrPM    = 0; 
   Double_t *arrMM    = 0;
   Double_t *arrSP    = 0; 
   Double_t *arrSM    = 0;
   Int_t    *arrXP    = 0; 
   Int_t    *arrXM    = 0;

// Calculate expression
   Int_t i, j, ij, x, y, start, end;
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) {err = errGetTree; break;}

   // Informing user
      if (XManager::fgVerbose) {
         cout << "      summarizing <" << datatree[k]->GetName() << ">..." << endl;
      }//if

   // Get tree info from datatree
      fDataFile->cd();
      TString name = datatree[k]->GetName();

      XDataTreeInfo *info = 0;
      info = (XDataTreeInfo*)datatree[k]->GetUserInfo()->FindObject(name);
      if (!info) {
         cerr << "Error: Could not get tree info for <" << name << ">." << endl;
         err = errGeneral;
         break;
      }//if

   // Get chip parameters from scheme file (also alternative CDFs)
      if (strcmp(fSchemeName.Data(), "") == 0) {
         fSchemeName = datatree[k]->GetTitle();
      } else if (!fSchemeName.Contains(datatree[k]->GetTitle())) {
         return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[k]->GetTitle());
      }//if
      if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

      XDNAChip *chip = (XDNAChip*)fSchemes->FindObject(fSchemeName, kTRUE);
      if (!chip) {
         err = fManager->HandleError(errGetScheme, fSchemeName);
         break;
      }//if
      Int_t numrows  = chip->GetNumRows();
      Int_t numcols  = chip->GetNumColumns();
      Int_t numctrls = chip->GetNumControls();
      Int_t numunits = chip->GetNumUnits();

   // Init min/max expression levels
      Double_t min = DBL_MAX;  //defined in float.h
      Double_t max = 0;

   // Get scheme tree for scheme
      XScheme *scheme = 0;
      TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
      if (scmtree == 0) {err = errGetTree; break;}
      scmtree->SetBranchAddress("ScmBranch", &scheme);

   // Get unit tree for scheme
      XGCUnit *unit = 0;
      TTree *idxtree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
      if (idxtree == 0) {err = errGetTree; break;}
      idxtree->SetBranchAddress("IdxBranch", &unit);

   // Get maximum number of pairs from tree info for unit tree
      XUnitTreeInfo *idxinfo = 0;
      idxinfo = (XUnitTreeInfo*)idxtree->GetUserInfo()->FindObject(idxtree->GetName());
      if (!idxinfo) {
         cerr << "Error: Could not get tree info for <" << idxtree->GetName() << ">." << endl;
         err = errGeneral;
         break;
      }//if
      Int_t maxnumpairs = (Int_t)idxinfo->GetValue("fMaxNPairs");

   // Create new tree exprtree
      if (!fFile->cd(fName)) {err = errGetDir; break;}

      TString dataname = Path2Name(datatree[k]->GetName(),"/",".");
      TString exprname = dataname + "." + fExpressor->GetTitle();
      TTree  *exprtree = new TTree(exprname, fSchemeName);
      if (exprtree == 0) {err = errCreateTree; break;}

      XGCExpression *expr = 0;
      expr = new XGCExpression();
      exprtree->Branch("ExprBranch", "XGCExpression", &expr, 64000, split);

   // Initialize memory for data arrays
      Int_t size = numrows*numcols;
      if (!(arrMask  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
      if (!(arrInten = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrStdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrNPix  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}

      for (i=0; i<size; i++) {
         arrMask[i]  = eINITMASK;  //init mask
         arrInten[i] = arrStdev[i] = 0.0;
         arrNPix[i]  = 0;
      }//for_i

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

   // Get mask for PM/MM from scheme tree and store in array 
      arrMask = this->FillMaskArray(chip, scmtree, scheme, 0, size, arrMask); //level=0?

   // Calculate mask for expression
      err = fExprSelector->Calculate(size, 0, 0, arrMask);
      if (err != errNoErr) goto cleanup;

   // Get data from datatree and store in arrays
      err = this->FillDataArrays(datatree[k], bgrdtree[k], doBg,
                       numrows, numcols, arrInten, arrStdev, arrNPix);
      if (err != errNoErr) goto cleanup;

   // Initialize maximum memory for PM/MM arrays (maxnumpairs+1 to avoid potential buffer overflow)
      if (!(arrPM = new (nothrow) Double_t[maxnumpairs+1])) {err = errInitMemory; goto cleanup;}
      if (!(arrMM = new (nothrow) Double_t[maxnumpairs+1])) {err = errInitMemory; goto cleanup;}
      if (!(arrSP = new (nothrow) Double_t[maxnumpairs+1])) {err = errInitMemory; goto cleanup;}
      if (!(arrSM = new (nothrow) Double_t[maxnumpairs+1])) {err = errInitMemory; goto cleanup;}
      if (!(arrXP = new (nothrow) Int_t[maxnumpairs+1]))    {err = errInitMemory; goto cleanup;}
      if (!(arrXM = new (nothrow) Int_t[maxnumpairs+1]))    {err = errInitMemory; goto cleanup;}

   // Calculate expression values
      start = 0;
      end   = 0;
      for (Int_t id=0; id<numunits; id++) { 
         idxtree->GetEntry(id);

         Int_t unitID   = id - numctrls;
         Int_t numcells = unit->GetNumCells();
         // skip negative unit entries (controls)
         if (unit->GetUnitID() < 0) {
            start += numcells;
            end = start;
            continue;
         }//if

         Int_t p = 0;
         Int_t m = 0;
         end += numcells;
         for (j=start; j<end; j++) {
            scmtree->GetEntry(j);

            if ((scheme->GetUnitID()) != unitID) {
               cerr << "Error: unitID is not equal to: " << unitID << endl;
               err = errAbort;
               goto cleanup;
            }//if

            x  = scheme->GetX();
            y  = scheme->GetY();
            ij = XY2Index(x, y, numcols);

            if (arrMask[ij] == 1) {
               arrPM[p] = arrInten[ij];
               arrSP[p] = arrStdev[ij];
               arrXP[p] = arrNPix[ij];
               p++;
            } else if (arrMask[ij] == 0) {
               arrMM[m] = arrInten[ij];
               arrSM[m] = arrStdev[ij];
               arrXM[m] = arrNPix[ij];
               m++;
            }//if

            if (p > maxnumpairs || m > maxnumpairs) {
               cerr << "Error: unitID <" << unitID << "> exceeds maximum number of pairs <"
                    << maxnumpairs << ">. Buffer overflow!" << endl;
               err = errAbort;
               goto cleanup;
            }//if
         }//for_j
         start += numcells;

         // continue if arrays arrPM etc are not filled
         if (p == 0) continue;
         if (p != m) {
            cout << "Warning: Skipping unitID <" << unitID
                 << "> with different numbers of PM and MM data." << endl;
//            continue;
            err = errAbort;
            goto cleanup;
         }//if

         // calculate mean expression level
         Int_t    arrlen = 0;
         Double_t mean   = 0;
         Double_t var    = 0;
         fExpressor->InitTreeInfo(info);
         fExpressor->InitArrays(p, arrPM, arrSP, arrXP, arrMM, arrSM, arrXM);
         if ((err = fExpressor->CreateArray(p)))               goto cleanup;
         if ((err = fExpressor->Calculate(mean, var, arrlen))) goto cleanup;
         fExpressor->DeleteArray();

         // get minimal/maximal expression levels
         if (mean < min) min = mean;
         if (mean > max) max = mean;

         // fill expression tree
         expr->SetUnitID(unitID);
         expr->SetLevel(mean);
         expr->SetStdev(TMath::Sqrt(var));
         expr->SetNumPairs(arrlen);
         exprtree->Fill();
      }//for_id

      if (XManager::fgVerbose) {
         cout << "      expression statistics: " << endl;
         cout << "         minimal expression level is <" << min << ">." << endl;
         cout << "         maximal expression level is <" << max << ">." << endl;
      }//if

   // Add tree info to tree
      AddExprTreeInfo(exprtree, exprtree->GetName(), fExpressor->GetOption(),
                      numunits-numctrls, min, max);

   // Write expression tree to file 
      if ((err = WriteTree(exprtree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(exprtree->GetName(), "Expr", 0, fExpressor->GetNumParameters(),
                       fExpressor->GetParameters());
      }//if

   cleanup:
      // delete arrays
      if (arrXM)    {delete [] arrXM;    arrXM    = 0;}
      if (arrXP)    {delete [] arrXP;    arrXP    = 0;}
      if (arrSM)    {delete [] arrSM;    arrSM    = 0;}
      if (arrSP)    {delete [] arrSP;    arrSP    = 0;}
      if (arrMM)    {delete [] arrMM;    arrMM    = 0;}
      if (arrPM)    {delete [] arrPM;    arrPM    = 0;}
      if (arrNPix)  {delete [] arrNPix;  arrNPix  = 0;}
      if (arrStdev) {delete [] arrStdev; arrStdev = 0;}
      if (arrInten) {delete [] arrInten; arrInten = 0;}
      if (arrMask)  {delete [] arrMask;  arrMask  = 0;}
//?      // delete scheme tree from RAM
//?      if (scmtree)  {scmtree->Delete(""); scmtree = 0;}
      // Note: do not remove exprtree and expr from RAM, needed later

      if (err != errNoErr) break;
   }//for_k

   return err;
}//DoExpress

//______________________________________________________________________________
Int_t XGCProcesSet::DoMedianPolish(Int_t numdata, TTree **datatree,
                    Int_t numbgrd, TTree **bgrdtree)
{
   // Compute expression values using mediapolish
   // Intensities from datatrees are stored in one large table in RAM for fast access
   // Note: all trees must have same number of entries (i.e. identical chip types)
   if(kCS) cout << "------XGCProcesSet::DoMedianPolish(table)------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "      Summarizing with medianpolish..." << endl;
   }//if

//TEST Benchmark
//gBenchmark->Reset(); 
//gBenchmark->Start("Bench_MedianPolish");

   Int_t x, y, ij, start, end, id;
   Int_t idx = 0;
   Int_t err = errNoErr;

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fExpressor, fExpressor->GetBgrdOption());

// Get chip parameters from scheme file (also alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XDNAChip *chip = (XDNAChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) {
      return fManager->HandleError(errGetScheme, fSchemeName);
   }//if
   Int_t numrows  = chip->GetNumRows();
   Int_t numcols  = chip->GetNumColumns();
   Int_t numctrls = chip->GetNumControls();
   Int_t numunits = chip->GetNumUnits();

// Get scheme tree for scheme
   XScheme *scheme = 0;
   TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
   if (scmtree == 0) return errGetTree;
   scmtree->SetBranchAddress("ScmBranch", &scheme);

// Get unit tree for scheme
   XGCUnit *unit = 0;
   TTree *idxtree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   if (idxtree == 0) return errGetTree;
   idxtree->SetBranchAddress("IdxBranch", &unit);

// Init size of arrays
   Int_t size     = numrows*numcols;
   Int_t numsels  = 0;  //number of selected entries

// Check for equal number of data tree entries
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;

      if ((Int_t)(datatree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, datatree[k]->GetName(), str);
      }//if
   }//for_k

// Check for equal number of background tree entries
   if (numbgrd > 0) {
      for (Int_t k=0; k<numdata; k++) {
         if (bgrdtree[k] == 0) return errGetTree;

         if ((Int_t)(bgrdtree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, bgrdtree[k]->GetName(), str);
         }//if
      }//for_k
   } else if (doBg == kTRUE) {
      cout << "Warning: No background trees available for background subtraction."
           << endl;
      doBg = kFALSE;
   } else {
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Init branch addresses
   XBgCell *bgcell[numdata];
   XGCCell *gccell[numdata];
   for (Int_t k=0; k<numdata; k++) {
      bgcell[k] = 0;
      gccell[k] = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell[k]);
      if (numbgrd > 0) bgrdtree[k]->SetBranchAddress("BgrdBranch", &bgcell[k]);
   }//for_k

// Init expression trees
   Int_t  split = 99;
   TTree *exprtree[numdata];
   XGCExpression *expr[numdata];

// Init min/max expression levels
   Double_t min = DBL_MAX;  //defined in float.h
   Double_t max = 0;

// Init local arrays
   Int_t     *arrMask = 0;
   Int_t     *arrIndx = 0;
   Double_t  *colmed  = 0;
   Double_t  *results = 0;
   Double_t **table   = 0;

// Create local arrays
   if (!(arrMask = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrIndx = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(colmed  = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}
   if (!(results = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<size; i++) {
      arrMask[i] = eINITMASK; 
      arrIndx[i] = 0;
   }//for_i
   for (Int_t i=0; i<numdata; i++) colmed[i] = results[i] = 0.0; 

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

// Get mask for PM from scheme tree and store in array 
   arrMask = this->FillMaskArray(chip, scmtree, scheme, 0, size, arrMask);

// Calculate mask for expression
   err = fExprSelector->Calculate(size, 0, 0, arrMask);
   if (err != errNoErr) goto cleanup;

// Get entry index from datatree
   idx = 0;
   for (Int_t i=0; i<size; i++) {
      datatree[0]->GetEntry(i);

      x  = gccell[0]->GetX();
      y  = gccell[0]->GetY();
      ij = XY2Index(x, y, numcols);

      if (arrMask[ij] == 1) {
         arrIndx[ij] = idx++;
      }//if
   }//for_i

// Get number of selected entries
   for (Int_t i=0; i<size; i++) {
      numsels = (arrMask[i] == 1) ? ++numsels : numsels;
   }//for_i

// Create table to store selected intensities of all datatrees
   if (!(table = new (nothrow) Double_t*[numdata])) {err = errInitMemory; goto cleanup;} 
   for (Int_t k=0; k<numdata; k++) {
      table[k] = 0;
      if (!(table[k] = new (nothrow) Double_t[numsels])) {err = errInitMemory; goto cleanup;} 
   }//for_i

//TEST
//gBenchmark->Start("Bench_table");
// Get data from datatrees (and bgrdtrees) and store in table
   if ((numbgrd > 0) && (doBg == kTRUE)) {
      for (Int_t k=0; k<numdata; k++) {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);
            bgrdtree[k]->GetEntry(i);
//unsafe?? use arrInten[ij] and arrBgrd[ij] ??

            // subtract background from intensity
            if (arrMask[i] == 1) {
               table[k][idx++] = this->AdjustIntensity(gccell[k]->GetIntensity(),
                                                       bgcell[k]->GetBackground(),
                                                       bgcell[k]->GetStdev());
            }//if
         }//for_i
      }//for_k
   } else {
      for (Int_t k=0; k<numdata; k++) {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);

            if (arrMask[i] == 1) {
               table[k][idx++] = gccell[k]->GetIntensity();
            }//if
         }//for_i
      }//for_k
   }//if
//TEST Benchmark
//gBenchmark->Show("Bench_table");

// Create new trees exprtree
   if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

   for (Int_t k=0; k<numdata; k++) {
      TString dataname = Path2Name(datatree[k]->GetName(),"/",".");
      TString exprname = dataname + "." + fExpressor->GetTitle();
      exprtree[k] = new TTree(exprname, fSchemeName);
      if (exprtree[k] == 0) {err = errCreateTree; goto cleanup;}

      expr[k] = new XGCExpression();
      exprtree[k]->Branch("ExprBranch", "XGCExpression", &expr[k], 64000, split);
   }//for_k

//TEST
//gBenchmark->Show("Bench_MedianPolish");
//gBenchmark->Reset();
//gBenchmark->Start("Bench_Loop");

// Calculate expression values
   start = 0;
   end   = 0;
   idx   = 0;
   for (id=0; id<numunits; id++) { 
      idxtree->GetEntry(id);

      Int_t unitID   = id - numctrls;
      Int_t numcells = unit->GetNumCells();
      // skip negative unit entries (controls)
      if (unit->GetUnitID() < 0) {
         start += numcells;
         end = start;
         continue;
      }//if

      // create array to store PM values for all probes with current unitID
      Int_t  numpairs = (Int_t)(numcells / 2);
      Double_t *arrPM = new Double_t[numpairs*numdata]; 
      for (Int_t i=0; i<numpairs*numdata;  i++) arrPM[i] = 0.0;

      // fill arrPM with PM values of current unitID
      Int_t p = 0;
      end += numcells;
      for (Int_t j=start; j<end; j++) {
         scmtree->GetEntry(j);

         if ((scheme->GetUnitID()) != unitID) {
            cerr << "Error: unitID is not equal to: " << unitID << endl;
            err = errAbort;
            goto cleanup;
         }//if

         x  = scheme->GetX();
         y  = scheme->GetY();
         ij = XY2Index(x, y, numcols);

         if (arrMask[ij] == 1) {
            if (p == 0) idx++;  //count number of units to be summarized

            for (Int_t k=0; k<numdata; k++) {
               arrPM[p] = table[k][arrIndx[ij]];
               p++;
            }//for_k
         }//if
      }//for_j
      start += numcells;

      if (XManager::fgVerbose && id%10000 == 0) {
         cout << "      calculating expression for <" << idx << "> of <"
              << numunits << "> units...\r" << flush;
      }//if

      // fill arrPM or continue if it is not filled
      if ((err = fExpressor->SetArray(p, arrPM)) != errNoErr) {
         delete [] arrPM;
         continue;
      }//if

      // calculate median polish for PMs of current unitID
      if ((err = fExpressor->Calculate(numdata, colmed, results, 0))) break;

//////////////////
// TO DO: return fResiduals for residual-plot!!!! (like affyPLM)
//      residuals = fExpressor->GetResiduals();
// store as residuals tree??
//////////////////

      // fill expression trees
      for (Int_t k=0; k<numdata; k++) {
         // get minimal/maximal expression levels
         if (results[k] < min) min = results[k];
         if (results[k] > max) max = results[k];

         expr[k]->SetUnitID(unitID);
         expr[k]->SetLevel(results[k]);
//??         expr[k]->SetStdev(colmed[k]);
         expr[k]->SetStdev(TMath::Abs(colmed[k]));
         expr[k]->SetNumPairs(numpairs);
         exprtree[k]->Fill();
      }//for_k

      delete [] arrPM;
   }//for_id
   if (XManager::fgVerbose) {
      cout << "      calculating expression for <" << idx << "> of <"
           << numunits << "> units...Finished." << endl;
   }//if

   if (XManager::fgVerbose) {
      cout << "      expression statistics: " << endl;
      cout << "         minimal expression level is <" << min << ">." << endl;
      cout << "         maximal expression level is <" << max << ">." << endl;
   }//if
//TEST
//gBenchmark->Show("Bench_Loop");

// Write expression trees to file 
   for (Int_t k=0; k<numdata; k++) {
   // Add tree info to tree
      AddExprTreeInfo(exprtree[k], exprtree[k]->GetName(), fExpressor->GetOption(),
                      idx, min, max);

      if ((err = WriteTree(exprtree[k], TObject::kOverwrite)) == errNoErr) {
//      if ((err = WriteTree(exprtree[k], 0)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(exprtree[k]->GetName(), "Expr", 0, fExpressor->GetNumParameters(),
                       fExpressor->GetParameters());
      } else {
         break;
      }//if
   }//for_k

// Cleanup
cleanup:
   // delete table
   for (Int_t k=0; k<numdata; k++) {
      if (table[k]) {delete [] table[k]; table[k] = 0;}
   }//for_k
   if (table) delete [] table;

   // delete arrays
   if (results) {delete [] results; results = 0;}
   if (colmed)  {delete [] colmed;  colmed  = 0;}
   if (arrIndx) {delete [] arrIndx; arrIndx = 0;}
   if (arrMask) {delete [] arrMask; arrMask = 0;}

   // delete scheme tree from RAM
//?   if (scmtree)  {scmtree->Delete(""); scmtree = 0;}

   return err;
}//DoMedianPolish

//______________________________________________________________________________
Int_t XGCProcesSet::DoMedianPolish(Int_t numdata, TTree **datatree,
                    Int_t numbgrd, TTree **bgrdtree, TFile *file)
{
   // Compute expression values using mediapolish
   // In order to reduce memory consumption intensities from datatrees are
   // stored in the order of the corresponding schemetree entries in temporary 
   // trees in a temporary root file
   // Note: all trees must have same number of entries (i.e. identical chip types)
   if(kCS) cout << "------XGCProcesSet::DoMedianPolish(file)------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "      Summarizing with medianpolish (using temporary file)..." << endl;
   }//if

//TEST
//gBenchmark->Reset(); 
//gBenchmark->Start("Bench_MedianPolish");

   Int_t x, y, ij, start, end, id, entry;
   Int_t idx = 0;
   Int_t err = errNoErr;

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fExpressor, fExpressor->GetBgrdOption());

// Get chip parameters from scheme file (also alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XDNAChip *chip = (XDNAChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) {
      return fManager->HandleError(errGetScheme, fSchemeName);
   }//if
   Int_t numrows  = chip->GetNumRows();
   Int_t numcols  = chip->GetNumColumns();
   Int_t numctrls = chip->GetNumControls();
   Int_t numunits = chip->GetNumUnits();

// Get scheme tree for scheme
   XScheme *scheme = 0;
   TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
   if (scmtree == 0) return errGetTree;
   scmtree->SetBranchAddress("ScmBranch", &scheme);

// Get unit tree for scheme
   XGCUnit *unit = 0;
   TTree *idxtree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   if (idxtree == 0) return errGetTree;
   idxtree->SetBranchAddress("IdxBranch", &unit);

// Init size of arrays
   Int_t size     = numrows*numcols;
   Int_t numsels  = 0;  //number of selected entries

// Check for equal number of data tree entries
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;

      if ((Int_t)(datatree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, datatree[k]->GetName(), str);
      }//if
   }//for_k

// Check for equal number of background tree entries
   if (numbgrd > 0) {
      for (Int_t k=0; k<numdata; k++) {
         if (bgrdtree[k] == 0) return errGetTree;

         if ((Int_t)(bgrdtree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, bgrdtree[k]->GetName(), str);
         }//if
      }//for_k
   } else if (doBg == kTRUE) {
      cout << "Warning: No background trees available for background subtraction."
           << endl;
      doBg = kFALSE;
   } else {
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Init branch addresses
   XBgCell *bgcell[numdata];
   XGCCell *gccell[numdata];
   for (Int_t k=0; k<numdata; k++) {
      bgcell[k] = 0;
      gccell[k] = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell[k]);
      if (numbgrd > 0) bgrdtree[k]->SetBranchAddress("BgrdBranch", &bgcell[k]);
   }//for_k

// Init temporary trees and expression trees
   Int_t    split = 99;
   Double_t sort  = 0.0;
   TTree *tmptree[numdata];
   TTree *exprtree[numdata];
   XGCExpression *expr[numdata];

// Init min/max expression levels
   Double_t min = DBL_MAX;  //defined in float.h
   Double_t max = 0;

// Init local arrays
   Int_t    *arrMask = 0;
   Int_t    *arrIndx = 0;
   Double_t *arrData = 0;
   Double_t *colmed  = 0;
   Double_t *results = 0;

// Create local arrays
   if (!(arrMask = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrIndx = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(colmed  = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}
   if (!(results = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<size; i++) {
      arrMask[i] = eINITMASK; 
      arrIndx[i] = 0;
   }//for_i 
   for (Int_t i=0; i<numdata; i++) colmed[i] = results[i] = 0.0; 

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

// Get mask for PM from scheme tree and store in array 
   arrMask = this->FillMaskArray(chip, scmtree, scheme, 0, size, arrMask);

// Calculate mask for expression
   err = fExprSelector->Calculate(size, 0, 0, arrMask);
   if (err != errNoErr) goto cleanup;

// Get entry index from datatree
   idx = 0;
   for (Int_t i=0; i<size; i++) {
      datatree[0]->GetEntry(i);

      x  = gccell[0]->GetX();
      y  = gccell[0]->GetY();
      ij = XY2Index(x, y, numcols);

      if (arrMask[ij] == 1) {
         arrIndx[ij] = idx++;
      }//if
   }//for_i

// Get number of selected entries
   for (Int_t i=0; i<size; i++) {
      numsels = (arrMask[i] == 1) ? ++numsels : numsels;
   }//for_i

// Create array to store selected intensities
   if (!(arrData = new (nothrow) Double_t[numsels])) {err = errInitMemory; goto cleanup;}
   for (Int_t i=0; i<numsels;  i++) arrData[i] = 0.0;

//TEST
//gBenchmark->Start("Bench_tmptree");
// Change directory to temporary file
   if (!file->cd()) {err = errGetDir; goto cleanup;}

// Get data from datatrees (and bgrdtrees) and store in temporary file
   for (Int_t k=0; k<numdata; k++) {
      // create temporary tree
      tmptree[k] = new TTree(datatree[k]->GetName(), "temporary tree");
      if (tmptree[k] == 0) {err = errCreateTree; goto cleanup;}
      tmptree[k]->Branch("sortBr", &sort, "sort/D");

      // fill array with (background corrected) intensities
      if ((numbgrd > 0) && (doBg == kTRUE)) {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);
            bgrdtree[k]->GetEntry(i);

            // subtract background from intensity
            if (arrMask[i] == 1) {
               arrData[idx++] = this->AdjustIntensity(gccell[k]->GetIntensity(),
                                                      bgcell[k]->GetBackground(),
                                                      bgcell[k]->GetStdev());
            }//if
         }//for_i
      } else {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);

            if (arrMask[i] == 1) {
               arrData[idx++] = gccell[k]->GetIntensity();
            }//if
         }//for_i
      }//if

      // fill tmptree with array in the order of scheme tree entries for (x,y)
      for (Int_t i=0; i<size; i++) {
         scmtree->GetEntry(i);

         x  = scheme->GetX();
         y  = scheme->GetY();
         ij = XY2Index(x, y, numcols);

         if (arrMask[ij] == 1) {
            sort = arrData[arrIndx[ij]];
            tmptree[k]->Fill();
         }//if
      }//for_i

      // write tmptree to temporary file
      tmptree[k]->Write();
//??      tmptree[k]->Write(TObject::kOverwrite);
   }//for_k
//TEST
//gBenchmark->Show("Bench_tmptree");

// Change directory to current directory for treeset in main file
   if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

// Create new trees exprtree
   for (Int_t k=0; k<numdata; k++) {
      TString dataname = Path2Name(datatree[k]->GetName(),"/",".");
      TString exprname = dataname + "." + fExpressor->GetTitle();
      exprtree[k] = new TTree(exprname, fSchemeName);
      if (exprtree[k] == 0) {err = errCreateTree; goto cleanup;}

      expr[k] = new XGCExpression();
      exprtree[k]->Branch("ExprBranch", "XGCExpression", &expr[k], 64000, split);
   }//for_k

//TEST
//gBenchmark->Show("Bench_MedianPolish");
//gBenchmark->Reset();
//gBenchmark->Start("Bench_Loop");

// Calculate expression values
   start = 0;
   end   = 0;
   idx   = 0;
   entry = 0;
   for (id=0; id<numunits; id++) { 
      idxtree->GetEntry(id);

      Int_t unitID   = id - numctrls;
      Int_t numcells = unit->GetNumCells();
      // skip negative unit entries (controls)
      if (unit->GetUnitID() < 0) {
         start += numcells;
         end = start;
         continue;
      }//if

      // create array to store PM values for all probes with current unitID
      Int_t  numpairs = (Int_t)(numcells / 2);
      Double_t *arrPM = new Double_t[numpairs*numdata]; 

      // fill arrPM with PM values of current unitID
      Int_t p = 0;
      end += numcells;
      for (Int_t j=start; j<end; j++) {
         scmtree->GetEntry(j);

         if ((scheme->GetUnitID()) != unitID) {
            cerr << "Error: unitID is not equal to: " << unitID << endl;
            err = errAbort;
            goto cleanup;
         }//if

         x  = scheme->GetX();
         y  = scheme->GetY();
         ij = XY2Index(x, y, numcols);

         if (arrMask[ij] == 1) {
            if (p == 0) idx++;  //count number of units to be summarized

            for (Int_t k=0; k<numdata; k++) {
               tmptree[k]->GetEntry(entry);
               arrPM[p] = sort;
               p++;
            }//for_k

            entry++;
         }//if
      }//for_j
      start += numcells;

      if (XManager::fgVerbose && id%10000 == 0) {
         cout << "      calculating expression for <" << idx << "> of <"
              << numunits << "> units...\r" << flush;
      }//if

      // fill arrPM or continue if it is not filled
      if ((err = fExpressor->SetArray(p, arrPM)) != errNoErr) {
         delete [] arrPM;
         continue;
      }//if

      // calculate median polish for PMs of current unitID
      if ((err = fExpressor->Calculate(numdata, colmed, results, 0))) break;

//////////////////
// TO DO: return fResiduals for residual-plot!!!! (like affyPLM)
//      residuals = fExpressor->GetResiduals();
// store as residuals tree??
//////////////////

      // fill expression trees
      for (Int_t k=0; k<numdata; k++) {
         // get minimal/maximal expression levels
         if (results[k] < min) min = results[k];
         if (results[k] > max) max = results[k];

         expr[k]->SetUnitID(unitID);
         expr[k]->SetLevel(results[k]);
//??         expr[k]->SetStdev(colmed[k]);
         expr[k]->SetStdev(TMath::Abs(colmed[k]));
         expr[k]->SetNumPairs(numpairs);
         exprtree[k]->Fill();
      }//for_k

      delete [] arrPM;
   }//for_id
   if (XManager::fgVerbose) {
      cout << "      calculating expression for <" << idx << "> of <"
           << numunits << "> units...Finished." << endl;
   }//if

   if (XManager::fgVerbose) {
      cout << "      expression statistics: " << endl;
      cout << "         minimal expression level is <" << min << ">." << endl;
      cout << "         maximal expression level is <" << max << ">." << endl;
   }//if
//TEST Benchmark
//gBenchmark->Show("Bench_Loop");

// Write expression trees to file 
   for (Int_t k=0; k<numdata; k++) {
   // Add tree info to tree
      AddExprTreeInfo(exprtree[k], exprtree[k]->GetName(), fExpressor->GetOption(),
                      idx, min, max);

      if ((err = WriteTree(exprtree[k], TObject::kOverwrite)) == errNoErr) {
//      if ((err = WriteTree(exprtree[k], 0)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(exprtree[k]->GetName(), "Expr", 0, fExpressor->GetNumParameters(),
                       fExpressor->GetParameters());
      } else {
         break;
      }//if
   }//for_k

// Cleanup
cleanup:
   // delete temporary trees
   for (Int_t k=0; k<numdata; k++) {
      tmptree[k]->Delete(""); tmptree[k] = 0;
   }//for_k

   // delete arrays
   if (arrData) {delete [] arrData; arrData = 0;}
   if (results) {delete [] results; results = 0;}
   if (colmed)  {delete [] colmed;  colmed  = 0;}
   if (arrIndx) {delete [] arrIndx; arrIndx = 0;}
   if (arrMask) {delete [] arrMask; arrMask = 0;}

   // delete scheme tree from RAM
//?   if (scmtree)  {scmtree->Delete(""); scmtree = 0;}

   return err;
}//DoMedianPolish


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeProcesSet                                                     //
//                                                                      //
// Class for GenomeChip oligonucleotide array pre-processing            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGenomeProcesSet::XGenomeProcesSet()
                 :XGCProcesSet()
{
   // Default GenomeProcesSet constructor
   if(kCS) cout << "---XGenomeProcesSet::XGenomeProcesSet(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XGenomeProcesSet::XGenomeProcesSet(const char *name, const char *type)
                 :XGCProcesSet(name, type)
{
   // Normal GenomeProcesSet constructor
   if(kCS) cout << "---XGenomeProcesSet::XGenomeProcesSet------" << endl;

}//Constructor

//______________________________________________________________________________
XGenomeProcesSet::~XGenomeProcesSet()
{
   // GenomeProcesSet destructor
   if(kCS) cout << "---XGenomeProcesSet::~XGenomeProcesSet------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XGenomeProcesSet::DetectCall(Int_t numdata, TTree **datatree,
                        Int_t &numbgrd, TTree **bgrdtree)
{
   // Detect presence call
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   // Note: each datatree is processed individually and thus can have different
   //       number of entries (i.e. belong to different chip types)
   if(kCS) cout << "------XGenomeProcesSet::DetectCall------" << endl;

   Int_t err = errNoErr;

   fFile->cd();

// Init local arrays to store data from trees
   Int_t    *arrUnit  = 0;
   Int_t    *mskUnit  = 0;
   Int_t    *arrMask  = 0;
   Double_t *arrInten = 0;
   Double_t *arrStdev = 0;  //not needed for current fCaller
   Int_t    *arrNPix  = 0;  //not needed for current fCaller
   Double_t *arrBgrd  = 0;
   Double_t *arrBgdev = 0;  //not needed for current fCaller
   Double_t *arrPM    = 0; 
   Double_t *arrMM    = 0;
   Double_t *arrSP    = 0;  //not needed for current fCaller
   Double_t *arrSM    = 0;  //not needed for current fCaller
   Int_t    *arrXP    = 0;  //not needed for current fCaller
   Int_t    *arrXM    = 0;  //not needed for current fCaller

   TTree  *calltree = 0;
   XPCall *call     = 0;
   Int_t   split    = 99;

// Calculate detection call
   Bool_t doGC = (Bool_t)(strcmp(fCaller->GetName(), "dabgcall") == 0);
   Int_t  ij, idx, x, y, start, end;
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;

   // Informing user
      TString name = datatree[k]->GetName();
      if (XManager::fgVerbose) {
         cout << "   Calculating present call for <" << name << ">..." << endl;
      }//if

   // Get tree info for datatree
      XDataTreeInfo *info = 0;
      info = (XDataTreeInfo*)datatree[k]->GetUserInfo()->FindObject(name);
      if (!info) {
         cerr << "Error: Could not get tree info for <" << name << ">." << endl;
         return errGeneral;
      }//if

   // Get chip parameters from scheme file (also alternative CDFs)
      if (strcmp(fSchemeName.Data(), "") == 0) {
         fSchemeName = datatree[k]->GetTitle();
      } else if (!fSchemeName.Contains(datatree[k]->GetTitle())) {
         return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[k]->GetTitle());
      }//if
      if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

      XGeneChip *chip = (XGeneChip*)fSchemes->FindObject(fSchemeName, kTRUE);
      if (!chip) {
         return fManager->HandleError(errGetScheme, fSchemeName);
      }//if
      Int_t numrows  = chip->GetNumRows();
      Int_t numcols  = chip->GetNumColumns();
      Int_t numctrls = chip->GetNumControls();
      Int_t numunits = chip->GetNumUnits();
      Int_t size     = numrows*numcols;

      Int_t numabsent  = 0;
      Int_t numarginal = 0;
      Int_t numpresent = 0;
      Double_t minpval = 1.0;
      Double_t maxpval = 0.0;

   // Get scheme tree for scheme
      XScheme *scheme = 0;
      TLeaf *scmleaf = 0;
      TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
      if (scmtree == 0) return errGetTree;
      scmtree->SetBranchAddress("ScmBranch", &scheme);

   // Get unit tree for scheme
      XGCUnit *unit = 0;
      TTree *idxtree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
      if (idxtree == 0) return errGetTree;
      idxtree->SetBranchAddress("IdxBranch", &unit);

   // Get maximum number of cells from tree info for unit tree
      XGenomeTreeInfo *idxinfo = 0;
      idxinfo = (XGenomeTreeInfo*)idxtree->GetUserInfo()->FindObject(idxtree->GetName());
      if (!idxinfo) {
         cerr << "Error: Could not get tree info for <" << idxtree->GetName() << ">." << endl;
         return errGeneral;
      }//if
      Int_t maxnumcells = (Int_t)idxinfo->GetValue("fMaxNCells");

   // Get exon level of annotation
      Int_t level = 0;
      if (fCallSelector->GetNumParameters() > 0) {
         level = (Int_t)(fCallSelector->GetParameters())[0];
      }//if

   // Init unit selector
      XUnitSelector *unitSelector = 0;

   // Initialize memory for unit arrays
      if (!(arrUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}
      if (!(mskUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 

   // Initialize memory for data arrays (uncomment when needed in fCaller)
      if (!(arrMask  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
      if (!(arrInten = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
//      if (!(arrStdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
//      if (!(arrNPix  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<size; i++) {
         arrMask[i]  = eINITMASK;  //init mask
         arrInten[i] = 0.0;
      }//for_i

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

   // Get mask from scheme tree and store in array 
      arrMask = this->FillMaskArray(chip, scmtree, scheme, level, size, arrMask);
      if (arrMask == 0) {err = errInitMemory; goto cleanup;}

   // Calculate units satisfying mask
      unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
      unitSelector->SetOption("genome");
      err = unitSelector->InitParameters(fCallSelector->GetNumParameters(),
                                         fCallSelector->GetParameters());
      if (err != errNoErr) goto cleanup;

   // Get mask from scheme tree and store in array 
      arrUnit = this->FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);

      err = unitSelector->Calculate(numunits, arrUnit, mskUnit);
      if (err != errNoErr) goto cleanup;

   // Calculate mask for detection call (set arrMask to 1 or 0)
      err = fCallSelector->Calculate(size, 0, 0, arrMask);
      if (err != errNoErr) goto cleanup;

   // Get data from datatree and store in arrays
      err = FillDataArrays(datatree[k], numrows, numcols, arrInten, arrStdev, arrNPix);
      if (err != errNoErr) goto cleanup;

   // For DABG only, fill arrMask with GC content (GC>=0 for PM)
   // for PM set GC >= 0 , i.e. GC = 0...kProbeLength,
   // for MM set GC < eINITMASK, i.e. GC = eINITMASK - (1...kProbeLength+1)
      if (doGC) {
      // Get probe tree for scheme
         XGCProbe *probe = 0;
         TTree *probetree = (TTree*)gDirectory->Get(chip->GetProbeTree()); 
         if (probetree == 0) {err = errGetTree; goto cleanup;}
         probetree->SetBranchAddress("PrbBranch", &probe);

      // Get GC content from probe tree and store in array
         for (Int_t i=0; i<size; i++) {
            probetree->GetEntry(i);

            x  = probe->GetX();
            y  = probe->GetY();
            ij = XY2Index(x, y, numcols);

            if (arrMask[ij] == 1) {
               arrMask[ij] = probe->GetNumberGC();
            } else if (arrMask[ij] == 0) {
               //need to use (numberGC + 1) to avoid setting arrMask=eINITMASK for GC=0!!
               arrMask[ij] = eINITMASK - (probe->GetNumberGC() + 1);
//no               arrMask[ij] = eINITMASK - probe->GetNumberGC();
            }//if
         }//for_i

         // get intensities and save in table sorted for GC-content
         if ((err = fCaller->Calculate(size, arrInten, arrMask))) goto cleanup;
      } else {
      // Need to get background data for MM values!
         if (!(arrBgrd  = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
//         if (!(arrBgdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}

         for (Int_t i=0; i<size; i++) arrBgrd[i] = 0.0;

         if (bgrdtree[k] != 0) {
           // get data from bgrdtree and store in arrays
//            err = FillBgrdArrays(bgrdtree[k], numrows, numcols, arrBgrd, arrBgdev);
            err = FillBgrdArrays(bgrdtree[k], numrows, numcols, arrBgrd, 0);
            if (err != errNoErr) goto cleanup;
         }//if
      }//if

   // Create new tree calltree
      if (!fFile->cd(fName)) return errGetDir;

      name = Path2Name(datatree[k]->GetName(),"/",".") + "." + fCaller->GetTitle();
      calltree = new TTree(name, fSchemeName);
      if (calltree == 0) return errCreateTree;
      call = new XPCall();
      calltree->Branch("CallBranch", "XPCall", &call, 64000, split);

   // Initialize memory for PM/MM arrays (uncomment when needed in fCaller)
      if (!(arrPM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrMM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
//      if (!(arrSP = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
//      if (!(arrSM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrXP = new (nothrow) Int_t[maxnumcells]))    {err = errInitMemory; goto cleanup;}
//      if (!(arrXM = new (nothrow) Int_t[maxnumcells]))    {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<maxnumcells; i++) {
         arrPM[i] = arrMM[i] = 0.0; 
         arrXP[i] = 0;
      }//for_i

   // Calculate detection call values
      start = 0;
      end   = 0;
      idx   = 0;
      for (Int_t id=0; id<numunits; id++) { 
         idxtree->GetEntry(id);

         Int_t unitID   = unit->GetUnitID();
         Int_t numcells = unit->GetNumCells();
         // skip masked unit entries
         if (mskUnit[id] <= 0) {
            start += numcells;
            end = start;
            continue;
         }//if

         Int_t p = 0;
         Int_t m = 0;
         end += numcells;
         for (Int_t j=start; j<end; j++) {
            scmtree->GetEntry(j);

            if ((scheme->GetUnitID()) != unitID) {
               cerr << "Error: unitID is not equal to: " << unitID << endl;
               err = errAbort;
               goto cleanup;
            }//if

            x  = scheme->GetX();
            y  = scheme->GetY();
            ij = XY2Index(x, y, numcols);

            if (doGC) {
               if (arrMask[ij] >= 0) {
                  if (p == 0) idx++;  //count number of units
                  arrPM[p] = arrInten[ij];
                  arrXP[p] = arrMask[ij];
                  p++;
               }//if
            } else {
               if (arrMask[ij] == 1) {
                  if (p == 0) idx++;  //count number of units
                  arrPM[p] = arrInten[ij];
//no                  arrSP[p] = arrStdev[ij];
//no                  arrXP[p] = arrNPix[ij];
                  p++;

               // set MM values to background data
                  arrMM[m] = arrBgrd[ij];
//no                  arrSM[m] = arrBgdev[ij];
//no                  arrXM[p] = arrNPix[ij];
                  m++;
               }//if
            }//if

            if (p > maxnumcells || m > maxnumcells) {
               cerr << "Error: unitID <" << unitID << "> exceeds maximum number of cells <"
                    << maxnumcells << ">. Buffer overflow!" << endl;
               err = errAbort;
               goto cleanup;
            }//if
         }//for_j
         start += numcells;

         // continue if arrays arrPM etc are not filled
         if (p == 0) continue;
         if (!doGC && (p != m)) {
            cout << "Warning: Skipping unitID <" << unitID
                 << "> with different numbers of PM and MM data." << endl;
//            continue;
            err = errAbort;
            goto cleanup;
         }//if

         if (XManager::fgVerbose && id%10000 == 0) {
            cout << "      <" << idx << "> of <" << numunits << "> calls processed...\r" << flush;
         }//if

         // calculate detection call
         Int_t    dummy    = 0;
         Double_t prescall = 0.0;
         Double_t pvalue   = 1.0;
         fCaller->InitTreeInfo(info);
         fCaller->InitArrays(p, arrPM, arrSP, arrXP, arrMM, arrSM, arrXM);
         if ((err = fCaller->Calculate(prescall, pvalue, dummy))) goto cleanup;

         // number of present/absent calls
         if      (prescall == 2.0) numpresent++;
         else if (prescall == 1.0) numarginal++;
         else if (prescall == 0.0) numabsent++;

         // minimal/maximal detection call p-value
         if      (pvalue < minpval) minpval = pvalue;
         else if (pvalue > maxpval) maxpval = pvalue;

         // fill call tree
         call->SetUnitID(unitID);
         call->SetCall((Short_t)prescall);
         call->SetPValue(pvalue);
         calltree->Fill();
      }//for_id

      if (XManager::fgVerbose) {
         cout << "      <" << idx << "> of <" << numunits << "> calls processed...Finished" << endl;
      }//if

      if (XManager::fgVerbose) {
         cout << "      detection call statistics: " << endl;
         cout << "         minimum detection p-value = " << minpval << endl;
         cout << "         maximum detection p-value = " << maxpval << endl;
         cout << "         P: <" << numpresent*100.0/idx << "> percent units present."  << endl;
         cout << "         M: <" << numarginal*100.0/idx << "> percent units marginal." << endl;
         cout << "         A: <" << numabsent*100.0/idx  << "> percent units absent."   << endl;
      }//if

   // Add tree info to tree
      AddCallTreeInfo(calltree, calltree->GetName(), fCaller->GetOption(),
                      numunits-numctrls, numabsent, numarginal, numpresent,
                      minpval, maxpval);

   // Write call tree to file 
      if ((err = WriteTree(calltree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(calltree->GetName(), "Call", 0, fCaller->GetNumParameters(),
                       fCaller->GetParameters());
      }//if

   cleanup:
      // delete arrays
      if (arrXM)    {delete [] arrXM;    arrXM    = 0;}
      if (arrXP)    {delete [] arrXP;    arrXP    = 0;}
      if (arrSM)    {delete [] arrSM;    arrSM    = 0;}
      if (arrSP)    {delete [] arrSP;    arrSP    = 0;}
      if (arrMM)    {delete [] arrMM;    arrMM    = 0;}
      if (arrPM)    {delete [] arrPM;    arrPM    = 0;}
      if (arrBgdev) {delete [] arrBgdev; arrBgdev = 0;}
      if (arrBgrd)  {delete [] arrBgrd;  arrBgrd  = 0;}
      if (arrNPix)  {delete [] arrNPix;  arrNPix  = 0;}
      if (arrStdev) {delete [] arrStdev; arrStdev = 0;}
      if (arrInten) {delete [] arrInten; arrInten = 0;}
      if (arrMask)  {delete [] arrMask;  arrMask  = 0;}
      if (mskUnit)  {delete [] mskUnit;  mskUnit  = 0;}
      if (arrUnit)  {delete [] arrUnit;  arrUnit  = 0;}

      if (err != errNoErr) break;
   }//for_k

   return err;
}//DetectCall

//______________________________________________________________________________
Int_t XGenomeProcesSet::DoExpress(Int_t numdata, TTree **datatree,
                        Int_t numbgrd, TTree **bgrdtree)
{
   // Condense data from  datatrees to expression values
   // Note: each datatree is processed individually and thus can have different
   //       number of entries (i.e. belong to different chip types)
   if(kCS) cout << "------XGenomeProcesSet::DoExpress------" << endl;

   Int_t err = errNoErr;
   Int_t split = 99;

// Init local arrays to store data from trees
   Int_t    *arrUnit  = 0;
   Int_t    *mskUnit  = 0;
   Int_t    *arrMask  = 0;
   Double_t *arrInten = 0;
   Double_t *arrStdev = 0;
   Int_t    *arrNPix  = 0;
   Double_t *arrBgrd  = 0;
   Double_t *arrBgdev = 0;
   Double_t *arrPM    = 0; 
   Double_t *arrMM    = 0;
   Double_t *arrSP    = 0; 
   Double_t *arrSM    = 0;
   Int_t    *arrXP    = 0; 
   Int_t    *arrXM    = 0;

// Calculate expression
   Int_t ij, idx, x, y, start, end;
   Int_t    arrlen = 0;
   Double_t mean   = 0;
   Double_t var    = 0;

   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;
      if (bgrdtree[k] == 0) {
         cerr << "Error: Could not get background tree for <"
              << fExpressor->GetName() << ">." << endl;
         return errGetTree;
      }//if

   // Informing user
      TString name = datatree[k]->GetName();
      if (XManager::fgVerbose) {
         cout << "      Summarizing <" << name << "> using <" << fExpressor->GetName()
              << ">..." << endl;
      }//if

   // Get tree info from datatree
      fDataFile->cd();

      XDataTreeInfo *info = 0;
      info = (XDataTreeInfo*)datatree[k]->GetUserInfo()->FindObject(name);
      if (!info) {
         cerr << "Error: Could not get tree info for <" << name << ">." << endl;
         return errGeneral;
      }//if

   // Get chip parameters from scheme file (also alternative CDFs)
      if (strcmp(fSchemeName.Data(), "") == 0) {
         fSchemeName = datatree[k]->GetTitle();
      } else if (!fSchemeName.Contains(datatree[k]->GetTitle())) {
         return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[k]->GetTitle());
      }//if
      if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

      XGeneChip *chip = (XGeneChip*)fSchemes->FindObject(fSchemeName, kTRUE);
      if (!chip) {
         return fManager->HandleError(errGetScheme, fSchemeName);
      }//if
      Int_t numrows  = chip->GetNumRows();
      Int_t numcols  = chip->GetNumColumns();
      Int_t numctrls = chip->GetNumControls();
      Int_t numunits = chip->GetNumUnits();
      Int_t size     = numrows*numcols;

   // Init min/max expression levels
      Double_t min = DBL_MAX;  //defined in float.h
      Double_t max = 0;

   // Get scheme tree for scheme
      XScheme *scheme = 0;
      TLeaf *scmleaf = 0;
      TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
      if (scmtree == 0) return errGetTree;
      scmtree->SetBranchAddress("ScmBranch", &scheme);
      scmleaf = scmtree->FindLeaf("fUnitID"); 

   // Get unit tree for scheme
      XGCUnit *unit = 0;
      TTree *idxtree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
      if (idxtree == 0) return errGetTree;
      idxtree->SetBranchAddress("IdxBranch", &unit);

   // Get maximum number of cells from tree info for unit tree
      XGenomeTreeInfo *idxinfo = 0;
      idxinfo = (XGenomeTreeInfo*)idxtree->GetUserInfo()->FindObject(idxtree->GetName());
      if (!idxinfo) {
         cerr << "Error: Could not get tree info for <" << idxtree->GetName() << ">." << endl;
         err = errGeneral;
         break;
      }//if
      Int_t maxnumcells = (Int_t)idxinfo->GetValue("fMaxNCells");

   // Create new tree exprtree
      if (!fFile->cd(fName)) return errGetDir;

      TString dataname = Path2Name(datatree[k]->GetName(),"/",".");
      TString exprname = dataname + "." + fExpressor->GetTitle();
      TTree  *exprtree = new TTree(exprname, fSchemeName);
      if (exprtree == 0) return errCreateTree;

      XGCExpression *expr = 0;
      expr = new XGCExpression();
      exprtree->Branch("ExprBranch", "XGCExpression", &expr, 64000, split);

   // Get exon level of annotation
      Int_t level = 0;
      if (fExprSelector->GetNumParameters() > 0) {
         level = (Int_t)(fExprSelector->GetParameters())[0];
      }//if

   // Init unit selector
      XUnitSelector *unitSelector = 0;

   // Initialize memory for unit arrays
      if (!(arrUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}
      if (!(mskUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 

   // Initialize memory for scheme arrays
      if (!(arrMask  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
      if (!(arrInten = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrStdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrNPix  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
      if (!(arrBgrd  = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrBgdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<size; i++) {
         arrMask[i]  = eINITMASK;  //init mask
         arrInten[i] = arrStdev[i] = 0.0; 
         arrNPix[i]  = 0; 
         arrBgrd[i]  = arrBgdev[i] = 0.0;
      }//for_i

   // Initialize memory for PM/MM arrays
      if (!(arrPM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrMM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrSP = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrSM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrXP = new (nothrow) Int_t[maxnumcells]))    {err = errInitMemory; goto cleanup;}
      if (!(arrXM = new (nothrow) Int_t[maxnumcells]))    {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<maxnumcells; i++) {
         arrPM[i] = arrMM[i] = 0.0; 
         arrSP[i] = arrSM[i] = 0.0; 
         arrXP[i] = arrXM[i] = 0; 
      }//for_i

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

   // Get mask from scheme tree and store in array 
      arrMask = this->FillMaskArray(chip, scmtree, scheme, level, size, arrMask);
      if (arrMask == 0) {err = errInitMemory; goto cleanup;}

   // Calculate units satisfying mask
      unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
      unitSelector->SetOption("genome");
      err = unitSelector->InitParameters(fExprSelector->GetNumParameters(),
                                         fExprSelector->GetParameters());
      if (err != errNoErr) goto cleanup;

   // Get mask from scheme tree and store in array 
      arrUnit = this->FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);

      err = unitSelector->Calculate(numunits, arrUnit, mskUnit);
      if (err != errNoErr) goto cleanup;

   // Calculate mask for expression
      err = fExprSelector->Calculate(size, 0, 0, arrMask);
      if (err != errNoErr) goto cleanup;

   // Get data from datatree and bgrdtree and store in arrays
      err = FillDataArrays(datatree[k], numrows, numcols, arrInten, arrStdev, arrNPix);
      if (err != errNoErr) goto cleanup;
      err = FillBgrdArrays(bgrdtree[k], numrows, numcols, arrBgrd, arrBgdev);
      if (err != errNoErr) goto cleanup;

   // Calculate expression values
      start = 0;
      end   = 0;
      idx   = 0;
      for (Int_t id=0; id<numunits; id++) { 
         idxtree->GetEntry(id);

         Int_t unitID   = unit->GetUnitID();
         Int_t numcells = unit->GetNumCells();
         // skip masked unit entries
         if (mskUnit[id] <= 0) {
            start += numcells;
            end = start;
            continue;
         }//if

         Int_t p = 0;
         Int_t m = 0;
         end += numcells;
         for (Int_t j=start; j<end; j++) {
            scmtree->GetEntry(j);

            if ((Int_t)(scmleaf->GetValue()) != unitID) {
               cerr << "Error: unitID is not equal to: " << unitID << endl;
               err = errAbort;
               goto cleanup;
            }//if

            x  = scheme->GetX();
            y  = scheme->GetY();
            ij = XY2Index(x, y, numcols);

            if (arrMask[ij] == 1) {
               if (p == 0) idx++;  //count number of units to be summarized
               arrPM[p] = arrInten[ij];
               arrSP[p] = arrStdev[ij];
               arrXP[p] = arrNPix[ij];
               p++;

               // set MM values to background data
               arrMM[m] = arrBgrd[ij];
               arrSM[m] = arrBgdev[ij];
               arrXM[m] = arrNPix[ij]; //use numpix from datatree
               m++;
            }//if
         }//for_j
         start += numcells;

         if (XManager::fgVerbose && id%10000 == 0) {
            cout << "      calculating expression for <" << idx << "> of <"
                 << numunits << "> units...\r" << flush;
         }//if

         // continue if arrays arrPM etc are not filled
         if (p == 0) continue;
         if (p != m) {
            cout << "Warning: Skipping unitID <" << unitID
                 << "> with different numbers of PM and MM data." << endl;
            continue;
//?            err = errAbort;
//?            goto cleanup;
         }//if

         // calculate mean expression level
         arrlen = 0;
         mean   = 0;
         var    = 0;
         fExpressor->InitTreeInfo(info);
         fExpressor->InitArrays(p, arrPM, arrSP, arrXP, arrMM, arrSM, arrXM);
         if ((err = fExpressor->CreateArray(p)))               goto cleanup;
         if ((err = fExpressor->Calculate(mean, var, arrlen))) goto cleanup;
         fExpressor->DeleteArray();

         // get minimal/maximal expression levels
         if (mean < min) min = mean;
         if (mean > max) max = mean;

         // fill expression tree
         expr->SetUnitID(unitID);
         expr->SetLevel(mean);
         expr->SetStdev(TMath::Sqrt(var));
         expr->SetNumPairs(arrlen);
         exprtree->Fill();
      }//for_id

      if (XManager::fgVerbose) {
         cout << "      calculating expression for <" << idx << "> of <" << numunits
              << "> units...Finished." << endl;
      }//if

      if (XManager::fgVerbose) {
         cout << "      expression statistics: " << endl;
         cout << "         minimal expression level is <" << min << ">." << endl;
         cout << "         maximal expression level is <" << max << ">." << endl;
      }//if

   // Add tree info to tree
      AddExprTreeInfo(exprtree, exprtree->GetName(), fExpressor->GetOption(),
                      numunits-numctrls, min, max);

   // Write expression tree to file 
      if ((err = WriteTree(exprtree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(exprtree->GetName(), "Expr", 0, fExpressor->GetNumParameters(),
                       fExpressor->GetParameters());
      }//if

   cleanup:
      // delete arrays
      if (arrXM)    {delete [] arrXM;    arrXM    = 0;}
      if (arrXP)    {delete [] arrXP;    arrXP    = 0;}
      if (arrSM)    {delete [] arrSM;    arrSM    = 0;}
      if (arrSP)    {delete [] arrSP;    arrSP    = 0;}
      if (arrMM)    {delete [] arrMM;    arrMM    = 0;}
      if (arrPM)    {delete [] arrPM;    arrPM    = 0;}
      if (arrBgdev) {delete [] arrBgdev; arrBgdev = 0;}
      if (arrBgrd)  {delete [] arrBgrd;  arrBgrd  = 0;}
      if (arrNPix)  {delete [] arrNPix;  arrNPix  = 0;}
      if (arrStdev) {delete [] arrStdev; arrStdev = 0;}
      if (arrInten) {delete [] arrInten; arrInten = 0;}
      if (arrMask)  {delete [] arrMask;  arrMask  = 0;}
      if (mskUnit)  {delete [] mskUnit;  mskUnit  = 0;}
      if (arrUnit)  {delete [] arrUnit;  arrUnit  = 0;}
      // Note: do not remove exprtree and expr from RAM, needed later!

      if (err != errNoErr) break;
   }//for_k

   return err;
}//DoExpress

//______________________________________________________________________________
Int_t XGenomeProcesSet::DoMedianPolish(Int_t numdata, TTree **datatree,
                        Int_t numbgrd, TTree **bgrdtree)
{
   // Compute expression values using mediapolish
   // Intensities from datatrees are stored in one large table in RAM for fast access
   // Note: all trees must have same number of entries (i.e. identical chip types)
   if(kCS) cout << "------XGenomeProcesSet::DoMedianPolish(table)------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "      Summarizing with medianpolish..." << endl;
   }//if

//TEST Benchmark
//gBenchmark->Reset(); 
//gBenchmark->Start("Bench_MedianPolish");

   Int_t ij, start, end, id;
   Int_t x   = 0;
   Int_t y   = 0;
   Int_t idx = 0;
   Int_t err = errNoErr;

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fExpressor, fExpressor->GetBgrdOption());

// Get chip parameters from scheme file (also alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XGeneChip *chip = (XGeneChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) {
      return fManager->HandleError(errGetScheme, fSchemeName);
   }//if
   Int_t numrows  = chip->GetNumRows();
   Int_t numcols  = chip->GetNumColumns();
   Int_t numunits = chip->GetNumUnits();

// Get scheme tree for scheme
   XScheme *scheme = 0;
   TLeaf *scmleaf = 0;
   TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
   if (scmtree == 0) return errGetTree;
   scmtree->SetBranchAddress("ScmBranch", &scheme);
   scmleaf = scmtree->FindLeaf("fUnitID");      

 // Get unit tree for scheme
   XGCUnit *unit = 0;
   TTree *idxtree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   if (idxtree == 0) return errGetTree;
   idxtree->SetBranchAddress("IdxBranch", &unit);

// Init size of arrays
   Int_t size    = numrows*numcols;
   Int_t numsels = 0;  //number of selected entries

// Check for equal number of data tree entries
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;

      if ((Int_t)(datatree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, datatree[k]->GetName(), str);
      }//if
   }//for_k

// Check for equal number of background tree entries
   if (numbgrd > 0) {
      for (Int_t k=0; k<numdata; k++) {
         if (bgrdtree[k] == 0) return errGetTree;

         if ((Int_t)(bgrdtree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, bgrdtree[k]->GetName(), str);
         }//if
      }//for_k
   } else if (doBg == kTRUE) {
      cout << "Warning: No background trees available for background subtraction."
           << endl;
      doBg = kFALSE;
   } else {
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Init branch addresses
   XBgCell *bgcell[numdata];
   XGCCell *gccell[numdata];
   for (Int_t k=0; k<numdata; k++) {
      bgcell[k] = 0;
      gccell[k] = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell[k]);
      if (numbgrd > 0) bgrdtree[k]->SetBranchAddress("BgrdBranch", &bgcell[k]);
   }//for_k

// Init expression trees
   Int_t  split = 99;
   TTree *exprtree[numdata];
   XGCExpression *expr[numdata];

// Init min/max expression levels
   Double_t min = DBL_MAX;  //defined in float.h
   Double_t max = 0;

// Get exon level of annotation
   Int_t level = 0;
   if (fExprSelector->GetNumParameters() > 0) {
      level = (Int_t)(fExprSelector->GetParameters())[0];
   }//if

// Init unit selector
   XUnitSelector *unitSelector = 0;

// Init 
   Int_t     *arrMask = 0;
   Int_t     *arrIndx = 0;
   Int_t     *arrUnit = 0;
   Int_t     *mskUnit = 0;
   Double_t  *colmed  = 0;
   Double_t  *results = 0;
   Double_t **table   = 0;

// Create local arrays
   if (!(arrMask = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrIndx = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrUnit = new (nothrow) Int_t[numunits]))   {err = errInitMemory; goto cleanup;}
   if (!(mskUnit = new (nothrow) Int_t[numunits]))   {err = errInitMemory; goto cleanup;}
   if (!(colmed  = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}
   if (!(results = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<size; i++)    {arrIndx[i] = 0; arrMask[i] = eINITMASK;}
   for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 
   for (Int_t i=0; i<numdata; i++)  colmed[i]  = results[i] = 0.0; 

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

// Get mask from scheme tree and store in array 
   arrMask = this->FillMaskArray(chip, scmtree, scheme, level, size, arrMask);
   if (arrMask == 0) {err = errInitMemory; goto cleanup;}

// Calculate units satisfying mask
   unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
   unitSelector->SetOption("genome");
   err = unitSelector->InitParameters(fExprSelector->GetNumParameters(),
                                      fExprSelector->GetParameters());
   if (err != errNoErr) goto cleanup;

// Get units and mask from unit tree and store in array 
   arrUnit = this->FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);

   err = unitSelector->Calculate(numunits, arrUnit, mskUnit);
   if (err != errNoErr) goto cleanup;

// Calculate mask for expression
   err = fExprSelector->Calculate(size, 0, 0, arrMask);
   if (err != errNoErr) goto cleanup;

// Get entry index from datatree
   idx = 0;
   for (Int_t i=0; i<size; i++) {
      datatree[0]->GetEntry(i);

      x  = gccell[0]->GetX();
      y  = gccell[0]->GetY();
      ij = XY2Index(x, y, numcols);

      if (arrMask[ij] == 1) {
         arrIndx[ij] = idx++;
      }//if
   }//for_i

// Get number of selected entries
   for (Int_t i=0; i<size; i++) {
      numsels = (arrMask[i] == 1) ? ++numsels : numsels;
   }//for_i

// Create table to store intensities of all datatrees
   if (!(table = new (nothrow) Double_t*[numdata])) {err = errInitMemory; goto cleanup;} 
   for (Int_t k=0; k<numdata; k++) {
      table[k] = 0;
      if (!(table[k] = new (nothrow) Double_t[numsels])) {err = errInitMemory; goto cleanup;} 
   }//for_k

//TEST Benchmark
//gBenchmark->Start("Bench_table");
// Get data from datatrees (and bgrdtrees) and store in table
   if ((numbgrd > 0) && (doBg == kTRUE)) {
      for (Int_t k=0; k<numdata; k++) {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);
            bgrdtree[k]->GetEntry(i);
//unsafe?? use arrInten[ij] and arrBgrd[ij] ??

            // subtract background from intensity
            if (arrMask[i] == 1) {
               table[k][idx++] = this->AdjustIntensity(gccell[k]->GetIntensity(),
                                                       bgcell[k]->GetBackground(),
                                                       bgcell[k]->GetStdev());
            }//if
         }//for_i
      }//for_k
   } else {
      for (Int_t k=0; k<numdata; k++) {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);

            if (arrMask[i] == 1) {
               table[k][idx++] = gccell[k]->GetIntensity();
//??               table[k][arrIndx[i]] = gccell[k]->GetIntensity();
            }//if
         }//for_i
      }//for_k
   }//if
//TEST Benchmark
//gBenchmark->Show("Bench_table");

// Create new trees exprtree
   if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

   for (Int_t k=0; k<numdata; k++) {
      TString dataname = Path2Name(datatree[k]->GetName(),"/",".");
      TString exprname = dataname + "." + fExpressor->GetTitle();
      exprtree[k] = new TTree(exprname, fSchemeName);
      if (exprtree[k] == 0) {err = errCreateTree; goto cleanup;}

      expr[k] = new XGCExpression();
      exprtree[k]->Branch("ExprBranch", "XGCExpression", &expr[k], 64000, split);
   }//for_k

//TEST Benchmark
//gBenchmark->Show("Bench_MedianPolish");
//gBenchmark->Reset();
//gBenchmark->Start("Bench_Loop");

// Calculate expression values
   start = 0;
   end   = 0;
   idx   = 0;
   for (id=0; id<numunits; id++) { 
      idxtree->GetEntry(id);

      Int_t unitID   = unit->GetUnitID();
      Int_t numcells = unit->GetNumCells();
      // skip masked unit entries
      if (mskUnit[id] <= 0) {
         start += numcells;
         end = start;
         continue;
      }//if

///////////////////////////
//TO DO
//Better above: arrPM = new Double_t[maxnumcells*numdata]; 
///////////////////////////
      // create array to store PM values for all probes with current unitID
      Int_t  numatoms = unit->GetNumAtoms();
      Double_t *arrPM = new Double_t[numatoms*numdata]; 

      // fill arrPM with PM values of current unitID
      Int_t p = 0;
      end += numcells;
      for (Int_t j=start; j<end; j++) {
         scmtree->GetEntry(j);

         if ((Int_t)(scmleaf->GetValue()) != unitID) {
            cerr << "Error: unitID is not equal to: " << unitID << endl;
            err = errAbort;
            goto cleanup;
         }//if

         x  = scheme->GetX();
         y  = scheme->GetY();
         ij = XY2Index(x, y, numcols);

         if (arrMask[ij] == 1) {
            if (p == 0) idx++;  //count number of units to be summarized

            for (Int_t k=0; k<numdata; k++) {
               arrPM[p] = table[k][arrIndx[ij]];
               p++;
            }//for_k
         }//if
      }//for_j
      start += numcells;

      if (XManager::fgVerbose && id%10000 == 0) {
         cout << "      calculating expression for <" << idx << "> of <" << numunits
              << "> units...\r" << flush;
      }//if

      // fill arrPM or continue if it is not filled
      if ((err = fExpressor->SetArray(p, arrPM)) != errNoErr) {
         delete [] arrPM;
         continue;
      }//if

      // calculate median polish for PMs of current unitID
      if ((err = fExpressor->Calculate(numdata, colmed, results, 0))) break;

//////////////////
// TO DO: return fResiduals for residual-plot!!!! (like affyPLM)
//      residuals = fExpressor->GetResiduals();
// store as residuals tree??
//////////////////

      // fill expression trees
      for (Int_t k=0; k<numdata; k++) {
         // get minimal/maximal expression levels
         if (results[k] < min) min = results[k];
         if (results[k] > max) max = results[k];

         expr[k]->SetUnitID(unitID);
         expr[k]->SetLevel(results[k]);
         expr[k]->SetStdev(TMath::Abs(colmed[k]));
         expr[k]->SetNumPairs((Int_t)(p / numdata));
         exprtree[k]->Fill();
      }//for_k

      delete [] arrPM;
   }//for_id

   if (XManager::fgVerbose) {
      cout << "      calculating expression for <" << idx << "> of <"
           << numunits << "> units...Finished." << endl;
   }//if

   if (XManager::fgVerbose) {
      cout << "      expression statistics: " << endl;
      cout << "         minimal expression level is <" << min << ">." << endl;
      cout << "         maximal expression level is <" << max << ">." << endl;
   }//if
//TEST Benchmark
//gBenchmark->Show("Bench_Loop");

// Write expression trees to file 
   for (Int_t k=0; k<numdata; k++) {
   // Add tree info to tree
      AddExprTreeInfo(exprtree[k], exprtree[k]->GetName(), fExpressor->GetOption(),
                      idx, min, max);

      if ((err = WriteTree(exprtree[k], TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(exprtree[k]->GetName(), "Expr", 0, fExpressor->GetNumParameters(),
                       fExpressor->GetParameters());
      } else {
         break;
      }//if
   }//for_k

// Cleanup
cleanup:
   SafeDelete(unitSelector);

   // delete table
   for (Int_t k=0; k<numdata; k++) {
      if (table[k]) {delete [] table[k]; table[k] = 0;}
   }//for_k
   if (table) delete [] table;

   // delete arrays
   if (results) {delete [] results; results = 0;}
   if (colmed)  {delete [] colmed;  colmed  = 0;}
   if (mskUnit) {delete [] mskUnit; mskUnit = 0;}
   if (arrUnit) {delete [] arrUnit; arrUnit = 0;}
   if (arrIndx) {delete [] arrIndx; arrIndx = 0;}
   if (arrMask) {delete [] arrMask; arrMask = 0;}

   return err;
}//DoMedianPolish

//______________________________________________________________________________
Int_t XGenomeProcesSet::DoMedianPolish(Int_t numdata, TTree **datatree,
                        Int_t numbgrd, TTree **bgrdtree, TFile *file)
{
   // Compute expression values using mediapolish
   // In order to reduce memory consumption intensities from datatrees are
   // stored in the order of the corresponding schemetree entries in temporary 
   // trees in a temporary root file
   // Note: all trees must have same number of entries (i.e. identical chip types)
   if(kCS) cout << "------XGenomeProcesSet::DoMedianPolish(file)------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "      Summarizing with medianpolish (using temporary file)..." << endl;
   }//if

//TEST Benchmark
//gBenchmark->Reset(); 
//gBenchmark->Start("Bench_MedianPolish");

   Int_t ij, start, end, id, entry;
   Int_t x   = 0;
   Int_t y   = 0;
   Int_t idx = 0;
   Int_t err = errNoErr;

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fExpressor, fExpressor->GetBgrdOption());

// Get chip parameters from scheme file (also alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XGeneChip *chip = (XGeneChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) {
      return fManager->HandleError(errGetScheme, fSchemeName);
   }//if
   Int_t numrows  = chip->GetNumRows();
   Int_t numcols  = chip->GetNumColumns();
   Int_t numunits = chip->GetNumUnits();

// Get scheme tree for scheme
   XScheme *scheme = 0;
   TLeaf *scmleaf = 0;
   TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
   if (scmtree == 0) return errGetTree;
   scmtree->SetBranchAddress("ScmBranch", &scheme);
   scmleaf = scmtree->FindLeaf("fUnitID");      

 // Get unit tree for scheme
   XGCUnit *unit = 0;
   TTree *idxtree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   if (idxtree == 0) return errGetTree;
   idxtree->SetBranchAddress("IdxBranch", &unit);

// Init size of arrays
   Int_t size    = numrows*numcols;
   Int_t numsels = 0;  //number of selected entries

// Check for equal number of data tree entries
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;

      if ((Int_t)(datatree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, datatree[k]->GetName(), str);
      }//if
   }//for_k

// Check for equal number of background tree entries
   if (numbgrd > 0) {
      for (Int_t k=0; k<numdata; k++) {
         if (bgrdtree[k] == 0) return errGetTree;

         if ((Int_t)(bgrdtree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, bgrdtree[k]->GetName(), str);
         }//if
      }//for_k
   } else if (doBg == kTRUE) {
      cout << "Warning: No background trees available for background subtraction."
           << endl;
      doBg = kFALSE;
   } else {
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Init branch addresses
   XBgCell *bgcell[numdata];
   XGCCell *gccell[numdata];
   for (Int_t k=0; k<numdata; k++) {
      bgcell[k] = 0;
      gccell[k] = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell[k]);
      if (numbgrd > 0) bgrdtree[k]->SetBranchAddress("BgrdBranch", &bgcell[k]);
   }//for_k

// Init temporary trees and expression trees
   Int_t    split = 99;
   Double_t sort  = 0.0;
   TTree *tmptree[numdata];
   TTree *exprtree[numdata];
   XGCExpression *expr[numdata];

// Init min/max expression levels
   Double_t min = DBL_MAX;  //defined in float.h
   Double_t max = 0;

// Get exon level of annotation
   Int_t level = 0;
   if (fExprSelector->GetNumParameters() > 0) {
      level = (Int_t)(fExprSelector->GetParameters())[0];
   }//if

// Init unit selector
   XUnitSelector *unitSelector = 0;

// Create local arrays
   Int_t    *arrMask = 0;
   Int_t    *arrIndx = 0;
   Int_t    *arrUnit = 0;
   Int_t    *mskUnit = 0;
   Double_t *colmed  = 0;
   Double_t *results = 0;
   Double_t *arrData = 0;

// Create local arrays
   if (!(arrMask = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrIndx = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrUnit = new (nothrow) Int_t[numunits]))   {err = errInitMemory; goto cleanup;}
   if (!(mskUnit = new (nothrow) Int_t[numunits]))   {err = errInitMemory; goto cleanup;}
   if (!(colmed  = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}
   if (!(results = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<size; i++)    {arrIndx[i] = 0; arrMask[i] = eINITMASK;}
   for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 
   for (Int_t i=0; i<numdata; i++)  colmed[i]  = results[i] = 0.0; 

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

// Get mask from scheme tree and store in array 
   arrMask = this->FillMaskArray(chip, scmtree, scheme, level, size, arrMask);
   if (arrMask == 0) {err = errInitMemory; goto cleanup;}

// Calculate units satisfying mask
   unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
   unitSelector->SetOption("genome");
   err = unitSelector->InitParameters(fExprSelector->GetNumParameters(),
                                      fExprSelector->GetParameters());
   if (err != errNoErr) goto cleanup;

// Get mask from scheme tree and store in array 
   arrUnit = this->FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);

   err = unitSelector->Calculate(numunits, arrUnit, mskUnit);
   if (err != errNoErr) goto cleanup;

// Calculate mask for expression
   err = fExprSelector->Calculate(size, 0, 0, arrMask);
   if (err != errNoErr) goto cleanup;

// Get entry index from datatree
   idx = 0;
   for (Int_t i=0; i<size; i++) {
      datatree[0]->GetEntry(i);

      x  = gccell[0]->GetX();
      y  = gccell[0]->GetY();
      ij = XY2Index(x, y, numcols);

      if (arrMask[ij] == 1) {
         arrIndx[ij] = idx++;
      }//if
   }//for_i

// Get number of selected entries
   for (Int_t i=0; i<size; i++) {
      numsels = (arrMask[i] == 1) ? ++numsels : numsels;
   }//for_i

// Create array to store selected intensities
   if (!(arrData = new Double_t[numsels])) {err = errInitMemory; goto cleanup;}
   for (Int_t i=0; i<numsels; i++) arrData[i] = 0.0;

//TEST Benchmark
//gBenchmark->Start("Bench_tmptree");
// Change directory to temporary file
   if (!file->cd()) {err = errGetDir; goto cleanup;}

// Get data from datatrees (and bgrdtrees) and store in temporary file
   for (Int_t k=0; k<numdata; k++) {
      // create temporary tree
      tmptree[k] = new TTree(datatree[k]->GetName(), "temporary tree");
      if (tmptree[k] == 0) {err = errCreateTree; goto cleanup;}
      tmptree[k]->Branch("sortBr", &sort, "sort/D");

      // fill array with (background corrected) intensities
      if ((numbgrd > 0) && (doBg == kTRUE)) {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);
            bgrdtree[k]->GetEntry(i);

            // subtract background from intensity
            if (arrMask[i] == 1) {
               arrData[idx++] = this->AdjustIntensity(gccell[k]->GetIntensity(),
                                                      bgcell[k]->GetBackground(),
                                                      bgcell[k]->GetStdev());
            }//if
         }//for_i
      } else {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);

            if (arrMask[i] == 1) {
               arrData[idx++] = gccell[k]->GetIntensity();
            }//if
         }//for_i
      }//if

      // fill tmptree with array in the order of scheme tree entries for (x,y)
      for (Int_t i=0; i<size; i++) {
         scmtree->GetEntry(i);

         x  = scheme->GetX();
         y  = scheme->GetY();
         ij = XY2Index(x, y, numcols);

         if (arrMask[ij] == 1) {
            sort = arrData[arrIndx[ij]];
            tmptree[k]->Fill();
         }//if
      }//for_i

      // write tmptree to temporary file
      tmptree[k]->Write();
//TEST??      tmptree[k]->Write(TObject::kOverwrite);
   }//for_k
//TEST Benchmark
//gBenchmark->Show("Bench_tmptree");

// Change directory to current directory for treeset in main file
   if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

// Create new trees exprtree
   for (Int_t k=0; k<numdata; k++) {
      TString dataname = Path2Name(datatree[k]->GetName(),"/",".");
      TString exprname = dataname + "." + fExpressor->GetTitle();
      exprtree[k] = new TTree(exprname, fSchemeName);
      if (exprtree[k] == 0) {err = errCreateTree; goto cleanup;}

      expr[k] = new XGCExpression();
      exprtree[k]->Branch("ExprBranch", "XGCExpression", &expr[k], 64000, split);
   }//for_k

//TEST Benchmark
//gBenchmark->Show("Bench_MedianPolish");
//gBenchmark->Reset();
//gBenchmark->Start("Bench_Loop");

// Calculate expression values
   start = 0;
   end   = 0;
   idx   = 0;
   entry = 0;
   for (id=0; id<numunits; id++) { 
      idxtree->GetEntry(id);

      Int_t unitID   = unit->GetUnitID();
      Int_t numcells = unit->GetNumCells();
      // skip masked unit entries
      if (mskUnit[id] <= 0) {
         start += numcells;
         end = start;
         continue;
      }//if

//Better above: arrPM = new Double_t[maxnumcells*numdata]; 
      // create array to store PM values for all probes with current unitID
      Int_t numatoms = unit->GetNumAtoms();
      Double_t *arrPM = new Double_t[numatoms*numdata]; 

      // fill arrPM with PM values of current unitID
      Int_t p = 0;
      end += numcells;
      for (Int_t j=start; j<end; j++) {
         scmtree->GetEntry(j);

         if ((Int_t)(scmleaf->GetValue()) != unitID) {
            cerr << "Error: unitID is not equal to: " << unitID << endl;
            err = errAbort;
            goto cleanup;
         }//if

         x  = scheme->GetX();
         y  = scheme->GetY();
         ij = XY2Index(x, y, numcols);

         if (arrMask[ij] == 1) {
            if (p == 0) idx++;  //count number of units to be summarized

            for (Int_t k=0; k<numdata; k++) {
               tmptree[k]->GetEntry(entry);
               arrPM[p] = sort;
               p++;
            }//for_k

            entry++;
         }//if
      }//for_j
      start += numcells;

      if (XManager::fgVerbose && id%100000 == 0) {
         cout << "      calculating expression for <" << idx << "> of <"
              << numunits << "> units...\r" << flush;
      }//if

      // fill arrPM or continue if it is not filled
      if ((err = fExpressor->SetArray(p, arrPM)) != errNoErr) {
         delete [] arrPM;
         continue;
      }//if

      // calculate median polish for PMs of current unitID
      if ((err = fExpressor->Calculate(numdata, colmed, results, 0))) break;

//////////////////
// TO DO: return fResiduals for residual-plot!!!! (like affyPLM)
//      residuals = fExpressor->GetResiduals();
// store as residuals tree??
//////////////////

      // fill expression trees
      for (Int_t k=0; k<numdata; k++) {
         // get minimal/maximal expression levels
         if (results[k] < min) min = results[k];
         if (results[k] > max) max = results[k];

         expr[k]->SetUnitID(unitID);
         expr[k]->SetLevel(results[k]);
         expr[k]->SetStdev(TMath::Abs(colmed[k]));
         expr[k]->SetNumPairs((Int_t)(p / numdata));
         exprtree[k]->Fill();
      }//for_k

      delete [] arrPM;
   }//for_id
   if (XManager::fgVerbose) {
      cout << "      calculating expression for <" << idx << "> of <"
           << numunits << "> units...Finished." << endl;
   }//if

   if (XManager::fgVerbose) {
      cout << "      expression statistics: " << endl;
      cout << "         minimal expression level is <" << min << ">." << endl;
      cout << "         maximal expression level is <" << max << ">." << endl;
   }//if
//TEST Benchmark
//gBenchmark->Show("Bench_Loop");

// Write expression trees to file 
   for (Int_t k=0; k<numdata; k++) {
   // Add tree info to tree
      AddExprTreeInfo(exprtree[k], exprtree[k]->GetName(), fExpressor->GetOption(),
                      idx, min, max);

      if ((err = WriteTree(exprtree[k], TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(exprtree[k]->GetName(), "Expr", 0, fExpressor->GetNumParameters(),
                       fExpressor->GetParameters());
      } else {
         break;
      }//if
   }//for_k

// Cleanup
cleanup:
   SafeDelete(unitSelector);

   // delete temporary trees
   for (Int_t k=0; k<numdata; k++) {
      tmptree[k]->Delete(""); tmptree[k] = 0;
   }//for_k

   // delete arrays
   if (arrData) {delete [] arrData; arrData = 0;}
   if (results) {delete [] results; results = 0;}
   if (colmed)  {delete [] colmed;  colmed  = 0;}
   if (mskUnit) {delete [] mskUnit; mskUnit = 0;}
   if (arrUnit) {delete [] arrUnit; arrUnit = 0;}
   if (arrIndx) {delete [] arrIndx; arrIndx = 0;}
   if (arrMask) {delete [] arrMask; arrMask = 0;}

   return err;
}//DoMedianPolish

//______________________________________________________________________________
Int_t XGenomeProcesSet::ExportExprTrees(Int_t n, TString *names, const char *varlist,
                        ofstream &output, const char *sep)
{
   // Export variables from varlist for expression tree(s) to file output
   if(kCS) cout << "------XGenomeProcesSet::ExportExprTrees------" << endl;

   Int_t err = errNoErr;

// Decompose varlist
   Bool_t hasUnit   = kFALSE;
   Bool_t hasName   = kFALSE;
   Bool_t hasSymbol = kFALSE;
   Bool_t hasCyto   = kFALSE;
   Bool_t hasAnnot  = kFALSE;
   Bool_t hasLevel  = kFALSE;
   Bool_t hasStdev  = kFALSE;
   Bool_t hasNAtoms = kFALSE;
////////////////
//TO DO
//fTranscriptID
//fAccession
//fEntrezID
//fChromosome
//fStart
//fStop
//fStrand
////////////////

   if (strcmp(varlist,"*")  == 0) {
      hasUnit   = kTRUE;
      hasName   = kTRUE;
      hasSymbol = kTRUE;
      hasCyto   = kTRUE;
      hasLevel  = kTRUE;
      hasStdev  = kTRUE;
      hasNAtoms = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName") == 0) {hasUnit   = kTRUE;}
         if (strcmp(name,"fName")     == 0) {hasName   = kTRUE;}
         if (strcmp(name,"fSymbol")   == 0) {hasSymbol = kTRUE;}
         if (strcmp(name,"fCytoBand") == 0) {hasCyto   = kTRUE;}
         if (strcmp(name,"fLevel")    == 0) {hasLevel  = kTRUE;}
         if (strcmp(name,"fStdev")    == 0) {hasStdev  = kTRUE;}
         if (strcmp(name,"fNPairs")   == 0) {hasNAtoms = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if
   hasAnnot = (hasName || hasSymbol || hasCyto);

// Get trees
   TTree       *tree[n];
   XExpression *expr[n];

   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         expr[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("ExprBranch", &expr[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         expr[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("ExprBranch", &expr[k]);
      }//for_k
   }//if

// Get scheme name (also for alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = tree[0]->GetTitle();
   } else if (!fSchemeName.Contains(tree[0]->GetTitle())) {
      cerr << "Error: Scheme <" << fSchemeName << "> is not derived from <"
           << tree[0]->GetTitle() << ">." << endl;
      return errAbort;
   }//if

// Get chip from scheme file
   if (fSchemeFile == 0) return errGetScheme;
   fSchemeFile->cd();
   XFolder *schemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!schemes) {
      return fManager->HandleError(errMissingContent, "Scheme", kContent);
   }//if

   XGeneChip *chip = (XGeneChip*)schemes->FindObject(fSchemeName, kTRUE);
   if (chip == 0) return errAbort;

   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

// Get unit tree for scheme
   XUnit *unit     = 0;
   TTree *unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   if (unittree == 0) return errGetTree;
   unittree->SetBranchAddress("IdxBranch", &unit);

// Get annotation tree for scheme
   XTransAnnotation *annot = 0;
   TTree *anntree = (TTree*)gDirectory->Get(chip->GetAnnotTree()); 
   if (anntree == 0) return errGetTree;
   anntree->SetBranchAddress("AnnBranch", &annot);

   Int_t numunits = (Int_t)(unittree->GetEntries());
   Int_t numannot = (Int_t)(anntree->GetEntries());

// Create hash table to store unit names from anntree
   THashTable *htable = 0;
   if (!(htable = new THashTable(2*numannot))) return errInitMemory;

   if (XManager::fgVerbose) {
      cout << "Reading entries from <" << anntree->GetName() << "> ...";
   }//if
   XIdxString *idxstr = 0;
   for (Int_t i=0; i<numannot; i++) {
      anntree->GetEntry(i);

      idxstr = new XIdxString(i, annot->GetTranscriptID());
      htable->Add(idxstr);
   }//for_i
   if (XManager::fgVerbose) {
      cout << "Finished" << endl;
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit)      output << sep << "UNIT_NAME";
   if (hasAnnot) {
      if (hasName)   output << sep << "GENE_NAME";
      if (hasSymbol) output << sep << "GENE_SYMBOL";
      if (hasCyto)   output << sep << "CYTOBAND";
   }//if
   if (n == 1) {
      if (hasLevel)  output << sep << "LEVEL";
      if (hasStdev)  output << sep << "STDEV";
      if (hasNAtoms) output << sep << "NUMBER_ATOMS";
   } else {
      for (Int_t k=0; k<n; k++) {
         if (hasLevel)  output << sep << (names[k] + "_LEVEL");
         if (hasStdev)  output << sep << (names[k] + "_STDEV");
         if (hasNAtoms) output << sep << (names[k] + "_NUMBER_ATOMS");
      }//for_k
   }//if
   output << endl;

// Loop over tree entries and trees
   Int_t idx = 0;
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);

         // export annotation
         if (k == 0) {
            Int_t unitID = expr[k]->GetUnitID();
            output << unitID;

            unittree->GetEntry(idx++);
            while (unitID != unit->GetUnitID()) {
               if (idx == numunits) {
                 cerr << "Error: UnitID <" << unitID << "> not found." << endl;
                 err = errAbort; goto cleanup;
               }//if

               unittree->GetEntry(idx++);
            }//while

            if (hasUnit) {
               output << sep << unit->GetUnitName();
            }//if

            if (hasAnnot) {
               idxstr = (XIdxString*)(htable->FindObject(unit->GetUnitName()));
               if (idxstr) {
                  anntree->GetEntry(idxstr->GetIndex());
                  if (hasName)   output << sep << annot->GetName();
                  if (hasSymbol) output << sep << annot->GetSymbol();
                  if (hasCyto)   output << sep << annot->GetCytoBand();
               } else {
                  if (hasName)   output << sep << "NA";
                  if (hasSymbol) output << sep << "NA";
                  if (hasCyto)   output << sep << "NA";
               }//if
            }//if
         }//if

         // export data
         if (hasLevel)  output << sep << expr[k]->GetLevel();
         if (hasStdev)  output << sep << ((XGCExpression*)expr[k])->GetStdev();
         if (hasNAtoms) output << sep << ((XGCExpression*)expr[k])->GetNumPairs();
      }//for_j
      output << endl;
   }//for_i

//Cleanup
cleanup:
   // remove trees from RAM
   if (anntree)  {anntree->Delete("");  anntree  = 0;}
   if (unittree) {unittree->Delete(""); unittree = 0;}
   if (htable)   {htable->Delete(); delete htable; htable = 0;}
   SafeDelete(schemes);

   return err;
}//ExportExprTrees

//______________________________________________________________________________
Int_t XGenomeProcesSet::ExportCallTrees(Int_t n, TString *names, const char *varlist,
                        ofstream &output, const char *sep)
{
   // Export variables from varlist for call tree(s) to file output
   if(kCS) cout << "------XGenomeProcesSet::ExportCallTrees------" << endl;

   Int_t err = errNoErr;

// Decompose varlist
   Bool_t hasUnit   = kFALSE;
   Bool_t hasName   = kFALSE;
   Bool_t hasSymbol = kFALSE;
   Bool_t hasCyto   = kFALSE;
   Bool_t hasAnnot  = kFALSE;
   Bool_t hasCall   = kFALSE;
   Bool_t hasPVal   = kFALSE;
////////////////
//TO DO
//fTranscriptID
//fAccession
//fEntrezID
//fChromosome
//fStart
//fStop
//fStrand
////////////////

   if (strcmp(varlist,"*")  == 0) {
      hasUnit   = kTRUE;
      hasName   = kTRUE;
      hasSymbol = kTRUE;
      hasCyto   = kTRUE;
      hasCall   = kTRUE;
      hasPVal   = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName") == 0) {hasUnit   = kTRUE;}
         if (strcmp(name,"fName")     == 0) {hasName   = kTRUE;}
         if (strcmp(name,"fSymbol")   == 0) {hasSymbol = kTRUE;}
         if (strcmp(name,"fCytoBand") == 0) {hasCyto   = kTRUE;}
         if (strcmp(name,"fCall")     == 0) {hasCall   = kTRUE;}
         if (strcmp(name,"fPValue")   == 0) {hasPVal   = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if
   hasAnnot = (hasName || hasSymbol || hasCyto);

// Get trees
   TTree  *tree[n];
   XPCall *call[n];

   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         call[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("CallBranch", &call[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         call[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("CallBranch", &call[k]);
      }//for_k
   }//if

// Get scheme name (also for alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = tree[0]->GetTitle();
   } else if (!fSchemeName.Contains(tree[0]->GetTitle())) {
      cerr << "Error: Scheme <" << fSchemeName << "> is not derived from <"
           << tree[0]->GetTitle() << ">." << endl;
      return errAbort;
   }//if

// Get chip from scheme file
   if (fSchemeFile == 0) return errGetScheme;
   fSchemeFile->cd();
   XFolder *schemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!schemes) {
      return fManager->HandleError(errMissingContent, "Scheme", kContent);
   }//if

   XGeneChip *chip = (XGeneChip*)schemes->FindObject(fSchemeName, kTRUE);
   if (chip == 0) return errAbort;

   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

// Get unit tree for scheme
   XUnit *unit     = 0;
   TTree *unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   if (unittree == 0) return errGetTree;
   unittree->SetBranchAddress("IdxBranch", &unit);

// Get annotation tree for scheme
   XTransAnnotation *annot = 0;
   TTree *anntree = (TTree*)gDirectory->Get(chip->GetAnnotTree()); 
   if (anntree == 0) return errGetTree;
   anntree->SetBranchAddress("AnnBranch", &annot);

   Int_t numunits = (Int_t)(unittree->GetEntries());
   Int_t numannot = (Int_t)(anntree->GetEntries());

// Create hash table to store unit names from anntree
   THashTable *htable = 0;
   if (!(htable = new THashTable(2*numannot))) return errInitMemory;

   if (XManager::fgVerbose) {
      cout << "Reading entries from <" << anntree->GetName() << "> ...";
   }//if
   XIdxString *idxstr = 0;
   for (Int_t i=0; i<numannot; i++) {
      anntree->GetEntry(i);

      idxstr = new XIdxString(i, annot->GetTranscriptID());
      htable->Add(idxstr);
   }//for_i
   if (XManager::fgVerbose) {
      cout << "Finished" << endl;
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit)      output << sep << "UNIT_NAME";
   if (hasAnnot) {
      if (hasName)   output << sep << "GENE_NAME";
      if (hasSymbol) output << sep << "GENE_SYMBOL";
      if (hasCyto)   output << sep << "CYTOBAND";
   }//if
   if (n == 1) {
      if (hasCall)   output << sep << "CALL";
      if (hasPVal)   output << sep << "PVALUE";
   } else {
      for (Int_t k=0; k<n; k++) {
         if (hasCall)  output << sep << (names[k] + "_CALL");
         if (hasPVal)  output << sep << (names[k] + "_PVALUE");
      }//for_k
   }//if
   output << endl;

// Loop over tree entries and trees
   Int_t idx = 0;
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);

         // export annotation
         if (k == 0) {
            Int_t unitID = call[k]->GetUnitID();
            output << unitID;

            unittree->GetEntry(idx++);
            while (unitID != unit->GetUnitID()) {
               if (idx == numunits) {
                 cerr << "Error: UnitID <" << unitID << "> not found." << endl;
                 err = errAbort; goto cleanup;
               }//if

               unittree->GetEntry(idx++);
            }//while

            if (hasUnit) {
               output << sep << unit->GetUnitName();
            }//if

            if (hasAnnot) {
               idxstr = (XIdxString*)(htable->FindObject(unit->GetUnitName()));
               if (idxstr) {
                  anntree->GetEntry(idxstr->GetIndex());
                  if (hasName)   output << sep << annot->GetName();
                  if (hasSymbol) output << sep << annot->GetSymbol();
                  if (hasCyto)   output << sep << annot->GetCytoBand();
               } else {
                  if (hasName)   output << sep << "NA";
                  if (hasSymbol) output << sep << "NA";
                  if (hasCyto)   output << sep << "NA";
               }//if
            }//if
         }//if

         // export data
         if (hasCall) {
            Int_t cl = call[k]->GetCall();
            char *ch = "NA";
            if      (cl == 2) ch = "P";
            else if (cl == 0) ch = "A";
            else if (cl == 1) ch = "M";
            output << sep << ch;
         }//if

         if (hasPVal) {
            output << sep << call[k]->GetPValue();
         }//if
      }//for_j
      output << endl;
   }//for_i

//Cleanup
cleanup:
   // remove trees from RAM
   if (anntree)  {anntree->Delete("");  anntree  = 0;}
   if (unittree) {unittree->Delete(""); unittree = 0;}
   if (htable)   {htable->Delete(); delete htable; htable = 0;}
   SafeDelete(schemes);

   return err;
}//ExportCallTrees

//______________________________________________________________________________
Int_t *XGenomeProcesSet::FillMaskArray(XDNAChip *chip, TTree *scmtree, XScheme *scheme,
                         Int_t level, Int_t n, Int_t *msk)
{
   // "level" contains currently used exon level(s) as combination of ELevel
   // Problem: scheme tree for whole genome arrays can contain oligos at same 
   // (x,y) position belonging to different genes and/or affx->controls
   // Solution: extract current exon levels and store first in XBitSet, then
   // fill array "msk" with resulting level
   if(kCS) cout << "------XGenomeProcesSet::FillMaskArray------" << endl;

   Int_t x, y, ij, numcols;
   Int_t mask = 0;
   Int_t res  = 0;

// Set bit mask to level
   XBitSet bitmsk;
   bitmsk.SetBit(level);

// Create array of XBitSet to save results of adding msk values satisfying level
   XBitSet *bitMask = 0;
   if (!(bitMask = new XBitSet[n])) return 0;

// Fill bitMask with bit positions satisfying level
   numcols = chip->GetNumColumns();
   for (Int_t i=0; i<scmtree->GetEntries(); i++) {
      scmtree->GetEntry(i);

      x  = scheme->GetX();
      y  = scheme->GetY();
      ij = XY2Index(x, y, numcols);

      mask = scheme->GetMask();
      mask = (mask >= 0) ? mask : (abs(mask) << 15); //convert negative mask to int
      if (bitmsk.TestBit(mask) == kTRUE) (bitMask[ij]).SetBit(mask);
   }//for_i

// Fill array msk with results from bitMask 
   for (Int_t i=0; i<scmtree->GetEntries(); i++) {
      scmtree->GetEntry(i);

      x  = scheme->GetX();
      y  = scheme->GetY();
      ij = XY2Index(x, y, numcols);

      mask = scheme->GetMask();
      mask = (mask >= 0) ? mask : (abs(mask) << 15); //convert negative mask to int
      res  = (bitMask[ij]).TestBits(XBitSet::kBitMask);
      msk[ij] = (bitmsk.TestBit(res)) ? res : mask;
   }//for_i

   if (bitMask) {delete [] bitMask; bitMask = 0;}

   return msk;
}//FillMaskArray

//______________________________________________________________________________
Int_t *XGenomeProcesSet::FillUnitArray(TTree *idxtree, XGCUnit *unit,
                         Int_t n, Int_t *arr, Int_t *msk)
{
   // Fill array arr with unitID and array msk with unittype
   // Note: Convert unittype < 0 to higher bytes of integer for use with XBitSet
   if(kCS) cout << "------XGenomeProcesSet::FillUnitArray------" << endl;

   Int_t mask = 0;
   for (Int_t i=0; i<n; i++) {
      idxtree->GetEntry(i);

      arr[i] = unit->GetUnitID();
      mask   = unit->GetUnitType();
      msk[i] = (mask >= 0) ? mask : (abs(mask) << 15); //convert negative mask to int
   }//for_i

   return arr;
}//FillUnitArray


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonProcesSet                                                       //
//                                                                      //
// Class for ExonChip oligonucleotide array pre-processing              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XExonProcesSet::XExonProcesSet()
               :XGenomeProcesSet()
{
   // Default ExonChipProcesSet constructor
   if(kCS) cout << "---XExonProcesSet::XExonProcesSet(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XExonProcesSet::XExonProcesSet(const char *name, const char *type)
               :XGenomeProcesSet(name, type)
{
   // Normal ExonChipProcesSet constructor
   if(kCS) cout << "---XExonProcesSet::XExonProcesSet------" << endl;

}//Constructor

//______________________________________________________________________________
XExonProcesSet::~XExonProcesSet()
{
   // ExonChipProcesSet destructor
   if(kCS) cout << "---XExonProcesSet::~XExonProcesSet------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XExonProcesSet::DetectCall(Int_t numdata, TTree **datatree,
                      Int_t &numbgrd, TTree **bgrdtree)
{
   // Detect presence call
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   // Note: each datatree is processed individually and thus can have different
   //       number of entries (i.e. belong to different chip types)
   if(kCS) cout << "------XExonProcesSet::DetectCall------" << endl;

   Int_t err = errNoErr;

   fFile->cd();

// Init local arrays to store data from trees
   Int_t    *arrUnit  = 0;
   Int_t    *mskUnit  = 0;
   Int_t    *arrMask  = 0;
   Double_t *arrInten = 0;
   Double_t *arrStdev = 0;  //not needed for current fCaller
   Int_t    *arrNPix  = 0;  //not needed for current fCaller
   Double_t *arrBgrd  = 0;
   Double_t *arrBgdev = 0;  //not needed for current fCaller
   Double_t *arrPM    = 0; 
   Double_t *arrMM    = 0;
   Double_t *arrSP    = 0;  //not needed for current fCaller
   Double_t *arrSM    = 0;  //not needed for current fCaller
   Int_t    *arrXP    = 0;  //not needed for current fCaller
   Int_t    *arrXM    = 0;  //not needed for current fCaller

   TTree  *calltree = 0;
   XPCall *call     = 0;
   Int_t   split    = 99;

// Calculate detection call
   Bool_t doGC = (Bool_t)(strcmp(fCaller->GetName(), "dabgcall") == 0);
   Int_t  ij, idx, x, y, start, end;
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;

   // Informing user
      TString name = datatree[k]->GetName();
      if (XManager::fgVerbose) {
         cout << "   Calculating present call for <" << name << ">..." << endl;
      }//if

   // Get tree info for datatree
      XDataTreeInfo *info = 0;
      info = (XDataTreeInfo*)datatree[k]->GetUserInfo()->FindObject(name);
      if (!info) {
         cerr << "Error: Could not get tree info for <" << name << ">." << endl;
         return errGeneral;
      }//if

   // Get chip parameters from scheme file (also alternative CDFs)
      if (strcmp(fSchemeName.Data(), "") == 0) {
         fSchemeName = datatree[k]->GetTitle();
      } else if (!fSchemeName.Contains(datatree[k]->GetTitle())) {
         return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[k]->GetTitle());
      }//if
      if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

      XExonChip *chip = (XExonChip*)fSchemes->FindObject(fSchemeName, kTRUE);
      if (!chip) {
         return fManager->HandleError(errGetScheme, fSchemeName);
      }//if
      Int_t numrows  = chip->GetNumRows();
      Int_t numcols  = chip->GetNumColumns();
      Int_t numctrls = chip->GetNumControls();
      Int_t size     = numrows*numcols;

      Int_t    numabsent  = 0;
      Int_t    numarginal = 0;
      Int_t    numpresent = 0;
      Double_t minpval    = 1.0;
      Double_t maxpval    = 0.0;

   // Get scheme tree for scheme
      XExonScheme *scheme = 0;
      TLeaf *scmleaf = 0;
      TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
      if (scmtree == 0) return errGetTree;
      scmtree->SetBranchAddress("ScmBranch", &scheme);

   // Get unit tree for scheme
      XGCUnit *unit    = 0;
      TTree   *idxtree = 0; 
      Int_t   numunits = 0;
      if (strcmp(fCaller->GetOption(), "exon") == 0) { 
         idxtree  = (TTree*)gDirectory->Get(chip->GetExonUnitTree()); //tree.exn
         if (idxtree == 0) return errGetTree;

         numunits = chip->GetNumExonUnits();
         scmleaf  = scmtree->FindLeaf("fExonID");      
      } else if (strcmp(fCaller->GetOption(), "probeset") == 0) {
         idxtree  = (TTree*)gDirectory->Get(chip->GetProbesetUnitTree()); //tree.pbs
         if (idxtree == 0) return errGetTree;

         numunits = idxtree->GetEntries();
         scmleaf  = scmtree->FindLeaf("fProbesetID");      
      } else {
         idxtree  = (TTree*)gDirectory->Get(chip->GetUnitTree());     //tree.idx
         if (idxtree == 0) return errGetTree;

         numunits = chip->GetNumUnits();
         scmleaf  = scmtree->FindLeaf("fUnitID");      
      }//if
      idxtree->SetBranchAddress("IdxBranch", &unit);

   // Get maximum number of cells from tree info for unit tree
      XGenomeTreeInfo *idxinfo = 0;
      idxinfo = (XGenomeTreeInfo*)idxtree->GetUserInfo()->FindObject(idxtree->GetName());
      if (!idxinfo) {
         cerr << "Error: Could not get tree info for <" << idxtree->GetName() << ">." << endl;
         return errGeneral;
      }//if
      Int_t maxnumcells = (Int_t)idxinfo->GetValue("fMaxNCells");

   // Get exon level of annotation
      Int_t level = 0;
      if (fCallSelector->GetNumParameters() > 0) {
         level = (Int_t)(fCallSelector->GetParameters())[0];
      }//if

   // Init unit selector
      XUnitSelector *unitSelector = 0;

   // Initialize memory for unit arrays
      if (!(arrUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}
      if (!(mskUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 

   // Initialize memory for data arrays (uncomment when needed in fCaller)
      if (!(arrMask  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
      if (!(arrInten = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
//      if (!(arrStdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
//      if (!(arrNPix  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<size; i++) {
         arrMask[i]  = eINITMASK;  //init mask
         arrInten[i] = 0.0; 
      }//for_i

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

   // Get mask for PM/MM from scheme tree and store in array 
      arrMask = this->FillMaskArray(chip, scmtree, scheme, 0, size, arrMask);
      if (arrMask == 0) {err = errInitMemory; goto cleanup;}

   // Calculate units satisfying mask
      unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
      unitSelector->SetOption("exon");
      err = unitSelector->InitParameters(fCallSelector->GetNumParameters(),
                                         fCallSelector->GetParameters());
      if (err != errNoErr) goto cleanup;

   // Get mask from scheme tree and store in array 
      arrUnit = this->FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);

      err = unitSelector->Calculate(numunits, arrUnit, mskUnit);
      if (err != errNoErr) goto cleanup;

   // Calculate mask for detection call
      err = fCallSelector->Calculate(size, 0, 0, arrMask);
      if (err != errNoErr) goto cleanup;

   // Get data from datatree and store in arrays
      err = FillDataArrays(datatree[k], numrows, numcols, arrInten, arrStdev, arrNPix);
      if (err != errNoErr) goto cleanup;

   // For DABG only, fill arrMask with GC content (GC>=0 for PM)
   // for PM set GC >= 0 , i.e. GC = 0...kProbeLength,
   // for MM set GC < eINITMASK, i.e. GC = eINITMASK - (1...kProbeLength+1)
      if (doGC) {
      // Get probe tree for scheme
         XGCProbe *probe = 0;
         TTree *probetree = (TTree*)gDirectory->Get(chip->GetProbeTree()); 
         if (probetree == 0) {err = errGetTree; goto cleanup;}
         probetree->SetBranchAddress("PrbBranch", &probe);

      // Get GC content from probe tree and store in array
         for (Int_t i=0; i<size; i++) {
            probetree->GetEntry(i);

            x  = probe->GetX();
            y  = probe->GetY();
            ij = XY2Index(x, y, numcols);

            if (arrMask[ij] == 1) {
               arrMask[ij] = probe->GetNumberGC();
            } else if (arrMask[ij] == 0) {
               //need to use (numberGC + 1) to avoid setting arrMask=eINITMASK for GC=0!!
               arrMask[ij] = eINITMASK - (probe->GetNumberGC() + 1);
//no               arrMask[ij] = eINITMASK - probe->GetNumberGC();
            }//if
         }//for_i

         // get intensities and save in table sorted for GC-content
         if ((err = fCaller->Calculate(size, arrInten, arrMask))) goto cleanup;
      } else {
      // Need to get background data for MM values!
         if (bgrdtree[k] == 0) {
            cerr << "Error: Could not get background tree for <" << fCaller->GetName() << ">."
                 << endl;
            err = errGetTree; goto cleanup;
         }//if

         if (!(arrBgrd  = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
//         if (!(arrBgdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}

        // get data from bgrdtree and store in arrays
//         err = FillBgrdArrays(bgrdtree[k], numrows, numcols, arrBgrd, arrBgdev);
         err = FillBgrdArrays(bgrdtree[k], numrows, numcols, arrBgrd, 0);
         if (err != errNoErr) goto cleanup;
      }//if

   // Create new tree calltree
      if (!fFile->cd(fName)) return errGetDir;

      name = Path2Name(datatree[k]->GetName(),"/",".") + "." + fCaller->GetTitle();
      calltree = new TTree(name, fSchemeName);
      if (calltree == 0) return errCreateTree;
      call = new XPCall();
      calltree->Branch("CallBranch", "XPCall", &call, 64000, split);

   // Initialize memory for PM/MM arrays (uncomment when needed in fCaller)
      if (!(arrPM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrMM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
//      if (!(arrSP = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
//      if (!(arrSM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrXP = new (nothrow) Int_t[maxnumcells]))    {err = errInitMemory; goto cleanup;}
//      if (!(arrXM = new (nothrow) Int_t[maxnumcells]))    {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<maxnumcells; i++) {
         arrPM[i] = arrMM[i] = 0.0; 
         arrXP[i] = 0; 
      }//for_i

   // Calculate detection call values
      start = 0;
      end   = 0;
      idx   = 0;
      for (Int_t id=0; id<numunits; id++) { 
         idxtree->GetEntry(id);

         Int_t unitID   = unit->GetUnitID();
         Int_t numcells = unit->GetNumCells();
         // skip masked unit entries
         if (mskUnit[id] <= 0) {
            start += numcells;
            end = start;
            continue;
         }//if

         Int_t p = 0;
         Int_t m = 0;
         end += numcells;
         for (Int_t j=start; j<end; j++) {
            scmtree->GetEntry(j);

            if ((Int_t)(scmleaf->GetValue()) != unitID) {
               cerr << "Error: unitID is not equal to: " << unitID << endl;
               err = errAbort;
               goto cleanup;
            }//if

            x  = scheme->GetX();
            y  = scheme->GetY();
            ij = XY2Index(x, y, numcols);

            if (doGC) {
               if (arrMask[ij] >= 0) {
                  if (p == 0) idx++;  //count number of units
                  arrPM[p] = arrInten[ij];
                  arrXP[p] = arrMask[ij];
                  p++;
               }//if
            } else {
               if (arrMask[ij] == 1) {
                  if (p == 0) idx++;  //count number of units
                  arrPM[p] = arrInten[ij];
//no                  arrSP[p] = arrStdev[ij];
//no                  arrXP[p] = arrNPix[ij];
                  p++;

               // set MM values to background data
                  arrMM[m] = arrBgrd[ij];
//no                  arrSM[m] = arrBgdev[ij];
//no                  arrXM[p] = arrNPix[ij];
                  m++;
               }//if
            }//if

            if (p > maxnumcells || m > maxnumcells) {
               cerr << "Error: unitID <" << unitID << "> exceeds maximum number of cells <"
                    << maxnumcells << ">. Buffer overflow!" << endl;
               err = errAbort;
               goto cleanup;
            }//if
         }//for_j
         start += numcells;

         // continue if arrays arrPM etc are not filled
         if (p == 0) continue;
         if (!doGC && (p != m)) {
            cout << "Warning: Skipping unitID <" << unitID
                 << "> with different numbers of PM and MM data." << endl;
//            continue;
            err = errAbort;
            goto cleanup;
         }//if

         if (XManager::fgVerbose && id%10000 == 0) {
            cout << "      <" << idx << "> of <" << numunits << "> calls processed...\r" << flush;
         }//if

         // calculate detection call
         Int_t    dummy    = 0;
         Double_t prescall = 0.0;
         Double_t pvalue   = 1.0;
         fCaller->InitTreeInfo(info);
         fCaller->InitArrays(p, arrPM, arrSP, arrXP, arrMM, arrSM, arrXM);
         if ((err = fCaller->Calculate(prescall, pvalue, dummy))) goto cleanup;

         // number of present/absent calls
         if      (prescall == 2.0) numpresent++;
         else if (prescall == 1.0) numarginal++;
         else if (prescall == 0.0) numabsent++;

         // minimal/maximal detection call p-value
         if      (pvalue < minpval) minpval = pvalue;
         else if (pvalue > maxpval) maxpval = pvalue;

         // fill call tree
         call->SetUnitID(unitID);
         call->SetCall((Short_t)prescall);
         call->SetPValue(pvalue);
         calltree->Fill();
      }//for_id

      if (XManager::fgVerbose) {
         cout << "      <" << idx << "> of <" << numunits << "> calls processed...Finished" << endl;
      }//if

      if (XManager::fgVerbose) {
         cout << "      detection call statistics: " << endl;
         cout << "         minimum detection p-value = " << minpval << endl;
         cout << "         maximum detection p-value = " << maxpval << endl;
         cout << "         P: <" << numpresent*100.0/idx << "> percent units present."  << endl;
         cout << "         M: <" << numarginal*100.0/idx << "> percent units marginal." << endl;
         cout << "         A: <" << numabsent*100.0/idx  << "> percent units absent."   << endl;
      }//if

   // Add tree info to tree
      AddCallTreeInfo(calltree, calltree->GetName(), fCaller->GetOption(),
                      numunits-numctrls, numabsent, numarginal, numpresent,
                      minpval, maxpval);

   // Write call tree to file 
      if ((err = WriteTree(calltree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(calltree->GetName(), "Call", 0, fCaller->GetNumParameters(),
                       fCaller->GetParameters());
      }//if

   cleanup:
      // delete arrays
      if (arrXM)    {delete [] arrXM;    arrXM    = 0;}
      if (arrXP)    {delete [] arrXP;    arrXP    = 0;}
      if (arrSM)    {delete [] arrSM;    arrSM    = 0;}
      if (arrSP)    {delete [] arrSP;    arrSP    = 0;}
      if (arrMM)    {delete [] arrMM;    arrMM    = 0;}
      if (arrPM)    {delete [] arrPM;    arrPM    = 0;}
      if (arrBgdev) {delete [] arrBgdev; arrBgdev = 0;}
      if (arrBgrd)  {delete [] arrBgrd;  arrBgrd  = 0;}
      if (arrNPix)  {delete [] arrNPix;  arrNPix  = 0;}
      if (arrStdev) {delete [] arrStdev; arrStdev = 0;}
      if (arrInten) {delete [] arrInten; arrInten = 0;}
      if (arrMask)  {delete [] arrMask;  arrMask  = 0;}
      if (mskUnit)  {delete [] mskUnit;  mskUnit  = 0;}
      if (arrUnit)  {delete [] arrUnit;  arrUnit  = 0;}

      if (err != errNoErr) break;
   }//for_k

   return err;
}//DetectCall

//______________________________________________________________________________
Int_t XExonProcesSet::DoExpress(Int_t numdata, TTree **datatree,
                      Int_t numbgrd, TTree **bgrdtree)
{
   // Condense data from  datatrees to expression values
   // Note: each datatree is processed individually and thus can have different
   //       number of entries (i.e. belong to different chip types)
   if(kCS) cout << "------XExonProcesSet::DoExpress------" << endl;

   Int_t err = errNoErr;
   Int_t split = 99;

// Init local arrays to store data from trees
   Int_t    *arrUnit  = 0;
   Int_t    *mskUnit  = 0;
   Int_t    *arrMask  = 0;
   Double_t *arrInten = 0;
   Double_t *arrStdev = 0;
   Int_t    *arrNPix  = 0;
   Double_t *arrBgrd  = 0;
   Double_t *arrBgdev = 0;
   Double_t *arrPM    = 0; 
   Double_t *arrMM    = 0;
   Double_t *arrSP    = 0; 
   Double_t *arrSM    = 0;
   Int_t    *arrXP    = 0; 
   Int_t    *arrXM    = 0;

// Calculate expression
   Int_t ij, idx, x, y, start, end;
   Int_t    arrlen = 0;
   Double_t mean   = 0;
   Double_t var    = 0;

   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;
      if (bgrdtree[k] == 0) {
         cerr << "Error: Could not get background tree for <"
              << fExpressor->GetName() << ">." << endl;
         return errGetTree;
      }//if

   // Informing user
      TString name = datatree[k]->GetName();
      if (XManager::fgVerbose) {
         cout << "      Summarizing <" << name << "> using <" << fExpressor->GetName()
              << ">..." << endl;
      }//if

   // Get tree info from datatree
      fDataFile->cd();

      XDataTreeInfo *info = 0;
      info = (XDataTreeInfo*)datatree[k]->GetUserInfo()->FindObject(name);
      if (!info) {
         cerr << "Error: Could not get tree info for <" << name << ">." << endl;
         return errGeneral;
      }//if

   // Get chip parameters from scheme file (also alternative CDFs)
      if (strcmp(fSchemeName.Data(), "") == 0) {
         fSchemeName = datatree[k]->GetTitle();
      } else if (!fSchemeName.Contains(datatree[k]->GetTitle())) {
         return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[k]->GetTitle());
      }//if
      if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

      XExonChip *chip = (XExonChip*)fSchemes->FindObject(fSchemeName, kTRUE);
      if (!chip) {
         return fManager->HandleError(errGetScheme, fSchemeName);
      }//if
      Int_t numrows  = chip->GetNumRows();
      Int_t numcols  = chip->GetNumColumns();
      Int_t numctrls = chip->GetNumControls();
      Int_t size     = numrows*numcols;

   // Init min/max expression levels
      Double_t min = DBL_MAX;  //defined in float.h
      Double_t max = 0;

   // Get scheme tree for scheme
      XExonScheme *scheme = 0;
      TLeaf *scmleaf = 0;
      TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
      if (scmtree == 0) return errGetTree;
      scmtree->SetBranchAddress("ScmBranch", &scheme);

   // Get unit tree for scheme
      XGCUnit *unit    = 0;
      TTree   *idxtree = 0; 
      Int_t   numunits = 0;
      if (strcmp(fExpressor->GetOption(), "exon") == 0) { 
         idxtree  = (TTree*)gDirectory->Get(chip->GetExonUnitTree()); //tree.exn
         if (idxtree == 0) return errGetTree;

         numunits = chip->GetNumExonUnits();
         scmleaf  = scmtree->FindLeaf("fExonID");      
      } else if (strcmp(fExpressor->GetOption(), "probeset") == 0) {
         idxtree  = (TTree*)gDirectory->Get(chip->GetProbesetUnitTree()); //tree.pbs
         if (idxtree == 0) return errGetTree;

         numunits = idxtree->GetEntries();
         scmleaf  = scmtree->FindLeaf("fProbesetID");      
      } else {
         idxtree  = (TTree*)gDirectory->Get(chip->GetUnitTree());     //tree.idx
         if (idxtree == 0) return errGetTree;

         numunits = chip->GetNumUnits();
         scmleaf  = scmtree->FindLeaf("fUnitID");      
      }//if
      idxtree->SetBranchAddress("IdxBranch", &unit);

   // Get maximum number of cells from tree info for unit tree
      XGenomeTreeInfo *idxinfo = 0;
      idxinfo = (XGenomeTreeInfo*)idxtree->GetUserInfo()->FindObject(idxtree->GetName());
      if (!idxinfo) {
         cerr << "Error: Could not get tree info for <" << idxtree->GetName() << ">." << endl;
         err = errGeneral;
         break;
      }//if
      Int_t maxnumcells = (Int_t)idxinfo->GetValue("fMaxNCells");

   // Create new tree exprtree
      if (!fFile->cd(fName)) return errGetDir;

      TString dataname = Path2Name(datatree[k]->GetName(),"/",".");
      TString exprname = dataname + "." + fExpressor->GetTitle();
      TTree  *exprtree = new TTree(exprname, fSchemeName);
      if (exprtree == 0) return errCreateTree;

      XGCExpression *expr = 0;
      expr = new XGCExpression();
      exprtree->Branch("ExprBranch", "XGCExpression", &expr, 64000, split);

   // Get exon level of annotation
      Int_t level = 0;
      if (fExprSelector->GetNumParameters() > 0) {
         level = (Int_t)(fExprSelector->GetParameters())[0];
      }//if

   // Init unit selector
      XUnitSelector *unitSelector = 0;

   // Initialize memory for unit arrays
      if (!(arrUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}
      if (!(mskUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 

   // Initialize memory for local arrays
      if (!(arrMask  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
      if (!(arrInten = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrStdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrNPix  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
      if (!(arrBgrd  = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
      if (!(arrBgdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<size; i++) {
         arrMask[i]  = eINITMASK;  //init mask
         arrInten[i] = arrStdev[i] = 0.0; 
         arrNPix[i]  = 0; 
         arrBgrd[i]  = arrBgdev[i] = 0.0;
      }//for_i

   // Initialize memory for PM/MM arrays
      if (!(arrPM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrMM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrSP = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrSM = new (nothrow) Double_t[maxnumcells])) {err = errInitMemory; goto cleanup;}
      if (!(arrXP = new (nothrow) Int_t[maxnumcells]))    {err = errInitMemory; goto cleanup;}
      if (!(arrXM = new (nothrow) Int_t[maxnumcells]))    {err = errInitMemory; goto cleanup;}

      for (Int_t i=0; i<maxnumcells; i++) {
         arrPM[i] = arrMM[i] = 0.0; 
         arrSP[i] = arrSM[i] = 0.0; 
         arrXP[i] = arrXM[i] = 0; 
      }//for_i

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

   // Get mask for PM/MM from scheme tree and store in array 
      arrMask = this->FillMaskArray(chip, scmtree, scheme, 0, size, arrMask);
      if (arrMask == 0) {err = errInitMemory; goto cleanup;}

   // Calculate units satisfying mask
      unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
      unitSelector->SetOption("exon");
      err = unitSelector->InitParameters(fExprSelector->GetNumParameters(),
                                         fExprSelector->GetParameters());
      if (err != errNoErr) goto cleanup;

   // Get mask from scheme tree and store in array 
      arrUnit = this->FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);

      err = unitSelector->Calculate(numunits, arrUnit, mskUnit);
      if (err != errNoErr) goto cleanup;

   // Calculate mask for expression
      err = fExprSelector->Calculate(size, 0, 0, arrMask);
      if (err != errNoErr) goto cleanup;

   // Get data from datatree and bgrdtree and store in arrays
      err = FillDataArrays(datatree[k], numrows, numcols, arrInten, arrStdev, arrNPix);
      if (err != errNoErr) goto cleanup;
      err = FillBgrdArrays(bgrdtree[k], numrows, numcols, arrBgrd, arrBgdev);
      if (err != errNoErr) goto cleanup;

   // Calculate expression values
      start = 0;
      end   = 0;
      idx   = 0;
      for (Int_t id=0; id<numunits; id++) { 
         idxtree->GetEntry(id);

         Int_t unitID   = unit->GetUnitID();
         Int_t numcells = unit->GetNumCells();
         // skip masked unit entries
         if (mskUnit[id] <= 0) {
            start += numcells;
            end = start;
            continue;
         }//if

         Int_t p = 0;
         Int_t m = 0;
         end += numcells;
         for (Int_t j=start; j<end; j++) {
            scmtree->GetEntry(j);

            if ((Int_t)(scmleaf->GetValue()) != unitID) {
               cerr << "Error: unitID is not equal to: " << unitID << endl;
               err = errAbort;
               goto cleanup;
            }//if

            x  = scheme->GetX();
            y  = scheme->GetY();
            ij = XY2Index(x, y, numcols);

            if (arrMask[ij] == 1) {
               if (p == 0) idx++;  //count number of units to be summarized
               arrPM[p] = arrInten[ij];
               arrSP[p] = arrStdev[ij];
               arrXP[p] = arrNPix[ij];
               p++;

               // set MM values to background data
               arrMM[m] = arrBgrd[ij];
               arrSM[m] = arrBgdev[ij];
               arrXM[m] = arrNPix[ij]; //use numpix from datatree
               m++;
            }//if
         }//for_j
         start += numcells;

         if (XManager::fgVerbose && id%100000 == 0) {
            cout << "      calculating expression for <" << idx << "> of <"
                 << numunits << "> units...\r" << flush;
         }//if

         // continue if arrays arrPM etc are not filled
         if (p == 0) continue;
         if (p != m) {
            cout << "Warning: Skipping unitID <" << unitID
                 << "> with different numbers of PM and MM data." << endl;
            continue;
//?            err = errAbort;
//?            goto cleanup;
         }//if

         // calculate mean expression level
         arrlen = 0;
         mean   = 0;
         var    = 0;
         fExpressor->InitTreeInfo(info);
         fExpressor->InitArrays(p, arrPM, arrSP, arrXP, arrMM, arrSM, arrXM);
         if ((err = fExpressor->CreateArray(p)))               goto cleanup;
         if ((err = fExpressor->Calculate(mean, var, arrlen))) goto cleanup;
         fExpressor->DeleteArray();

         // get minimal/maximal expression levels
         if (mean < min) min = mean;
         if (mean > max) max = mean;

         // fill expression tree
         expr->SetUnitID(unitID);
         expr->SetLevel(mean);
         expr->SetStdev(TMath::Sqrt(var));
         expr->SetNumPairs(arrlen);
         exprtree->Fill();
      }//for_id

      if (XManager::fgVerbose) {
         cout << "      calculating expression for <" << idx << "> of <" << numunits
              << "> units...Finished." << endl;
      }//if

      if (XManager::fgVerbose) {
         cout << "      expression statistics: " << endl;
         cout << "         minimal expression level is <" << min << ">." << endl;
         cout << "         maximal expression level is <" << max << ">." << endl;
      }//if

   // Add tree info to tree
      AddExprTreeInfo(exprtree, exprtree->GetName(), fExpressor->GetOption(),
                      numunits-numctrls, min, max);

   // Write expression tree to file 
      if ((err = WriteTree(exprtree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(exprtree->GetName(), "Expr", 0, fExpressor->GetNumParameters(),
                       fExpressor->GetParameters());
      }//if

   cleanup:
      // delete arrays
      if (arrXM)    {delete [] arrXM;    arrXM    = 0;}
      if (arrXP)    {delete [] arrXP;    arrXP    = 0;}
      if (arrSM)    {delete [] arrSM;    arrSM    = 0;}
      if (arrSP)    {delete [] arrSP;    arrSP    = 0;}
      if (arrMM)    {delete [] arrMM;    arrMM    = 0;}
      if (arrPM)    {delete [] arrPM;    arrPM    = 0;}
      if (arrBgdev) {delete [] arrBgdev; arrBgdev = 0;}
      if (arrBgrd)  {delete [] arrBgrd;  arrBgrd  = 0;}
      if (arrNPix)  {delete [] arrNPix;  arrNPix  = 0;}
      if (arrStdev) {delete [] arrStdev; arrStdev = 0;}
      if (arrInten) {delete [] arrInten; arrInten = 0;}
      if (arrMask)  {delete [] arrMask;  arrMask  = 0;}
      if (mskUnit)  {delete [] mskUnit;  mskUnit  = 0;}
      if (arrUnit)  {delete [] arrUnit;  arrUnit  = 0;}
      // Note: do not remove exprtree and expr from RAM, needed later!

      if (err != errNoErr) break;
   }//for_k

   return err;
}//DoExpress

//______________________________________________________________________________
Int_t XExonProcesSet::DoMedianPolish(Int_t numdata, TTree **datatree,
                      Int_t numbgrd, TTree **bgrdtree)
{
   // Compute expression values using mediapolish
   // Intensities from datatrees are stored in one large table in RAM for fast access
   // Note: all trees must have same number of entries (i.e. identical chip types)
   if(kCS) cout << "------XExonProcesSet::DoMedianPolish(table)------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "      Summarizing with medianpolish..." << endl;
   }//if

//TEST Benchmark
//gBenchmark->Reset(); 
//gBenchmark->Start("Bench_MedianPolish");

   Int_t ij, start, end, id;
   Int_t x   = 0;
   Int_t y   = 0;
   Int_t idx = 0;
   Int_t err = errNoErr;

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fExpressor, fExpressor->GetBgrdOption());

// Get chip parameters from scheme file (also alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XExonChip *chip = (XExonChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) {
      return fManager->HandleError(errGetScheme, fSchemeName);
   }//if
   Int_t numrows  = chip->GetNumRows();
   Int_t numcols  = chip->GetNumColumns();

// Get scheme tree for scheme
   XExonScheme *scheme = 0;
   TLeaf *scmleaf = 0;
   TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
   if (scmtree == 0) return errGetTree;
   scmtree->SetBranchAddress("ScmBranch", &scheme);

// Get unit tree for scheme
   XGCUnit *unit    = 0;
   TTree   *idxtree = 0; 
   Int_t   numunits = 0;
   if (strcmp(fExpressor->GetOption(), "exon") == 0) { 
      idxtree  = (TTree*)gDirectory->Get(chip->GetExonUnitTree()); //tree.exn
      if (idxtree == 0) return errGetTree;

      numunits = chip->GetNumExonUnits();
      scmleaf  = scmtree->FindLeaf("fExonID");      
   } else if (strcmp(fExpressor->GetOption(), "probeset") == 0) {
      idxtree  = (TTree*)gDirectory->Get(chip->GetProbesetUnitTree()); //tree.pbs
      if (idxtree == 0) return errGetTree;

      numunits = idxtree->GetEntries();
      scmleaf  = scmtree->FindLeaf("fProbesetID");      
   } else {
      idxtree  = (TTree*)gDirectory->Get(chip->GetUnitTree());     //tree.idx
      if (idxtree == 0) return errGetTree;

      numunits = chip->GetNumUnits();
      scmleaf  = scmtree->FindLeaf("fUnitID");      
   }//if
   idxtree->SetBranchAddress("IdxBranch", &unit);

// Init size of arrays
   Int_t size    = numrows*numcols;
   Int_t numsels = 0;  //number of selected entries

// Check for equal number of data tree entries
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;

      if ((Int_t)(datatree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, datatree[k]->GetName(), str);
      }//if
   }//for_k

// Check for equal number of background tree entries
   if (numbgrd > 0) {
      for (Int_t k=0; k<numdata; k++) {
         if (bgrdtree[k] == 0) return errGetTree;

         if ((Int_t)(bgrdtree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, bgrdtree[k]->GetName(), str);
         }//if
      }//for_k
//   } else if (subbgrd || corbgrd) {
   } else if (doBg == kTRUE) {
      cout << "Warning: No background trees available for background subtraction."
           << endl;
      doBg = kFALSE;
   } else {
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Init branch addresses
   XBgCell *bgcell[numdata];
   XGCCell *gccell[numdata];
   for (Int_t k=0; k<numdata; k++) {
      bgcell[k] = 0;
      gccell[k] = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell[k]);
      if (numbgrd > 0) bgrdtree[k]->SetBranchAddress("BgrdBranch", &bgcell[k]);
   }//for_k

// Init expression trees
   Int_t  split = 99;
   TTree *exprtree[numdata];
   XGCExpression *expr[numdata];

// Init min/max expression levels
   Double_t min = DBL_MAX;  //defined in float.h
   Double_t max = 0;

// Get exon level of annotation
   Int_t level = 0;
   if (fExprSelector->GetNumParameters() > 0) {
      level = (Int_t)(fExprSelector->GetParameters())[0];
   }//if

// Init unit selector
   XUnitSelector *unitSelector = 0;

// Init 
   Int_t     *arrMask = 0;
   Int_t     *arrIndx = 0;
   Int_t     *arrUnit = 0;
   Int_t     *mskUnit = 0;
   Double_t  *colmed  = 0;
   Double_t  *results = 0;
   Double_t **table   = 0;

// Create local arrays
   if (!(arrMask = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrIndx = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrUnit = new (nothrow) Int_t[numunits]))   {err = errInitMemory; goto cleanup;}
   if (!(mskUnit = new (nothrow) Int_t[numunits]))   {err = errInitMemory; goto cleanup;}
   if (!(colmed  = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}
   if (!(results = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<size; i++)    {arrIndx[i] = 0; arrMask[i] = eINITMASK;}
   for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 
   for (Int_t i=0; i<numdata; i++)  colmed[i]  = results[i] = 0.0; 

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

// Get mask for PM from scheme tree and store in array 
   arrMask = this->FillMaskArray(chip, scmtree, scheme, level, size, arrMask);
   if (arrMask == 0) {err = errInitMemory; goto cleanup;}

// Calculate units satisfying mask
   unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
   unitSelector->SetOption("exon");
   err = unitSelector->InitParameters(fExprSelector->GetNumParameters(),
                                      fExprSelector->GetParameters());
   if (err != errNoErr) goto cleanup;

// Get units and mask from unit tree and store in array 
   arrUnit = FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);

   err = unitSelector->Calculate(numunits, arrUnit, mskUnit);
   if (err != errNoErr) goto cleanup;

// Calculate mask for expression
   err = fExprSelector->Calculate(size, 0, 0, arrMask);
   if (err != errNoErr) goto cleanup;

// Get entry index from datatree
   idx = 0;
   for (Int_t i=0; i<size; i++) {
      datatree[0]->GetEntry(i);

      x  = gccell[0]->GetX();
      y  = gccell[0]->GetY();
      ij = XY2Index(x, y, numcols);

      if (arrMask[ij] == 1) {
         arrIndx[ij] = idx++;
      }//if
   }//for_i

// Get number of selected entries
   for (Int_t i=0; i<size; i++) {
      numsels = (arrMask[i] == 1) ? ++numsels : numsels;
   }//for_i

// Create table to store intensities of all datatrees
   if (!(table = new (nothrow) Double_t*[numdata])) {err = errInitMemory; goto cleanup;} 
   for (Int_t k=0; k<numdata; k++) {
      table[k] = 0;
      if (!(table[k] = new (nothrow) Double_t[numsels])) {err = errInitMemory; goto cleanup;} 
   }//for_k

//TEST Benchmark
//gBenchmark->Start("Bench_table");
// Get data from datatrees (and bgrdtrees) and store in table
   if ((numbgrd > 0) && (doBg == kTRUE)) {
      for (Int_t k=0; k<numdata; k++) {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);
            bgrdtree[k]->GetEntry(i);
//unsafe?? use arrInten[ij] and arrBgrd[ij] ??

            // subtract background from intensity
            if (arrMask[i] == 1) {
               table[k][idx++] = this->AdjustIntensity(gccell[k]->GetIntensity(),
                                                       bgcell[k]->GetBackground(),
                                                       bgcell[k]->GetStdev());
            }//if
         }//for_i
      }//for_k
   } else {
      for (Int_t k=0; k<numdata; k++) {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);

            if (arrMask[i] == 1) {
               table[k][idx++] = gccell[k]->GetIntensity();
            }//if
         }//for_i
      }//for_k
   }//if
//TEST Benchmark
//gBenchmark->Show("Bench_table");

// Create new trees exprtree
   if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

   for (Int_t k=0; k<numdata; k++) {
      TString dataname = Path2Name(datatree[k]->GetName(),"/",".");
      TString exprname = dataname + "." + fExpressor->GetTitle();
      exprtree[k] = new TTree(exprname, fSchemeName);
      if (exprtree[k] == 0) {err = errCreateTree; goto cleanup;}

      expr[k] = new XGCExpression();
      exprtree[k]->Branch("ExprBranch", "XGCExpression", &expr[k], 64000, split);
   }//for_k

//TEST Benchmark
//gBenchmark->Show("Bench_MedianPolish");
//gBenchmark->Reset();
//gBenchmark->Start("Bench_Loop");

// Calculate expression values
   start = 0;
   end   = 0;
   idx   = 0;
   for (id=0; id<numunits; id++) { 
      idxtree->GetEntry(id);

      Int_t unitID   = unit->GetUnitID();
      Int_t numcells = unit->GetNumCells();
      // skip masked unit entries
      if (mskUnit[id] <= 0) {
         start += numcells;
         end = start;
         continue;
      }//if

//TO DO
//Better above: arrPM = new Double_t[maxnumcells*numdata]; 
      // create array to store PM values for all probes with current unitID
      Int_t  numatoms = unit->GetNumAtoms();
      Double_t *arrPM = new Double_t[numatoms*numdata]; 

      // fill arrPM with PM values of current unitID
      Int_t p = 0;
      end += numcells;
      for (Int_t j=start; j<end; j++) {
         scmtree->GetEntry(j);

         if ((Int_t)(scmleaf->GetValue()) != unitID) {
            cerr << "Error: unitID is not equal to: " << unitID << endl;
            err = errAbort;
            goto cleanup;
         }//if

         x  = scheme->GetX();
         y  = scheme->GetY();
         ij = XY2Index(x, y, numcols);

         if (arrMask[ij] == 1) {
            if (p == 0) idx++;  //count number of units to be summarized

            for (Int_t k=0; k<numdata; k++) {
               arrPM[p] = table[k][arrIndx[ij]];
               p++;
            }//for_k
         }//if
      }//for_j
      start += numcells;

      if (XManager::fgVerbose && id%100000 == 0) {
         cout << "      calculating expression for <" << idx << "> of <" << numunits
              << "> units...\r" << flush;
      }//if

      // fill arrPM or continue if it is not filled
      if ((err = fExpressor->SetArray(p, arrPM)) != errNoErr) {
         delete [] arrPM;
         continue;
      }//if

      // calculate median polish for PMs of current unitID
      if ((err = fExpressor->Calculate(numdata, colmed, results, 0))) break;

//////////////////
// TO DO: return fResiduals for residual-plot!!!! (like affyPLM)
//      residuals = fExpressor->GetResiduals();
// store as residuals tree??
//////////////////

      // fill expression trees
      for (Int_t k=0; k<numdata; k++) {
         // get minimal/maximal expression levels
         if (results[k] < min) min = results[k];
         if (results[k] > max) max = results[k];

         expr[k]->SetUnitID(unitID);
         expr[k]->SetLevel(results[k]);
         expr[k]->SetStdev(TMath::Abs(colmed[k]));
         expr[k]->SetNumPairs((Int_t)(p / numdata));
         exprtree[k]->Fill();
      }//for_k

      delete [] arrPM;
   }//for_id
   if (XManager::fgVerbose) {
      cout << "      calculating expression for <" << idx << "> of <"
           << numunits << "> units...Finished." << endl;
   }//if

   if (XManager::fgVerbose) {
      cout << "      expression statistics: " << endl;
      cout << "         minimal expression level is <" << min << ">." << endl;
      cout << "         maximal expression level is <" << max << ">." << endl;
   }//if
//TEST Benchmark
//gBenchmark->Show("Bench_Loop");

// Write expression trees to file 
   for (Int_t k=0; k<numdata; k++) {
   // Add tree info to tree
      AddExprTreeInfo(exprtree[k], exprtree[k]->GetName(), fExpressor->GetOption(),
                      idx, min, max);

      if ((err = WriteTree(exprtree[k], TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(exprtree[k]->GetName(), "Expr", 0, fExpressor->GetNumParameters(),
                       fExpressor->GetParameters());
      } else {
         break;
      }//if
   }//for_k

// Cleanup
cleanup:
   // delete table
   for (Int_t k=0; k<numdata; k++) {
      if (table[k]) {delete [] table[k]; table[k] = 0;}
   }//for_k
   if (table) delete [] table;

   // delete arrays
   if (results) {delete [] results; results = 0;}
   if (colmed)  {delete [] colmed;  colmed  = 0;}
   if (mskUnit) {delete [] mskUnit; mskUnit = 0;}
   if (arrUnit) {delete [] arrUnit; arrUnit = 0;}
   if (arrMask) {delete [] arrMask; arrMask = 0;}
   if (arrIndx) {delete [] arrIndx; arrIndx = 0;}

   return err;
}//DoMedianPolish

//______________________________________________________________________________
Int_t XExonProcesSet::DoMedianPolish(Int_t numdata, TTree **datatree,
                      Int_t numbgrd, TTree **bgrdtree, TFile *file)
{
   // Compute expression values using mediapolish
   // In order to reduce memory consumption intensities from datatrees are
   // stored in the order of the corresponding schemetree entries in temporary 
   // trees in a temporary root file
   // Note: all trees must have same number of entries (i.e. identical chip types)
   if(kCS) cout << "------XExonProcesSet::DoMedianPolish(file)------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "      Summarizing with medianpolish (using temporary file)..." << endl;
   }//if

//TEST Benchmark
//gBenchmark->Reset(); 
//gBenchmark->Start("Bench_MedianPolish");

   Int_t ij, start, end, id, entry;
   Int_t x   = 0;
   Int_t y   = 0;
   Int_t idx = 0;
   Int_t err = errNoErr;

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fExpressor, fExpressor->GetBgrdOption());

// Get chip parameters from scheme file (also alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XExonChip *chip = (XExonChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) {
      return fManager->HandleError(errGetScheme, fSchemeName);
   }//if
   Int_t numrows  = chip->GetNumRows();
   Int_t numcols  = chip->GetNumColumns();

// Get scheme tree for scheme
   XExonScheme *scheme = 0;
   TLeaf *scmleaf = 0;
   TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
   if (scmtree == 0) return errGetTree;
   scmtree->SetBranchAddress("ScmBranch", &scheme);

// Get unit tree for scheme
   XGCUnit *unit  = 0;
   TTree *idxtree = 0; 
   Int_t numunits = 0;
   if (strcmp(fExpressor->GetOption(), "exon") == 0) { 
      idxtree  = (TTree*)gDirectory->Get(chip->GetExonUnitTree()); //tree.exn
      if (idxtree == 0) return errGetTree;

      numunits = chip->GetNumExonUnits();
      scmleaf  = scmtree->FindLeaf("fExonID");      
   } else if (strcmp(fExpressor->GetOption(), "probeset") == 0) {
      idxtree  = (TTree*)gDirectory->Get(chip->GetProbesetUnitTree()); //tree.pbs
      if (idxtree == 0) return errGetTree;

      numunits = idxtree->GetEntries();
      scmleaf  = scmtree->FindLeaf("fProbesetID");      
   } else {
      idxtree  = (TTree*)gDirectory->Get(chip->GetUnitTree());     //tree.idx
      if (idxtree == 0) return errGetTree;

      numunits = chip->GetNumUnits();
      scmleaf  = scmtree->FindLeaf("fUnitID");      
   }//if
   idxtree->SetBranchAddress("IdxBranch", &unit);

// Init size of arrays
   Int_t size    = numrows*numcols;
   Int_t numsels = 0;  //number of selected entries

// Check for equal number of data tree entries
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) return errGetTree;

      if ((Int_t)(datatree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, datatree[k]->GetName(), str);
      }//if
   }//for_k

// Check for equal number of background tree entries
   if (numbgrd > 0) {
      for (Int_t k=0; k<numdata; k++) {
         if (bgrdtree[k] == 0) return errGetTree;

         if ((Int_t)(bgrdtree[k]->GetEntries()) != size) {
            TString str = ""; str += size;
            return fManager->HandleError(errNumTreeEntries, bgrdtree[k]->GetName(), str);
         }//if
      }//for_k
//   } else if (subbgrd || corbgrd) {
   } else if (doBg == kTRUE) {
      cout << "Warning: No background trees available for background subtraction."
           << endl;
      doBg = kFALSE;
   } else {
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Init branch addresses
   XBgCell *bgcell[numdata];
   XGCCell *gccell[numdata];
   for (Int_t k=0; k<numdata; k++) {
      bgcell[k] = 0;
      gccell[k] = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell[k]);
      if (numbgrd > 0) bgrdtree[k]->SetBranchAddress("BgrdBranch", &bgcell[k]);
   }//for_k

// Init temporary trees and expression trees
   Int_t    split = 99;
   Double_t sort  = 0.0;
   TTree *tmptree[numdata];
   TTree *exprtree[numdata];
   XGCExpression *expr[numdata];

// Init min/max expression levels
   Double_t min = DBL_MAX;  //defined in float.h
   Double_t max = 0;

// Get exon level of annotation
   Int_t level = 0;
   if (fExprSelector->GetNumParameters() > 0) {
      level = (Int_t)(fExprSelector->GetParameters())[0];
   }//if

// Init unit selector
   XUnitSelector *unitSelector = 0;

// Create local arrays
   Int_t    *arrMask = 0;
   Int_t    *arrIndx = 0;
   Int_t     *arrUnit = 0;
   Int_t     *mskUnit = 0;
   Double_t *colmed  = 0;
   Double_t *results = 0;
   Double_t *arrData = 0;

// Create local arrays
   if (!(arrMask = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrIndx = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrUnit = new (nothrow) Int_t[numunits]))   {err = errInitMemory; goto cleanup;}
   if (!(mskUnit = new (nothrow) Int_t[numunits]))   {err = errInitMemory; goto cleanup;}
   if (!(colmed  = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}
   if (!(results = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<size; i++)    {arrIndx[i] = 0; arrMask[i] = eINITMASK;}
   for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 
   for (Int_t i=0; i<numdata; i++)  colmed[i]  = results[i] = 0.0; 

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

// Get mask for PM from scheme tree and store in array 
   arrMask = this->FillMaskArray(chip, scmtree, scheme, level, size, arrMask);
   if (arrMask == 0) {err = errInitMemory; goto cleanup;}

// Calculate units satisfying mask
   unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
   unitSelector->SetOption("exon");
   err = unitSelector->InitParameters(fExprSelector->GetNumParameters(),
                                      fExprSelector->GetParameters());
   if (err != errNoErr) goto cleanup;

// Get units and mask from unit tree and store in array 
   arrUnit = FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);

   err = unitSelector->Calculate(numunits, arrUnit, mskUnit);
   if (err != errNoErr) goto cleanup;

// Calculate mask for expression
   err = fExprSelector->Calculate(size, 0, 0, arrMask);
   if (err != errNoErr) goto cleanup;

// Get entry index from datatree
   idx = 0;
   for (Int_t i=0; i<size; i++) {
      datatree[0]->GetEntry(i);

      x  = gccell[0]->GetX();
      y  = gccell[0]->GetY();
      ij = XY2Index(x, y, numcols);

      if (arrMask[ij] == 1) {
         arrIndx[ij] = idx++;
      }//if
   }//for_i

// Get number of selected entries
   for (Int_t i=0; i<size; i++) {
      numsels = (arrMask[i] == 1) ? ++numsels : numsels;
   }//for_i

// Create array to store selected intensities
   if (!(arrData = new Double_t[numsels])) {err = errInitMemory; goto cleanup;}
   for (Int_t i=0; i<numsels; i++) arrData[i] = 0.0;

//TEST Benchmark
//gBenchmark->Start("Bench_tmptree");
// Change directory to temporary file
   if (!file->cd()) {err = errGetDir; goto cleanup;}

// Get data from datatrees (and bgrdtrees) and store in temporary file
   for (Int_t k=0; k<numdata; k++) {
      // create temporary tree
      tmptree[k] = new TTree(datatree[k]->GetName(), "temporary tree");
      if (tmptree[k] == 0) {err = errCreateTree; goto cleanup;}
      tmptree[k]->Branch("sortBr", &sort, "sort/D");

      // fill array with (background corrected) intensities
      if ((numbgrd > 0) && (doBg == kTRUE)) {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);
            bgrdtree[k]->GetEntry(i);

            // subtract background from intensity
            if (arrMask[i] == 1) {
               arrData[idx++] = this->AdjustIntensity(gccell[k]->GetIntensity(),
                                                      bgcell[k]->GetBackground(),
                                                      bgcell[k]->GetStdev());
            }//if
         }//for_i
      } else {
         idx = 0;
         for (Int_t i=0; i<size; i++) {
            datatree[k]->GetEntry(i);

            if (arrMask[i] == 1) {
               arrData[idx++] = gccell[k]->GetIntensity();
            }//if
         }//for_i
      }//if

      // fill tmptree with array in the order of scheme tree entries for (x,y)
      for (Int_t i=0; i<size; i++) {
         scmtree->GetEntry(i);

         x  = scheme->GetX();
         y  = scheme->GetY();
         ij = XY2Index(x, y, numcols);

         if (arrMask[ij] == 1) {
            sort = arrData[arrIndx[ij]];
            tmptree[k]->Fill();
         }//if
      }//for_i

      // write tmptree to temporary file
      tmptree[k]->Write();
//TEST??      tmptree[k]->Write(TObject::kOverwrite);
   }//for_k
//TEST Benchmark
//gBenchmark->Show("Bench_tmptree");

// Change directory to current directory for treeset in main file
   if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

// Create new trees exprtree
   for (Int_t k=0; k<numdata; k++) {
      TString dataname = Path2Name(datatree[k]->GetName(),"/",".");
      TString exprname = dataname + "." + fExpressor->GetTitle();
      exprtree[k] = new TTree(exprname, fSchemeName);
      if (exprtree[k] == 0) {err = errCreateTree; goto cleanup;}

      expr[k] = new XGCExpression();
      exprtree[k]->Branch("ExprBranch", "XGCExpression", &expr[k], 64000, split);
   }//for_k

//TEST Benchmark
//gBenchmark->Show("Bench_MedianPolish");
//gBenchmark->Reset();
//gBenchmark->Start("Bench_Loop");

// Calculate expression values
   start = 0;
   end   = 0;
   idx   = 0;
   entry = 0;
   for (id=0; id<numunits; id++) { 
      idxtree->GetEntry(id);

      Int_t unitID   = unit->GetUnitID();
      Int_t numcells = unit->GetNumCells();
      // skip masked unit entries
      if (mskUnit[id] <= 0) {
         start += numcells;
         end = start;
         continue;
      }//if

//Better above: arrPM = new Double_t[maxnumcells*numdata]; 
      // create array to store PM values for all probes with current unitID
      Int_t numatoms = unit->GetNumAtoms();
      Double_t *arrPM = new Double_t[numatoms*numdata]; 

      // fill arrPM with PM values of current unitID
      Int_t p = 0;
      end += numcells;
      for (Int_t j=start; j<end; j++) {
         scmtree->GetEntry(j);

         if ((Int_t)(scmleaf->GetValue()) != unitID) {
            cerr << "Error: unitID is not equal to: " << unitID << endl;
            err = errAbort;
            goto cleanup;
         }//if

         x  = scheme->GetX();
         y  = scheme->GetY();
         ij = XY2Index(x, y, numcols);

         if (arrMask[ij] == 1) {
            if (p == 0) idx++;  //count number of units to be summarized

            for (Int_t k=0; k<numdata; k++) {
               tmptree[k]->GetEntry(entry);
               arrPM[p] = sort;
               p++;
            }//for_k

            entry++;
         }//if
      }//for_j
      start += numcells;

      if (XManager::fgVerbose && id%100000 == 0) {
         cout << "      calculating expression for <" << idx << "> of <" << numunits
              << "> units...\r" << flush;
      }//if

      // fill arrPM or continue if it is not filled
      if ((err = fExpressor->SetArray(p, arrPM)) != errNoErr) {
         delete [] arrPM;
         continue;
      }//if

      // calculate median polish for PMs of current unitID
      if ((err = fExpressor->Calculate(numdata, colmed, results, 0))) break;

//////////////////
// TO DO: return fResiduals for residual-plot!!!! (like affyPLM)
//      residuals = fExpressor->GetResiduals();
// store as residuals tree??
//////////////////

      // fill expression trees
      for (Int_t k=0; k<numdata; k++) {
         // get minimal/maximal expression levels
         if (results[k] < min) min = results[k];
         if (results[k] > max) max = results[k];

         expr[k]->SetUnitID(unitID);
         expr[k]->SetLevel(results[k]);
         expr[k]->SetStdev(TMath::Abs(colmed[k]));
         expr[k]->SetNumPairs((Int_t)(p / numdata));
         exprtree[k]->Fill();
      }//for_k

      delete [] arrPM;
   }//for_id
   if (XManager::fgVerbose) {
      cout << "      calculating expression for <" << idx << "> of <"
           << numunits << "> units...Finished." << endl;
   }//if

   if (XManager::fgVerbose) {
      cout << "      expression statistics: " << endl;
      cout << "         minimal expression level is <" << min << ">." << endl;
      cout << "         maximal expression level is <" << max << ">." << endl;
   }//if
//TEST Benchmark
//gBenchmark->Show("Bench_Loop");

// Write expression trees to file 
   for (Int_t k=0; k<numdata; k++) {
   // Add tree info to tree
      AddExprTreeInfo(exprtree[k], exprtree[k]->GetName(), fExpressor->GetOption(),
                      idx, min, max);

      if ((err = WriteTree(exprtree[k], TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(exprtree[k]->GetName(), "Expr", 0, fExpressor->GetNumParameters(),
                       fExpressor->GetParameters());
      } else {
         break;
      }//if
   }//for_k

// Cleanup
cleanup:
   // delete temporary trees
   for (Int_t k=0; k<numdata; k++) {
      tmptree[k]->Delete(""); tmptree[k] = 0;
   }//for_k

   // delete arrays
   if (arrData) {delete [] arrData; arrData = 0;}
   if (results) {delete [] results; results = 0;}
   if (colmed)  {delete [] colmed;  colmed  = 0;}
   if (mskUnit) {delete [] mskUnit; mskUnit = 0;}
   if (arrUnit) {delete [] arrUnit; arrUnit = 0;}
   if (arrIndx) {delete [] arrIndx; arrIndx = 0;}
   if (arrMask) {delete [] arrMask; arrMask = 0;}

   return err;
}//DoMedianPolish

//______________________________________________________________________________
Int_t XExonProcesSet::ExportExprTrees(Int_t n, TString *names, const char *varlist,
                      ofstream &output, const char *sep)
{
   // Export variables from varlist for expression tree(s) to file output
   if(kCS) cout << "------XExonProcesSet::ExportExprTrees------" << endl;

   Int_t err = errNoErr;

// Decompose varlist
   Bool_t hasUnit   = kFALSE;
   Bool_t hasName   = kFALSE;
   Bool_t hasSymbol = kFALSE;
   Bool_t hasCyto   = kFALSE;
   Bool_t hasAnnot  = kFALSE;
   Bool_t hasLevel  = kFALSE;
   Bool_t hasStdev  = kFALSE;
   Bool_t hasNAtoms = kFALSE;
////////////////
//TO DO
//fTranscriptID
//fAccession
//fEntrezID
//fChromosome
//fStart
//fStop
//fStrand
//for exon:
//fExonID
//for probeset
//??
//??
////////////////

   if (strcmp(varlist,"*")  == 0) {
      hasUnit   = kTRUE;
      hasName   = kTRUE;
      hasSymbol = kTRUE;
      hasCyto   = kTRUE;
      hasLevel  = kTRUE;
      hasStdev  = kTRUE;
      hasNAtoms = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName") == 0) {hasUnit   = kTRUE;}
         if (strcmp(name,"fName")     == 0) {hasName   = kTRUE;}
         if (strcmp(name,"fSymbol")   == 0) {hasSymbol = kTRUE;}
         if (strcmp(name,"fCytoBand") == 0) {hasCyto   = kTRUE;}
         if (strcmp(name,"fLevel")    == 0) {hasLevel  = kTRUE;}
         if (strcmp(name,"fStdev")    == 0) {hasStdev  = kTRUE;}
         if (strcmp(name,"fNPairs")   == 0) {hasNAtoms = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if
   hasAnnot = (hasName || hasSymbol || hasCyto);

// Get trees
   TTree       *tree[n];
   XExpression *expr[n];

   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         expr[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("ExprBranch", &expr[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         expr[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("ExprBranch", &expr[k]);
      }//for_k
   }//if

// Get treeinfo and its option for selection of unittree and anntree
   XTreeInfo *info   = (XTreeInfo*)tree[0]->GetUserInfo()->At(0);
   Option_t  *option = info->GetOption();

// Get scheme name (also for alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = tree[0]->GetTitle();
   } else if (!fSchemeName.Contains(tree[0]->GetTitle())) {
      cerr << "Error: Scheme <" << fSchemeName << "> is not derived from <"
           << tree[0]->GetTitle() << ">." << endl;
      return errAbort;
   }//if

// Get chip from scheme file
   if (fSchemeFile == 0) return errGetScheme;
   fSchemeFile->cd();
   XFolder *schemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!schemes) {
      return fManager->HandleError(errMissingContent, "Scheme", kContent);
   }//if

   XExonChip *chip = (XExonChip*)schemes->FindObject(fSchemeName, kTRUE);
   if (chip == 0) return errAbort;

   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

// Get unit tree for scheme (needed for annotation)
   XExonUnit *unit = 0;
   TTree *unittree = 0; 
   if (strcmp(option, "exon") == 0) { 
      unittree = (TTree*)gDirectory->Get(chip->GetExonUnitTree()); 
   } else if (strcmp(option, "probeset") == 0) {
      unittree = (TTree*)gDirectory->Get(chip->GetProbesetUnitTree()); 
   } else {
      unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   }//if
   if (unittree == 0) return errGetTree;
   unittree->SetBranchAddress("IdxBranch", &unit);

// Get annotation tree for scheme
   TTree *anntree  = 0; 
   XTransAnnotation *annot = 0;
   if (strcmp(option, "exon") == 0) { 
      anntree = (TTree*)gDirectory->Get(chip->GetExonAnnotTree()); 
   } else if (strcmp(option, "probeset") == 0) {
      anntree = (TTree*)gDirectory->Get(chip->GetProbesetAnnotTree()); 
   } else {
      anntree = (TTree*)gDirectory->Get(chip->GetAnnotTree()); 
   }//if
   if (anntree == 0) return errGetTree;
   anntree->SetBranchAddress("AnnBranch", &annot);

   Int_t numunits = (Int_t)(unittree->GetEntries());
   Int_t numannot = (Int_t)(anntree->GetEntries());

// Create hash table to store unit names from anntree
   THashTable *htable = 0;
   if (!(htable = new THashTable(2*numannot))) return errInitMemory;

   if (XManager::fgVerbose) {
      cout << "Reading entries from <" << anntree->GetName() << "> ...";
   }//if
   TString str;
   XIdxString *idxstr = 0;
   for (Int_t i=0; i<numannot; i++) {
      anntree->GetEntry(i);

      if (strcmp(option, "exon") == 0) { 
         str.Form("%d", ((XExonAnnotation*)annot)->GetExonID());
         idxstr = new XIdxString(i, str.Data());
      } else if (strcmp(option, "probeset") == 0) {
         str.Form("%d", ((XProbesetAnnotation*)annot)->GetProbesetID());
         idxstr = new XIdxString(i, str.Data());
      } else {
         str = ((XGenomeAnnotation*)annot)->GetTranscriptID();
         idxstr = new XIdxString(i, str.Data());
      }//if

      htable->Add(idxstr);
   }//for_i
   if (XManager::fgVerbose) {
      cout << "Finished" << endl;
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit)      output << sep << "UNIT_NAME";
   if (hasAnnot) {
      if (hasName)   output << sep << "GENE_NAME";
      if (hasSymbol) output << sep << "GENE_SYMBOL";
      if (hasCyto)   output << sep << "CYTOBAND";
   }//if
   if (n == 1) {
      if (hasLevel)  output << sep << "LEVEL";
      if (hasStdev)  output << sep << "STDEV";
      if (hasNAtoms) output << sep << "NUMBER_ATOMS";
   } else {
      for (Int_t k=0; k<n; k++) {
         if (hasLevel)  output << sep << (names[k] + "_LEVEL");
         if (hasStdev)  output << sep << (names[k] + "_STDEV");
         if (hasNAtoms) output << sep << (names[k] + "_NUMBER_ATOMS");
      }//for_k
   }//if
   output << endl;

// Loop over tree entries and trees
   Int_t idx = 0;
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);

         // export annotation
         if (k == 0) {
            Int_t unitID = expr[k]->GetUnitID();
            output << unitID;

            unittree->GetEntry(idx++);
            while (unitID != unit->GetUnitID()) {
               if (idx == numunits) {
                 cerr << "Error: UnitID <" << unitID << "> not found." << endl;
                 err = errAbort; goto cleanup;
               }//if

               unittree->GetEntry(idx++);
            }//while

            if (hasUnit) {
               output << sep << unit->GetSubUnitID();
            }//if

            if (hasAnnot) {
               str.Form("%d", unit->GetSubUnitID());
               idxstr = (XIdxString*)(htable->FindObject(str.Data()));
               if (idxstr) {
                  anntree->GetEntry(idxstr->GetIndex());
                  if (hasName)   output << sep << annot->GetName();
                  if (hasSymbol) output << sep << annot->GetSymbol();
                  if (hasCyto)   output << sep << annot->GetCytoBand();
               } else {
                  if (hasName)   output << sep << "NA";
                  if (hasSymbol) output << sep << "NA";
                  if (hasCyto)   output << sep << "NA";
               }//if
            }//if
         }//if

         // export data
         if (hasLevel)  output << sep << expr[k]->GetLevel();
         if (hasStdev)  output << sep << ((XGCExpression*)expr[k])->GetStdev();
         if (hasNAtoms) output << sep << ((XGCExpression*)expr[k])->GetNumPairs();
      }//for_k
      output << endl;
   }//for_i

//Cleanup
cleanup:
   // remove trees from RAM
   if (anntree)  {anntree->Delete("");  anntree  = 0;}
   if (unittree) {unittree->Delete(""); unittree = 0;}
   if (htable)   {htable->Delete(); delete htable; htable = 0;}
   SafeDelete(schemes);

   return err;
}//ExportExprTrees

//______________________________________________________________________________
Int_t XExonProcesSet::ExportCallTrees(Int_t n, TString *names, const char *varlist,
                      ofstream &output, const char *sep)
{
   // Export variables from varlist for call tree(s) to file output
   if(kCS) cout << "------XExonProcesSet::ExportCallTrees------" << endl;

   Int_t err = errNoErr;

// Decompose varlist
   Bool_t hasUnit   = kFALSE;
   Bool_t hasName   = kFALSE;
   Bool_t hasSymbol = kFALSE;
   Bool_t hasCyto   = kFALSE;
   Bool_t hasAnnot  = kFALSE;
   Bool_t hasCall   = kFALSE;
   Bool_t hasPVal   = kFALSE;
////////////////
//TO DO
//fTranscriptID
//fAccession
//fEntrezID
//fChromosome
//fStart
//fStop
//fStrand
//for exon:
//fExonID
//for probeset
//??
//??
////////////////

   if (strcmp(varlist,"*")  == 0) {
      hasUnit   = kTRUE;
      hasName   = kTRUE;
      hasSymbol = kTRUE;
      hasCyto   = kTRUE;
      hasCall   = kTRUE;
      hasPVal   = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName") == 0) {hasUnit   = kTRUE;}
         if (strcmp(name,"fName")     == 0) {hasName   = kTRUE;}
         if (strcmp(name,"fSymbol")   == 0) {hasSymbol = kTRUE;}
         if (strcmp(name,"fCytoBand") == 0) {hasCyto   = kTRUE;}
         if (strcmp(name,"fCall")     == 0) {hasCall   = kTRUE;}
         if (strcmp(name,"fPValue")   == 0) {hasPVal   = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if
   hasAnnot = (hasName || hasSymbol || hasCyto);

// Get trees
   TTree  *tree[n];
   XPCall *call[n];

   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         call[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("CallBranch", &call[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         call[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("CallBranch", &call[k]);
      }//for_k
   }//if

// Get treeinfo and its option for selection of unittree and anntree
   XTreeInfo *info   = (XTreeInfo*)tree[0]->GetUserInfo()->At(0);
   Option_t  *option = info->GetOption();

// Get scheme name (also for alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = tree[0]->GetTitle();
   } else if (!fSchemeName.Contains(tree[0]->GetTitle())) {
      cerr << "Error: Scheme <" << fSchemeName << "> is not derived from <"
           << tree[0]->GetTitle() << ">." << endl;
      return errAbort;
   }//if

// Get chip from scheme file
   if (fSchemeFile == 0) return errGetScheme;
   fSchemeFile->cd();
   XFolder *schemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!schemes) {
      return fManager->HandleError(errMissingContent, "Scheme", kContent);
   }//if

   XExonChip *chip = (XExonChip*)schemes->FindObject(fSchemeName, kTRUE);
   if (chip == 0) return errAbort;

   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

// Get unit tree for scheme (needed for annotation)
   XExonUnit *unit = 0;
   TTree *unittree = 0; 
   if (strcmp(option, "exon") == 0) { 
      unittree = (TTree*)gDirectory->Get(chip->GetExonUnitTree()); 
   } else if (strcmp(option, "probeset") == 0) {
      unittree = (TTree*)gDirectory->Get(chip->GetProbesetUnitTree()); 
   } else {
      unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   }//if
   if (unittree == 0) return errGetTree;
   unittree->SetBranchAddress("IdxBranch", &unit);

// Get annotation tree for scheme
   TTree *anntree  = 0; 
   XTransAnnotation *annot = 0;
   if (strcmp(option, "exon") == 0) { 
      anntree = (TTree*)gDirectory->Get(chip->GetExonAnnotTree()); 
   } else if (strcmp(option, "probeset") == 0) {
      anntree = (TTree*)gDirectory->Get(chip->GetProbesetAnnotTree()); 
   } else {
      anntree = (TTree*)gDirectory->Get(chip->GetAnnotTree()); 
   }//if
   if (anntree == 0) return errGetTree;
   anntree->SetBranchAddress("AnnBranch", &annot);

   Int_t numunits = (Int_t)(unittree->GetEntries());
   Int_t numannot = (Int_t)(anntree->GetEntries());

// Create hash table to store unit names from anntree
   THashTable *htable = 0;
   if (!(htable = new THashTable(2*numannot))) return errInitMemory;

   if (XManager::fgVerbose) {
      cout << "Reading entries from <" << anntree->GetName() << "> ...";
   }//if
   TString str;
   XIdxString *idxstr = 0;
   for (Int_t i=0; i<numannot; i++) {
      anntree->GetEntry(i);

      if (strcmp(option, "exon") == 0) { 
         str.Form("%d", ((XExonAnnotation*)annot)->GetExonID());
         idxstr = new XIdxString(i, str.Data());
      } else if (strcmp(option, "probeset") == 0) {
         str.Form("%d", ((XProbesetAnnotation*)annot)->GetProbesetID());
         idxstr = new XIdxString(i, str.Data());
      } else {
         str = ((XGenomeAnnotation*)annot)->GetTranscriptID();
         idxstr = new XIdxString(i, str.Data());
      }//if

      htable->Add(idxstr);
   }//for_i
   if (XManager::fgVerbose) {
      cout << "Finished" << endl;
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit)      output << sep << "UNIT_NAME";
   if (hasAnnot) {
      if (hasName)   output << sep << "GENE_NAME";
      if (hasSymbol) output << sep << "GENE_SYMBOL";
      if (hasCyto)   output << sep << "CYTOBAND";
   }//if
   if (n == 1) {
      if (hasCall)   output << sep << "CALL";
      if (hasPVal)   output << sep << "PVALUE";
   } else {
      for (Int_t k=0; k<n; k++) {
         if (hasCall)  output << sep << (names[k] + "_CALL");
         if (hasPVal)  output << sep << (names[k] + "_PVALUE");
      }//for_k
   }//if
   output << endl;

// Loop over tree entries and trees
   Int_t idx = 0;
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);

         // export annotation
         if (k == 0) {
            Int_t unitID = call[k]->GetUnitID();
            output << unitID;

            unittree->GetEntry(idx++);
            while (unitID != unit->GetUnitID()) {
               if (idx == numunits) {
                 cerr << "Error: UnitID <" << unitID << "> not found." << endl;
                 err = errAbort; goto cleanup;
               }//if

               unittree->GetEntry(idx++);
            }//while

            if (hasUnit) {
               output << sep << unit->GetSubUnitID();
            }//if

            if (hasAnnot) {
               str.Form("%d", unit->GetSubUnitID());
               idxstr = (XIdxString*)(htable->FindObject(str.Data()));
               if (idxstr) {
                  anntree->GetEntry(idxstr->GetIndex());
                  if (hasName)   output << sep << annot->GetName();
                  if (hasSymbol) output << sep << annot->GetSymbol();
                  if (hasCyto)   output << sep << annot->GetCytoBand();
               } else {
                  if (hasName)   output << sep << "NA";
                  if (hasSymbol) output << sep << "NA";
                  if (hasCyto)   output << sep << "NA";
               }//if
            }//if
         }//if

         // export data
         if (hasCall) {
            Int_t cl = call[k]->GetCall();
            char *ch = "NA";
            if      (cl == 2) ch = "P";
            else if (cl == 0) ch = "A";
            else if (cl == 1) ch = "M";
            output << sep << ch;
         }//if

         if (hasPVal) {
            output << sep << call[k]->GetPValue();
         }//if
      }//for_k
      output << endl;
   }//for_i

//Cleanup
cleanup:
   // remove trees from RAM
   if (anntree)  {anntree->Delete("");  anntree  = 0;}
   if (unittree) {unittree->Delete(""); unittree = 0;}
   if (htable)   {htable->Delete(); delete htable; htable = 0;}
   SafeDelete(schemes);

   return err;
}//ExportCallTrees

//______________________________________________________________________________
Int_t *XExonProcesSet::FillMaskArray(XDNAChip *chip, TTree *scmtree, XScheme *scheme,
                       Int_t level, Int_t n, Int_t *msk)
{
   // "level=0" is cutoff to prevent PM to be redefined as MM for alternative CDFs
   // Fill array "msk" with mask from scheme tree
   if(kCS) cout << "------XExonProcesSet::FillMaskArray------" << endl;

//?   return msk = XGCProcesSet::FillMaskArray(chip, scmtree, scheme, level, n, msk);
   return msk = XGenomeProcesSet::FillMaskArray(chip, scmtree, scheme, level, n, msk);
}//FillMaskArray

