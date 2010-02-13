// File created: 08/05/2002                          last modified: 01/03/2010
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2010 Dr. Christian Stratowa
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
* Aug 2008 - Add summarization algorithms FARMS and DFW, classes XFARMS, XDFW
* Oct 2008 - Add I/NI-call algorithm class XINICall
* Jul 2009 - Allow to change bufsize of tree branches
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
      if (!err) err = InitAlgorithm("normalizer","quantile","transcript:together:none:0", 0, 1,0.0);
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
   //    - numpars: number of other parameters as integer, i.e. numpars = 1-5:
   //    - trim:    trim value for mean, in range [0, 0.5]
   //    - delta:   for robust ties, default 1.0 or 0.4
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
   // type = "MedianPolish": median polish (multichip algorithm), with parameters:
   //    parameters are: numpars, maxiter, eps, neglog, (nfrac, l, h)
   //    - numpars: number of other parameters as integer, i.e. numpars = 3-5:
   //    - maxiter: maximal number of iterations, default is 10
   //    - eps:     epsilon of test for convergence, default is 0.01
   //    - neglog:  substitution for logarithm of negative values
   //    - nfrac:   noise fraction for bgrd option "correctbg", or
   //    - l:       optional tunable parameter, 0<=l<=1 (default is 0.005), and     
   //    - h:       optional parameter (default is -1) for "attenuatebg"
   //
   // type = "FARMS": FARMS (multichip algorithm), with parameters:
   //    parameters are: numpars,version, weight, mu, scale, tol, cyc, weighted, neglog
   //    - numpars: number of other parameters as integer, i.e. numpars = 6-8:
   //    - version: version of farms package, 131 (farms_1.3.1); 130 (farms_1.3)
   //    - weight:  hyperparameter, default is 0.5
   //    - mu:      hyperparameter, default is 0.0
   //    - scale:   scaling parameter, default is 1.0
   //    - tol:     termination tolerance, default is 0.00001
   //    - cyc:     maximum number of cycles, default is 0.0
   //    - weighted:weighted mean (for 131 only)
   //    - neglog:  substitution for logarithm of negative values
   //
   // type = "DFW": DFW (multichip algorithm), with parameters:
   //    parameters are: numpars, m, n, c, neglog,
   //    - numpars: number of other parameters as integer, i.e. numpars = 3-4:
   //    - m:       exponent for range WR, default is 3
   //    - n:       exponent for stdev WSD, default is 1
   //    - c:       scale parameter, default is 0.01
   //    - neglog:  substitution for logarithm of negative values
   //
   //    filename = "": data for "multichip" algorithms will be stored as table in RAM  
   //             = "tmp": optional filename to create temporary file "tmp_exten",
   //               where data will be stored temporarily
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
   } else if (strcmp(exten.Data(), kExtenExpr[8]) == 0) {
      fExpressor = new XFARMS(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenExpr[9]) == 0) {
      fExpressor = new XDFW(stype.Data(), exten.Data());
   } else {
      cerr << "Error: Expressor <" << type << "> is not known." << endl;
      return errInitSetting;
   }//if
   if (fExpressor == 0) return errInitMemory;

   fExpressor->SetOptions(options);
   fExpressor->NewFile(filename, exten.Data());

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
   //
   // type = "INICall": Informative call of FARMS (multichip algorithm), with parameters:
   //    parameters are: numpars,version, weight, mu, scale, tol, cyc, alpha1, alpha2
   //    - numpars: number of other parameters as integer, i.e. numpars = 6-8:
   //    - version: version of farms package, 131 (farms_1.3.1); 130 (farms_1.3)
   //    - weight:  hyperparameter, default is 0.5
   //    - mu:      hyperparameter, default is 0.0
   //    - scale:   scaling parameter, default is 1.0
   //    - tol:     termination tolerance, default is 0.00001
   //    - cyc:     maximum number of cycles, default is 0.0
   //    - alpha1:  significance level between P and M call, default is 0.4
   //    - alpha2:  significance level between M and A call, default is 0.6
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
   } else if (strcmp(exten.Data(), kExtenCall[4]) == 0) {
      fCaller = new XINICall(stype.Data(), exten.Data());
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

// Initialize data trees and background trees
   TTree **datatree = new TTree*[numdata+1]; //for possible Reference tree for normalization
   TTree **bgrdtree = new TTree*[numdata];
   for (Int_t k=0; k<numdata; k++) datatree[k] = bgrdtree[k] = 0;

   err = this->InitTrees(numdata, datatree, numbgrd, bgrdtree);
   if (err != errNoErr) goto cleanup;

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

// Qualify

// Condense
   if (fExpressor && fExprSelector && (doExpr || doAll)) {
//? if () cerr << "Error: At least two trees need to be selected." << endl;
      err = this->Express(numdata, datatree, numbgrd, bgrdtree);
      if (err != errNoErr) goto cleanup;
   } else if (doExpr && !fExpressor) {
      cerr << "Error: Expressor algorithm is not initialized!" << endl;
      err = errAbort; goto cleanup;
   }//if

// Informing user
   if (XManager::fgVerbose) {
      cout << "   preprocessing finished." << endl;
   }//if

// Cleanup
cleanup:
   delete [] bgrdtree;
   delete [] datatree;

   SafeDelete(fData);
   SafeDelete(fSchemes);

   return err;
}//Preprocess

//______________________________________________________________________________
Int_t XGCProcesSet::AdjustBackground(Int_t numdata, TTree **datatree,
                    Int_t &numbgrd, TTree **bgrdtree)
{
   // Adjust background, i.e. calculate background and subtact from intensity.
   // All datatrees will be replaced with bg-corrected datatrees
   // Note: Parameter numbgrd will be set to numbgrd=0 to prevent further use
   //       of externally added bgrdtrees, since datatrees are already corrected!
   if(kCS) cout << "------XGCProcesSet::AdjustBackground------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "   Background correcting raw data..." << endl;
   }//if

   Int_t err   = errNoErr;
   Int_t split = 99;
   Int_t i, j, ij, x, y;

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
   Double_t *dummy    = 0;  //to prevent compilation error

// Get chip parameters from scheme file (also for alternative CDFs)
   if (datatree[0] == 0) return errGetTree;
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XDNAChip *chip = (XDNAChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) return fManager->HandleError(errGetScheme, fSchemeName);
   Int_t numrows = chip->GetNumRows();
   Int_t numcols = chip->GetNumColumns();
   Int_t size    = numrows*numcols;

// Get exon level of annotation
   Int_t level = 0;
   if (fBgrdSelector->GetNumParameters() > 0) {
      level = (Int_t)(fBgrdSelector->GetParameters())[0];
   }//if

// Initialize memory for mask array
   if (!(arrMask  = new (nothrow) Int_t[size])) return errInitMemory;
   //not =0,since msk=1 and msk=0 must be determined!
   for (i=0; i<size; i++) arrMask[i] = eINITMASK;

// Get mask for PM/MM from scheme tree and store in array 
   err = this->SchemeMask(chip, level, size, arrMask);
   if (err != errNoErr) goto cleanup;

// Calculate mask for background
   err = fBgrdSelector->Calculate(size, dummy, dummy, arrMask);
   if (err != errNoErr) goto cleanup;

// For GCBackground only, fill arrMask with GC content:
   if (strcmp(fBackgrounder->GetName(), "gccontent") == 0) {
      err = this->ProbeMask(chip, size, arrMask);
      if (err != errNoErr) goto cleanup;
   }//if

// Change directory
   if (tmpfile != 0) tmpfile->cd();
   else if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

// Initialize memory for data arrays
   if (!(arrInten = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
   if (!(arrStdev = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
   if (!(arrBgrd  = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
   if (!(arrNoise = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
   for (i=0; i<size; i++) {
      arrInten[i] = arrStdev[i] = arrBgrd[i] = arrNoise[i] = 0.0;
   }//for_i

// Calculate background
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) {err = errGetTree; break;}

   // Informing user
      if (XManager::fgVerbose) {
         cout << "      calculating background for <" << datatree[k]->GetName() << ">..."
              << endl;
      }//if

   // Init data tree
      XGCCell *gccell = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell);

   // Init counters for min/max values
      Double_t min    = DBL_MAX;  //defined in float.h
      Double_t max    = 0;
      Int_t    nummin = 0;
      Int_t    nummax = 0;

   // Create new tree bgrdtree
      TString dataname = Path2Name(datatree[k]->GetName(), dSEP, ".");
      TString bgrdname = dataname + "." + fBackgrounder->GetTitle();
      bgrdtree[k] = new TTree(bgrdname, fSchemeName);
      if (bgrdtree[k] == 0) return errCreateTree;
      XBgCell *bgcell  = new XBgCell();
      Int_t    bufsize = XManager::GetBufSize();
      bgrdtree[k]->Branch("BgrdBranch", "XBgCell", &bgcell, bufsize, split);

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
      if (tmpfile == 0) {
         AddTreeHeader(bgrdtree[k]->GetName(), "Bgrd", 0, fBackgrounder->GetNumParameters(),
                       fBackgrounder->GetParameters());
      }//if

   // Write background-corrected data tree to file 
      arrStdev = fBackgrounder->AdjustError(size, arrStdev, arrNoise);
      datatree[k] = this->FillDataTree(datatree[k], exten, fBackgrounder,
                          numrows, numcols, arrInten, arrStdev);
      if (datatree[k] == 0) {err = errCreateTree; goto cleanup;}

      // reset branches
      SafeDelete(bgcell);
      bgrdtree[k]->ResetBranchAddress(bgrdtree[k]->GetBranch("BgrdBranch"));

      SafeDelete(gccell);
      datatree[k]->DropBaskets();  //to remove baskets from memory
      datatree[k]->ResetBranchAddress(datatree[k]->GetBranch("DataBranch"));

      if (err != errNoErr) break;
   }//for_k

cleanup:
   // delete arrays
   if (arrNoise) {delete [] arrNoise; arrNoise = 0;}
   if (arrBgrd)  {delete [] arrBgrd;  arrBgrd  = 0;}
   if (arrStdev) {delete [] arrStdev; arrStdev = 0;}
   if (arrInten) {delete [] arrInten; arrInten = 0;}
   if (arrMask)  {delete [] arrMask;  arrMask  = 0;}
//?   // delete scheme tree from RAM
//?   if (scmtree)  {scmtree->Delete(""); scmtree = 0;}
   // note: keep datatree[k] and bgrdtree[k] for later use!

//not allowed   SafeDelete(chip);

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

// Informing user
   if (XManager::fgVerbose) {
      cout << "   Normalizing raw data..." << endl;
   }//if

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
//TO DO: case of pars[npars-1] == -100 in BackgroundParameters()
   } else {
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Initialize reference trees
   Int_t  numsels = fSelections ? fSelections->GetSize() : 0;
   Int_t  numrefs = fReferences ? fReferences->GetSize() : 0;
   TTree **reftree = new TTree*[numrefs];
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
   Double_t *dummy   = 0;  //to prevent compilation error

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
      cout << "      normalizing data using method <" << fNormalizer->GetName() << ">..." << endl;
   }//if

// Compute reference arrays
   if (strcmp(fNormalizer->GetName(), "quantile") == 0) {

      // get mask for data to be used for normalization
      err = fNormSelector->Calculate(size, dummy, dummy, arrMask);
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
            cout << "         filling array " << k+1 << " of " << numdata
                 << ": <" << treename.Data() << ">...        \r" << flush;
         }//if

         err = this->FillDataArrays(datatree[k], bgrdtree[k], doBg,
                     numrows, numcols, arrIntx, 0, 0);
         if (err != errNoErr) goto cleanup;

         err = fNormalizer->AddArray(size, arrIntx, arrMask, treename);
         if (err != errNoErr) goto cleanup;
      }//for_k
      if (XManager::fgVerbose) {
         cout << "         finished filling <" << numdata << "> arrays."
              << "                   " << endl;
      }//if
   } else if (numrefs == 1) {
//////////
//TO DO: if reftree is not one of datatrees, ev. from different root file
//       datatree[refid] not possible, need to use reftree[0] (w/o bgrd??)
//////////
      // informing user
      if (XManager::fgVerbose) {
         cout << "         filling array <" << datatree[refid]->GetName() << ">..." << endl;
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
         cout << "         filling array <" << kReference << ">..." << endl;
      }//if

      TTree **rbgtree = new TTree*[numrefs];
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
            cerr << "Error: <" << (numrefs - numrbgs) 
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

      delete [] rbgtree;

      numsels++;
   } else {
         cerr << "Error: No reference tree is selected." << endl;
         err = errGeneral;
         goto cleanup;
   }//if

// Normalization
   if (strcmp(fNormalizer->GetName(), "quantile") == 0) {
      // quantile normalization
      if (XManager::fgVerbose) {
         cout << "         computing common mean..." << endl;
      }//if
      if ((err = fNormalizer->Calculate(size, arrIntx, arrInty, arrMask))) goto cleanup;

      // get normalized arrays and fill data trees
      for (Int_t k=0; k<numdata; k++) {
         treename = Path2Name(datatree[k]->GetName(),"",".");
         TString exten = fNormalizer->GetTitle();

         // informing user
         if (XManager::fgVerbose) {
            cout << "         filling tree " << k+1 << " of " << numdata
                 << ": <" << (treename + "." + exten).Data() << ">..."
                 << "              \r" << flush;
         }//if

         arrInty = fNormalizer->GetArray(size, arrInty, arrMask, treename);
         if (arrInty == 0) {err = errAbort; goto cleanup;}

         datatree[k] = this->FillDataTree(datatree[k], exten, fNormalizer,
                             numrows, numcols, arrInty, 0);
         if (datatree[k] == 0) {err = errCreateTree; goto cleanup;}
      }//for_k

      if (XManager::fgVerbose) {
         cout << "         finished filling <" << numdata << "> trees."
              << "                " << endl;
      }//if
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

   delete [] reftree;
   // delete scheme tree from RAM
//   if (scmtree) {scmtree->Delete(""); scmtree = 0;}

   savedir->cd();

   return err;
}//Normalize

//______________________________________________________________________________
Int_t XGCProcesSet::DetectCall(Int_t numdata, TTree **datatree,
                    Int_t &numbgrd, TTree **bgrdtree)
{
   // Detect presence call
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCS) cout << "------XGCProcesSet::DetectCall------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "   Calculating detection calls..." << endl;
   }//if

   Int_t err = errNoErr;

   if (fCaller->IsMultichip()) {
      err = this->DoMultichipCall(numdata, datatree, numbgrd, bgrdtree,
                                  fCaller->GetFile());
   } else {
      err = this->DoCall(numdata, datatree, numbgrd, bgrdtree);
   }//if

   return err;
}//DetectCall

//______________________________________________________________________________
Int_t XGCProcesSet::DoCall(Int_t numdata, TTree **datatree,
                    Int_t &numbgrd, TTree **bgrdtree)
{
   // Detect presence call
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCS) cout << "------XGCProcesSet::DoCall------" << endl;

   Int_t err   = errNoErr;
   Int_t level = 0;
   Int_t idx   = 0;
   Int_t ij, start, end;

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
   Double_t *dummy    = 0;  //to prevent compilation error

   TTree   *scmtree = 0; 
   TLeaf   *scmleaf = 0;
   XScheme *scheme  = 0;

   TTree   *idxtree = 0; 
   XGCUnit *unit    = 0;

   TTree  *calltree = 0;
   XPCall *call     = 0;
   Int_t   bufsize  = XManager::GetBufSize(numdata, 10000);
   Int_t   split    = 99;

   fFile->cd();

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fCaller, fCaller->GetBgrdOption());
   Bool_t doGC = (Bool_t)(strcmp(fCaller->GetName(), "dabgcall") == 0);

// Check for presence of background trees
   if (numbgrd == 0) {
      if (doBg == kTRUE) {
         cout << "Warning: No background trees available for background subtraction."
              << endl;
      }//if
//TO DO: case of pars[npars-1] == -100 in BackgroundParameters()
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Get chip parameters from scheme file (also alternative CDFs)
   if (datatree[0] == 0) return errGetTree;
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if

// Change directory to scheme file
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XGeneChip *chip = (XGeneChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) return fManager->HandleError(errGetScheme, fSchemeName);
   Int_t numrows  = chip->GetNumRows();
   Int_t numcols  = chip->GetNumColumns();
   Int_t numunits = chip->GetNumUnits();
   Int_t size     = numrows*numcols;

// Get scheme tree  and unit tree for scheme
   idxtree = this->UnitTree(fCaller, &unit, numunits); 
   if (idxtree == 0) return errGetTree;

   scmtree = this->SchemeTree(fCaller, &scheme, &scmleaf); 
   if (scmtree == 0) return errGetTree;

// Get maximum number of pairs/cells
   Int_t maxnumcells = this->MaxNumberCells(idxtree);
   if (maxnumcells <= 0) return errGeneral;

// Init unit selector
   XUnitSelector *unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
   unitSelector->SetOption((this->ChipType(fTitle)).Data());
   err = unitSelector->InitParameters(fCallSelector->GetNumParameters(),
                                      fCallSelector->GetParameters());
   if (err != errNoErr) goto cleanup;

// Initialize memory for unit arrays
   if (!(arrUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}
   if (!(mskUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 

// Get mask from scheme tree and store in array 
   arrUnit = this->FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);
   if (arrUnit == 0) {err = errInitMemory; goto cleanup;}

// Calculate units satisfying mask
   err = unitSelector->Calculate(numunits, arrUnit, mskUnit);
   if (err != errNoErr) goto cleanup;

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

// Get exon level of annotation
   if (fCallSelector->GetNumParameters() > 0) {
      level = (Int_t)(fCallSelector->GetParameters())[0];
   }//if

// Get mask for PM/MM from scheme tree and store in array 
   arrMask = this->FillMaskArray(chip, scmtree, scheme, level, size, arrMask);
   if (arrMask == 0) {err = errInitMemory; goto cleanup;}

// Calculate mask for detection call (set arrMask to 1 or 0)
   err = fCallSelector->Calculate(size, dummy, dummy, arrMask);
   if (err != errNoErr) goto cleanup;

// For DABG only, fill arrMask with GC content (GC>=0 for PM)
   if (doGC) if ((err = this->MaskArray2GC(chip, arrMask))) goto cleanup;

// Initialize maximum memory for PM/MM arrays (maxnumcells+1 to avoid potential buffer overflow)
   if (!(arrPM = new (nothrow) Double_t[maxnumcells+1])) {err = errInitMemory; goto cleanup;}
   if (!(arrMM = new (nothrow) Double_t[maxnumcells+1])) {err = errInitMemory; goto cleanup;}
   if (!(arrSP = new (nothrow) Double_t[maxnumcells+1])) {err = errInitMemory; goto cleanup;}
   if (!(arrSM = new (nothrow) Double_t[maxnumcells+1])) {err = errInitMemory; goto cleanup;}
   if (!(arrXP = new (nothrow) Int_t[maxnumcells+1]))    {err = errInitMemory; goto cleanup;}
   if (!(arrXM = new (nothrow) Int_t[maxnumcells+1]))    {err = errInitMemory; goto cleanup;}

// Change directory to data file
   if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

// Calculate detection call
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) {err = errGetTree; goto cleanup;}

   // Informing user
      TString name = datatree[k]->GetName();
      if (XManager::fgVerbose) {
         cout << "      calculating present call for <" << name.Data() << ">..." << endl;
      }//if

   // Get tree info for datatree name.exten
      XDataTreeInfo *info = 0;
      info = (XDataTreeInfo*)datatree[k]->GetUserInfo()->FindObject(name);
      if (!info) {
         cerr << "Error: Could not get tree info for <" << name.Data() << ">." << endl;
         return errGeneral;
      }//if

      Int_t numabsent  = 0;
      Int_t numarginal = 0;
      Int_t numpresent = 0;
      Double_t minpval = 1.0;
      Double_t maxpval = 0.0;

   // Get data from datatree and bgrdtree and store in arrays
      err = this->FillDataArrays(datatree[k], bgrdtree[k], doBg, numrows, numcols,
                                 arrInten, arrStdev, arrNPix, arrBgrd, arrBgdev);
      if (err != errNoErr) goto cleanup;

   // For DABG only get intensities and save in table sorted for GC-content
      if (doGC) if ((err = fCaller->Calculate(size, arrInten, arrMask))) goto cleanup;

   // Create new tree calltree
      name = Path2Name(datatree[k]->GetName(), dSEP, ".") + "." + fCaller->GetTitle();
      calltree = new TTree(name, fSchemeName);
      if (calltree == 0) {err = errCreateTree; goto cleanup;}
      call = new XPCall();
      calltree->Branch("CallBranch", "XPCall", &call, bufsize, split);

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

            ij = XY2Index(scheme->GetX(), scheme->GetY(), numcols);

            if (doGC) {
               this->FillProbeSets(p, idx,
                                   arrPM, arrSP, arrXP,
                                   arrMask[ij], arrInten[ij], arrStdev[ij]);
            } else {
               this->FillProbeSets(p, m, idx,
                                   arrPM, arrMM, arrSP, arrSM, arrXP, arrXM, 
                                   arrMask[ij], arrInten[ij], arrStdev[ij],
                                   arrBgrd[ij], arrBgdev[ij], arrNPix[ij]);
            }//if

            if (p > maxnumcells || m > maxnumcells) {
               cerr << "Error: unitID <" << unitID << "> exceeds maximum number of cells <"
                    << maxnumcells << ">. Buffer overflow!" << endl;
               err = errAbort;
               goto cleanup;
            }//if
         }//for_j
//         start += numcells;

         // continue if arrays arrPM etc are not filled
         if (doGC || p == m) {
            start += numcells;
         } else if (p == 0) {
            start += numcells;
            continue;
         } else if (!doGC && m == 0) {
            // for m=0 refill arrMM with background values
            for (Int_t j=start; j<end; j++) {
               scmtree->GetEntry(j);
               ij = XY2Index(scheme->GetX(), scheme->GetY(), numcols);
               this->FillBgrdProbeSets(m, arrMM, arrSM, arrXM,
                                       arrMask[ij], arrBgrd[ij], arrBgdev[ij], arrNPix[ij]);
            }//for_j
            start += numcells;
         } else if (!doGC && p != m) {
            cerr << "Error: UnitID <" << unitID << "> has different numbers of PM <"
                 << p << "> and MM <" << m << "> data." << endl;
            start += numcells;
            continue;
//x            err = errAbort;
//x            goto cleanup;
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

         if (XManager::fgVerbose && idx%1000 == 0) {
            cout << "      <" << idx << "> of <" << numunits << "> calls processed...\r" << flush;
         }//if
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
                      idx, numabsent, numarginal, numpresent, minpval, maxpval);

   // Write call tree to file 
      if ((err = WriteTree(calltree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(calltree->GetName(), "Call", 0, fCaller->GetNumParameters(),
                       fCaller->GetParameters());
      }//if

      SafeDelete(call);
      calltree->ResetBranchAddress(calltree->GetBranch("CallBranch"));
//?? calltree needed later?
      SafeDelete(calltree);

      if (err != errNoErr) break;
   }//for_k

// Cleanup
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

   SafeDelete(unitSelector);
   SafeDelete(unit);
   idxtree->ResetBranchAddress(idxtree->GetBranch("IdxBranch"));
   SafeDelete(scheme);
   scmtree->ResetBranchAddress(scmtree->GetBranch("ScmBranch"));
//?   // delete scheme tree from RAM
//?   SafeDelete(scmtree);

   return err;
}//DoCall

//______________________________________________________________________________
Int_t XGCProcesSet::DoMultichipCall(Int_t numdata, TTree **datatree,
                    Int_t &numbgrd, TTree **bgrdtree, TFile *file)
{
   // Detect presence call
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCS) cout << "------XGCProcesSet::DoMultichipCall------" << endl;

   Int_t err = errNoErr;
   Int_t idx = 0;
   Int_t x   = 0;
   Int_t y   = 0;
   Int_t ij, start, end, id, entry;

   Int_t numsels = 0;  //number of selected entries
   Int_t exlevel = 0;
   Int_t stepout = (Int_t)((100000.0 + 10.0*numdata)/(Float_t)numdata); //step size for verbose output

// Init min/max p-values
   Int_t numabsent  = 0;
   Int_t numarginal = 0;
   Int_t numpresent = 0;
   Double_t minpval = 1.0;
   Double_t maxpval = 0.0;

// Init 
   Int_t     *arrMask  = 0;
   Int_t     *arrIndx  = 0;
   Int_t     *arrUnit  = 0;
   Int_t     *mskUnit  = 0;
   Double_t  *prescall = 0;
   Double_t  *pvalue   = 0;

   Double_t **table   = 0;
   Double_t  *arrData = 0;
   Double_t  *arrPM   = 0; 
   Double_t  *dummy   = 0;  //to prevent compilation error

   TTree   *scmtree = 0; 
   TLeaf   *scmleaf = 0;
   XScheme *scheme  = 0;

   TTree   *idxtree = 0; 
   XGCUnit *unit    = 0;

   XGCCell **gccell = 0;
   XBgCell **bgcell = 0;

   TTree  **tmptree  = 0;
   TTree  **calltree = 0;
   XPCall **call     = 0;
   Int_t    bufsize  = XManager::GetBufSize(numdata, 10000);
   Int_t    split    = 99;
   Double_t sort;

   fFile->cd();

// Get chip parameters from scheme file (also alternative CDFs)
   if (datatree[0] == 0) return errGetTree;
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if

// Change directory to scheme file
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XGeneChip *chip = (XGeneChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) return fManager->HandleError(errGetScheme, fSchemeName);
   Int_t numrows  = chip->GetNumRows();
   Int_t numcols  = chip->GetNumColumns();
   Int_t numunits = chip->GetNumUnits();
   Int_t size     = numrows*numcols;

// Get scheme tree  and unit tree for scheme
   idxtree = this->UnitTree(fCaller, &unit, numunits); 
   if (idxtree == 0) return errGetTree;

   scmtree = this->SchemeTree(fCaller, &scheme, &scmleaf); 
   if (scmtree == 0) return errGetTree;

// Get maximum number of pairs/cells
   Int_t maxnumcells = this->MaxNumberCells(idxtree);
   if (maxnumcells <= 0) return errGeneral;

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fCaller, fCaller->GetBgrdOption());

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
//?   } else if (subbgrd || corbgrd) {
   } else if (doBg == kTRUE) {
      cout << "Warning: No background trees available for background subtraction."
           << endl;
      doBg = kFALSE;
//TO DO: case of pars[npars-1] == -100 in BackgroundParameters()
   } else {
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Init unit selector
   XUnitSelector *unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
   unitSelector->SetOption((this->ChipType(fTitle)).Data());
   err = unitSelector->InitParameters(fCallSelector->GetNumParameters(),
                                      fCallSelector->GetParameters());
   if (err != errNoErr) goto cleanup;

// Initialize memory for unit arrays
   if (!(arrUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}
   if (!(mskUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 

// Get mask from scheme tree and store in array 
   arrUnit = this->FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);
   if (arrUnit == 0) {err = errInitMemory; goto cleanup;}

// Calculate units satisfying mask
   err = unitSelector->Calculate(numunits, arrUnit, mskUnit); //arrUnit to prevent compile error
   if (err != errNoErr) goto cleanup;

// Create local arrays
   if (!(arrMask  = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrIndx  = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(prescall = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}
   if (!(pvalue   = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<size; i++) {
      arrMask[i] = eINITMASK; 
      arrIndx[i] = 0;
   }//for_i
   for (Int_t i=0; i<numdata; i++) {
      prescall[i] = 0.0; 
      pvalue[i]   = 1.0; 
   }//for_i

// Get exon level of annotation
   if (fCallSelector->GetNumParameters() > 0) {
      exlevel = (Int_t)(fCallSelector->GetParameters())[0];
   }//if

// Get mask for PM from scheme tree and store in array 
   arrMask = this->FillMaskArray(chip, scmtree, scheme, exlevel, size, arrMask);
   if (arrMask == 0) {err = errInitMemory; goto cleanup;}

// Calculate mask for expression
   err = fCallSelector->Calculate(size, dummy, dummy, arrMask);
   if (err != errNoErr) goto cleanup;

// Init branch addresses
   gccell = new XGCCell*[numdata];
   bgcell = new XBgCell*[numdata];
   for (Int_t k=0; k<numdata; k++) {
      gccell[k] = 0;
      bgcell[k] = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell[k]);
      if (numbgrd > 0) bgrdtree[k]->SetBranchAddress("BgrdBranch", &bgcell[k]);
   }//for_k

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
//DB   datatree[0]->DropBaskets();  //to remove baskets from memory

// Get number of selected entries
   for (Int_t i=0; i<size; i++) {
      numsels = (arrMask[i] == 1) ? ++numsels : numsels;
   }//for_i

// Get data from datatrees (and bgrdtrees) and store in table or array
   if (file == 0) {
      // create table
      if (!(table = CreateTable(numdata, numsels))) {err = errInitMemory; goto cleanup;}

      // fill table with selected intensities of all datatrees
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
   } else {
      // create array to store selected intensities
      if (!(arrData = new (nothrow) Double_t[numsels])) {err = errInitMemory; goto cleanup;}
      for (Int_t i=0; i<numsels;  i++) arrData[i] = 0.0;

      // change directory to temporary file
      if (!file->cd()) {err = errGetDir; goto cleanup;}

      // get data from datatrees (and bgrdtrees) and store in temporary file
      tmptree  = new TTree*[numdata];
      for (Int_t k=0; k<numdata; k++) {
         // create temporary tree
         tmptree[k] = new TTree(datatree[k]->GetName(), "temporary tree");
         if (tmptree[k] == 0) {err = errCreateTree; goto cleanup;}
         tmptree[k]->Branch("sort", &sort, "sort/D", bufsize);

         // informing user
         if (XManager::fgVerbose) {
            cout << "         filling temporary tree <" << tmptree[k]->GetName()
                 << ">...              \r" << flush;
         }//if

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

            datatree[k]->DropBaskets();  //to remove baskets from memory
            bgrdtree[k]->DropBaskets();  //to remove baskets from memory
         } else {
            idx = 0;
            for (Int_t i=0; i<size; i++) {
               datatree[k]->GetEntry(i);

               if (arrMask[i] == 1) {
                  arrData[idx++] = gccell[k]->GetIntensity();
               }//if
            }//for_i

            datatree[k]->DropBaskets();  //to remove baskets from memory
         }//if

         // fill tmptree with array in the order of scheme tree entries for (x,y)
         start = 0;
         end   = 0;
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
                  sort = arrData[arrIndx[ij]];
                  tmptree[k]->Fill();
               }//if
            }//for_j
            start += numcells;
         }//for_id

         // write tmptree to temporary file
         tmptree[k]->Write();
//??         tmptree[k]->Write(TObject::kOverwrite);
         tmptree[k]->DropBaskets();  //to remove baskets from memory
      }//for_k

      if (XManager::fgVerbose) {
         cout << "         finished filling <" << numdata << "> temporary trees.          " << endl;
      }//if
   }//if

// Change directory to current directory for treeset in main file
   if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

// Create new trees calltree
   calltree = new TTree*[numdata];
   call     = new XPCall*[numdata];

   for (Int_t k=0; k<numdata; k++) {
      TString dataname = Path2Name(datatree[k]->GetName(), dSEP, ".");
      TString callname = dataname + "." + fCaller->GetTitle();
      calltree[k] = new TTree(callname, fSchemeName);
      if (calltree[k] == 0) {err = errCreateTree; goto cleanup;}

      call[k] = new XPCall();
      calltree[k]->Branch("CallBranch", "XPCall", &call[k], bufsize, split);

      // to reduce number of baskets in memory when reading trees
      if (file) tmptree[k]->SetMaxVirtualSize(bufsize);
   }//for_k

// Create array to store PM values for all probes with current unitID
   if (!(arrPM = new (nothrow) Double_t[(maxnumcells+1)*numdata])) {
      err = errInitMemory; goto cleanup;
   }//if

// Calculate detection call values
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

            if (file == 0) {
               for (Int_t k=0; k<numdata; k++) {
                  arrPM[p] = table[k][arrIndx[ij]];
                  p++;
               }//for_k
            } else {
               for (Int_t k=0; k<numdata; k++) {
                  tmptree[k]->GetEntry(entry);
                  arrPM[p] = sort;
                  p++;
               }//for_k
            }//if

            entry++;
         }//if
      }//for_j
      start += numcells;

      // fill arrPM or continue if it is not filled
      if ((err = fCaller->SetArray(p, arrPM)) != errNoErr) continue;

      // calculate detection call for PMs of current unitID
      if ((err = fCaller->Calculate(numdata, prescall, pvalue, 0))) break;

      // fill call trees
      for (Int_t k=0; k<numdata; k++) {
         // number of present/absent calls
         if      (prescall[k] == 2.0) numpresent++;
         else if (prescall[k] == 1.0) numarginal++;
         else if (prescall[k] == 0.0) numabsent++;

         // minimal/maximal detection call p-value
         if      (pvalue[k] < minpval) minpval = pvalue[k];
         else if (pvalue[k] > maxpval) maxpval = pvalue[k];

         call[k]->SetUnitID(unitID);
         call[k]->SetCall((Short_t)prescall[k]);
         call[k]->SetPValue(pvalue[k]);
         calltree[k]->Fill();
      }//for_k

      if (XManager::fgVerbose && (idx == 1 || id%stepout == 0)) {
         cout << "      calculating detection call for <" << idx << "> of <"
              << numunits << "> units...\r" << flush;
      }//if
   }//for_id

   if (XManager::fgVerbose) {
      cout << "      calculating detection call for <" << idx << "> of <"
           << numunits << "> units...Finished." << endl;
   }//if

   if (XManager::fgVerbose) {
      cout << "      detection call statistics: " << endl;
      cout << "         minimum detection p-value = " << minpval << endl;
      cout << "         maximum detection p-value = " << maxpval << endl;
      cout << "         P: <" << numpresent*100.0/idx << "> percent units present."  << endl;
      cout << "         M: <" << numarginal*100.0/idx << "> percent units marginal." << endl;
      cout << "         A: <" << numabsent*100.0/idx  << "> percent units absent."   << endl;
   }//if

// Write call trees to file 
   for (Int_t k=0; k<numdata; k++) {
   // Add tree info to tree
      AddCallTreeInfo(calltree[k], calltree[k]->GetName(), fCaller->GetOption(),
                      idx, numabsent, numarginal, numpresent, minpval, maxpval);
//need to do:
//                      idx, numabsent[k], numarginal[k], numpresent[k], minpval[k], maxpval[k]);

      if ((err = WriteTree(calltree[k], TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(calltree[k]->GetName(), "Call", 0, fCaller->GetNumParameters(),
                       fCaller->GetParameters());
      } else {
         break;
      }//if
   }//for_k

// Cleanup
cleanup:
   // delete table
   DeleteTable(table, numdata);

   // delete temporary trees
   if (tmptree) {
      for (Int_t k=0; k<numdata; k++) {
         tmptree[k]->Delete(""); tmptree[k] = 0;
      }//for_k
      delete [] tmptree;
   }//if
   // delete arrays
   if (arrPM)    {delete [] arrPM;    arrPM    = 0;}
   if (arrData)  {delete [] arrData;  arrData  = 0;}
   if (pvalue)   {delete [] pvalue;   pvalue   = 0;}
   if (prescall) {delete [] prescall; prescall = 0;}
   if (arrIndx)  {delete [] arrIndx;  arrIndx  = 0;}
   if (arrMask)  {delete [] arrMask;  arrMask  = 0;}
   if (mskUnit)  {delete [] mskUnit;  mskUnit  = 0;}
   if (arrUnit)  {delete [] arrUnit;  arrUnit  = 0;}

   for (Int_t k=0; k<numdata; k++) {
      SafeDelete(gccell[k]);
      datatree[k]->DropBaskets();  //to remove baskets from memory
      datatree[k]->ResetBranchAddress(datatree[k]->GetBranch("DataBranch"));

      if (numbgrd > 0) {
         SafeDelete(bgcell[k]);
         bgrdtree[k]->DropBaskets();  //to remove baskets from memory
         bgrdtree[k]->ResetBranchAddress(bgrdtree[k]->GetBranch("BgrdBranch"));
      }//if

      SafeDelete(call[k]);
      calltree[k]->ResetBranchAddress(calltree[k]->GetBranch("CallBranch"));
      SafeDelete(calltree[k]);
   }//for_k

   delete [] gccell;
   delete [] bgcell;
   delete [] call;
   delete [] calltree;

   SafeDelete(unitSelector);
   SafeDelete(unit);
   idxtree->ResetBranchAddress(idxtree->GetBranch("IdxBranch"));
   SafeDelete(scheme);
   scmtree->ResetBranchAddress(scmtree->GetBranch("ScmBranch"));
//?   // delete scheme tree from RAM
//?   SafeDelete(scmtree);

   return err;
}//DoMultichipCall

//______________________________________________________________________________
Int_t XGCProcesSet::Express(Int_t numdata, TTree **datatree,
                    Int_t &numbgrd, TTree **bgrdtree)
{
   // Convert intensity data to expression levels
   if(kCS) cout << "------XGCProcesSet::Express------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "   Converting raw data to expression levels..." << endl;
      cout << "      summarizing with <" << fExpressor->GetName() << ">..." << endl;
   }//if

   Int_t err = errNoErr;

   if (fExpressor->IsMultichip()) {
      err = this->DoMultichipExpress(numdata, datatree, numbgrd, bgrdtree,
                                     fExpressor->GetFile());
   } else {
      err = this->DoExpress(numdata, datatree, numbgrd, bgrdtree);
   }//if

   return err;
}//Express

//______________________________________________________________________________
Int_t XGCProcesSet::DoExpress(Int_t numdata, TTree **datatree,
                    Int_t numbgrd, TTree **bgrdtree)
{
   // Condense data from  datatrees to expression values
   if(kCS) cout << "------XGCProcesSet::DoExpress------" << endl;

   Int_t err   = errNoErr;
   Int_t level = 0;
   Int_t idx   = 0;
   Int_t ij, start, end;

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
   Double_t *dummy    = 0;  //to prevent compilation error

   TTree   *scmtree = 0; 
   TLeaf   *scmleaf = 0;
   XScheme *scheme  = 0;

   TTree   *idxtree = 0; 
   XGCUnit *unit    = 0;

   TTree  *exprtree    = 0;
   XGCExpression *expr = 0;
   Int_t   bufsize     = XManager::GetBufSize(numdata, 10000);
   Int_t   split       = 99;

   fFile->cd();

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fExpressor, fExpressor->GetBgrdOption());

// Check for presence of background trees
   if (numbgrd == 0) {
      if (doBg == kTRUE) {
         cout << "Warning: No background trees available for background subtraction."
              << endl;
      }//if
//TO DO: case of pars[npars-1] == -100 in BackgroundParameters()
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Get chip parameters from scheme file (also alternative CDFs)
   if (datatree[0] == 0) return errGetTree;
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if

// Change directory to scheme file
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XGeneChip *chip = (XGeneChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) return fManager->HandleError(errGetScheme, fSchemeName);
   Int_t numrows  = chip->GetNumRows();
   Int_t numcols  = chip->GetNumColumns();
   Int_t numunits = chip->GetNumUnits();
   Int_t size     = numrows*numcols;

// Get scheme tree  and unit tree for scheme
   idxtree = this->UnitTree(fExpressor, &unit, numunits); 
   if (idxtree == 0) return errGetTree;

   scmtree = this->SchemeTree(fExpressor, &scheme, &scmleaf); 
   if (scmtree == 0) return errGetTree;

// Get maximum number of pairs/cells
   Int_t maxnumcells = this->MaxNumberCells(idxtree);
   if (maxnumcells <= 0) return errGeneral;

// Init unit selector
   XUnitSelector *unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
   unitSelector->SetOption((this->ChipType(fTitle)).Data());
   err = unitSelector->InitParameters(fExprSelector->GetNumParameters(),
                                      fExprSelector->GetParameters());
   if (err != errNoErr) goto cleanup;

// Initialize memory for unit arrays
   if (!(arrUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}
   if (!(mskUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 

// Get mask from scheme tree and store in array 
   arrUnit = this->FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);
   if (arrUnit == 0) {err = errInitMemory; goto cleanup;}

// Calculate units satisfying mask
   err = unitSelector->Calculate(numunits, arrUnit, mskUnit);
   if (err != errNoErr) goto cleanup;

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

// Get exon level of annotation
   if (fExprSelector->GetNumParameters() > 0) {
      level = (Int_t)(fExprSelector->GetParameters())[0];
   }//if

// Get mask for PM/MM from scheme tree and store in array 
   arrMask = this->FillMaskArray(chip, scmtree, scheme, level, size, arrMask);
   if (arrMask == 0) {err = errInitMemory; goto cleanup;}

// Calculate mask for expression
   err = fExprSelector->Calculate(size, dummy, dummy, arrMask);
   if (err != errNoErr) goto cleanup;

// Initialize maximum memory for PM/MM arrays (maxnumpairs+1 to avoid potential buffer overflow)
   if (!(arrPM = new (nothrow) Double_t[maxnumcells+1])) {err = errInitMemory; goto cleanup;}
   if (!(arrMM = new (nothrow) Double_t[maxnumcells+1])) {err = errInitMemory; goto cleanup;}
   if (!(arrSP = new (nothrow) Double_t[maxnumcells+1])) {err = errInitMemory; goto cleanup;}
   if (!(arrSM = new (nothrow) Double_t[maxnumcells+1])) {err = errInitMemory; goto cleanup;}
   if (!(arrXP = new (nothrow) Int_t[maxnumcells+1]))    {err = errInitMemory; goto cleanup;}
   if (!(arrXM = new (nothrow) Int_t[maxnumcells+1]))    {err = errInitMemory; goto cleanup;}

// Change directory to data file
   if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

// Calculate expression
   for (Int_t k=0; k<numdata; k++) {
      if (datatree[k] == 0) {err = errGetTree; goto cleanup;}

   // Informing user
      TString name = datatree[k]->GetName();
      if (XManager::fgVerbose) {
         cout << "      summarizing <" << name.Data() << ">..." << endl;
      }//if

   // Get tree info for datatree name.exten
      XDataTreeInfo *info = 0;
      info = (XDataTreeInfo*)datatree[k]->GetUserInfo()->FindObject(name);
      if (!info) {
         cerr << "Error: Could not get tree info for <" << name.Data() << ">." << endl;
         err = errGeneral;
         break;
      }//if

   // Init min/max expression levels
      Double_t min = DBL_MAX;  //defined in float.h
      Double_t max = 0;

   // Get data from datatree and bgrdtree and store in arrays
      err = this->FillDataArrays(datatree[k], bgrdtree[k], doBg, numrows, numcols,
                                 arrInten, arrStdev, arrNPix, arrBgrd, arrBgdev);
      if (err != errNoErr) goto cleanup;

   // Create new tree exprtree
      name = Path2Name(datatree[k]->GetName(), dSEP, ".") + "." + fExpressor->GetTitle();
      exprtree = new TTree(name, fSchemeName);
      if (exprtree == 0) {err = errCreateTree; goto cleanup;}
      expr = new XGCExpression();
      exprtree->Branch("ExprBranch", "XGCExpression", &expr, bufsize, split);

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

            ij = XY2Index(scheme->GetX(), scheme->GetY(), numcols);

            this->FillProbeSets(p, m, idx,
                                arrPM, arrMM, arrSP, arrSM, arrXP, arrXM, 
                                arrMask[ij], arrInten[ij], arrStdev[ij],
                                arrBgrd[ij], arrBgdev[ij], arrNPix[ij]);

            if (p > maxnumcells || m > maxnumcells) {
               cerr << "Error: unitID <" << unitID << "> exceeds maximum number of cells <"
                    << maxnumcells << ">.  (p/m = " << p << "/" << m << ")  Buffer overflow!"
                    << endl;
               err = errAbort;
               goto cleanup;
            }//if
         }//for_j
//x         start += numcells;

         // continue if arrays arrPM etc are not filled
         if (p == m) {
            start += numcells;
         } else if (p == 0) {
            start += numcells;
            continue;
         } else if (m == 0) {
            // for m=0 refill arrMM with background values
            for (Int_t j=start; j<end; j++) {
               scmtree->GetEntry(j);
               ij = XY2Index(scheme->GetX(), scheme->GetY(), numcols);
               this->FillBgrdProbeSets(m, arrMM, arrSM, arrXM,
                                       arrMask[ij], arrBgrd[ij], arrBgdev[ij], arrNPix[ij]);
            }//for_j
            start += numcells;
         } else if (p != m) {
            cerr << "Error: UnitID <" << unitID << "> has different numbers of PM <"
                 << p << "> and MM <" << m << "> data." << endl;
            start += numcells;
            continue;
//x            err = errAbort;
//x            goto cleanup;
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

         if (XManager::fgVerbose && idx%1000 == 0) {
            cout << "      <" << idx << "> of <" << numunits << "> calls processed...\r" << flush;
         }//if
      }//for_id

      if (XManager::fgVerbose) {
         cout << "      <" << idx << "> of <" << numunits << "> calls processed...Finished" << endl;
      }//if

      if (XManager::fgVerbose) {
         cout << "      expression statistics: " << endl;
         cout << "         minimal expression level is <" << min << ">." << endl;
         cout << "         maximal expression level is <" << max << ">." << endl;
      }//if

   // Add tree info to tree
      AddExprTreeInfo(exprtree, exprtree->GetName(), fExpressor->GetOption(), idx, min, max);

   // Write expression tree to file 
      if ((err = WriteTree(exprtree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(exprtree->GetName(), "Expr", 0, fExpressor->GetNumParameters(),
                       fExpressor->GetParameters());
      }//if

// ??? do not remove exprtree and expr from RAM, needed later
      SafeDelete(expr);
      exprtree->ResetBranchAddress(exprtree->GetBranch("ExprBranch"));
      SafeDelete(exprtree);

      if (err != errNoErr) break;
   }//for_k

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

   SafeDelete(unitSelector);
   SafeDelete(unit);
   idxtree->ResetBranchAddress(idxtree->GetBranch("IdxBranch"));
   SafeDelete(scheme);
   scmtree->ResetBranchAddress(scmtree->GetBranch("ScmBranch"));
//?   // delete scheme tree from RAM
//?   SafeDelete(scmtree);

   return err;
}//DoExpress

//______________________________________________________________________________
Int_t XGCProcesSet::DoMultichipExpress(Int_t numdata, TTree **datatree,
                    Int_t numbgrd, TTree **bgrdtree, TFile *file)
{
   // Compute expression values using multichip algorithms, e.g. medpol, farms, DFW
   // Intensities from datatrees are stored:
   // either in one large table in RAM for fast access
   // or in temporary trees in a temporary root file in order to reduce memory consumption
   // Note: all trees must have same number of entries (i.e. identical chip types)
   if(kCS) cout << "------XGCProcesSet::DoMultichipExpress------" << endl;

//TEST Benchmark
//gBenchmark->Reset(); 
//gBenchmark->Start("Bench_MedianPolish");

   Int_t err = errNoErr;
   Int_t idx = 0;
   Int_t x   = 0;
   Int_t y   = 0;
   Int_t ij, start, end, id, entry;

   Int_t numsels = 0;  //number of selected entries
   Int_t exlevel = 0;
   Int_t stepout = (Int_t)((100000.0 + 10.0*numdata)/(Float_t)numdata); //step size for verbose output

// Init min/max expression levels
   Double_t min = DBL_MAX;  //defined in float.h
   Double_t max = 0;

// Init 
   Int_t     *arrMask = 0;
   Int_t     *arrIndx = 0;
   Int_t     *arrUnit = 0;
   Int_t     *mskUnit = 0;
   Double_t  *level   = 0;
   Double_t  *stdev   = 0;

   Double_t **table   = 0;
   Double_t  *arrData = 0;
   Double_t  *arrPM   = 0; 
   Double_t  *dummy   = 0;  //to prevent compilation error

   TTree   *scmtree = 0; 
   TLeaf   *scmleaf = 0;
   XScheme *scheme  = 0;

   TTree   *idxtree = 0; 
   XGCUnit *unit    = 0;

   XGCCell **gccell = 0;
   XBgCell **bgcell = 0;

   TTree     **tmptree  = 0;
   TTree     **exprtree = 0;
   XGCExpression **expr = 0;

   Int_t    bufsize = XManager::GetBufSize(numdata, 10000);
   Int_t    split   = 99;
   Double_t sort    = 0;

   fFile->cd();

// Get chip parameters from scheme file (also alternative CDFs)
   if (datatree[0] == 0) return errGetTree;
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if

// Change directory to scheme file
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

   XGeneChip *chip = (XGeneChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) return fManager->HandleError(errGetScheme, fSchemeName);
   Int_t numrows  = chip->GetNumRows();
   Int_t numcols  = chip->GetNumColumns();
   Int_t numunits = chip->GetNumUnits();
   Int_t size     = numrows*numcols;

// Get scheme tree  and unit tree for scheme
   idxtree = this->UnitTree(fExpressor, &unit, numunits); 
   if (idxtree == 0) return errGetTree;

   scmtree = this->SchemeTree(fExpressor, &scheme, &scmleaf); 
   if (scmtree == 0) return errGetTree;

// Get maximum number of pairs/cells
   Int_t maxnumcells = this->MaxNumberCells(idxtree);
   if (maxnumcells <= 0) return errGeneral;

// Get parameters for background subtraction
   Bool_t doBg = this->BackgroundParameters(fExpressor, fExpressor->GetBgrdOption());

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
//?   } else if (subbgrd || corbgrd) {
   } else if (doBg == kTRUE) {
      cout << "Warning: No background trees available for background subtraction."
           << endl;
      doBg = kFALSE;
//TO DO: case of pars[npars-1] == -100 in BackgroundParameters()
   } else {
      // to prevent subtraction of background in FillDataArrays()
      // e.g. AdjustBackground() creates bgrdtrees but sets numbgrd=0 !!
      doBg = kFALSE;
   }//if

// Init unit selector
   XUnitSelector *unitSelector = new XUnitSelector(kTypeSlct[3], kExtenSlct[3]);
   unitSelector->SetOption((this->ChipType(fTitle)).Data());
   err = unitSelector->InitParameters(fExprSelector->GetNumParameters(),
                                      fExprSelector->GetParameters());
   if (err != errNoErr) goto cleanup;

// Initialize memory for unit arrays
   if (!(arrUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}
   if (!(mskUnit = new (nothrow) Int_t[numunits])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<numunits; i++) arrUnit[i] = mskUnit[i] = 0; 

// Get mask from scheme tree and store in array 
   arrUnit = this->FillUnitArray(idxtree, unit, numunits, arrUnit, mskUnit);
   if (arrUnit == 0) {err = errInitMemory; goto cleanup;}

// Calculate units satisfying mask
   err = unitSelector->Calculate(numunits, arrUnit, mskUnit); //arrUnit to prevent compile error
   if (err != errNoErr) goto cleanup;

// Create local arrays
   if (!(arrMask = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(arrIndx = new (nothrow) Int_t[size]))       {err = errInitMemory; goto cleanup;}
   if (!(level   = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}
   if (!(stdev   = new (nothrow) Double_t[numdata])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<size; i++) {
      arrMask[i] = eINITMASK; 
      arrIndx[i] = 0;
   }//for_i
   for (Int_t i=0; i<numdata; i++) level[i] = stdev[i] = 0.0; 

// Get exon level of annotation
   if (fExprSelector->GetNumParameters() > 0) {
      exlevel = (Int_t)(fExprSelector->GetParameters())[0];
   }//if

// Get mask for PM from scheme tree and store in array 
   arrMask = this->FillMaskArray(chip, scmtree, scheme, exlevel, size, arrMask);
   if (arrMask == 0) {err = errInitMemory; goto cleanup;}

// Calculate mask for expression
   err = fExprSelector->Calculate(size, dummy, dummy, arrMask);
   if (err != errNoErr) goto cleanup;

// Init branch addresses
   gccell = new XGCCell*[numdata];
   bgcell = new XBgCell*[numdata];
   for (Int_t k=0; k<numdata; k++) {
      gccell[k] = 0;
      bgcell[k] = 0;
      datatree[k]->SetBranchAddress("DataBranch", &gccell[k]);
      if (numbgrd > 0) bgrdtree[k]->SetBranchAddress("BgrdBranch", &bgcell[k]);
   }//for_k

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
//DB   datatree[0]->DropBaskets();  //to remove baskets from memory

// Get number of selected entries
   for (Int_t i=0; i<size; i++) {
      numsels = (arrMask[i] == 1) ? ++numsels : numsels;
   }//for_i

//TEST
//gBenchmark->Start("Bench_fill");
// Get data from datatrees (and bgrdtrees) and store in table or array
   if (file == 0) {
      // create table
      if (!(table = CreateTable(numdata, numsels))) {err = errInitMemory; goto cleanup;}

      // fill table with selected intensities of all datatrees
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
   } else {
      // create array to store selected intensities
      if (!(arrData = new (nothrow) Double_t[numsels])) {err = errInitMemory; goto cleanup;}
      for (Int_t i=0; i<numsels;  i++) arrData[i] = 0.0;

      // change directory to temporary file
      if (!file->cd()) {err = errGetDir; goto cleanup;}

      // get data from datatrees (and bgrdtrees) and store in temporary file
      tmptree  = new TTree*[numdata];
      for (Int_t k=0; k<numdata; k++) {
         // create temporary tree
         tmptree[k] = new TTree(datatree[k]->GetName(), "temporary tree");
         if (tmptree[k] == 0) {err = errCreateTree; goto cleanup;}
         tmptree[k]->Branch("sort", &sort, "sort/D", bufsize);

         // informing user
         if (XManager::fgVerbose) {
            cout << "         filling temporary tree <" << tmptree[k]->GetName()
                 << ">...              \r" << flush;
         }//if

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

            datatree[k]->DropBaskets();  //to remove baskets from memory
            bgrdtree[k]->DropBaskets();  //to remove baskets from memory
         } else {
            idx = 0;
            for (Int_t i=0; i<size; i++) {
               datatree[k]->GetEntry(i);

               if (arrMask[i] == 1) {
                  arrData[idx++] = gccell[k]->GetIntensity();
               }//if
            }//for_i

            datatree[k]->DropBaskets();  //to remove baskets from memory
         }//if

         // fill tmptree with array in the order of scheme tree entries for (x,y)
         start = 0;
         end   = 0;
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
                  sort = arrData[arrIndx[ij]];
                  tmptree[k]->Fill();
               }//if
            }//for_j
            start += numcells;
         }//for_id

         // write tmptree to temporary file
         tmptree[k]->Write();
//??         tmptree[k]->Write(TObject::kOverwrite);
         tmptree[k]->DropBaskets();  //to remove baskets from memory
      }//for_k

      if (XManager::fgVerbose) {
         cout << "         finished filling <" << numdata << "> temporary trees.          " << endl;
      }//if
   }//if
//TEST Benchmark
//gBenchmark->Show("Bench_table");

// Change directory to current directory for treeset in main file
   if (!fFile->cd(fName)) {err = errGetDir; goto cleanup;}

// Create new trees exprtree
   exprtree = new TTree*[numdata];
   expr     = new XGCExpression*[numdata];

   for (Int_t k=0; k<numdata; k++) {
      TString dataname = Path2Name(datatree[k]->GetName(), dSEP, ".");
      TString exprname = dataname + "." + fExpressor->GetTitle();
      exprtree[k] = new TTree(exprname, fSchemeName);
      if (exprtree[k] == 0) {err = errCreateTree; goto cleanup;}

      expr[k] = new XGCExpression();
      exprtree[k]->Branch("ExprBranch", "XGCExpression", &expr[k], bufsize, split);

      // to reduce number of baskets in memory when reading trees
      if (file) tmptree[k]->SetMaxVirtualSize(bufsize);
   }//for_k

//TEST
//gBenchmark->Show("Bench_MedianPolish");
//gBenchmark->Reset();
//gBenchmark->Start("Bench_Loop");

// Create array to store PM values for all probes with current unitID
   if (!(arrPM = new (nothrow) Double_t[(maxnumcells+1)*numdata])) {
      err = errInitMemory; goto cleanup;
   }//if

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

//TO DO
//Better above: arrPM = new Double_t[maxnumcells*numdata]; 
      // create array to store PM values for all probes with current unitID
      Int_t  numatoms = unit->GetNumAtoms();
//x      Double_t *arrPM = new Double_t[numcells*numdata]; 

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

            if (file == 0) {
               for (Int_t k=0; k<numdata; k++) {
                  arrPM[p] = table[k][arrIndx[ij]];
                  p++;
               }//for_k
            } else {
               for (Int_t k=0; k<numdata; k++) {
                  tmptree[k]->GetEntry(entry);
                  arrPM[p] = sort;
                  p++;
               }//for_k
            }//if

            entry++;
         }//if
      }//for_j
      start += numcells;

      // fill arrPM or continue if it is not filled
      if ((err = fExpressor->SetArray(p, arrPM)) != errNoErr) {
//x         delete [] arrPM;
         continue;
      }//if

      // calculate expression level for PMs of current unitID
      if ((err = fExpressor->Calculate(numdata, level, stdev, 0))) break;

//////////////////
// TO DO: return fResiduals for residual-plot!!!! (like affyPLM)
//      residuals = fExpressor->GetResiduals();
// store as residuals tree??
//////////////////

      // fill expression trees
      for (Int_t k=0; k<numdata; k++) {
         // get minimal/maximal expression levels
         if (level[k] < min) min = level[k];
         if (level[k] > max) max = level[k];

         expr[k]->SetUnitID(unitID);
         expr[k]->SetLevel(level[k]);
//??         expr[k]->SetStdev(stdev[k]);
         expr[k]->SetStdev(TMath::Abs(stdev[k]));
         expr[k]->SetNumPairs(numatoms);
         exprtree[k]->Fill();
      }//for_k

//x      delete [] arrPM;

/////////////////////////////////
// TO DO??: if (id%(10000/numdata)??) for k: exprtree[k]->DropBaskets??? ev also tmptree[k]->DropBaskets???
/////////////////////////////////

//      if (XManager::fgVerbose && id%10000 == 0) {
      if (XManager::fgVerbose && (idx == 1 || id%stepout == 0)) {
         cout << "      calculating expression for <" << idx << "> of <"
              << numunits << "> units...\r" << flush;
      }//if
   }//for_id

   if (XManager::fgVerbose) {
      cout << "      calculating expression for <" << idx << "> of <"
           << numunits << "> units...Finished." << endl;
   }//if

   if (XManager::fgVerbose) {
      cout << "      expression statistics: " << endl;
      cout << "         minimal expression level is <" << min << ">" << endl;
      cout << "         maximal expression level is <" << max << ">" << endl;
   }//if
//TEST
//gBenchmark->Show("Bench_Loop");

// Write expression trees to file 
   for (Int_t k=0; k<numdata; k++) {
   // Add tree info to tree
      AddExprTreeInfo(exprtree[k], exprtree[k]->GetName(), fExpressor->GetOption(),
                      idx, min, max);
//need to do:
//                      idx, min[k], max[k]);

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
   DeleteTable(table, numdata);

   // delete temporary trees
   if (tmptree) {
      for (Int_t k=0; k<numdata; k++) {
         tmptree[k]->Delete(""); tmptree[k] = 0;
      }//for_k
      delete [] tmptree;
   }//if
   // delete arrays
   if (arrPM)   {delete [] arrPM;   arrPM   = 0;}
   if (arrData) {delete [] arrData; arrData = 0;}
   if (stdev)   {delete [] stdev;   stdev   = 0;}
   if (level)   {delete [] level;   level   = 0;}
   if (arrIndx) {delete [] arrIndx; arrIndx = 0;}
   if (arrMask) {delete [] arrMask; arrMask = 0;}
   if (mskUnit) {delete [] mskUnit; mskUnit = 0;}
   if (arrUnit) {delete [] arrUnit; arrUnit = 0;}

   for (Int_t k=0; k<numdata; k++) {
      SafeDelete(gccell[k]);
      datatree[k]->DropBaskets();  //to remove baskets from memory
      datatree[k]->ResetBranchAddress(datatree[k]->GetBranch("DataBranch"));

      if (numbgrd > 0) {
         SafeDelete(bgcell[k]);
         bgrdtree[k]->DropBaskets();  //to remove baskets from memory
         bgrdtree[k]->ResetBranchAddress(bgrdtree[k]->GetBranch("BgrdBranch"));
      }//if

      SafeDelete(expr[k]);
      exprtree[k]->ResetBranchAddress(exprtree[k]->GetBranch("ExprBranch"));
      SafeDelete(exprtree[k]);
   }//for_k

   delete [] gccell;
   delete [] bgcell;
   delete [] expr;
   delete [] exprtree;

   SafeDelete(unitSelector);
   SafeDelete(unit);
   idxtree->ResetBranchAddress(idxtree->GetBranch("IdxBranch"));
   SafeDelete(scheme);
   scmtree->ResetBranchAddress(scmtree->GetBranch("ScmBranch"));
//?   // delete scheme tree from RAM
//?   SafeDelete(scmtree);

   return err;
}//DoMultichipExpress

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
      name = strtok(strcpy(name, varlist), ":");
      while(name) {
         if (strcmp(name,"fBg")    == 0) {hasBgrd = kTRUE;}
         if (strcmp(name,"fStdev") == 0) {hasStdv = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get trees
   TTree   **tree = new TTree*[n];
   XBgCell **cell = new XBgCell*[n];
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
         if (hasBgrd) output << sep << (names[i] + "_BGRD").Data();
         if (hasStdv) output << sep << (names[i] + "_STDV").Data();
      }//for_i
   } else {
      if (hasBgrd) output << sep << "BGRD";
      if (hasStdv) output << sep << "STDV";
   }//if
   output << endl;

// Loop over tree entries and trees
   Int_t nentries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);
         if (k == 0)  output << cell[k]->GetX() << sep << cell[k]->GetY();
         if (hasBgrd) output << sep << cell[k]->GetBackground();
         if (hasStdv) output << sep << cell[k]->GetStdev();
      }//for_k
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported." << endl;
   }//if

// Cleanup
   for (Int_t k=0; k<n; k++) {
      SafeDelete(cell[k]);
      tree[k]->ResetBranchAddress(tree[k]->GetBranch("BgrdBranch"));
//?? keep tree for later???
      SafeDelete(tree[k]);
   }//for_k

   delete [] cell;
   delete [] tree;

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
   TTree   **tree = new TTree*[n];
   XGCCell **cell = new XGCCell*[n];
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
         if (hasMean) output << sep << (names[i] + "_MEAN").Data();
         if (hasStdv) output << sep << (names[i] + "_STDV").Data();
         if (hasNPix) output << sep << (names[i] + "_NPIXELS").Data();
      }//for_i
   } else {
      if (hasMean) output << sep << "MEAN";
      if (hasStdv) output << sep << "STDV";
      if (hasNPix) output << sep << "NPIXELS";
   }//if
   output << endl;

// Loop over tree entries and tree branches
   Int_t nentries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);
         if (k == 0)  output << cell[k]->GetX() << sep << cell[k]->GetY();
         if (hasMean) output << sep << cell[k]->GetIntensity();
         if (hasStdv) output << sep << cell[k]->GetStdev();
         if (hasNPix) output << sep << cell[k]->GetNumPixels();
      }//for_k
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported." << endl;
   }//if

// Cleanup
   for (Int_t k=0; k<n; k++) {
      SafeDelete(cell[k]);
      tree[k]->ResetBranchAddress(tree[k]->GetBranch("DataBranch"));
//?? keep tree for later???
      SafeDelete(tree[k]);
   }//for_k

   delete [] cell;
   delete [] tree;

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

   Int_t err = errNoErr;

// Decompose varlist
   Short_t hasUnit   = 0;  //unit name
   Short_t hasTrans  = 0;  //transcript_id
   Short_t hasName   = 0;  //gene name
   Short_t hasSymbol = 0;  //gens symbol
   Short_t hasAccess = 0;  //mRNA accession
   Short_t hasEntrez = 0;  //entrez ID
   Short_t hasChromo = 0;  //chromosome
   Short_t hasStart  = 0;  //start position
   Short_t hasStop   = 0;  //stop position
   Short_t hasStrand = 0;  //strand
   Short_t hasCyto   = 0;  //cytoband
   Short_t hasLevel  = 0;  //expression level
   Short_t hasStdev  = 0;  //standard deviation
   Short_t hasNPairs = 0;  //number of pairs; number of atoms

   Short_t hasAnnot  = 0;  //annotation
   Short_t hasData   = 0;  //data

   Short_t idx = 0;
   if (strcmp(varlist,"*")  == 0) {
      hasUnit   = ++idx;
      hasTrans  = ++idx;
      hasName   = ++idx;
      hasSymbol = ++idx;
      hasAccess = ++idx;
      hasEntrez = ++idx;
      hasChromo = ++idx;
      hasStart  = ++idx;
      hasStop   = ++idx;
      hasStrand = ++idx;
      hasCyto   = ++idx;
      hasLevel  = ++idx;
      hasStdev  = ++idx;
      hasNPairs = ++idx;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName")     == 0) {hasUnit   = ++idx;}
         if (strcmp(name,"fTranscriptID") == 0) {hasTrans  = ++idx;}
         if (strcmp(name,"fName")         == 0) {hasName   = ++idx;}
         if (strcmp(name,"fSymbol")       == 0) {hasSymbol = ++idx;}
         if (strcmp(name,"fAccession")    == 0) {hasAccess = ++idx;}
         if (strcmp(name,"fEntrezID")     == 0) {hasEntrez = ++idx;}
         if (strcmp(name,"fChromosome")   == 0) {hasChromo = ++idx;}
         if (strcmp(name,"fStart")        == 0) {hasStart  = ++idx;}
         if (strcmp(name,"fStop")         == 0) {hasStop   = ++idx;}
         if (strcmp(name,"fStrand")       == 0) {hasStrand = ++idx;}
         if (strcmp(name,"fCytoBand")     == 0) {hasCyto   = ++idx;}
         if (strcmp(name,"fLevel")        == 0) {hasLevel  = ++idx;}
         if (strcmp(name,"fStdev")        == 0) {hasStdev  = ++idx;}
         if (strcmp(name,"fNPairs")       == 0) {hasNPairs = ++idx;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

   // check for presence of at least one annotation variable
   hasAnnot = (hasTrans + hasName  + hasSymbol + hasAccess + hasEntrez
            + hasChromo + hasStart + hasStop   + hasStrand + hasCyto);
   // check for presence of at least one of fLevel, fStdev, fNPairs
   hasData = (hasStdev > 0 
           ? (hasStdev <= hasNPairs ? hasStdev : (hasNPairs > 0 ? hasNPairs : hasStdev))
           : hasNPairs);
   hasData = (hasLevel > 0 
           ? (hasLevel <= hasData ? hasLevel : (hasData > 0 ? hasData : hasLevel))
           : hasData);

// Get trees
   TTree         **tree = new TTree*[n];
   XGCExpression **expr = new XGCExpression*[n];
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

   Int_t type = eTRANSCRIPT;
   if      (strcmp(option, "exon")     == 0) type = eEXONTYPE; 
   else if (strcmp(option, "probeset") == 0) type = ePROBESET; 

// Get scheme name (also for alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = tree[0]->GetTitle();
   } else if (!fSchemeName.Contains(tree[0]->GetTitle())) {
      cerr << "Error: Scheme <" << fSchemeName.Data() << "> is not derived from <"
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

// Get unit tree for scheme (needed for annotation)
   XGCUnit *unit   = 0;
//?   XExonUnit *unit = 0;
   TTree *unittree = this->GetUnitTree(chip, type);
   if (unittree == 0) return errGetTree;
   unittree->SetBranchAddress("IdxBranch", &unit);

   Int_t numunits = (Int_t)(unittree->GetEntries());

// Get annotation tree for scheme
   XTransAnnotation *annot = 0;
   TTree *anntree = this->GetAnnotationTree(chip, type);
   Int_t numannot = 0;

// Create hash table to store unit names from anntree
   THashTable *htable = 0;
   if (anntree) {
      anntree->SetBranchAddress("AnnBranch", &annot);
      numannot = (Int_t)(anntree->GetEntries());

      if (!(htable = new THashTable(2*numannot))) return errInitMemory;
      htable = this->FillHashTable(htable, anntree, annot, type);
   } else {
      if (hasAnnot) {
         cout << "Warning: Missing annotation, gene info not exported." << endl;
      }//if
      hasAnnot = kFALSE;
   }//if

// Output header
   output << "UNIT_ID";
   for (Short_t j=1; j<=idx; j++) {
      if (hasUnit == j) output << sep << "UnitName";

      if (hasAnnot > 0) {
         if (hasTrans  == j) output << sep << "ProbesetID";
         if (hasName   == j) output << sep << "GeneName";
         if (hasSymbol == j) output << sep << "GeneSymbol";
         if (hasAccess == j) output << sep << "GeneAccession";
         if (hasEntrez == j) output << sep << "EntrezID";
         if (hasChromo == j) output << sep << "Chromosome";
         if (hasStart  == j) output << sep << "Start";
         if (hasStop   == j) output << sep << "Stop";
         if (hasStrand == j) output << sep << "Strand";
         if (hasCyto   == j) output << sep << "Cytoband";
      }//if

      if (hasData == j) {
         if (n == 1) {
            if (hasLevel)  output << sep << "LEVEL";
            if (hasStdev)  output << sep << "STDEV";
            if (hasNPairs) output << sep << "NPAIRS";
         } else {
            for (Int_t k=0; k<n; k++) {
               if (hasLevel)  output << sep << (names[k] + "_LEVEL").Data();
               if (hasStdev)  output << sep << (names[k] + "_STDEV").Data();
               if (hasNPairs) output << sep << (names[k] + "_NPAIRS").Data();
            }//for_k
         }//if
      }//if
   }//for_j
   output << endl;

// Loop over tree entries and trees
   XIdxString *idxstr = 0;
   Int_t index = 0;
   Int_t cnt   = 0;
   Int_t nentries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      tree[0]->GetEntry(i);

      Int_t unitID = expr[0]->GetUnitID();
      output << unitID;

      unittree->GetEntry(index++);
      while (unitID != unit->GetUnitID()) {
         if (index == numunits) {
           cerr << "Error: UnitID <" << unitID << "> not found." << endl;
           err = errAbort; goto cleanup;
         }//if

         unittree->GetEntry(index++);
      }//while

      for (Short_t j=1; j<=idx; j++) {
         // export unitname
         if (hasUnit == j) {
            output << sep << unit->GetUnitName();
         }//if

         // export annotation
         if (hasAnnot > 0) {
            idxstr = FindUnitID(htable, unit);
            if (idxstr) {
               anntree->GetEntry(idxstr->GetIndex());
               if (hasTrans  == j) output << sep << this->GetTranscriptID(unit, annot, type);
               if (hasName   == j) output << sep << annot->GetName();
               if (hasSymbol == j) output << sep << annot->GetSymbol();
               if (hasAccess == j) output << sep << annot->GetAccession();
               if (hasEntrez == j) output << sep << annot->GetEntrezID();
               if (hasChromo == j) output << sep << annot->GetChromosome();
               if (hasStart  == j) output << sep << annot->GetStart();
               if (hasStop   == j) output << sep << annot->GetStop();
               if (hasStrand == j) output << sep << annot->GetStrand();
               if (hasCyto   == j) output << sep << annot->GetCytoBand();
            } else {
               if (hasTrans  == j) output << sep << "NA";
               if (hasName   == j) output << sep << "NA";
               if (hasSymbol == j) output << sep << "NA";
               if (hasAccess == j) output << sep << "NA";
               if (hasEntrez == j) output << sep << "NA";
               if (hasChromo == j) output << sep << "NA";
               if (hasStart  == j) output << sep << "-1";
               if (hasStop   == j) output << sep << "-1";
               if (hasStrand == j) output << sep << "?";
               if (hasCyto   == j) output << sep << "NA";
            }//if
         }//if

         // export data
         if (hasData == j) {
            for (Int_t k=0; k<n; k++) {
               tree[k]->GetEntry(i);
               if (hasLevel)  output << sep << expr[k]->GetLevel();
               if (hasStdev)  output << sep << expr[k]->GetStdev();
               if (hasNPairs) output << sep << expr[k]->GetNumPairs();
            }//for_k
         }//if
      }//for_j
      output << endl;

/////////////////////////////////
// TO DO??: if (id%10000??) for k: tree[k]->DropBaskets???
// OR???:   tree[k]->SetMaxVirtualSize(bufsize);
/////////////////////////////////

      cnt++;
      if (XManager::fgVerbose && cnt%1000 == 0) {
         cout << "<" << cnt << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << cnt << "> of " << "<" << nentries << "> records exported." << endl;
   }//if

//Cleanup
cleanup:
   if (htable)   {htable->Delete(); delete htable; htable = 0;}
   SafeDelete(schemes);

   // remove trees from RAM
   for (Int_t k=0; k<n; k++) {
      SafeDelete(expr[k]);
      tree[k]->ResetBranchAddress(tree[k]->GetBranch("ExprBranch"));
      SafeDelete(tree[k]);
   }//for_k

   delete [] expr;
   delete [] tree;

   SafeDelete(annot);
   SafeDelete(anntree);

   SafeDelete(unit);
   SafeDelete(unittree);

   return err;
}//ExportExprTrees

//______________________________________________________________________________
Int_t XGCProcesSet::ExportCallTrees(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export data stored in call tree to file output
   if(kCS) cout << "------XGCProcesSet::ExportCallTrees------" << endl;

   Int_t err = errNoErr;

// Decompose varlist
   Short_t hasUnit   = 0;  //unit name
   Short_t hasTrans  = 0;  //transcript_id
   Short_t hasName   = 0;  //gene name
   Short_t hasSymbol = 0;  //gens symbol
   Short_t hasAccess = 0;  //mRNA accession
   Short_t hasEntrez = 0;  //entrez ID
   Short_t hasChromo = 0;  //chromosome
   Short_t hasStart  = 0;  //start position
   Short_t hasStop   = 0;  //stop position
   Short_t hasStrand = 0;  //strand
   Short_t hasCyto   = 0;  //cytoband
   Short_t hasCall   = 0;  //detection call
   Short_t hasPVal   = 0;  //detection p-value

   Short_t hasAnnot  = 0;  //annotation
   Short_t hasData   = 0;  //data

   Short_t idx = 0;
   if (strcmp(varlist,"*")  == 0) {
      hasUnit   = ++idx;
      hasTrans  = ++idx;
      hasName   = ++idx;
      hasSymbol = ++idx;
      hasAccess = ++idx;
      hasEntrez = ++idx;
      hasChromo = ++idx;
      hasStart  = ++idx;
      hasStop   = ++idx;
      hasStrand = ++idx;
      hasCyto   = ++idx;
      hasCall   = ++idx;
      hasPVal   = ++idx;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while (name) {
         if (strcmp(name,"fUnitName")     == 0) {hasUnit   = ++idx;}
         if (strcmp(name,"fTranscriptID") == 0) {hasTrans  = ++idx;}
         if (strcmp(name,"fName")         == 0) {hasName   = ++idx;}
         if (strcmp(name,"fSymbol")       == 0) {hasSymbol = ++idx;}
         if (strcmp(name,"fAccession")    == 0) {hasAccess = ++idx;}
         if (strcmp(name,"fEntrezID")     == 0) {hasEntrez = ++idx;}
         if (strcmp(name,"fChromosome")   == 0) {hasChromo = ++idx;}
         if (strcmp(name,"fStart")        == 0) {hasStart  = ++idx;}
         if (strcmp(name,"fStop")         == 0) {hasStop   = ++idx;}
         if (strcmp(name,"fStrand")       == 0) {hasStrand = ++idx;}
         if (strcmp(name,"fCytoBand")     == 0) {hasCyto   = ++idx;}
         if (strcmp(name,"fCall")         == 0) {hasCall   = ++idx;}
         if (strcmp(name,"fPValue")       == 0) {hasPVal   = ++idx;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

   // check for presence of at least one annotation variable
   hasAnnot = (hasTrans + hasName  + hasSymbol + hasAccess + hasEntrez
            + hasChromo + hasStart + hasStop   + hasStrand + hasCyto);
   // check for presence of at least one of fCall, fPValue
   hasData = (hasCall > 0 
           ? (hasCall <= hasPVal ? hasCall : (hasPVal > 0 ? hasPVal : hasCall))
           : hasPVal);

// Get trees
   TTree  **tree = new TTree*[n];
   XPCall **call = new XPCall*[n];
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

   Int_t type = eTRANSCRIPT;
   if      (strcmp(option, "exon")     == 0) type = eEXONTYPE; 
   else if (strcmp(option, "probeset") == 0) type = ePROBESET; 

// Get scheme name (also for alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = tree[0]->GetTitle();
   } else if (!fSchemeName.Contains(tree[0]->GetTitle())) {
      cerr << "Error: Scheme <" << fSchemeName.Data() << "> is not derived from <"
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

// Get unit tree for scheme (needed for annotation)
   XGCUnit *unit   = 0;
//?   XExonUnit *unit = 0;
   TTree *unittree = this->GetUnitTree(chip, type);
   if (unittree == 0) return errGetTree;
   unittree->SetBranchAddress("IdxBranch", &unit);

   Int_t numunits = (Int_t)(unittree->GetEntries());

// Get annotation tree for scheme
   XTransAnnotation *annot = 0;
   TTree *anntree = this->GetAnnotationTree(chip, type);
   Int_t numannot = 0;

// Create hash table to store unit names from anntree
   THashTable *htable = 0;
   if (anntree) {
      anntree->SetBranchAddress("AnnBranch", &annot);
      numannot = (Int_t)(anntree->GetEntries());

      if (!(htable = new THashTable(2*numannot))) return errInitMemory;
      htable = this->FillHashTable(htable, anntree, annot, type);
   } else {
      if (hasAnnot) {
         cout << "Warning: Missing annotation, gene info not exported." << endl;
      }//if
      hasAnnot = kFALSE;
   }//if

// Output header
   output << "UNIT_ID";
   for (Short_t j=1; j<=idx; j++) {
      if (hasUnit == j) output << sep << "UnitName";

      if (hasAnnot > 0) {
         if (hasTrans  == j) output << sep << "ProbesetID";
         if (hasName   == j) output << sep << "GeneName";
         if (hasSymbol == j) output << sep << "GeneSymbol";
         if (hasAccess == j) output << sep << "GeneAccession";
         if (hasEntrez == j) output << sep << "EntrezID";
         if (hasChromo == j) output << sep << "Chromosome";
         if (hasStart  == j) output << sep << "Start";
         if (hasStop   == j) output << sep << "Stop";
         if (hasStrand == j) output << sep << "Strand";
         if (hasCyto   == j) output << sep << "Cytoband";
      }//if

      if (hasData == j) {
         if (n == 1) {
            if (hasCall)  output << sep << "CALL";
            if (hasPVal)  output << sep << "PVALUE";
         } else {
            for (Int_t k=0; k<n; k++) {
               if (hasCall)  output << sep << (names[k] + "_CALL").Data();
               if (hasPVal)  output << sep << (names[k] + "_PVALUE").Data();
            }//for_k
         }//if
      }//if
   }//for_j
   output << endl;

// Loop over tree entries and tree branches
   XIdxString *idxstr = 0;
   Int_t index = 0;
   Int_t cnt   = 0;
   Int_t nentries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      tree[0]->GetEntry(i);

      Int_t unitID = call[0]->GetUnitID();
      output << unitID;

      unittree->GetEntry(index++);
      while (unitID != unit->GetUnitID()) {
         if (index == numunits) {
           cerr << "Error: UnitID <" << unitID << "> not found." << endl;
           err = errAbort; goto cleanup;
         }//if

         unittree->GetEntry(index++);
      }//while

      for (Short_t j=1; j<=idx; j++) {
         // export unitname
         if (hasUnit == j) {
            output << sep << unit->GetUnitName();
         }//if

         // export annotation
         if (hasAnnot > 0) {
            idxstr = FindUnitID(htable, unit);
            if (idxstr) {
               anntree->GetEntry(idxstr->GetIndex());
               if (hasTrans  == j) output << sep << this->GetTranscriptID(unit, annot, type);
               if (hasName   == j) output << sep << annot->GetName();
               if (hasSymbol == j) output << sep << annot->GetSymbol();
               if (hasAccess == j) output << sep << annot->GetAccession();
               if (hasEntrez == j) output << sep << annot->GetEntrezID();
               if (hasChromo == j) output << sep << annot->GetChromosome();
               if (hasStart  == j) output << sep << annot->GetStart();
               if (hasStop   == j) output << sep << annot->GetStop();
               if (hasStrand == j) output << sep << annot->GetStrand();
               if (hasCyto   == j) output << sep << annot->GetCytoBand();
            } else {
               if (hasTrans  == j) output << sep << "NA";
               if (hasName   == j) output << sep << "NA";
               if (hasSymbol == j) output << sep << "NA";
               if (hasAccess == j) output << sep << "NA";
               if (hasEntrez == j) output << sep << "NA";
               if (hasChromo == j) output << sep << "NA";
               if (hasStart  == j) output << sep << "-1";
               if (hasStop   == j) output << sep << "-1";
               if (hasStrand == j) output << sep << "?";
               if (hasCyto   == j) output << sep << "NA";
            }//if
         }//if

         // export data
         if (hasData == j) {
            for (Int_t k=0; k<n; k++) {
               tree[k]->GetEntry(i);

               if (hasCall) {
                  Int_t cl = call[k]->GetCall();
/*                  char *ch = "NA";
                  if      (cl == 2) ch = "P";
                  else if (cl == 0) ch = "A";
                  else if (cl == 1) ch = "M";
                  output << sep << ch;
*/
                  const char *ch[1];
                  ch[0] = "NA";
                  if      (cl == 2) ch[0] = "P";
                  else if (cl == 0) ch[0] = "A";
                  else if (cl == 1) ch[0] = "M";
                  output << sep << ch[0];
               }//if

               if (hasPVal) {
                  output << sep << call[k]->GetPValue();
               }//if
            }//for_k
         }//if
      }//for_j
      output << endl;

/////////////////////////////////
// TO DO??: if (id%10000??) for k: tree[k]->DropBaskets???
// OR???:   tree[k]->SetMaxVirtualSize(bufsize);
/////////////////////////////////

      cnt++;
      if (XManager::fgVerbose && cnt%1000 == 0) {
         cout << "<" << cnt << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << cnt << "> of " << "<" << nentries << "> records exported." << endl;
   }//if

//Cleanup
cleanup:
   if (htable)   {htable->Delete(); delete htable; htable = 0;}
   SafeDelete(schemes);

   // remove trees from RAM
   for (Int_t k=0; k<n; k++) {
      SafeDelete(call[k]);
      tree[k]->ResetBranchAddress(tree[k]->GetBranch("CallBranch"));
      SafeDelete(tree[k]);
   }//for_k

   delete [] call;
   delete [] tree;

   SafeDelete(annot);
   SafeDelete(anntree);

   SafeDelete(unit);
   SafeDelete(unittree);

   return err;
}//ExportCallTrees

//______________________________________________________________________________
Int_t XGCProcesSet::ProbeMask(XDNAChip *chip, Int_t n, Int_t *msk)
{
   // Get name of probe tree from "chip"
   // For GCBackground fill array msk with GC content:
   // for PM set GC >= 0 , i.e. GC = 0...kProbeLength,
   // for MM set GC < eINITMASK, i.e. GC = eINITMASK - (1...kProbeLength+1)
   if(kCS) cout << "------XGCProcesSet::ProbeMask------" << endl;

// Change dir to scheme file
   TDirectory *savedir = gDirectory;
   if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

// Get probe tree for scheme
   XGCProbe *probe = 0;
   TTree *probetree = (TTree*)gDirectory->Get(chip->GetProbeTree()); 
   if (probetree == 0) return errGetTree;
   probetree->SetBranchAddress("PrbBranch", &probe);

// Get GC content from probe tree and store in array
   Int_t x, y, ij;
   Int_t numcols = chip->GetNumColumns();
   for (Int_t i=0; i<n; i++) {
      probetree->GetEntry(i);

      x  = probe->GetX();
      y  = probe->GetY();
      ij = XY2Index(x, y, numcols);

// IMPORTANT NOTE: do not use the following code although it is faster:
// for (i=0; i<size; i++) {xxtree->GetEntry(i); arrXX[i] = xx->GetXX();}
// Every tree contains the (x,y) coordinates as unique identifier, thus
// it is safer to get (x,y) coordinates and to use ij = x + y*numcols

      if (msk[ij] == 1) {
         msk[ij] = probe->GetNumberGC();
      } else if (msk[ij] == 0) {
         //need to use (numberGC + 1) to avoid setting arrMask=eINITMASK for GC=0!!
         msk[ij] = eINITMASK - (probe->GetNumberGC() + 1);
//no!!         msk[ij] = eINITMASK - probe->GetNumberGC();
      }//if
   }//for_i

   // delete probetree
   SafeDelete(probe);
   probetree->ResetBranchAddress(probetree->GetBranch("PrbBranch"));
   SafeDelete(probetree);

   savedir->cd();

   return errNoErr;
}//ProbeMask

//______________________________________________________________________________
Int_t XGCProcesSet::SchemeMask(XDNAChip *chip, Int_t level, Int_t n, Int_t *msk)
{
   // Get name of scheme tree from "chip"
   // "level=0" is cutoff to prevent PM to be redefined as MM for alternative CDFs
   // Fill array "msk" with mask from scheme tree
   if(kCS) cout << "------XGCProcesSet::SchemeMask------" << endl;

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
   SafeDelete(scheme);
   scmtree->ResetBranchAddress(scmtree->GetBranch("ScmBranch"));
   SafeDelete(scmtree);

   savedir->cd();

   return errNoErr;
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
Int_t *XGCProcesSet::FillUnitArray(TTree *idxtree, XGCUnit *unit,
                     Int_t n, Int_t *arr, Int_t *msk)
{
   // Fill array arr with unitID and array msk with unittype
   if(kCS) cout << "------XGCProcesSet::FillUnitArray------" << endl;

   for (Int_t i=0; i<n; i++) {
      idxtree->GetEntry(i);

      arr[i] = unit->GetUnitID();
      msk[i] = unit->GetUnitType();
   }//for_i

   return arr;
}//FillUnitArray

//______________________________________________________________________________
void XGCProcesSet::FillProbeSets(Int_t &p, Int_t &m, Int_t &idx, Double_t *pm, 
                   Double_t *mm, Double_t *sp, Double_t *sm, Int_t *xp, 
                   Int_t *xm, Int_t msk, Double_t inten, Double_t stdev,
                   Double_t bgrd, Double_t bgdev, Int_t npix)
{
   // Fill PM and MM probesets
   if(kCSa) cout << "------XGCProcesSet::FillProbeSets------" << endl;

   if (msk == 1) {
      if (p == 0) idx++;  //count number of units
      pm[p] = inten;
      sp[p] = stdev;
      xp[p] = npix;
      p++;
   } else if (msk == 0) {
      mm[m] = inten;
      sm[m] = stdev;
      xm[m] = npix;
      m++;
   }//if
}//FillProbeSets

//______________________________________________________________________________
void XGCProcesSet::FillBgrdProbeSets(Int_t &m, Double_t *mm, Double_t *sm, Int_t *xm,
                   Int_t msk, Double_t bgrd, Double_t bgdev, Int_t npix)
{
   // Fill MM probesets with background values
   if(kCSa) cout << "------XGCProcesSet::FillBgrdProbeSets------" << endl;

   if (msk == 1) {
      mm[m] = bgrd;
      sm[m] = bgdev;
      xm[m] = npix; //use numpix from datatree
      m++;
   }//if
}//FillBgrdProbeSets

//______________________________________________________________________________
Int_t XGCProcesSet::MaxNumberCells(TTree *idxtree)
{
   // Get maximum number of cells from tree info for unit tree
   // Note: need to return maximum number of cells, i.e. 2*fMaxNPairs, not pairs 
   //       since for miRNA-1_0.CDF unit tree has: NumCells = NumAtoms
   if(kCS) cout << "------XGCProcesSet::MaxNumberCells------" << endl;

   XUnitTreeInfo *idxinfo = 0;
   idxinfo = (XUnitTreeInfo*)idxtree->GetUserInfo()->FindObject(idxtree->GetName());
   if (!idxinfo) {
      cerr << "Error: Could not get tree info for <" << idxtree->GetName() << ">."
           << endl;
      return errGeneral;
   }//if

   return 2*(Int_t)idxinfo->GetValue("fMaxNPairs");
}//MaxNumberCells

//______________________________________________________________________________
TTree *XGCProcesSet::SchemeTree(XAlgorithm *algorithm, void *scheme, TLeaf **scmleaf)
{
   // Get scheme tree for scheme
   if(kCS) cout << "------XGCProcesSet::SchemeTree------" << endl;


   XGeneChip *chip = (XGeneChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) return 0;

   TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
   if (scmtree == 0) return 0;
   scmtree->SetBranchAddress("ScmBranch", scheme);

   TLeaf *leaf = scmtree->FindLeaf("fUnitID");
   *scmleaf = leaf;

   return scmtree;
}//SchemeTree

//______________________________________________________________________________
TTree *XGCProcesSet::UnitTree(XAlgorithm *algorithm, void *unit, Int_t &numunits)
{
   // Get unit tree for scheme
   if(kCS) cout << "------XGCProcesSet::UnitTree------" << endl;


   XGeneChip *chip = (XGeneChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) return 0;
   numunits = chip->GetNumUnits();

   TTree *idxtree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   if (idxtree == 0) return 0;
   idxtree->SetBranchAddress("IdxBranch", unit);

   return idxtree;
}//UnitTree

//______________________________________________________________________________
Int_t XGCProcesSet::InitTrees(Int_t &numdata, TTree **datatree,
                    Int_t &numbgrd, TTree **bgrdtree)
{
   // Adjust intensity based on fOption
   if(kCS) cout << "------XGCProcesSet::InitTrees------" << endl;

// Create hash table to store datatree names
   TString     tmpstr;
   XIdxString *idxstr = 0;
   THashTable *htable = 0;
   if (!(htable = new THashTable(2*numdata))) return errInitMemory;

// Initialize temporary trees for sorting background trees
   TTree **temptree = new TTree*[numdata];
   for (Int_t k=0; k<numdata; k++) temptree[k] = 0;

   numdata = 0;
   numbgrd = 0;
   for (Int_t k=0; k<fTrees->GetSize(); k++) {
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

   delete [] temptree;

   // delete idxstr before deleting htable
   if (htable) {htable->Delete(); delete htable; htable = 0;}

// Check again for equal number of bgrd trees and corresponding data trees
   if ((numbgrd > 0) && (numbgrd != numdata)) {
      cerr << "Error: <" << (numdata - numbgrd) 
           << "> data trees have no corresponding background tree!" << endl;
      return errAbort;
   }//if

   return errNoErr;
}//InitTrees

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
   Int_t npar1 = npars -1;
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

//TO DO      if (pars[npars-1] == -100) npar1 = npars -2;

      fBgPars[0] = 2; //"correctbg"
      fBgPars[1] = (npars > 0) ? pars[npar1] : 0.5;  //nfrac
   } else if (strcmp(option, "attenuatebg") == 0) {
      fNBgPar = 3;
      fBgPars = new (nothrow) Double_t[fNBgPar];

//TO DO      if (pars[npars-1] == -100) npar1 = npars -2;

      fBgPars[0] = 3; //"attenuatebg"
      fBgPars[1] = (npars > 1) ? pars[npar1-1] : 0.005;  //l
      fBgPars[2] = (npars > 1) ? pars[npar1]   : -1.0;   //h
   }//if

// Do bgrd subtraction?
   Bool_t doBg = kTRUE;
   if (npars > 0 && pars[npars-1] == -100) doBg = kFALSE;

   return doBg;
}//BackgroundParameters

//______________________________________________________________________________
TString XGCProcesSet::ChipType(const char *type)
{
   // Convert ChipType to option type for e.g. XSelector
   if(kCS) cout << "------XGCProcesSet::ChipType------" << endl;

   TString chiptype = "none";

   if (strcmp(type,"GeneChip") == 0) {
      chiptype = "gene";
   } else if (strcmp(type,"SNPChip") == 0) {
      chiptype = "snp";
   } else if (strcmp(type,"GenomeChip") == 0) {
      chiptype = "genome";
   } else if (strcmp(type,"ExonChip") == 0) {
      chiptype = "exon";
   } else if (strcmp(type,"GenePix") == 0) {
      chiptype = "gene";
   }//if

   return chiptype;
}//ChipType

//______________________________________________________________________________
Int_t XGCProcesSet::FillDataArrays(TTree *datatree, TTree *bgrdtree, Bool_t doBg,
                    Int_t nrow, Int_t ncol, Double_t *inten, Double_t *stdev,
                    Int_t *npix, Double_t *bgrd, Double_t *bgdev)
{
   // Fill arrays inten, stdev, npix with data from datatree.
   // If bgrdtree != 0 fill arrays bgrd, bgdev with data from bgrdtree
   // If bgrdtree != 0 and doBg == kTRUE then subtract background from intensity
   if(kCS) cout << "------XGCProcesSet::FillDataArrays------" << endl;

// Init datatree
   XGCCell *gccell = 0;
   datatree->SetBranchAddress("DataBranch", &gccell);

// Get data from datatree and store in arrays
   Int_t size = datatree->GetEntries();
   // compare entries with number of entries in data tree
   if (size != nrow*ncol) {
      TString str = ""; str += nrow*ncol;
      return fManager->HandleError(errNumTreeEntries, datatree->GetName(), str);
   }//if

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

// Get background from bgrdtree and subtract optionally from intensity
   if (bgrdtree) {
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

         if (bgrd)  bgrd[ij]  = bgcell->GetBackground();
         if (bgdev) bgdev[ij] = bgcell->GetStdev();

         // subtract background from intensity
         if (inten && doBg) {
            inten[ij] = this->AdjustIntensity(inten[ij], bgcell->GetBackground(), bgcell->GetStdev());
         }//if
      }//for_i

      SafeDelete(bgcell);
      bgrdtree->DropBaskets();
      bgrdtree->ResetBranchAddress(bgrdtree->GetBranch("BgrdBranch"));
   }//if

   SafeDelete(gccell);
   datatree->DropBaskets();  //to remove baskets from memory
   datatree->ResetBranchAddress(datatree->GetBranch("DataBranch"));

   return errNoErr;
}//FillDataArrays

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
   // compare entries with number of entries in probe tree
   if (size != nrow*ncol) {
      TString str = ""; str += nrow*ncol;
      return fManager->HandleError(errNumTreeEntries, datatree->GetName(), str);
   }//if

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

      SafeDelete(bgcell);
      bgrdtree->DropBaskets();
      bgrdtree->ResetBranchAddress(bgrdtree->GetBranch("BgrdBranch"));
   }//if

   SafeDelete(gccell);
   datatree->DropBaskets();  //to remove baskets from memory
   datatree->ResetBranchAddress(datatree->GetBranch("DataBranch"));

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
   Int_t size = datatree->GetEntries();
   // compare entries with number of entries in probe tree
   if (size != nrow*ncol) {
      TString str = ""; str += nrow*ncol;
      return fManager->HandleError(errNumTreeEntries, datatree->GetName(), str);
   }//if

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

   SafeDelete(gccell);
   datatree->DropBaskets();  //to remove baskets from memory
   datatree->ResetBranchAddress(datatree->GetBranch("DataBranch"));

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

   SafeDelete(bgcell);
   bgrdtree->DropBaskets();
   bgrdtree->ResetBranchAddress(bgrdtree->GetBranch("BgrdBranch"));

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

   XGCCell *newcell = new XGCCell();
   Int_t    bufsize = XManager::GetBufSize();
   Int_t    split   = 99;
   newtree->Branch("DataBranch", "XGCCell", &newcell, bufsize, split);

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

   SafeDelete(newcell);
   newtree->DropBaskets();  //to remove baskets from memory
   newtree->ResetBranchAddress(newtree->GetBranch("DataBranch"));

   SafeDelete(oldcell);
   oldtree->ResetBranchAddress(oldtree->GetBranch("DataBranch"));
//?? delete old tree??
   SafeDelete(oldtree);

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

   XGCCell *cell    = new XGCCell();
   Int_t    bufsize = XManager::GetBufSize();
   Int_t    split   = 99;
   tree->Branch("DataBranch", "XGCCell", &cell, bufsize, split);

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
      if (algorithm->GetFile() == 0) {
         AddTreeHeader(tree->GetName(), algorithm->GetName(), 0,
                       algorithm->GetNumParameters(), algorithm->GetParameters());
      }//if
   }//if

   SafeDelete(cell);
   tree->ResetBranchAddress(tree->GetBranch("DataBranch"));

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

   XCellMask *mask    = new XCellMask();
   Int_t      bufsize = XManager::GetBufSize();
   Int_t      split   = 99;
   tree->Branch("MaskBranch", "XCellMask", &mask, bufsize, split);

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
      if (algorithm->GetFile() == 0) {
         AddTreeHeader(tree->GetName(), algorithm->GetName(), 0,
                       algorithm->GetNumParameters(), algorithm->GetParameters());
      }//if
   }//if

// Cleanup
   SafeDelete(mask);
   tree->ResetBranchAddress(tree->GetBranch("MaskBranch"));
   SafeDelete(tree);

   return err;
}//FillMaskTree

//______________________________________________________________________________
void XGCProcesSet::FillProbeSets(Int_t &p, Int_t &idx, Double_t *pm, Double_t *sp,
                   Int_t *xp, Int_t msk, Double_t inten, Double_t stdev)
{
   // Fill probesets (DABG calls only)
   if(kCSa) cout << "------XGCProcesSet::FillProbeSets------" << endl;

   if (msk >= 0) {
      if (p == 0) idx++;  //count number of units
      pm[p] = inten;
      sp[p] = stdev;
      xp[p] = msk;
      p++;
   }//if
}//FillProbeSets

//______________________________________________________________________________
Int_t XGCProcesSet::MaskArray2GC(XDNAChip *chip, Int_t *msk)
{
   // For DABG only, fill arrMask with GC content (GC>=0 for PM)
   // for PM set GC >= 0 , i.e. GC = 0...kProbeLength,
   // for MM set GC < eINITMASK, i.e. GC = eINITMASK - (1...kProbeLength+1)
   if(kCS) cout << "------XGCProcesSet::MaskArray2GC------" << endl;

// Get mask for PM/MM from scheme tree and store in array 
   Int_t x, y, ij;
   Int_t nrow = chip->GetNumRows();
   Int_t ncol = chip->GetNumColumns();
   Int_t size = nrow*ncol;

// Get probe tree for scheme
   XGCProbe *probe = 0;
   TTree *probetree = (TTree*)gDirectory->Get(chip->GetProbeTree()); 
   if (probetree == 0) return errGetTree;
   probetree->SetBranchAddress("PrbBranch", &probe);

// Get GC content from probe tree and store in array
   for (Int_t i=0; i<size; i++) {
      probetree->GetEntry(i);

      x  = probe->GetX();
      y  = probe->GetY();
      ij = XY2Index(x, y, ncol);

      if (msk[ij] == 1) {
         msk[ij] = probe->GetNumberGC();
      } else if (msk[ij] == 0) {
         //need to use (numberGC + 1) to avoid setting arrMask=eINITMASK for GC=0!!
         msk[ij] = eINITMASK - (probe->GetNumberGC() + 1);
//no         msk[ij] = eINITMASK - probe->GetNumberGC();
      }//if
   }//for_i

   SafeDelete(probe);
   probetree->ResetBranchAddress(probetree->GetBranch("PrbBranch"));

   return errNoErr;
}//MaskArray2GC

//______________________________________________________________________________
Int_t XGCProcesSet::MeanReference(Int_t numdata, TTree **datatree, Int_t numbgrd,
                    TTree **bgrdtree, Int_t nrow, Int_t ncol, Double_t *arr, Bool_t doBg)
{
   // Fill array with trimmed mean values from numdata datatrees
   if(kCS) cout << "------XGCProcesSet::MeanReference------" << endl;

// Init branch addresses
   XBgCell **bgcell = new XBgCell*[numdata];
   XGCCell **gccell = new XGCCell*[numdata];
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

// Cleanup
   for (Int_t k=0; k<numdata; k++) {
      SafeDelete(gccell[k]);
      datatree[k]->DropBaskets();
      datatree[k]->ResetBranchAddress(datatree[k]->GetBranch("DataBranch"));

      if (numbgrd > 0) {
         SafeDelete(bgcell[k]);
         bgrdtree[k]->DropBaskets();
         bgrdtree[k]->ResetBranchAddress(bgrdtree[k]->GetBranch("BgrdBranch"));
      }//if
   }//for_k

   delete [] gccell;
   delete [] bgcell;
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
void XGenomeProcesSet::FillProbeSets(Int_t &p, Int_t &m, Int_t &idx, Double_t *pm, 
                       Double_t *mm, Double_t *sp, Double_t *sm, Int_t *xp, 
                       Int_t *xm, Int_t msk, Double_t inten, Double_t stdev,
                       Double_t bgrd, Double_t bgdev, Int_t npix)
{
   // Fill PM probesets and set MM probesets to background
   if(kCSa) cout << "------XGenomeProcesSet::FillProbeSets------" << endl;

   if (msk == 1) {
      if (p == 0) idx++;  //count number of units
      pm[p] = inten;
      sp[p] = stdev;
      xp[p] = npix;
      p++;

      // set MM values to background data
      mm[m] = bgrd;
      sm[m] = bgdev;
      xm[m] = npix; //use numpix from datatree
      m++;
   }//if
}//FillProbeSets

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

//______________________________________________________________________________
Int_t XGenomeProcesSet::MaxNumberCells(TTree *idxtree)
{
   // Get maximum number of cells from tree info for unit tree
   if(kCS) cout << "------XGenomeProcesSet::MaxNumberCells------" << endl;

   XGenomeTreeInfo *idxinfo = 0;
   idxinfo = (XGenomeTreeInfo*)idxtree->GetUserInfo()->FindObject(idxtree->GetName());
   if (!idxinfo) {
      cerr << "Error: Could not get tree info for <" << idxtree->GetName() << ">."
           << endl;
      return errGeneral;
   }//if

   return (Int_t)idxinfo->GetValue("fMaxNCells");
}//MaxNumberCells


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
Int_t *XExonProcesSet::FillMaskArray(XDNAChip *chip, TTree *scmtree, XScheme *scheme,
                       Int_t level, Int_t n, Int_t *msk)
{
   // "level=0" is cutoff to prevent PM to be redefined as MM for alternative CDFs
   // Fill array "msk" with mask from scheme tree
   if(kCS) cout << "------XExonProcesSet::FillMaskArray------" << endl;

//?   return msk = XGCProcesSet::FillMaskArray(chip, scmtree, scheme, level, n, msk);
   return msk = XGenomeProcesSet::FillMaskArray(chip, scmtree, scheme, level, n, msk);
}//FillMaskArray

//______________________________________________________________________________
TTree *XExonProcesSet::SchemeTree(XAlgorithm *algorithm, void *scheme, TLeaf **scmleaf)
{
   // Get scheme tree for scheme
   if(kCS) cout << "------XExonProcesSet::SchemeTree------" << endl;


   XExonChip *chip = (XExonChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) return 0;

   TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
   if (scmtree == 0) return 0;

//   scheme = ((XExonScheme*)scheme);
   XExonScheme **scheme_ptr_ptr = (XExonScheme**)scheme;
   scmtree->SetBranchAddress("ScmBranch", scheme);

   TLeaf *leaf = 0;
   if (strcmp(algorithm->GetOption(), "exon") == 0) { 
      leaf = scmtree->FindLeaf("fExonID");      
   } else if (strcmp(algorithm->GetOption(), "probeset") == 0) {
      leaf = scmtree->FindLeaf("fProbesetID");      
   } else {
      leaf = scmtree->FindLeaf("fUnitID");      
   }//if
   *scmleaf = leaf;

   return scmtree;
}//SchemeTree

//______________________________________________________________________________
TTree *XExonProcesSet::UnitTree(XAlgorithm *algorithm, void *unit, Int_t &numunits)
{
   // Get unit tree for scheme
   if(kCS) cout << "------XExonProcesSet::UnitTree------" << endl;


   XExonChip *chip = (XExonChip*)fSchemes->FindObject(fSchemeName, kTRUE);
   if (!chip) return 0;

   TTree *idxtree = 0; 
   if (strcmp(algorithm->GetOption(), "exon") == 0) { 
      idxtree  = (TTree*)gDirectory->Get(chip->GetExonUnitTree()); //tree.exn
      if (idxtree == 0) return 0;
      numunits = chip->GetNumExonUnits();
   } else if (strcmp(algorithm->GetOption(), "probeset") == 0) {
      idxtree  = (TTree*)gDirectory->Get(chip->GetProbesetUnitTree()); //tree.pbs
      if (idxtree == 0) return 0;
      numunits = idxtree->GetEntries();
   } else {
      idxtree  = (TTree*)gDirectory->Get(chip->GetUnitTree());     //tree.idx
      if (idxtree == 0) return 0;
      numunits = chip->GetNumUnits();
   }//if
   idxtree->SetBranchAddress("IdxBranch", unit);

   return idxtree;
}//UnitTree

//______________________________________________________________________________
const char *XExonProcesSet::GetTranscriptID(XUnit *unit, XTransAnnotation *annot, Int_t type)
{
   // Get transcript IDs from unit tree

   Int_t id = (type == eTRANSCRIPT) 
            ? ((XExonUnit*)unit)->GetSubUnitID()
            : ((XProbesetAnnotation*)annot)->GetTranscriptIX();

   return Form("%d", id);
}//GetTranscriptID
