// File created: 10/27/2001                          last modified: 12/11/2005
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
* Oct 2001 - Initial versions finished
* Jun 2005 - Derive classes from base class XDataTypeInfo.
* Oct 2005 - Add member fDataType to base class XDataTypeInfo.
*
******************************************************************************/

#include <Riostream.h>
#include "TError.h"

#include "XPSDataTypes.h"

const Bool_t  kCS  = 0; //debug: print function names

ClassImp(XDataTypeInfo);
ClassImp(XDatabaseInfo);
ClassImp(XProjectInfo);
ClassImp(XAuthorInfo);
ClassImp(XLoginInfo);
ClassImp(XDatasetInfo);
ClassImp(XSourceInfo);
ClassImp(XArrayInfo);
ClassImp(XHybInfo);
ClassImp(XSampleInfo);
ClassImp(XCellLineInfo);
ClassImp(XPrimaryCellInfo);
ClassImp(XTissueInfo);
ClassImp(XBiopsyInfo);
ClassImp(XTreatment);
ClassImp(XTreatmentInfo);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataTypeInfo                                                        //
//                                                                      //
// Base class for data types                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDataTypeInfo::XDataTypeInfo()
              :TNamed()
{
   // Default datatype information constructor
   if(kCS) cout << "---XDataTypeInfo::XDataTypeInfo(default)------" << endl;

   fDataType = "";
   fHasInfo  = kFALSE;
}//Constructor

//______________________________________________________________________________
XDataTypeInfo::XDataTypeInfo(const char *name, const char *title, const char *type)
              :TNamed(name, title)
{
   // Datatype information constructor
   if(kCS) cout << "---XDataTypeInfo::XDataTypeInfo------" << endl;

   fDataType = type;
   fHasInfo  = kFALSE;
}//Constructor

//______________________________________________________________________________
XDataTypeInfo::XDataTypeInfo(const XDataTypeInfo &info) 
              :TNamed(info),fHasInfo(info.fHasInfo)
{
   // Datatype information copy constructor
   if(kCS) cout << "---XDataTypeInfo::XDataTypeInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XDataTypeInfo& XDataTypeInfo::operator=(const XDataTypeInfo& rhs)
{
   // Datatype information assignment operator.
   if(kCS) cout << "---XDataTypeInfo::operator=------" << endl;

   if (this != &rhs) {
      TNamed::operator=(rhs);
      fHasInfo = rhs.fHasInfo;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XDataTypeInfo::~XDataTypeInfo()
{
   // Datatype information destructor
   if(kCS) cout << "---XDataTypeInfo::~XDataTypeInfo------" << endl;

}//Destructor



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDatabaseInfo                                                        //
//                                                                      //
// Class describing XPS database                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDatabaseInfo::XDatabaseInfo()
              :XDataTypeInfo()
{
   // Default database information constructor
   if(kCS) cout << "---XDatabaseInfo::XDatabaseInfo(default)------" << endl;

   fDirectory = "";
   fAdminID   = "";
   fLoginList = new TList();
}//Constructor

//______________________________________________________________________________
XDatabaseInfo::XDatabaseInfo(const char *name, const char *directory)
              :XDataTypeInfo("Database", name)
{
   // Database information constructor
   if(kCS) cout << "---XDatabaseInfo::XDatabaseInfo------" << endl;

   fDirectory = directory;
   fAdminID   = "";
   fLoginList = new TList();
   fHasInfo   = kTRUE;
}//Constructor

//______________________________________________________________________________
XDatabaseInfo::XDatabaseInfo(const XDatabaseInfo &info) 
              :XDataTypeInfo(info),fAdminID(info.fAdminID),fDirectory(info.fDirectory)
{
   // Database information copy constructor
   if(kCS) cout << "---XDatabaseInfo::XDatabaseInfo(copy)------" << endl;

//to do: test if correct?
   fLoginList = 0;
   if (info.fLoginList != 0) {
      fLoginList = new TList();
      for (Int_t i=0; i<(info.fLoginList)->GetSize(); i++) {
         fLoginList->AddAt((info.fLoginList)->At(i), i);
      }//for_i
   }//if
}//CopyConstructor

//______________________________________________________________________________
XDatabaseInfo& XDatabaseInfo::operator=(const XDatabaseInfo& rhs)
{
   // Database information assignment operator.
   if(kCS) cout << "---XDatabaseInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fAdminID   = rhs.fAdminID;
      fDirectory = rhs.fDirectory;
//to do: need to copy fLoginList!!!
   cout << "Error: Copy of fLoginList not yet implemented." << endl;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XDatabaseInfo::~XDatabaseInfo()
{
   // Database information destructor
   if(kCS) cout << "---XDatabaseInfo::~XDatabaseInfo" << endl;

   if (fLoginList) {fLoginList->Delete(); delete fLoginList; fLoginList = 0;}
}//Destructor

//______________________________________________________________________________
void XDatabaseInfo::AddLoginInfo(XLoginInfo *login)
{
   // Add Login information to Login list
   if(kCS) cout << "------XDatabaseInfo::AddLoginInfo------" << endl;

   fLoginList->Add(login); 
}//AddLoginInfo

//______________________________________________________________________________
void XDatabaseInfo::RemoveLoginInfo(const char *userID)
{
   // Remove Login information to Login list
   if(kCS) cout << "------XDatabaseInfo::RemoveLoginInfo------" << endl;

   XLoginInfo *login = 0;
   login = (XLoginInfo*)fLoginList->FindObject(userID); 
   fLoginList->Remove(login); 
}//RemoveLoginInfo

//______________________________________________________________________________
TString XDatabaseInfo::GetPassword(const char *userID)
{
   // Return password for userID
   if(kCS) cout << "------XDatabaseInfo::GetPassword------" << endl;

   TString password = "";

   TIter next(fLoginList);
   XLoginInfo *login = 0;
   while ((login = (XLoginInfo*)next())) {
      if (strcmp(login->GetName(),userID) == 0) {
         password = login->GetPassword();
         break;
      }//if
   }//while

   return password;
}//GetPassword

//______________________________________________________________________________
Bool_t XDatabaseInfo::IsPresentID(const char *userID)
{
   // Return true if userID is in login list
   if(kCS) cout << "------XDatabaseInfo::IsPresentID------" << endl;

   Bool_t isPresent = kFALSE;

   TIter next(fLoginList);
   XLoginInfo *login = 0;
   while ((login = (XLoginInfo*)next())) {
      if (strcmp(login->GetName(),userID) == 0) {
         isPresent = kTRUE;
         break;
      }//if
   }//while

   return isPresent;
}//IsPresentID


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProjectInfo                                                         //
//                                                                      //
// Class describing XPS project                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XProjectInfo::XProjectInfo()
             :XDataTypeInfo()
{
   // Default Project information constructor
   if(kCS) cout << "---XProjectInfo::XProjectInfo(default)------" << endl;

   fDate        = 0;
   fDescription = "";
}//Constructor

//______________________________________________________________________________
XProjectInfo::XProjectInfo(const char *name, Long_t date, const char *description)
             :XDataTypeInfo("Project", name)
{
   // Project information constructor
   if(kCS) cout << "---XProjectInfo::XProjectInfo------" << endl;

   fDate        = date;
   fDescription = description;
//?? check if necessary infos name and date are set!
//also in other datatype infos!!!
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XProjectInfo::XProjectInfo(const XProjectInfo &info) 
             :XDataTypeInfo(info),fDescription(info.fDescription)
{
   // Project information copy constructor
   if(kCS) cout << "---XProjectInfo::XProjectInfo(copy)------" << endl;

   fDate = info.fDate;
}//CopyConstructor

//______________________________________________________________________________
XProjectInfo& XProjectInfo::operator=(const XProjectInfo& rhs)
{
   // Project information assignment operator.
   if(kCS) cout << "---XProjectInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fDate        = rhs.fDate;
      fDescription = rhs.fDescription;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XProjectInfo::~XProjectInfo()
{
   // Project information destructor
   if(kCS) cout << "---XProjectInfo::~XProjectInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAuthorInfo                                                          //
//                                                                      //
// Class describing author, e.g. administrator, submitter               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XAuthorInfo::XAuthorInfo()
            :XDataTypeInfo()
{
   // Default author information constructor
   if(kCS) cout << "---XAuthorInfo::XAuthorInfo(default)------" << endl;

   fFirstName = "";
   fCompany   = "";
   fLabName   = "";
   fMail      = "";
   fPhone     = "";
}//Constructor

//______________________________________________________________________________
XAuthorInfo::XAuthorInfo(const char *last, const char *first, const char *company,
             const char *lab, const char *mail, const char *phone)
            :XDataTypeInfo("Author", last)
{
   // Author information constructor
   if(kCS) cout << "---XAuthorInfo::XAuthorInfo------" << endl;

   fFirstName = first;
   fCompany   = company;
   fLabName   = lab;
   fMail      = mail;
   fPhone     = phone;
   fHasInfo   = kTRUE;
}//Constructor

//______________________________________________________________________________
XAuthorInfo::XAuthorInfo(const XAuthorInfo &info) 
            :XDataTypeInfo(info),fFirstName(info.fFirstName),fCompany(info.fCompany),
             fLabName(info.fLabName),fMail(info.fMail),fPhone(info.fPhone)
{
   // Author information copy constructor
   if(kCS) cout << "---XAuthorInfo::XAuthorInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XAuthorInfo& XAuthorInfo::operator=(const XAuthorInfo& rhs)
{
   // Author information assignment operator.
   if(kCS) cout << "---XAuthorInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fFirstName = rhs.fFirstName;
      fCompany   = rhs.fCompany;
      fLabName   = rhs.fLabName;
      fMail      = rhs.fMail;
      fPhone     = rhs.fPhone;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XAuthorInfo::~XAuthorInfo()
{
   // Author information destructor
   if(kCS) cout << "---XAuthorInfo::~XAuthorInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XLoginInfo                                                           //
//                                                                      //
// Class to store userID and password (and encrpyt/decrypt password)    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XLoginInfo::XLoginInfo()
           :XDataTypeInfo()
{
   //Default Login information constructor
   if(kCS) cout << "---XLoginInfo::XLoginInfo(default)------" << endl;

   fPassword = "";
}//Constructor

//______________________________________________________________________________
XLoginInfo::XLoginInfo(const char *userID, const char *password)
           :XDataTypeInfo("Login", userID)
{
   // Login information constructor
   if(kCS) cout << "---XLoginInfo::XLoginInfo------" << endl;

   SetName(userID);
   this->SetPassword(password);

   fHasInfo = kTRUE;
}//Constructor

//______________________________________________________________________________
XLoginInfo::XLoginInfo(const XLoginInfo &info) 
           :XDataTypeInfo(info),fPassword(info.fPassword)
{
   // Login information copy constructor
   if(kCS) cout << "---XLoginInfo::XLoginInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XLoginInfo& XLoginInfo::operator=(const XLoginInfo& rhs)
{
   // Login information assignment operator.
   if(kCS) cout << "---XLoginInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fPassword = rhs.fPassword;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XLoginInfo::~XLoginInfo()
{
   // Login information destructor
   if(kCS) cout << "---XLoginInfo::~XLoginInfo------" << endl;

   fPassword = "";
}//Destructor

//______________________________________________________________________________
void XLoginInfo::SetPassword(const char *name)
{
   // Encrypt and store encrypted password
   if(kCS) cout << "------XLoginInfo::SetPassword------" << endl;

   TString password = name; 

// encrypt password
   if (password != "") {
      for (Int_t i = 0; i < password.Length(); i++) {
         char inv = ~password(i);
         password.Replace(i, 1, inv);
      }//for
   }//if

   fPassword = password; 
}//SetPassword

//______________________________________________________________________________
TString XLoginInfo::GetPassword()
{
   // Decrypt and return decrypted password 
   if(kCS) cout << "------XLoginInfo::GetPassword------" << endl;

   TString password = fPassword; 

// decrypt password
   if (password != "") {
      for (Int_t i = 0; i < password.Length(); i++) {
         char inv = ~password(i);
         password.Replace(i, 1, inv);
      }//for
   }//if

   return password;
}//GetPassword


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDatasetInfo                                                         //
//                                                                      //
// Class describing XPS dataset                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDatasetInfo::XDatasetInfo()
             :XDataTypeInfo()
{
   // Default dataset information constructor
   if(kCS) cout << "---XDatasetInfo::XDatasetInfo(default)------" << endl;

   fType        = "";
   fSample      = "";
   fSubmitter   = "";
   fDate        = 0;
   fDescription = "";
}//Constructor

//______________________________________________________________________________
XDatasetInfo::XDatasetInfo(const char *name, const char *type, const char *sample,
              const char *submitter, Long_t date, const char *description)
             :XDataTypeInfo("Dataset", name)
{
   // Dataset information constructor
   if(kCS) cout << "---XDatasetInfo::XDatasetInfo------" << endl;

   fType        = type;
   fSample      = sample;
   fSubmitter   = submitter;
   fDate        = date;
   fDescription = description;
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XDatasetInfo::XDatasetInfo(const XDatasetInfo &info) 
             :XDataTypeInfo(info),fType(info.fType),fSample(info.fSample),
              fSubmitter(info.fSubmitter),fDescription(info.fDescription)
{
   // Dataset information copy constructor
   if(kCS) cout << "---XDatasetInfo::XDatasetInfo(copy)------" << endl;

   fDate = info.fDate;
}//CopyConstructor

//______________________________________________________________________________
XDatasetInfo& XDatasetInfo::operator=(const XDatasetInfo& rhs)
{
   // Dataset information assignment operator.
   if(kCS) cout << "---XDatasetInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fType        = rhs.fType;
      fSample      = rhs.fSample;
      fSubmitter   = rhs.fSubmitter;
      fDate        = rhs.fDate;
      fDescription = rhs.fDescription;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XDatasetInfo::~XDatasetInfo()
{
   // Dataset information destructor
   if(kCS) cout << "---XDatasetInfo::~XDatasetInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSourceInfo                                                          //
//                                                                      //
// Class describing sample source                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSourceInfo::XSourceInfo()
            :XDataTypeInfo()
{
   // Default source information constructor
   if(kCS) cout << "---XSourceInfo::XSourceInfo(default)------" << endl;

   fSpecies     = "";
   fSubSpecies  = "";
   fDescription = "";
}//Constructor

//______________________________________________________________________________
XSourceInfo::XSourceInfo(const char *source, const char *species,
             const char *subspecies, const char *description)
            :XDataTypeInfo("Source", source)
{
   // Source information constructor
   if(kCS) cout << "---XSourceInfo::XSourceInfo------" << endl;

   fSpecies     = species;
   fSubSpecies  = subspecies;
   fDescription = description;
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XSourceInfo::XSourceInfo(const XSourceInfo &info) 
            :XDataTypeInfo(info),fSpecies(info.fSpecies),fSubSpecies(info.fSubSpecies),
             fDescription(info.fDescription)
{
   // Source information copy constructor
   if(kCS) cout << "---XSourceInfo::XSourceInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XSourceInfo& XSourceInfo::operator=(const XSourceInfo& rhs)
{
   // Source information assignment operator.
   if(kCS) cout << "---XSourceInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fSpecies     = rhs.fSpecies;
      fSubSpecies  = rhs.fSubSpecies;
      fDescription = rhs.fDescription;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XSourceInfo::~XSourceInfo()
{
   // Source information destructor
   if(kCS) cout << "---XSourceInfo::~XSourceInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XArrayInfo                                                           //
//                                                                      //
// Class describing microarray used for hybridization                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XArrayInfo::XArrayInfo()
           :XDataTypeInfo()
{
   // Default array information constructor
   if(kCS) cout << "---XArrayInfo::XArrayInfo(default)------" << endl;

   fType        = "";
   fDescription = "";
}//Constructor

//______________________________________________________________________________
XArrayInfo::XArrayInfo(const char *name, const char *type, const char *description)
           :XDataTypeInfo("Array", name)
{
   // Array information constructor
   if(kCS) cout << "---XArrayInfo::XArrayInfo------" << endl;

   fType        = type;
   fDescription = description;
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XArrayInfo::XArrayInfo(const XArrayInfo &info) 
           :XDataTypeInfo(info),fType(info.fType),fDescription(info.fDescription)
{
   // Array information copy constructor
   if(kCS) cout << "---XArrayInfo::XArrayInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XArrayInfo& XArrayInfo::operator=(const XArrayInfo& rhs)
{
   // Array information assignment operator.
   if(kCS) cout << "---XArrayInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fType        = rhs.fType;
      fDescription = rhs.fDescription;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XArrayInfo::~XArrayInfo()
{
   // Array information destructor
   if(kCS) cout << "---XArrayInfo::~XArrayInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHybInfo                                                             //
//                                                                      //
// Class describing hybridization parameters                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XHybInfo::XHybInfo()
         :XDataTypeInfo()
{
   // Default hybridization information constructor
   if(kCS) cout << "---XHybInfo::XHybInfo(default)------" << endl;

   fInput       = "";
   fPreparation = "";
   fProtocol    = "";
   fReplicaName = "";
   fComment     = "";
   fDate        = 0;
   fReplica     = 0;
}//Constructor

//______________________________________________________________________________
XHybInfo::XHybInfo(const char *hybname, const char *input, Long_t date,
          Int_t replica, const char *prep, const char *protocol, const char *comment)
         :XDataTypeInfo("Hybridization", hybname)
{
   // Hybridization information constructor
   if(kCS) cout << "---XHybInfo::XHybInfo------" << endl;

//?? in title full input?? then fInput not needed!!
//   TString name = Path2Name(input,"/",".");
//   SetTitle(name);
   fInput       = input;
   fPreparation = prep;
   fProtocol    = protocol;
   fReplicaName = "";
   fComment     = comment;
   fDate        = date;
   fReplica     = replica;
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XHybInfo::XHybInfo(const XHybInfo &info) 
         :XDataTypeInfo(info),fInput(info.fInput),fPreparation(info.fPreparation),
          fProtocol(info.fProtocol),fReplicaName(info.fReplicaName),
          fComment(info.fComment)
{
   // Hybridization information copy constructor
   if(kCS) cout << "---XHybInfo::XHybInfo(copy)------" << endl;

   fDate    = info.fDate;
   fReplica = info.fReplica;
}//CopyConstructor

//______________________________________________________________________________
XHybInfo& XHybInfo::operator=(const XHybInfo& rhs)
{
   // Hybridization information assignment operator.
   if(kCS) cout << "---XHybInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fInput       = rhs.fInput;
      fPreparation = rhs.fPreparation;
      fProtocol    = rhs.fProtocol;
      fReplicaName = rhs.fReplicaName;
      fComment     = rhs.fComment;
      fDate        = rhs.fDate;
      fReplica     = rhs.fReplica;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XHybInfo::~XHybInfo()
{
   // Hybridization information destructor
   if(kCS) cout << "---XHybInfo::~XHybInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSampleInfo                                                          //
//                                                                      //
// Class describing sample used for hybridization                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSampleInfo::XSampleInfo()
            :XDataTypeInfo()
{
   // Default sample information constructor
   if(kCS) cout << "---XSampleInfo::XSampleInfo(default)------" << endl;

   fSex         = "";
   fPhenotype   = "";
   fGenotype    = "";
   fExtraction  = "";
   fComment     = "";
   fXenoStrain  = "";
   fXenoSex     = "";
   fAgeUnits    = "";
   fXenoAge     = 0;
   fIsXenograft = kFALSE;
}//Constructor

//______________________________________________________________________________
XSampleInfo::XSampleInfo(const char *name, const char *sex, const char *pheno, 
             const char *geno, const char *extract, const char *comment,
             Bool_t isXeno, const char *xenostrain, const char *xenosex,
             Double_t xenoage, const char *ageunits)
            :XDataTypeInfo("Sample", name)
{
   // Sample information constructor
   if(kCS) cout << "---XSampleInfo::XSampleInfo------" << endl;

   fSex         = sex;
   fPhenotype   = pheno;
   fGenotype    = geno;
   fExtraction  = extract;
   fComment     = comment;
   fXenoStrain  = xenostrain;
   fXenoSex     = xenosex;
   fAgeUnits    = ageunits;
   fXenoAge     = xenoage;
   fIsXenograft = isXeno;
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XSampleInfo::XSampleInfo(const XSampleInfo &info) 
            :XDataTypeInfo(info),fSex(info.fSex),fPhenotype(info.fPhenotype),
             fGenotype(info.fGenotype),fExtraction(info.fExtraction),
             fComment(info.fComment),fXenoStrain(info.fXenoStrain),
             fXenoSex(info.fXenoSex),fAgeUnits(info.fAgeUnits)
{
   // Sample information copy constructor
   if(kCS) cout << "---XSampleInfo::XSampleInfo(copy)------" << endl;

   fXenoAge     = info.fXenoAge;
   fIsXenograft = info.fIsXenograft;
}//CopyConstructor

//______________________________________________________________________________
XSampleInfo& XSampleInfo::operator=(const XSampleInfo& rhs)
{
   // Sample information assignment operator.
   if(kCS) cout << "---XSampleInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fSex         = rhs.fSex;
      fPhenotype   = rhs.fPhenotype;
      fGenotype    = rhs.fGenotype;
      fExtraction  = rhs.fExtraction;
      fComment     = rhs.fComment;
      fXenoStrain  = rhs.fXenoStrain;
      fXenoSex     = rhs.fXenoSex;
      fAgeUnits    = rhs.fAgeUnits;
      fXenoAge     = rhs.fXenoAge;
      fIsXenograft = rhs.fIsXenograft;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XSampleInfo::~XSampleInfo()
{
   // Sample information destructor
   if(kCS) cout << "---XSampleInfo::~XSampleInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCellLineInfo                                                        //
//                                                                      //
// Class describing sample of type cell line                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XCellLineInfo::XCellLineInfo()
              :XSampleInfo()
{
   // Default cell line information constructor
   if(kCS) cout << "---XCellLineInfo::XCellLineInfo(default)------" << endl;

   fCellType     = "";
   fParentCell   = "";
   fATCC         = "";
   fModification = "";
}//Constructor

//______________________________________________________________________________
XCellLineInfo::XCellLineInfo(const char *name, const char *type, const char *parent, 
               const char *atcc, const char *mod, const char *sex,
               const char *pheno, const char *geno, const char *extract,
               const char *comment, Bool_t isXeno, const char *xenostrain,
               const char *xenosex, Double_t xenoage, const char *ageunits)
              :XSampleInfo(name, sex, pheno, geno, extract, comment,
               isXeno, xenostrain, xenosex, xenoage, ageunits)
{
   // Cell line information constructor
   if(kCS) cout << "---XCellLineInfo::XCellLineInfo------" << endl;

   fCellType     = type;
   fParentCell   = parent;
   fATCC         = atcc;
   fModification = mod;
}//Constructor

//______________________________________________________________________________
XCellLineInfo::XCellLineInfo(const XCellLineInfo &info) 
              :XSampleInfo(info),fCellType(info.fCellType),fParentCell(info.fParentCell),
               fATCC(info.fATCC),fModification(info.fModification)
{
   // Cell line information copy constructor
   if(kCS) cout << "---XCellLineInfo::XCellLineInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XCellLineInfo& XCellLineInfo::operator=(const XCellLineInfo& rhs)
{
   // Cell line information assignment operator.
   if(kCS) cout << "---XProjectInfo::operator=------" << endl;

   if (this != &rhs) {
      XSampleInfo::operator=(rhs);
      fCellType     = rhs.fCellType;
      fParentCell   = rhs.fParentCell;
      fATCC         = rhs.fATCC;
      fModification = rhs.fModification;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XCellLineInfo::~XCellLineInfo()
{
   // Cell line information destructor
   if(kCS) cout << "---XCellLineInfo::~XCellLineInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPrimaryCellInfo                                                     //
//                                                                      //
// Class describing sample of type primary cell                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XPrimaryCellInfo::XPrimaryCellInfo()
                 :XSampleInfo()
{
   // Default primary cell information constructor
   if(kCS) cout << "---XPrimaryCellInfo::XPrimaryCellInfo(default)------" << endl;

   fCellType      = "";
   fIsolationDate = 0;
   fDescription   = "";
}//Constructor

//______________________________________________________________________________
XPrimaryCellInfo::XPrimaryCellInfo(const char *name, const char *type,  
                  Long_t date, const char *description, const char *sex,
                  const char *pheno, const char *geno, const char *extract,
                  const char *comment, Bool_t isXeno, const char *xenostrain,
                  const char *xenosex, Double_t xenoage, const char *xageunits)
                 :XSampleInfo(name, sex, pheno, geno, extract, comment,
                  isXeno, xenostrain, xenosex, xenoage, xageunits)
{
   // Primary cell information constructor
   if(kCS) cout << "---XPrimaryCellInfo::XPrimaryCellInfo------" << endl;

   fCellType      = type;
   fIsolationDate = date;
   fDescription   = description;
}//Constructor

//______________________________________________________________________________
XPrimaryCellInfo::XPrimaryCellInfo(const XPrimaryCellInfo &info) 
                 :XSampleInfo(info),fCellType(info.fCellType),
                  fDescription(info.fDescription)
{
   // Primary cell information copy constructor
   if(kCS) cout << "---XPrimaryCellInfo::XPrimaryCellInfo(copy)------" << endl;

   fIsolationDate = info.fIsolationDate;
}//CopyConstructor

//______________________________________________________________________________
XPrimaryCellInfo& XPrimaryCellInfo::operator=(const XPrimaryCellInfo& rhs)
{
   // Primary cell information assignment operator.
   if(kCS) cout << "---XPrimaryCellInfo::operator=------" << endl;

   if (this != &rhs) {
      XSampleInfo::operator=(rhs);
      fCellType      = rhs.fCellType;
      fIsolationDate = rhs.fIsolationDate;
      fDescription   = rhs.fDescription;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XPrimaryCellInfo::~XPrimaryCellInfo()
{
   // Primary cell information destructor
   if(kCS) cout << "---XPrimaryCellInfo::~XPrimaryCellInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTissueInfo                                                          //
//                                                                      //
// Class describing sample of type tissue                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XTissueInfo::XTissueInfo()
            :XSampleInfo()
{
   // Default tissue information constructor
   if(kCS) cout << "---XTissueInfo::XTissueInfo(default)------" << endl;

   fTissueType   = "";
   fDevelopment  = "";
   fMorphology   = "";
   fDisease      = "";
   fDiseaseStage = "";
   fDonorAge     = 0;
   fAgeUnits     = "";
   fDonorStatus  = "";
}//Constructor

//______________________________________________________________________________
XTissueInfo::XTissueInfo(const char *name, const char *type, const char *development, 
             const char *morphology, const char *disease, const char *stage,
             Double_t donorage, const char *ageunits, const char *status,
             const char *sex, const char *pheno, const char *geno, const char *extract,
             const char *comment, Bool_t isXeno, const char *xenostrain,
             const char *xenosex, Double_t xenoage, const char *xageunits)
            :XSampleInfo(name, sex, pheno, geno, extract, comment,
             isXeno, xenostrain, xenosex, xenoage, xageunits)
{
   // Tissue information constructor
   if(kCS) cout << "---XTissueInfo::XTissueInfo------" << endl;

   fTissueType   = type;
   fDevelopment  = development;
   fMorphology   = morphology;
   fDisease      = disease;
   fDiseaseStage = stage;
   fDonorAge     = donorage;
   fAgeUnits     = ageunits;
   fDonorStatus  = status;
}//Constructor

//______________________________________________________________________________
XTissueInfo::XTissueInfo(const XTissueInfo &info) 
            :XSampleInfo(info),fTissueType(info.fTissueType),
             fDevelopment(info.fDevelopment),fMorphology(info.fMorphology),
             fDisease(info.fDisease),fDiseaseStage(info.fDiseaseStage),
             fAgeUnits(info.fAgeUnits),fDonorStatus(info.fDonorStatus)
{
   // Tissue information copy constructor
   if(kCS) cout << "---XTissueInfo::XTissueInfo(copy)------" << endl;

   fDonorAge = info.fDonorAge;
}//CopyConstructor

//______________________________________________________________________________
XTissueInfo& XTissueInfo::operator=(const XTissueInfo& rhs)
{
   // Tissue information assignment operator.
   if(kCS) cout << "---XTissueInfo::operator=------" << endl;

   if (this != &rhs) {
      XSampleInfo::operator=(rhs);
      fTissueType   = rhs.fTissueType;
      fDevelopment  = rhs.fDevelopment;
      fMorphology   = rhs.fMorphology;
      fDisease      = rhs.fDisease;
      fDiseaseStage = rhs.fDiseaseStage;
      fDonorAge     = rhs.fDonorAge;
      fAgeUnits     = rhs.fAgeUnits;
      fDonorStatus  = rhs.fDonorStatus;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XTissueInfo::~XTissueInfo()
{
   // Tissue information destructor
   if(kCS) cout << "---XTissueInfo::~XTissueInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBiopsyInfo                                                          //
//                                                                      //
// Class describing sample of type biopsy                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XBiopsyInfo::XBiopsyInfo()
            :XSampleInfo()
{
   // Default biopsy information constructor
   if(kCS) cout << "---XBiopsyInfo::XBiopsyInfo(default)------" << endl;

   fBiopsyType   = "";
   fMorphology   = "";
   fDisease      = "";
   fDiseaseStage = "";
   fDonorAge     = 0;
   fAgeUnits     = "";
   fDonorStatus  = "";
}//Constructor

//______________________________________________________________________________
XBiopsyInfo::XBiopsyInfo(const char *name, const char *type,
             const char *morphology, const char *disease, const char *stage,
             Double_t donorage, const char *ageunits, const char *status,
             const char *sex, const char *pheno, const char *geno, const char *extract,
             const char *comment, Bool_t isXeno, const char *xenostrain,
             const char *xenosex, Double_t xenoage, const char *xageunits)
            :XSampleInfo(name, sex, pheno, geno, extract, comment,
             isXeno, xenostrain, xenosex, xenoage, xageunits)
{
   // Biopsy information constructor
   if(kCS) cout << "---XBiopsyInfo::XBiopsyInfo------" << endl;

   fBiopsyType   = type;
   fMorphology   = morphology;
   fDisease      = disease;
   fDiseaseStage = stage;
   fDonorAge     = donorage;
   fAgeUnits     = ageunits;
   fDonorStatus  = status;
}//Constructor

//______________________________________________________________________________
XBiopsyInfo::XBiopsyInfo(const XBiopsyInfo &info) 
            :XSampleInfo(info),fBiopsyType(info.fBiopsyType),
             fMorphology(info.fMorphology),fDisease(info.fDisease),
             fDiseaseStage(info.fDiseaseStage),fAgeUnits(info.fAgeUnits),
             fDonorStatus(info.fDonorStatus)
{
   // Biopsy information copy constructor
   if(kCS) cout << "---XBiopsyInfo::XBiopsyInfo(copy)------" << endl;

   fDonorAge = info.fDonorAge;
}//CopyConstructor

//______________________________________________________________________________
XBiopsyInfo& XBiopsyInfo::operator=(const XBiopsyInfo& rhs)
{
   // Biopsy information assignment operator.
   if(kCS) cout << "---XBiopsyInfo::operator=------" << endl;

   if (this != &rhs) {
      XSampleInfo::operator=(rhs);
      fBiopsyType   = rhs.fBiopsyType;
      fMorphology   = rhs.fMorphology;
      fDisease      = rhs.fDisease;
      fDiseaseStage = rhs.fDiseaseStage;
      fDonorAge     = rhs.fDonorAge;
      fAgeUnits     = rhs.fAgeUnits;
      fDonorStatus  = rhs.fDonorStatus;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XBiopsyInfo::~XBiopsyInfo()
{
   // Biopsy information destructor
   if(kCS) cout << "---XBiopsyInfo::~XBiopsyInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreatment                                                           //
//                                                                      //
// Class describing a single treatment                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XTreatment::XTreatment()
           :TNamed()
{
   // Default treatment constructor
   if(kCS) cout << "---XTreatment::XTreatment(default)------" << endl;

   fType           = "";
   fConcentration  = 0;
   fConcUnit       = "";
   fTime           = 0;
   fTimeUnit       = "";
   fAdministration = "";
}//Constructor

//______________________________________________________________________________
XTreatment::XTreatment(const char *name, const char *type, Double_t conc,
            const char *concunit, Double_t time, const char *timeunit,
            const char *admin)
           :TNamed("Treatment", name)
{
   // Treatment constructor
   if(kCS) cout << "---XTreatment::XTreatment------" << endl;

   fType           = type;
   fConcentration  = conc;
   fConcUnit       = concunit;
   fTime           = time;
   fTimeUnit       = timeunit;
   fAdministration = admin;
}//Constructor

//______________________________________________________________________________
XTreatment::XTreatment(const XTreatment &info) 
           :TNamed(info),fType(info.fType),fConcUnit(info.fConcUnit),
            fTimeUnit(info.fTimeUnit),fAdministration(info.fAdministration)
{
   // Treatment copy constructor
   if(kCS) cout << "---XTreatment::XTreatment(copy)------" << endl;

   fConcentration = info.fConcentration;
   fTime          = info.fTime;
}//CopyConstructor

//______________________________________________________________________________
XTreatment& XTreatment::operator=(const XTreatment& rhs)
{
   // Treatment assignment operator.
   if(kCS) cout << "---XTreatment::operator=------" << endl;

   if (this != &rhs) {
      TNamed::operator=(rhs);
      fType           = rhs.fType;
      fConcentration  = rhs.fConcentration;
      fConcUnit       = rhs.fConcUnit;
      fTime           = rhs.fTime;
      fTimeUnit       = rhs.fTimeUnit;
      fAdministration = rhs.fAdministration;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XTreatment::~XTreatment()
{
   // Treatment destructor
   if(kCS) cout << "---XTreatment::~XTreatment------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreatmentInfo                                                       //
//                                                                      //
// Class describing sample treatments                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XTreatmentInfo::XTreatmentInfo()
               :XDataTypeInfo()
{
   // Default treatment information constructor
   if(kCS) cout << "---XTreatmentInfo::XTreatmentInfo(default)------" << endl;

   fTreatments = 0;
   fComment    = "";
}//Constructor

//______________________________________________________________________________
XTreatmentInfo::XTreatmentInfo(const char *type, const char *comment)
               :XDataTypeInfo("Treatments", type)
{
   // Treatment information constructor
   if(kCS) cout << "---XTreatmentInfo::XTreatmentInfo------" << endl;

   fTreatments = new TList();
   fComment    = comment;
   fHasInfo    = kTRUE;
}//Constructor

//______________________________________________________________________________
XTreatmentInfo::XTreatmentInfo(const XTreatmentInfo &info) 
               :XDataTypeInfo(info),fComment(info.fComment)
{
   // Treatment information copy constructor
   if(kCS) cout << "---XTreatmentInfo::XTreatmentInfo(copy)------" << endl;

//to do: test if correct?
   fTreatments = 0;
   if (info.fTreatments != 0) {
      fTreatments = new TList();
      for (Int_t i=0; i<(info.fTreatments)->GetSize(); i++) {
         fTreatments->AddAt((info.fTreatments)->At(i), i);
      }//for_i
   }//if
}//CopyConstructor

//______________________________________________________________________________
XTreatmentInfo& XTreatmentInfo::operator=(const XTreatmentInfo& rhs)
{
   // Treatment information assignment operator.
   if(kCS) cout << "---XTreatmentInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fComment = rhs.fComment;
//to do: need to copy fTreatments!!!
   cout << "Error: Copy of fTreatments not yet implemented." << endl;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XTreatmentInfo::~XTreatmentInfo()
{
   // Treatment information destructor
   if(kCS) cout << "---XTreatmentInfo::~XTreatmentInfo------" << endl;

   if(fTreatments) {fTreatments->Delete(); delete fTreatments; fTreatments = 0;}
}//Destructor

//______________________________________________________________________________
void XTreatmentInfo::AddTreatment(XTreatment *treat)
{
   // Add treatment to treatment list
   if(kCS) cout << "------XTreatmentInfo::AddTreatment------" << endl;

//??
   fTreatments->Add(treat);
}//AddTreatment

//______________________________________________________________________________
Int_t XTreatmentInfo::RemoveTreatment(const char *name)
{
   // Remove treatment name from treatment list 
   if(kCS) cout << "------XTreatmentInfo::RemoveTreatment------" << endl;

   Int_t size = fTreatments->GetSize();
   if (size == 0) return 0;

// Loop over treatment list
   TIter next(fTreatments);
   XTreatment *treat = 0;
   while ((treat = (XTreatment*)next())) {
      TString tname = treat->GetTreatment();
      if (strcmp(name,  tname.Data())  == 0) {
         fTreatments->Remove(treat);
         size--;
      }//if
   }//while

   return size;
}//RemoveTreatment

//______________________________________________________________________________
XTreatment *XTreatmentInfo::GetTreatment(const char *name)
{
   // Find treatment with name stored in list fTreatments
   if(kCS) cout << "------XTreatmentInfo::GetTreatment------" << endl;

   XTreatment *treat = 0;
   treat = (XTreatment*)fTreatments->FindObject(name);
   if (treat) {
      return treat;
   } else {
      return 0;
   }//if
}//GetTreatment

