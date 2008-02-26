// Author: Christian Stratowa 11/25/2002             last modified: 11/05/2007

/*
 *******************************************************************************
 ***********************  Statistics Package for ROOT  *************************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2008 Dr. Christian Stratowa
 *
 *  Written by: Christian Stratowa, Vienna, Austria <cstrato@aon.at>
 *
 *******************************************************************************
 * Algorithms for statistical functions partially adapted from:                *
 * R: A Computer Language for Statistical Data Analysis                        *
 * Copyright (C) 1995, 1996 Robert Gentleman and Ross Ihaka                    *
 * Copyright (C) 1999-2001 The R Development Core Team                         *
 * R is free software, for licensing see the GNU General Public License        *
 * http:/cran.r-project.org/                                                   *
 *******************************************************************************
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

#ifndef __TStat__
#define __TStat__

#include <Riostream.h>
#include <cmath>

// ROOT includes
#include "TMath.h"


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStat                                                                //
//                                                                      //
// Class containing statistical methods                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TStat {

   private:
      
   public :
      TStat() {}
      virtual ~TStat() {}

      static Bool_t    Ident(Int_t n, const Double_t *arr);
      static Double_t  Max(Int_t n, const Double_t *arr);
      static Double_t  Min(Int_t n, const Double_t *arr);
      static Double_t *PMax(Int_t n, const Double_t *arr1, const Double_t *arr2,
                          Double_t *max);
      static Double_t *PMax(Int_t n, Int_t m, const Double_t **table,
                          Double_t *max);
      static Double_t *PMin(Int_t n, const Double_t *arr1, const Double_t *arr2,
                          Double_t *min);
      static Double_t *PMin(Int_t n, Int_t m, const Double_t **table, 
                          Double_t *min);
      static Double_t *CumMax(Int_t n, const Double_t *arr, Double_t *max);
      static Double_t *CumMin(Int_t n, const Double_t *arr, Double_t *min);
      static Double_t  CumProd(Int_t n, const Double_t *arr);
      static Double_t  CumSum(Int_t n, const Double_t *arr);
      static Double_t  LnFact(Int_t n);
      static Double_t  BinomCoeff(Int_t n, Int_t k);
      static void      Rank(Int_t n, Double_t *arr, Int_t *index, Int_t *rank,
                          Bool_t down = kTRUE);
      static void      Rank(Int_t n, Double_t *arr, Double_t *rank);
      static Int_t    *TrueRank(Int_t n, Double_t *arr, Int_t *rank,
                          Bool_t down = kTRUE);
      static Double_t *TrueRank1(Int_t n, Double_t *arr, Int_t *index,
                          Double_t *rank);

      static Double_t GeoMean(Int_t n, const Double_t *arr);
      static Double_t GeoMean(Int_t n, const Double_t *arr, Int_t &len,
                         const Double_t na);
      static Double_t GeoMean(Int_t n, const Double_t *arr, const Double_t *w);
      static Double_t GeoMean(Int_t n, const Double_t *arr, const Double_t trim,
                         Double_t &var, Int_t &len);
      static Double_t Mean(Int_t n, const Double_t *arr);
      static Double_t Mean(Int_t n, const Double_t *arr, Int_t &len,
                         const Double_t na);
      static Double_t Mean(Int_t n, const Double_t *arr, const Double_t *w);
      static Double_t Mean(Int_t n, const Double_t *arr, const Double_t trim);
      static Double_t Mean(Int_t n, const Double_t *arr, const Double_t trim,
                         Double_t &var, Int_t &len);
      static Double_t Mean(Int_t n, const Double_t *arr, const Double_t *w,
                         const Double_t trim, Double_t &var, Int_t &len);
      static Double_t Median(Int_t n, const Double_t *arr);
      static Double_t MAD(Int_t n, const Double_t *arr, Float_t constant = 1.4826);
      static Double_t MAD(Int_t n, const Double_t *arr, const Double_t trim,
                         Float_t constant = 1.4826);
      static Double_t StDev(Int_t n, const Double_t *arr, const Double_t mean);
      static Double_t Var(Int_t n, const Double_t *arr, const Double_t mean);
      static Double_t Var(Int_t n, const Double_t *arr, const Double_t mean,
                         Int_t &len, const Double_t na);

      static Double_t MedianPolish(Int_t nrow, Int_t ncol, Double_t *x,
                         Double_t *rowmed, Double_t *colmed, Double_t *residu,
                         Int_t maxiter = 10, Double_t eps = 0.01,
                         Bool_t verbose = kFALSE);
      static Double_t TukeyBiweight(Int_t n, const Double_t *x, Double_t &var,
                         Double_t c = 5.0, Double_t eps = 0.0001);

      static Double_t Quantile(Int_t n, const Double_t *arr, const Double_t q);
      static Double_t IQR(Int_t n, const Double_t *arr, const Double_t qlo = 0.25,
                         const Double_t qhi = 0.75);

      static void     NextPerm(Int_t n, Int_t k, Int_t *grp);
      static void     NextPerm(Int_t n, Int_t k, Int_t *grp1, Int_t nk, Int_t *grp2);
      static void     Sample(Int_t n, Int_t k, Int_t *grp);
      static void     Sample(Int_t n, Int_t k, Int_t *grp1, Int_t nk, Int_t *grp2);

      static Double_t CorPearson(Int_t n, Double_t *x, Double_t *y);
      static Double_t CorSpearman(Int_t n, Double_t *x, Double_t *y);

      static void     MassDist(Int_t nx, Double_t *x, Double_t *w, Double_t xmin,
                         Double_t xmax, Int_t ny, Double_t *y);
      static void     TwiddleFactor4FFT(Int_t n, Int_t i, Double_t &tf_re, Double_t &tf_im);
      static void     TwiddleFactor4IFFT(Int_t n, Int_t i, Double_t &tf_re, Double_t &tf_im);
      static void     FFT(Int_t n, Double_t *f_re, Double_t *f_im);
      static void     IFFT(Int_t n, Double_t *f_re, Double_t *f_im);
      static Int_t    FFTDensityConvolve(Int_t n, Double_t *x_re, Double_t *y_re);
      static void     Kernelize(Int_t n, Double_t *x, Double_t bw, const char *kernel);
      static Double_t Bandwidth(Int_t n, Double_t *x, Double_t iqr);
      static void     LinearInterpolate(Double_t *xin, Double_t *yin,
                         Int_t nout, Double_t *xout, Double_t *yout);
      static Int_t    Density(Int_t n, Double_t *x, Double_t *w, Int_t nout,
                         Double_t *xout, Double_t *yout, const char *kernel = "gaussian");
      static Double_t MaxDensity(Int_t n, Double_t *x, Double_t *w, Int_t npts = 512,
                         const char *kernel = "gaussian");

      static Double_t PNormApprox(Double_t x);

      ClassDef(TStat,0) //Stat
};

//______________________________________________________________________________
inline Bool_t TStat::Ident(Int_t n, const Double_t *arr)
{
   // Return true if all data in array are equal

   Bool_t ident = kTRUE;
   for (Int_t i=0; i<n-1; i++) {
      if (arr[i] != arr[i+1]) {ident = kFALSE; break;}
   }//for

   return ident;
}//Max

//______________________________________________________________________________
inline Double_t TStat::Max(Int_t n, const Double_t *arr)
{
   // Maximum of array arr.

   return arr[TMath::LocMax(n,arr)];
}//Max

//______________________________________________________________________________
inline Double_t TStat::Min(Int_t n, const Double_t *arr)
{
   // Minimum of array arr.

   return arr[TMath::LocMin(n,arr)];
}//Min

//______________________________________________________________________________
inline Double_t *TStat::PMax(Int_t n, const Double_t *arr1, const Double_t *arr2,
                        Double_t *max)
{
  // Returns maximum of arr1 and arr2 at each index in array max.

   for (Int_t i=0; i<n; i++) max[i] = (arr1[i] > arr2[i]) ? arr1[i] : arr2[i];
   return max;
}//PMax

//______________________________________________________________________________
inline Double_t *TStat::PMin(Int_t n, const Double_t *arr1, const Double_t *arr2,
                        Double_t *min)
{
  // Returns minimum of arr1 and arr2 at each index in array min.

   for (Int_t i=0; i<n; i++) min[i] = (arr1[i] < arr2[i]) ? arr1[i] : arr2[i];
   return min;
}//PMin

//______________________________________________________________________________
inline Double_t *TStat::CumMax(Int_t n, const Double_t *arr, Double_t *max)
{
  // Returns cumulative maxima of arr in array max.

   max[0] = arr[0];
   for (Int_t i=1; i<n; i++) max[i] = (arr[i] > max[i-1]) ? arr[i] : max[i-1];
   return max;
}//CumMax

//______________________________________________________________________________
inline Double_t *TStat::CumMin(Int_t n, const Double_t *arr, Double_t *min)
{
  // Returns cumulative minima of arr in array min.

   min[0] = arr[0];
   for (Int_t i=1; i<n; i++) min[i] = (arr[i] < min[i-1]) ? arr[i] : min[i-1];
   return min;
}//CumMin

//______________________________________________________________________________
inline Double_t TStat::CumProd(Int_t n, const Double_t *arr)
{
   // Cumulative sum.

   Double_t prod = 0;
   for (Int_t i=0; i<n; i++) prod *= arr[i];
   return prod;
}//CumProd

//______________________________________________________________________________
inline Double_t TStat::CumSum(Int_t n, const Double_t *arr)
{
   // Cumulative sum.

   Double_t sum = 0;
   for (Int_t i=0; i<n; i++) sum += arr[i];
   return sum;
}//CumSum

//______________________________________________________________________________
inline Double_t TStat::LnFact(Int_t n)
{
   // ln(n!)

   if (n > 1) return TMath::LnGamma(n+1.0);
   else       return 0.0;
}//LnFact

//______________________________________________________________________________
inline Double_t TStat::BinomCoeff(Int_t n, Int_t k)
{
  // Binomial coefficient "n chooose k" as double.

  return floor(0.5 + exp(LnFact(n)-LnFact(k)-LnFact(n-k)));
}//CumProd

//______________________________________________________________________________
inline Double_t TStat::StDev(Int_t n, const Double_t *arr, const Double_t mean)
{
   // Calculate standard deviation for arithmetic mean

   return TMath::Sqrt(Var(n, arr, mean));
}//StDev


#endif
