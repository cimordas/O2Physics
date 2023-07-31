// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \brief  Calculation methods for the AC analysis.
/// \author Cindy Mordasini (cindy.mordasini@cern.ch)
/// \since  July 2023

// Includes.
#include <TMath.h>

#include "Framework/StaticFor.h"
#include "PWGCF/JCorran/Core/JAC2hAnalysis.h"

using namespace o2;
using namespace o2::framework;
namespace o2::analysis::PWGCF
{
/// \brief Initialize to defaults values the arrays used as data members.
/// .
void JAC2hAnalysis::InitArrays()
{
  for (int iH = 0; iH < mNharmos; iH++) {
    for (int iP = 0; iP < mNpowers; iP++) {mQvectors[iH][iP] = TComplex(0.,0.);}
  }

  int powers2h[14][2] = {{1,1}, {2,2}, {2,0}, {2,1}, {3,0}, {3,1}, {4,0}, {4,1},
                                      {0,2}, {1,2}, {0,3}, {1,3}, {0,4}, {1,4}};
  for (int iP = 0; iP < 14; iP++) {
    for (int iH = 0; iH < 2; iH++) {m2hPowers[iP][iH] = powers2h[iP][iH];}
  }

  int harmo2h[2] = {2,3};
  for (int iH = 0; iH < 2; iH++) {
    m2hPairs[iH] = harmo2h[iH];
  }
}

// ------------------------------------------------------------------------- //
/// \brief Calculate the weighted Q-vectors for the given azimuthal angles.
/// .
void JAC2hAnalysis::CalculateQvectors(int nPart,
                    double myAngles[], double myWeights[])
{
  // Ensure first all Q-vectors are initially zero for this event.
  for (int iH = 0; iH < mNharmos; iH++) {
    for (int iP = 0; iP < mNpowers; iP++) {mQvectors[iH][iP] = 0.;}
  }

  // Fill then Q(n,p) for all pairs of n and p.
  for (int iPart = 0; iPart < nPart; iPart++) {
    double iAngle = myAngles[iPart];
    double iWeight = myWeights[iPart];

    for (int iH = 0; iH < mNharmos; iH++) {
      for (int iP = 0; iP < mNpowers; iP++) {
        double weightToP = TMath::Power(iWeight, iP);
        mQvectors[iH][iP] += TComplex(weightToP*TMath::Cos(iH*iAngle),
                                      weightToP*TMath::Sin(iH*iAngle));
      }
    }
  }
  if (mDebugLevel > 1) {printf("All Q-vectors calculated.\n");}
}

/// \brief Simplify the use of Q-vectors with Q(-n,p)=Q(n,p)*.
/// . 
TComplex JAC2hAnalysis::Q(int harmoN, int power)
{
// Simplify the usage of the Q-vectors using that Q(-n,p)=Q(n,p)*
  if (harmoN >= 0) {return mQvectors[harmoN][power];}
  return TComplex::Conjugate(mQvectors[-harmoN][power]);
}

/// \brief Calculate multi-particle correlators using recursion.
/// Calculate multi-particle correlators using recursion (an improved faster
/// version) originally developed by Kristjan Gulbrandsen (gulbrand@nbi.dk).
TComplex JAC2hAnalysis::CalculateRecursion(int n, int *harmonic,
                        int mult, int skip)
{
  Int_t nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= CalculateRecursion(nm1, harmonic);
  if (nm1 == skip) return c;

  Int_t multp1 = mult+1;
  Int_t nm2 = n-2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(CalculateRecursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n-3;

  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += CalculateRecursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;
}

// ------------------------------------------------------------------------- //
/// \brief Compute and fill the profiles with all the needed correlation terms.
void JAC2hAnalysis::ComputeAllTerms(const int myCent, const int mySample)
{
  if (mDebugLevel > 1) {printf("myCent: %d mySample: %d\n", myCent, mySample);}

  // Start first by calculating the denominators (= event weights) for each
  // number of particles of interest.
  double eventWeights[5] = {0.};  // Values of each event weight in this event.
  int denom2h[2] = {0};
  int denom4h[4] = {0};
  int denom6h[6] = {0};
  int denom8h[8] = {0};
  int denom10h[10] = {0};

  eventWeights[0] = (CalculateRecursion(2, denom2h)).Re();
  eventWeights[1] = (CalculateRecursion(4, denom4h)).Re();
  eventWeights[2] = (CalculateRecursion(6, denom6h)).Re();
  eventWeights[3] = (CalculateRecursion(8, denom8h)).Re();
  eventWeights[4] = (CalculateRecursion(10, denom10h)).Re();

  // Calculate the correlators of interest by first calculating the numerators
  // and then the real part of the correlators using the denominators above.
  // We start first with the two-particle case.
  int nPartCorrel = 2;  // Number of particles in the correlators.
  int iEW = (nPartCorrel/2)-1;  // Index of the event weight to use.
  int list_harmos[7] = {0};   // We work with symmetric correlators of max 14
                              // particles --> max 7 harmonics.
  double correl2p[6] = {0.};     // Value of the correlation terms to use in profiles.

  for (int iVn = 1; iVn <= 6; iVn++) {  // We go to v6 maximum.
    list_harmos[0] = iVn;   // Only the first element is filled in this case.
    correl2p[iVn-1] = ComputeCorrel(nPartCorrel, eventWeights[iEW], list_harmos);
    if (mDebugLevel > 1) {printf("correl2p: %e.\n", correl2p[iVn-1]);}
  } // Go to the next harmonic with 2-particles.

  // All 2-particle terms have been obtained, we pass to the 2-harmonic cases,
  // following 'm2hPowers'. Each element is a bin of the 2hProfile.
  // LOKI: Add here the gestion of the different pairs of harmonics.
  double correl2h[14] = {0.}; // Values of the correlators to fill in profiles.
  double eWeight2h[14] = {0.}; // Corresponding event weights.
  for (int iPow = 0; iPow < 14; iPow++) {   // Dimensions should match f2hPowers.
    // Reset the number of correlated particles and index in harmonic array (0-6)
    // for this correlator.
    nPartCorrel = 0; iEW = 0;
    int cHarmo = 0;   //

    // Update the power list for the current values only if they are non-zero.
    for (int iHarmo = 0; iHarmo < 2; iHarmo++) {
      if (m2hPowers[iPow][iHarmo] == 0) {continue;}   // Skip the unneeded harmonics.
      for (int jP = 0; jP < m2hPowers[iPow][iHarmo]; jP++) {
        list_harmos[cHarmo] = m2hPairs[iHarmo]; // Write the harmonic as many times as its power is.
        cHarmo++;
      }

      // Determine the number of particles in the correlator. 2-harmonic terms in AC have
      // twice as many particles as their cumulant order.
      nPartCorrel += 2*m2hPowers[iPow][iHarmo];
    } // Go to the second harmonic of the pair.

    iEW = (nPartCorrel/2)-1;  // Recalculate the index for the event weight.

    // Compute the correlator for this harmonic array and save it as before
    // in the corresponding profiles and branches, as well as in the full profiles.
    eWeight2h[iPow] = eventWeights[iEW];
    correl2h[iPow] = ComputeCorrel(nPartCorrel, eventWeights[iEW], list_harmos);
    //myHistManager->prof2hCorrel[myCent][mySample]->Fill((float)iPow+0.5, correlTerm, eventWeights[iEW]);
    //myHistManager->prof2hCorrel[myCent][mySample]->GetXaxis()->SetBinLabel(iPow+1, Form("{%d,%d}", list_powers[0], list_powers[1]));
    //myHistManager->prof2hCorrel[myCent][nSamples]->Fill((float)iPow+0.5, correlTerm, eventWeights[iEW]);
    //myHistManager->prof2hCorrel[myCent][nSamples]->GetXaxis()->SetBinLabel(iPow+1, Form("{%d,%d}", list_powers[0], list_powers[1]));
  } // Go to the next power for this specific pair of harmonics.

  // Save the correlators and associated event weights in the given outputs.
  // First the profile for the sample and then for the full dataset.
  // The static_for loop assumes we have maximum 20 samples.
  // TODO: Find a way to quit the loop once we finished with the sample.
  static_for<0, 19>([&](auto iS) {
    constexpr int indexS = iS.value;
    if (indexS == mySample) {
      printf("We have a match at indexS = %d\n", indexS);

      switch (myCent) {
      case 0:
        mAC2hHistManager.Fill2pProf<0,indexS>(correl2p, eventWeights[iEW]);
        mAC2hHistManager.Fill2hProf<0,indexS,0>(correl2h, eWeight2h);
        break;
      case 1:
        mAC2hHistManager.Fill2pProf<1,indexS>(correl2p, eventWeights[iEW]);
        mAC2hHistManager.Fill2hProf<1,indexS,0>(correl2h, eWeight2h);
        break;
      case 2:
        mAC2hHistManager.Fill2pProf<2,indexS>(correl2p, eventWeights[iEW]);
        mAC2hHistManager.Fill2hProf<2,indexS,0>(correl2h, eWeight2h);
        break;
      case 3:
        mAC2hHistManager.Fill2pProf<3,indexS>(correl2p, eventWeights[iEW]);
        mAC2hHistManager.Fill2hProf<3,indexS,0>(correl2h, eWeight2h);
        break;
      case 4:
        mAC2hHistManager.Fill2pProf<4,indexS>(correl2p, eventWeights[iEW]);
        mAC2hHistManager.Fill2hProf<4,indexS,0>(correl2h, eWeight2h);
        break;
      case 5:
        mAC2hHistManager.Fill2pProf<5,indexS>(correl2p, eventWeights[iEW]);
        mAC2hHistManager.Fill2hProf<5,indexS,0>(correl2h, eWeight2h);
        break;
      case 6:
        mAC2hHistManager.Fill2pProf<6,indexS>(correl2p, eventWeights[iEW]);
        mAC2hHistManager.Fill2hProf<6,indexS,0>(correl2h, eWeight2h);
        break;
      case 7:
        mAC2hHistManager.Fill2pProf<7,indexS>(correl2p, eventWeights[iEW]);
        mAC2hHistManager.Fill2hProf<7,indexS,0>(correl2h, eWeight2h);
        break;
      case 8:
        mAC2hHistManager.Fill2pProf<8,indexS>(correl2p, eventWeights[iEW]);
        mAC2hHistManager.Fill2hProf<8,indexS,0>(correl2h, eWeight2h);
        break;
      }
    }
  });

  if (mDebugLevel > 1) {printf("Terms obtained for all pairs of interest.\n");}
}

/// \brief Compute the multiparticle correlator term.
/// \param nPC: Order of the correlator (i.e number of particles).
/// \return The correlator since the event weight is already calculated.
double JAC2hAnalysis::ComputeCorrel(int nPC, const double myEventWeight,
                                    const int (&myHarmos)[7])
{
  if (mDebugLevel > 1) {
    printf("Computing correlation term with %d particles.\n", nPC);
  }

  // Initialize all possible cases of numerator arrays with the correct signs.
  int numer2h[2] = {myHarmos[0], -1*myHarmos[0]};
  int numer4h[4] = {myHarmos[0], myHarmos[1],
    -1*myHarmos[0], -1*myHarmos[1]};
  int numer6h[6] = {myHarmos[0], myHarmos[1], myHarmos[2],
    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2]};
  int numer8h[8] = {myHarmos[0], myHarmos[1], myHarmos[2], myHarmos[3],
    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2], -1*myHarmos[3]};
  int numer10h[10] = {myHarmos[0], myHarmos[1], myHarmos[2], myHarmos[3],
    myHarmos[4], -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2],
    -1*myHarmos[3], -1*myHarmos[4]};

  // Calculate the numerator and final correlator for each possible case of nPC.
  TComplex cCorrel = TComplex(0.,0.);   // Complex form of the correlator.
  switch(nPC) {
    case 2:
      cCorrel = (CalculateRecursion(2, numer2h))/myEventWeight;
      break;
    case 4:
      cCorrel = (CalculateRecursion(4, numer4h))/myEventWeight;
      break;
    case 6:
      cCorrel = (CalculateRecursion(6, numer6h))/myEventWeight;
      break;
    case 8:
      cCorrel = (CalculateRecursion(8, numer8h))/myEventWeight;
      break;
    case 10:
      cCorrel = (CalculateRecursion(10, numer10h))/myEventWeight;
      break;
    default :
      printf("Error: invalid number of particles. Doing nothing...\n");
      break;
  }
  return cCorrel.Re();  // Real part of the correlator of interest to be used.
}


} // namespace o2::analysis::PWGCF