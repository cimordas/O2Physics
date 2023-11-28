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

/* Header files. */

// O2 headers. //

// O2 Physics headers.
#include "PWGCF/Flow/Core/FlowAC2hAnalysis.h"

/* Namespaces. */
using namespace o2;
using namespace o2::framework;

namespace o2::analysis::PWGCF
{
/// \brief Decompose the provided integer into his two components.
/// \param inPair Integer representing the input pair of harmonics.
/// \return 2D array with the decomposed pair.
std::array<int, 2> FlowAC2hAnalysis::GetHarmoPair(int inPair)
{
  std::array<int, 2> pair;
  /* Security check to ensure proper behaviour. */
  if (inPair == 0) {
    pair.at(0) = 0; pair.at(1) = 0;
    return pair;
  }

  /* Extract each digit starting by the unit and save them in reverse. */
  pair.at(1) = inPair % 10;
  pair.at(0) = inPair/10;   // This removes the units as we deal with integers.
  return pair;
}

/// \brief Calculate the weighted Q-vectors for the vector of azimuthal angles.
/// \param myAngles Azimuthal angles of the tracks in the collision.
/// \param myWeights Associated particle weights.
void FlowAC2hAnalysis::CalculateQvectors(const std::vector<float>& myAngles,
                                         const std::vector<float>& myWeights)
{
  /* Security checks: all Q-vectors must initially be zero. */
  for (auto& Qn : mQvectors) {
    for (auto& Qnp : Qn) {
      Qnp = TComplex(0.,0.);
    }
  }

  /* Calculate the Q-vectors for the set of provided tracks. */
  for (std::size_t iT = 0; iT < myAngles.size(); iT++) {
    double iAngle = myAngles.at(iT);
    double iWeight = myWeights.at(iT);

    for (int iH = 0; iH < mNqHarmos; iH++) {
      for (int iP = 0; iP < mNqPowers; iP++) {
        double weightToP = TMath::Power(iWeight, iP);
        mQvectors[iH][iP] += TComplex(weightToP*TMath::Cos(iH*iAngle),
                                      weightToP*TMath::Sin(iH*iAngle));
      }
    }
  }

  if (mDebugPrint) {LOGF(info, "All Q-vectors obtained for this collision.");}
}

/// \brief Simplification of the use of Q-vectors from Q(-n,p) = Q(n,p)*.
/// \param harmoN Harmonic n of Q(n,p).
/// \param power Power p of Q(n,p).
/// \return Corresponding Q-vector value.
TComplex FlowAC2hAnalysis::Q(const int harmoN, const int power)
{
  if (harmoN >= 0) {return mQvectors[harmoN][power];}
  return TComplex::Conjugate(mQvectors[-harmoN][power]);
}

/// \brief Calculate multi-particle correlators using recursion.
/// \param n Number of particles in the correlator.
/// \param harmonic Array of length n of harmonics.
/// \return Complex value of the multiparticle correlator.
/// \note Improved faster version) originally developed by Kristjan Gulbrandsen
/// (gulbrand@nbi.dk).
TComplex FlowAC2hAnalysis::CalculateRecursion(int n, int *harmonic,
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

/// \brief Compute and fill the profiles with all the needed correlators.
//
void FlowAC2hAnalysis::ComputeAllCorrel(const int myCentBin, const int mySample)
{
  if (mDebugPrint) {LOGF(info, "myCentBin: %d mySample: %d", myCentBin, mySample);}

  /* Compute the denominators (= event weights) for each case of interest. */
  std::array<double, 5> eventWeights; // Event weights for all even-particle correlators.
  int denom2h[2] = {0};
  int denom4h[4] = {0};
  int denom6h[6] = {0};
  int denom8h[8] = {0};
  int denom10h[10] = {0};
  eventWeights.at(0) = (CalculateRecursion(2, denom2h)).Re();
  eventWeights.at(1) = (CalculateRecursion(4, denom4h)).Re();
  eventWeights.at(2) = (CalculateRecursion(6, denom6h)).Re();
  eventWeights.at(3) = (CalculateRecursion(8, denom8h)).Re();
  eventWeights.at(4) = (CalculateRecursion(10, denom10h)).Re();

  /* Compute the 2-particle correlators for v_1 to v_6 without eta gap. */
  int list_harmos[7] = {0};   // We work with symmetric correlators of max 14
                              // particles --> max 7 harmonics.
  std::array<double, 6> correl2p;   // Value of the correlation terms to use in profiles.

  for (int iVn = 1; iVn <= 6; iVn++) {  // We go to v6 maximum.
    list_harmos[0] = iVn;   // Only the first element is filled in this case.
    correl2p.at(iVn-1) = ComputeThisCorrel(2, eventWeights.at(0), list_harmos);
    ///printf("correl2p: %e.\n", correl2p.at(iVn-1));  // For prelim debug.
  } // Go to the next harmonic with 2-particles.

  /* All 2-particle terms have been obtained, we pass to the 2-harmonic case. */
  std::vector<std::array<double, 14>> correl2h;   // Correlations for all combinations.
  std::vector<std::array<double, 14>> eWeight2h;  // Corresponding event weights.
  int hIndex = 0;   // Index for the i-th harmonic pair.  // NOTE: Update to init-statement of c++20.
  for (auto hPair : m2hPairs) {
    /* Decompose first the pair of harmonics into its components. */
    std::array<int, 2> harmoPair = GetHarmoPair(hPair);
    mAC2hHistManager.FillPairProf(hIndex, harmoPair);
    LOGF(info, "Current harmoPair: [%d,%d]", harmoPair.at(0), harmoPair.at(1));

    /* For each power in m2hPowers, get the corresponding correlator. */
    for (int iPow = 0; iPow < 14; iPow++) { // Dimension should match m2hPowers.
      /* Reset the variables for this new "power" term. */
      int cHarmo = 0;   // Counter for the filling of list_harmos.

      /* Update the list of harmonics based on the current power. */
      // For that, we write each harmonic as many times as its power is.
      for (int iHarmo = 0; iHarmo < 2; iHarmo++) {
        if (m2hPowers.at(iPow).at(iHarmo) == 0) {continue;} // Skip the powers set to zero.

        for (int jPow = 0; jPow < m2hPowers.at(iPow).at(iHarmo); jPow++) {
          list_harmos[cHarmo] = harmoPair.at(iHarmo);
          cHarmo++;
        }

        printf("m2hPowers.at(%d).at(%d): %d\n", iPow, iHarmo, m2hPowers.at(iPow).at(iHarmo));


      }

    }   // Go to the next "power" term.


    /* Save the value of the correlator for this harmonic pair. */
    hIndex++;   // Raise the index for this harmonic pair.
  }   // Go to the next pair of harmonics.

  if (mDebugPrint) {LOGF(info, "Calculation of all terms done for this collision.");}
}

/// \brief Compute the multiparticle correlator for the set of harmonics.

double FlowAC2hAnalysis::ComputeThisCorrel(const int nPC, const double myEventWeight,
                                           const int (&myHarmos)[7])
{
  if (mDebugPrint) {
    LOGF(info, "Computing %d-particle correlator without eta gap", nPC);
  }

  /* Initialise all possible cases of numerator arrays with the proper signs. */
  int numer2h[2] = {myHarmos[0], -1*myHarmos[0]};
  int numer4h[4] = {myHarmos[0], myHarmos[1],
    -1*myHarmos[0], -1*myHarmos[1]};
  int numer6h[6] = {myHarmos[0], myHarmos[1], myHarmos[2],
    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2]};
  int numer8h[8] = {myHarmos[0], myHarmos[1], myHarmos[2], myHarmos[3],
    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2], -1*myHarmos[3]};
  int numer10h[10] = {myHarmos[0], myHarmos[1], myHarmos[2], myHarmos[3], myHarmos[4],
    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2],-1*myHarmos[3], -1*myHarmos[4]};

  /* Compute the complex form of the correlator term. */
  TComplex cCorrel = TComplex(0.,0.);
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
  default:
    printf("Error: invalid number of particles. Doing nothing...\n");
    break;
  }

  return cCorrel.Re();  // Real part of the correlator of interest.
}

} // namespace o2::analysis::PWGCF
