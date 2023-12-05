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
    float iAngle = myAngles.at(iT);
    float iWeight = myWeights.at(iT);

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
/// \param myCentBin Bin index for the centrality class of the collision.
/// \param mySample Sample index attributed to this collision.
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

  for (int iVn = 0; iVn < 6; iVn++) {  // We go to v6 maximum.
    list_harmos[0] = iVn+1;   // Only the first element is filled in this case.
    correl2p.at(iVn) = ComputeThisCorrel(2, eventWeights.at(0), list_harmos);
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
    std::array<double, 14> thisCorrel;  // Values of the correlators for each power.
    std::array<double, 14> thisWeight;  // Associated event weights.

    for (int iPow = 0; iPow < 14; iPow++) { // Dimension should match m2hPowers.
      /* Reset the variables for this new "power" term. */
      int cHarmo = 0;   // Counter for the filling of list_harmos.
      int nPartCorrel = 0;  // Number of particles in the current correlator.
      int iEW = 0;  // Associated event weight in the array.

      /* Update the list of harmonics based on the current power. */
      // For that, we write each harmonic as many times as its power is.
      for (int iHarmo = 0; iHarmo < 2; iHarmo++) {
        if (m2hPowers.at(iPow).at(iHarmo) == 0) {continue;} // Skip the powers set to zero.

        for (int jPow = 0; jPow < m2hPowers.at(iPow).at(iHarmo); jPow++) {
          list_harmos[cHarmo] = harmoPair.at(iHarmo);
          cHarmo++;
        }
        ///printf("m2hPowers.at(%d).at(%d): %d\n", iPow, iHarmo, m2hPowers.at(iPow).at(iHarmo));

        /* Determine the number of particles required for the correlator. */
        nPartCorrel += 2*m2hPowers.at(iPow).at(iHarmo);

      }
      /* Get the right event weight using 'iEW'. */
      iEW = (nPartCorrel/2)-1;
      thisWeight.at(iPow) = eventWeights.at(iEW);

      /* Compute the correlator for these pairs of harmonics and powers. */
      thisCorrel.at(iPow) = ComputeThisCorrel(nPartCorrel, eventWeights.at(iEW),
                                              list_harmos);
    }   // Go to the next "power" term.

      /* Save the correlators and weights combined for this pair of harmonics. */
      correl2h.push_back(thisCorrel);
      eWeight2h.push_back(thisWeight);

      hIndex++; // Raise the index for the next combination of harmonics.
      thisCorrel.fill(0); // Reset.
      thisWeight.fill(0); // Reset.
  }   // Go to the next pair of harmonics.

  /* Fill the sample distribution and 2D profiles for this collision. */
  // TODO: Find a less ugly way to go that.
  switch (myCentBin) {
  case 0:
    mAC2hHistManager.FillHistSamples<0>(mySample);
    mAC2hHistManager.FillProf2p<0>(correl2p, eventWeights.at(0), mySample, m2pCorrelEtaGap);
    mAC2hHistManager.FillProf2h<0,0>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<0,1>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<0,2>(correl2h, eWeight2h, mySample);
    break;
  case 1:
    mAC2hHistManager.FillHistSamples<1>(mySample);
    mAC2hHistManager.FillProf2p<1>(correl2p, eventWeights.at(0), mySample, m2pCorrelEtaGap);
    mAC2hHistManager.FillProf2h<1,0>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<1,1>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<1,2>(correl2h, eWeight2h, mySample);
    break;
  case 2:
    mAC2hHistManager.FillHistSamples<2>(mySample);
    mAC2hHistManager.FillProf2p<2>(correl2p, eventWeights.at(0), mySample, m2pCorrelEtaGap);
    mAC2hHistManager.FillProf2h<2,0>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<2,1>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<2,2>(correl2h, eWeight2h, mySample);
    break;
  case 3:
    mAC2hHistManager.FillHistSamples<3>(mySample);
    mAC2hHistManager.FillProf2p<3>(correl2p, eventWeights.at(0), mySample, m2pCorrelEtaGap);
    mAC2hHistManager.FillProf2h<3,0>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<3,1>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<3,2>(correl2h, eWeight2h, mySample);
    break;
  case 4:
    mAC2hHistManager.FillHistSamples<4>(mySample);
    mAC2hHistManager.FillProf2p<4>(correl2p, eventWeights.at(0), mySample, m2pCorrelEtaGap);
    mAC2hHistManager.FillProf2h<4,0>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<4,1>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<4,2>(correl2h, eWeight2h, mySample);
    break;
  case 5:
    mAC2hHistManager.FillHistSamples<5>(mySample);
    mAC2hHistManager.FillProf2p<5>(correl2p, eventWeights.at(0), mySample, m2pCorrelEtaGap);
    mAC2hHistManager.FillProf2h<5,0>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<5,1>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<5,2>(correl2h, eWeight2h, mySample);
    break;
  case 6:
    mAC2hHistManager.FillHistSamples<6>(mySample);
    mAC2hHistManager.FillProf2p<6>(correl2p, eventWeights.at(0), mySample, m2pCorrelEtaGap);
    mAC2hHistManager.FillProf2h<6,0>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<6,1>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<6,2>(correl2h, eWeight2h, mySample);
    break;
  case 7:
    mAC2hHistManager.FillHistSamples<7>(mySample);
    mAC2hHistManager.FillProf2p<7>(correl2p, eventWeights.at(0), mySample, m2pCorrelEtaGap);
    mAC2hHistManager.FillProf2h<7,0>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<7,1>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<7,2>(correl2h, eWeight2h, mySample);
    break;
  case 8:
    mAC2hHistManager.FillHistSamples<8>(mySample);
    mAC2hHistManager.FillProf2p<8>(correl2p, eventWeights.at(0), mySample, m2pCorrelEtaGap);
    mAC2hHistManager.FillProf2h<8,0>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<8,1>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<8,2>(correl2h, eWeight2h, mySample);
    break;
  case 9:
    mAC2hHistManager.FillHistSamples<9>(mySample);
    mAC2hHistManager.FillProf2p<9>(correl2p, eventWeights.at(0), mySample, m2pCorrelEtaGap);
    mAC2hHistManager.FillProf2h<9,0>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<9,1>(correl2h, eWeight2h, mySample);
    mAC2hHistManager.FillProf2h<9,2>(correl2h, eWeight2h, mySample);
    break;
  }

  /* Reset the arrays for the next collision. */
  correl2p.fill(0);
  eventWeights.fill(0);
  if (mDebugPrint) {LOGF(info, "Calculation of all terms done for this collision.");}
}

/// \brief Compute the multiparticle correlator for the set of harmonics.
/// \param nPC Number of particles in the correlator.
/// \param myEventWeight Value of the event weight (= denominator).
/// \param myHarmos Array of harmonics for the numerator.
/// \return Real value of the multiparticle correlator.
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


void FlowAC2hAnalysis::ComputeCorrelEtaGap(const std::vector<float>& myAngles,
                                             const std::vector<float>& myWeights,
                                             const std::vector<float>& myEtas)
{
  if (mDebugPrint) {
    LOGF(info, "Computing 2-particle correlators with eta gap.");
  }

  /* Define all the arrays for the quantities of interest. */
  std::array<float, 6> multiNeg;  // Number of tracks in the negative subevent.
  std::array<float, 6> multiPos;  // Nymber of tracks in the positive subevent.
  std::array<TComplex, 6> QnNeg;  // Q_n in the negative subevent (n = 1..6).
  std::array<TComplex, 6> QnPos;  // Q_n in the positive subevent (n = 1..6).

  for (std::size_t iT = 0; iT < myAngles.size(); iT++) {
    /* Get the information on the track itself. */
    // NOTE: myWeights contain NUE*NUA values already, and is a unit vector 
    // if the analysis does not use particle weight.
    float iAngle = myAngles.at(iT);
    float iWeight = myWeights.at(iT);
    float iEta = myEtas.at(iT);
    float weightToP = iWeight;  // Particle weight raised to the power 'p' = 1.
    
    /* Calculate Q_n and add it to the correct subevent. */
    if (iEta < (-0.5*mEtaGap)) {   // Track is in the negative/left subevent.
      for (int iH = 0; iH < 6; iH++) {
        QnNeg.at(iH) += TComplex(weightToP*TMath::Cos((iH+1)*iAngle),
                                weightToP*TMath::Sin((iH+1)*iAngle));
        multiNeg.at(iH) += weightToP;
      }
    }
    else if (iEta > (0.5*mEtaGap)) {  // Track is in the positive/right subevent.
      for (int iH = 0; iH < 6; iH++) {
        QnPos.at(iH) += TComplex(weightToP*TMath::Cos((iH+1)*iAngle),
                                weightToP*TMath::Sin((iH+1)*iAngle));
        multiPos.at(iH) += weightToP;
      }
    }
    else {continue;}  // Track is between the two subevents.
  }   // Go to the next track.

  /* Calculate the 2-particle correlators. */
  // The calculation is done only if both positive and negative Q_n and numbers
  // of tracks are strictly positive.
  for (int iH = 0; iH < 6; iH++) {
    if (!((QnNeg[iH].TComplex::Rho() > 0.) && (QnPos[iH].TComplex::Rho() > 0.))) {
      continue;
    }
    if (!((multiNeg[iH] > 0.) && (multiPos[iH] > 0.))) {continue;}

    TComplex cCorrel = QnNeg[iH]*TComplex::Conjugate(QnPos[iH]);
    m2pCorrelEtaGap[iH][0] = (cCorrel.Re())/(multiNeg[iH]*multiPos[iH]);  // 2-particle correlator.
    m2pCorrelEtaGap[iH][1] = multiNeg[iH]*multiPos[iH];   // Associated weight for the profile2D.
  }   // Go to the next harmonic.

}

} // namespace o2::analysis::PWGCF
