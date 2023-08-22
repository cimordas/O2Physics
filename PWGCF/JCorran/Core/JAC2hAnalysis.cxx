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

/* Header files. */
#include <vector>
#include <TMath.h>

// O2 headers.
#include "Framework/StaticFor.h"

// O2Physics headers.

// JCorran headers.
#include "PWGCF/JCorran/Core/JAC2hAnalysis.h"

/* Namespaces. */
using namespace o2;
using namespace o2::framework;
namespace o2::analysis::PWGCF
{
/// \brief Calculate the weighted Q-vectors for the given azimuthal angles.
/// \param myAngles Azimuthal angles of all the tracks in the collision.
/// \param myWeights Particle weights corresponding to the azimuthal angles.
void JAC2hAnalysis::CalculateQvectors(const std::vector<double>& myAngles,
        const std::vector<double>& myWeights)
{
    /* Ensure first the Q-vectors are initially set to 0. */
    for (int iH = 0; iH < mNqHarmos; iH++) {
        for (int iP = 0; iP < mNqPowers; iP++) {
            mQvectors[iH][iP] = TComplex(0.,0.);
        }
    }

    /* Calculate then the Q-vectors for the provided tracks. */
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
    if (mDebugLevel > 1) {printf("All Q-vectors calculated.\n");}
}

/// \brief Simplify the use of Q-vectors with Q(-n,p)=Q(n,p)*.
/// \param harmoN Harmonic value of Q(n,p).
/// \param power Power of Q(n,p).
/// \return The corresponding Q-vector value.
TComplex JAC2hAnalysis::Q(const int harmoN, const int power)
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

/// \brief Decompose the provided integer into his two components.
/// \param inPair Integer representing the input pair of harmonics.
/// \return 2D array for the decomposed pair.
std::array<int, 2> JAC2hAnalysis::GetHarmoPair(int inPair)
{
    std::array<int, 2> pair;
    /* Security check to ensure proper behaviour. */
    if (inPair == 0) {
        pair.at(0) = 0; pair.at(1) = 0;
        return pair;
    }

    /* Extract each digit starting by the unit and save them in reverse. */
    pair.at(1) = inPair % 10;
    pair.at(0) = inPair/10;    // This removes the units as we deal with integers.
    return pair;
}

/// \brief Compute and fill the profiles with all the needed correlation terms.
/// \param myCent Centrality class where the collision belongs to.
/// \param mySample Bootstrap sample associated with the collision.
void JAC2hAnalysis::ComputeAllTerms(const int myCent, const int mySample)
{
    if (mDebugLevel > 1) {printf("myCent: %d mySample: %d\n", myCent, mySample);}

    /* Compute the denominators (= event weights) for each case of interest. */
    std::array<double, 5> eventWeights;  // All needed event weights for this event.
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

    /* Compute the 2-particle correlators for v1 to v6 without eta gap. */
    int list_harmos[7] = {0};   // We work with symmetric correlators of max 14
                                // particles --> max 7 harmonics.
    std::array<double, 6> correl2p;     // Value of the correlation terms to use in profiles.

    for (int iVn = 1; iVn <= 6; iVn++) {  // We go to v6 maximum.
        list_harmos[0] = iVn;   // Only the first element is filled in this case.
        correl2p.at(iVn-1) = ComputeCorrel(2, eventWeights.at(0), list_harmos);
        if (mDebugLevel >= 2) {printf("correl2p: %e.\n", correl2p.at(iVn-1));}
    } // Go to the next harmonic with 2-particles.

    /* All 2-particle terms have been obtained, we pass to the 2-harmonic case. */
    std::vector<std::array<double, 14>> correl2h;   // Values of the correlators.
    std::vector<std::array<double, 14>> eWeight2h;  // Corresponding event weights.
    int hIndex = 0;   // Index for the i-th harmonic pair.
    for (auto hPair : m2hPairs) {
        std::array<double, 14> thisCorrel;  // Values of the correlators.
        std::array<double, 14> thisWeight;  // Corresponding event weights.

        /* Decompose first the current integer into its individual harmonics. */
        std::array<int, 2> harmoPair = GetHarmoPair(hPair);
        mAC2hHistManager.FillPairProf(hIndex, harmoPair);
        if (mDebugLevel >= 1) {printf("Current harmoPair: [%d,%d]\n", harmoPair.at(0), harmoPair.at(1));}

        /* Loop over each power and get the corresponding multiparticle correlator. */
        for (int iPow = 0; iPow < 14; iPow++) { // Dimension should match 'm2hPowers'.
            /* Reset the variables for this new power. */
            int nPartCorrel = 0; int iEW = 0;
            int cHarmo = 0;   // Counter for the list_harmos filling.

            /* Update the power list. */
            // Skip the harmonics that are zero.
            for (int iHarmo = 0; iHarmo < 2; iHarmo++) {
                if (m2hPowers.at(iPow).at(iHarmo) == 0) {continue;}

                /* Write the harmonic as many times as its power is. */
                for (int jP = 0; jP < m2hPowers.at(iPow).at(iHarmo); jP++) {
                    list_harmos[cHarmo] = harmoPair.at(iHarmo);
                    cHarmo++;
                }

                /* Determine the number of particles in the correlator. */
                // 2-harmonic terms in AC have twice as many particles as their cumulant order.
                nPartCorrel += 2*m2hPowers.at(iPow).at(iHarmo);
                if (mDebugLevel >= 2) {
                    printf("m2hPowers.at(iPow).at(iHarmo): %d\n", m2hPowers.at(iPow).at(iHarmo));
                }
            } // Go to the second harmonic of the pair.
            
            /* Get the right event weight using 'iEW'. */
            iEW = (nPartCorrel/2)-1;
            printf("nPartCorrel: %d iEW: %d\n", nPartCorrel, iEW);
            thisWeight.at(iPow) = eventWeights.at(iEW);

            /* Compute the correlator for this harmonic array. */
            thisCorrel.at(iPow) = ComputeCorrel(nPartCorrel, eventWeights.at(iEW),
                list_harmos);
        } // Go to the next power.

        /* Save the correlators and weights obtained for this pair of harmonics. */
        correl2h.push_back(thisCorrel);
        eWeight2h.push_back(thisWeight);
        hIndex++;   // Raise the index for the next harmonic pair.
    }   // Go to the next harmonic pair in the list.

    /* Fill the profile2D with the correlators. */
    // TODO: Find a less ugly way to implement the filling of the 2hProf.
    switch (myCent) {
    case 0:
        mAC2hHistManager.Fill2pProf<0>(correl2p, mySample, eventWeights.at(0));
        mAC2hHistManager.Fill2hProf<0,0>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<0,1>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<0,2>(correl2h, eWeight2h, mySample);
        break;
    case 1:
        mAC2hHistManager.Fill2pProf<1>(correl2p, mySample, eventWeights.at(0));
        mAC2hHistManager.Fill2hProf<1,0>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<1,1>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<1,2>(correl2h, eWeight2h, mySample);
        break;
    case 2:
        mAC2hHistManager.Fill2pProf<2>(correl2p, mySample, eventWeights.at(0));
        mAC2hHistManager.Fill2hProf<2,0>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<2,1>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<2,2>(correl2h, eWeight2h, mySample);
        break;
    case 3:
        mAC2hHistManager.Fill2pProf<3>(correl2p, mySample, eventWeights.at(0));
        mAC2hHistManager.Fill2hProf<3,0>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<3,1>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<3,2>(correl2h, eWeight2h, mySample);
        break;
    case 4:
        mAC2hHistManager.Fill2pProf<4>(correl2p, mySample, eventWeights.at(0));
        mAC2hHistManager.Fill2hProf<4,0>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<4,1>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<4,2>(correl2h, eWeight2h, mySample);
        break;
    case 5:
        mAC2hHistManager.Fill2pProf<5>(correl2p, mySample, eventWeights.at(0));
        mAC2hHistManager.Fill2hProf<5,0>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<5,1>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<5,2>(correl2h, eWeight2h, mySample);
        break;
    case 6:
        mAC2hHistManager.Fill2pProf<6>(correl2p, mySample, eventWeights.at(0));
        mAC2hHistManager.Fill2hProf<6,0>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<6,1>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<6,2>(correl2h, eWeight2h, mySample);
        break;
    case 7:
        mAC2hHistManager.Fill2pProf<7>(correl2p, mySample, eventWeights.at(0));
        mAC2hHistManager.Fill2hProf<7,0>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<7,1>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<7,2>(correl2h, eWeight2h, mySample);
        break;
    case 8:
        mAC2hHistManager.Fill2pProf<8>(correl2p, mySample, eventWeights.at(0));
        mAC2hHistManager.Fill2hProf<8,0>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<8,1>(correl2h, eWeight2h, mySample);
        mAC2hHistManager.Fill2hProf<8,2>(correl2h, eWeight2h, mySample);
        break;
    }

  if (mDebugLevel > 1) {printf("Terms obtained for all pairs of interest.\n");}
  return;
}

/// \brief Compute the multiparticle correlator term.
/// \param nPC: Order of the correlator (i.e number of particles).
/// \return The correlator since the event weight is already calculated.
double JAC2hAnalysis::ComputeCorrel(const int nPC, const double myEventWeight,
    const int (&myHarmos)[7])
{
    if (mDebugLevel > 1) {
        printf("Computing correlation term with %d particles.\n", nPC);
    }

    /* Initialize all possible cases of numerator arrays with the correct signs. */
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

    /* Compute the complex for of the multiparticle correlator term. */
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
    default:
        printf("Error: invalid number of particles. Doing nothing...\n");
        break;
    }

    return cCorrel.Re();  // Real part of the correlator of interest.
}

} // namespace o2::analysis::PWGCF