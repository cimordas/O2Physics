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

#ifndef JAC2HANALYSIS_H
#define JAC2HANALYSIS_H

/* Header files. */
#include <array>
#include <vector>
#include <TComplex.h>

// O2 headers.
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"

// O2Physics headers.

// JCorran headers.
#include "PWGCF/JCorran/DataModel/JCatalyst.h"
#include "PWGCF/JCorran/Core/JAC2hHistManager.h"

/* Namespaces. */
using namespace o2;
using namespace o2::framework;

// ----------------------------------------------------------------------------
// Analysis class with the methods to compute the correlator terms needed in
// AC studies with experimental data.
// ----------------------------------------------------------------------------
namespace o2::analysis::PWGCF
{
class JAC2hAnalysis
{
public:
    JAC2hAnalysis() = default;  // Constructor.

    /* Setters and getters. */
    void SetAC2hHistManager(JAC2hHistManager thisManager)
    {
        mAC2hHistManager = thisManager;
        LOGF(info, "Histogram registry successfully set.");
    }
    JAC2hHistManager GetAC2hHistManager() const {return mAC2hHistManager;}

    void SetDebugLevel(int debug)
    {
        mDebugLevel = debug;
        LOGF(info, "Debug level: %d.", mDebugLevel);
    }
    int GetDebugLevel() const {return mDebugLevel;}

    // No setter for the dimensions as many things rely on their fixed values.
    int GetHarmos() const {return mNqHarmos;}
    int GetPowers() const {return mNqPowers;}

    void Set2hPairs(const std::vector<int>& harmos)
    {
        m2hPairs = harmos;
        LOGF(info, "Pairs of 2h provided.");
    }

    /* Methods related to this class. */
    void CalculateQvectors(const std::vector<double>& myAngles,
        const std::vector<double>& myWeights);
    TComplex Q(const int harmoN, const int power);
    TComplex CalculateRecursion(int n, int *harmonic, int mult=1, int skip=0);
    std::array<int, 2> GetHarmoPair(int inPair);
    void ComputeAllTerms(const int myCent, const int mySample);
    double ComputeCorrel(const int nPC, const double myEventWeight,
        const int (&myHarmos)[7]);

private:
    JAC2hHistManager mAC2hHistManager;  ///< Predefined histogram manager.

    int mDebugLevel = 0;    ///< Class verbosity. 0: none, 1: low, 2: high/debug.
    const int mNqHarmos = 61;   ///< Max number of harmonics for Qn: (v6*10part)+1.
    const int mNqPowers = 11;   ///< Max number of powers for Qn: 10part+1.
    std::vector<int> m2hPairs;  ///< All pairs of harmonics of interest.
                                //   Dimensions must be at least 'cfgNcombis2h'.
    std::array<std::array<int, 2>, 14> m2hPowers = {{
        {1,1}, {2,2}, {2,0}, {2,1}, {3,0}, {3,1}, {4,0}, {4,1},
        {0,2}, {1,2}, {0,3}, {1,3}, {0,4}, {1,4}
        }};
        ///< All combinations of powers for the 2h AC.
    std::array<std::array<TComplex, 11>, 61> mQvectors;
        ///< All needed combinations of Q-vectors.
        //   Dimensions must match [mNqHarmos][mNqPowers].

    ClassDefNV(JAC2hAnalysis, 1);
};
} // namespace o2::analysis::PWGCF

#endif // JAC2HANALYSIS_H