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

// \brief   Calculation class for the AC-related analyses.
// \author  Cindy Mordasini (cindy.mordasini@cern.ch)

#ifndef PWGCF_FLOW_CORE_FLOWAC2HANALYSIS_H
#define PWGCF_FLOW_CORE_FLOWAC2HANALYSIS_H

/* Header files. */
#include <array>
#include <vector>
#include <TComplex.h>

// O2 headers. //
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"

// O2 Physics headers.
#include "PWGCF/Flow/Core/FlowAC2hHistManager.h"

/* Namespaces. */
using namespace o2;
using namespace o2::framework;

/* Class for the calculations of the correlator terms needed in AC studies. */
namespace o2::analysis::PWGCF
{
class FlowAC2hAnalysis
{
public:
  FlowAC2hAnalysis() = default;   // Constructor.

  /* Setters and getters. */
  void SetAC2hHistManager(FlowAC2hHistManager thisManager)
  {
    mAC2hHistManager = thisManager;
    LOGF(info, "Histogram manager successfully set.");
  }
  FlowAC2hHistManager GetAC2hHistManager() const {return mAC2hHistManager;}

  void SetDebugPrint(bool debug)
  {
    mDebugPrint = debug;
    LOGF(info, "Debug level: %d", mDebugPrint);
  }
  bool GetDebugPrint() const {return mDebugPrint;}

  void Set2hPairs(const std::vector<int>& harmos)
  {
    m2hPairs = harmos;
    LOGF(info, "Pairs of 2h successfully provided.");
  }

  // No setter for the dimensions as many things rely on their fixed values.
  int GetHarmos() const {return mNqHarmos;}
  int GetPowers() const {return mNqPowers;}

  /* Methods specific to this class. */
  std::array<int, 2> GetHarmoPair(int inPair);

  void CalculateQvectors(const std::vector<float>& myAngles,
                         const std::vector<float>& myWeights);
  TComplex Q(const int harmoN, const int power);
  TComplex CalculateRecursion(int n, int *harmonic, int mult=1, int skip=0);
  void ComputeAllCorrel(const int myCentBin, const int mySample);
  double ComputeThisCorrel(const int nPC, const double myEventWeight,
                           const int (&myHarmos)[7]);

private:
  FlowAC2hHistManager mAC2hHistManager;  ///< Link to the histogram manager.

  bool mDebugPrint = true;    ///< Class verbosity. 0: none, 1: debug.
  const int mNqHarmos = 61;   ///< Highest harmo for Q(n,p): (v6*10part)+1.
  const int mNqPowers = 11;   ///< Max power for Q(n,p): 10part+1.
  std::vector<int> m2hPairs;  ///< All pairs of harmonics of interest.
                                // Dimensions must be at least 'cfgNcombis2h'.

  std::array<std::array<int, 2>, 14> m2hPowers = {{
    {1,1}, {2,2}, {2,0}, {2,1}, {3,0}, {3,1}, {4,0}, {4,1},
                  {0,2}, {1,2}, {0,3}, {1,3}, {0,4}, {1,4}
  }};   ///< All combinations of powers for the 2h AC.
          // Must correspond to the one in FlowAC2hHistManager.h.

  std::array<std::array<TComplex, 11>, 61> mQvectors;   ///< Q(n,p) for all (n,p).
                                // Dimensions must match [mNqHarmos][mNqPowers].

  ClassDefNV(FlowAC2hAnalysis, 1);  
};
} // namespace o2::analysis::PWGCF

#endif  // PWGCF_FLOW_CORE_FLOWAC2HANALYSIS_H
