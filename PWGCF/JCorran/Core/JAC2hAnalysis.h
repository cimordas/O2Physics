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

// Includes.
#include <TComplex.h>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"

#include "PWGCF/JCorran/DataModel/JCatalyst.h"
#include "PWGCF/JCorran/Core/JAC2hHistManager.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis::PWGCF
{
class JAC2hAnalysis
{
public:
  JAC2hAnalysis() = default;  // Constructor.

  // Setters and getters.
  void SetAC2hHistManager(JAC2hHistManager thisManager)
  {
    mAC2hHistManager = thisManager;
  }
  JAC2hHistManager GetAC2hHistManager() const {return mAC2hHistManager;}

  void SetDebugLevel(int debug)
  {
    mDebugLevel = debug;
    LOGF(info, "Debug level set to %d.", mDebugLevel);
  }
  int GetDebugLevel() const {return mDebugLevel;}

  void SetHP(int harmos = 61, int powers = 11)
  {
    mNharmos = harmos; mNpowers = powers;
    LOGF(info, "Number of harmos set to %d.", mNharmos);
    LOGF(info, "Number of powers set to %d.", mNpowers);
  }
  int GetHarmos() const {return mNharmos;}
  int GetPowers() const {return mNpowers;}

  // Class-specific methods.
  void InitArrays();  // Initialize to default values the array data members.

  void CalculateQvectors(int nPart, double myAngles[], double myWeights[]);
  TComplex Q(int harmoN, int power);
  TComplex CalculateRecursion(int n, int *harmonic, int mult=1, int skip=0);

  void ComputeAllTerms(const int myCent, const int mySample);
  double ComputeCorrel(int nPC, const double myEventWeight,
                      const int (&myHarmos)[7]);

private:
  JAC2hHistManager mAC2hHistManager;  ///< Predefined histogram manager.
  int mDebugLevel = 0;    ///< Class verbosity. 0: minimal, 1: end printing,
                          //   2: high (debug).

  int m2hPowers[14][2];   ///< All combinations of powers for the 2h AC.
  int m2hPairs[2];        ///< All pairs of harmonics of interest.
  int mNharmos = 61;      ///< Max number of harmonics: (v6*10part)+1.
  int mNpowers = 11;      ///< Max number of powers: 10part+1.
  TComplex mQvectors[61][11];   ///< All needed combinations of Q-vectors.
                                //   Dimensions must match [mNharmos][mNpowers].


  ClassDefNV(JAC2hAnalysis, 1);
};
} // namespace o2::analysis::PWGCF

#endif // JAC2HANALYSIS_H