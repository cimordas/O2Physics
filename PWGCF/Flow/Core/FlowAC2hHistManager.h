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

// \brief   Histogram manager for the AC-related analyses.
// \author  Cindy Mordasini (cindy.mordasini@cern.ch)

#ifndef PWGCF_FLOW_CORE_FLOWAC2HHISTMANAGER_H
#define PWGCF_FLOW_CORE_FLOWAC2HHISTMANAGER_H

/* Header files. */
#include <array>
#include <string>
#include <string_view>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>

// O2 headers. //
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"

// O2 Physics headers.

/* Namespaces. */
using namespace o2;
using namespace o2::framework;

// ----------------------------------------------------------------------------
// Histogram manager to fill both QA and AN registries at the same time.
// ----------------------------------------------------------------------------
namespace o2::analysis::PWGCF
{
class FlowAC2hHistManager
{
public:
  FlowAC2hHistManager() = default;

  /* Setters and getters. */
  void SetHistRegistryQA(HistogramRegistry *myRegistry)
  {
    mQAHistRegistry = myRegistry;
    LOGF(info, "QA histogram registry successfully set.");
  }
  HistogramRegistry *GetQAHistRegistry() const {return mQAHistRegistry;}

  void SetNcombis2h(int myCombis)
  {
    mNcombis2h = myCombis;
    LOGF(info, "Number of pairs of harmonics: %d.", mNcombis2h);
  }
  int GetNcombis2h() const {return mNcombis2h;}

  void SetNsamples(int mySamples)
  {
      mNsamples = mySamples;
      LOGF(info, "Number of samples: %d.", mNsamples);
  }
  int GetNsamples() const {return mNsamples;}

  /* Methods specific to this class. */
  // The template functions are defined here to prevent compilation errors.
  void CreateHistsQA();

  int GetCentBin(float cValue);

  /// \brief Fill the event QA histograms in the centrality class.
  /// \param cBin Centrality bin of the collision.
  /// \param mode 0: Before cuts, 1: After cuts.
  /// \param coll Collision entry of the table.
  /// \param multi Collision multiplicity at this step.
  template <int cBin, int mode, typename T>
  void FillEventQA(const T& coll, float cent, int multi)
  {
    if (!mQAHistRegistry) {
      LOGF(fatal, "QA histogram registry missing. Quitting...");
      return;
    }
    static constexpr std::string_view subDir[] = {"Before/", "After/"};

    mQAHistRegistry->fill(HIST(mCentClasses[cBin])+HIST(subDir[mode])+HIST("histCent"), cent);
    mQAHistRegistry->fill(HIST(mCentClasses[cBin])+HIST(subDir[mode])+HIST("histMulti"), multi);
    mQAHistRegistry->fill(HIST(mCentClasses[cBin])+HIST(subDir[mode])+HIST("histZvtx"), coll.posZ());
  }

private:
  HistogramRegistry *mQAHistRegistry = nullptr;   ///< For the QA output.

  int mNcombis2h = 3;   ///< Number of pairs of harmonics.
  int mNsamples = 20;   ///< Number of samples for the bootstrap.

  static const int mNcentBins = 10;   ///< Number of centrality classes.
  static constexpr std::string_view mCentClasses[] = {  ///< Centrality classes.
    "Centrality_00-01/", "Centrality_01-02/", "Centrality_02-05/",
    "Centrality_05-10/", "Centrality_10-20/", "Centrality_20-30/",
    "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/",
    "Centrality_60-70/"
  };

  ClassDefNV(FlowAC2hHistManager, 1);
};
} // namespace o2::analysis::PWGCF

#endif  // PWGCF_FLOW_CORE_FLOWAC2HHISTMANAGER_H
