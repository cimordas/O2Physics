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
    mHistRegistryQA = myRegistry;
    LOGF(info, "QA histogram registry successfully set.");
  }
  HistogramRegistry *GetHistRegistryQA() const {return mHistRegistryQA;}

  void SetHistRegistryAN(HistogramRegistry *myRegistry)
  {
    mHistRegistryAN = myRegistry;
    LOGF(info, "AN histogram registry successfully set.");
  }
  HistogramRegistry *GetHistRegistryAN() const {return mHistRegistryAN;}

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

  void SetEtaGap(float myGap)
  {
    mEtaGap = myGap;
    LOGF(info, "Eta gap value: %.2f.", mEtaGap);
  }
  int GetEtaGap() const {return mEtaGap;}

  /* Methods specific to this class. */
  // The template functions are defined here to prevent compilation errors.
  void CreateHistQA();
  void CreateHistAN();
  int GetCentBin(float cValue);
  void FillPairProf(const int hBin, const std::array<int, 2>& pairHarmo);

  /// \brief Fill the event QA histograms in the centrality class.
  /// \param cBin Centrality bin of the collision.
  /// \param coll Collision entry of the table.
  /// \param multi Collision multiplicity at this step.
  template <int cBin, typename T>
  void FillEventQA(const T& coll, float cent, int multi)
  {
    if (!mHistRegistryQA) {
      LOGF(fatal, "QA histogram registry missing. Quitting...");
      return;
    }
    mHistRegistryQA->fill(HIST(mCentClasses[cBin])+HIST("histCent"), cent);
    mHistRegistryQA->fill(HIST(mCentClasses[cBin])+HIST("histMulti"), multi);
    mHistRegistryQA->fill(HIST(mCentClasses[cBin])+HIST("histZvtx"), coll.posZ());
  }

  /// \brief Fill the track QA histograms in the centrality class.
  /// \param cBin Centrality bin of the collision.
  /// \param track Track entry of the table.
  // TODO: Add filling of the weight histograms.
  template <int cBin, typename T>
  void FillTrackQA(const T& track)
  {
    if (!mHistRegistryQA) {
      LOGF(fatal, "QA histogram registry missing. Quitting...");
      return;
    }

    mHistRegistryQA->fill(HIST(mCentClasses[cBin])+HIST("histPt"), track.pt());
    mHistRegistryQA->fill(HIST(mCentClasses[cBin])+HIST("histEta"), track.eta());
    mHistRegistryQA->fill(HIST(mCentClasses[cBin])+HIST("histPhi"), track.phi());
    mHistRegistryQA->fill(HIST(mCentClasses[cBin])+HIST("histCharge"), track.sign());
  }

private:
  HistogramRegistry *mHistRegistryQA = nullptr;   ///< For the QA output.
  HistogramRegistry *mHistRegistryAN = nullptr;   ///< For the analysis output.

  bool mSaveQABefore = true;  ///< Save the QA output before the selection cuts?
  int mNcombis2h = 3;   ///< Number of pairs of harmonics.
  int mNsamples = 20;   ///< Number of samples for the bootstrap.
  float mEtaGap = 1.0;  ///< Value of the pseudorapidity gap.

  static const int mNcentBins = 10;   ///< Number of centrality classes.
  static constexpr std::string_view mCentClasses[] = {  ///< Centrality classes.
    "Centrality_00-01/", "Centrality_01-02/", "Centrality_02-05/",
    "Centrality_05-10/", "Centrality_10-20/", "Centrality_20-30/",
    "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/",
    "Centrality_60-70/"
  };

  static constexpr std::string_view mPowers[] = {   ///< Labels for 2h profile.
    "{1,1}", "{2,2}", "{2,0}", "{2,1}", "{3,0}", "{3,1}", "{4,0}", "{4,1}",
                      "{0,2}", "{1,2}", "{0,3}", "{1,3}", "{0,4}", "{1,4}"
  };

  static constexpr std::string_view mCombi[] = {     ///< 
      "Combi0", "Combi1", "Combi2"
  };

  ClassDefNV(FlowAC2hHistManager, 1);
};
} // namespace o2::analysis::PWGCF

#endif  // PWGCF_FLOW_CORE_FLOWAC2HHISTMANAGER_H
