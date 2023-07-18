
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

/// \brief  Histogram manager for the AC analysis.
/// \author Cindy Mordasini (cindy.mordasini@cern.ch)
/// \since  July 2023

#ifndef JAC2HHISTMANAGER_H
#define JAC2HHISTMANAGER_H

// Includes.
#include <string>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
//#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "PWGCF/JCorran/DataModel/JCatalyst.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis::PWGCF
{
class JAC2hHistManager
{
public:
  JAC2hHistManager() = default;

  void SetHistManager(HistogramRegistry *myRegistry)
  {
    mHistoRegistry = myRegistry;
    LOGF(info, "Histogram registry has been set.");
  }
  HistogramRegistry *GetHistManager() const {return mHistoRegistry;}

  void SetNsamples(int mySamples)
  {
    mNsamples = mySamples;
    LOGF(info, "Number of samples set to %d.", mNsamples);
  }
  int GetNsamples() const {return mNsamples;}

  void CreateHistosQA();
  void CreateHistosAN();

  // The template functions are defined in the header to prevent compilation issues.
  template <int cBin, typename T>
  void FillEventQA(const T& coll, int multi)
  {
    if (!mHistoRegistry) {
      LOGF(error, "No histogram manager provided. Quit.");
      return;
    }

    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histCent"), coll.cent());
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histZvtx"), coll.posZ());
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histMulti"), multi);
  }

  template <int cBin, typename T>
  void FillTrackQA(const T& track)
  {
    if (!mHistoRegistry) {
      LOGF(error, "No histogram manager provided. Quit.");
      return;
    }

    // NOTE: Crosscheck again that the corrected histograms are filled correctly!!!
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histPtUncorrected"),
                        track.pt());
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histPtCorrected"),
                        track.pt(), 1./(track.weightEff()));
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histEta"),
                        track.eta());
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histPhiUncorrected"),
                        track.phi());
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histPhiCorrected"),
                        track.phi(), 1./(track.weightNUA()));
  }

private:
  HistogramRegistry *mHistoRegistry = nullptr;  ///< QA + analysis output.
  int mNsamples = 10; ///< Number of samples for the bootstrap.

  static constexpr std::string_view mCentClasses[] = {  // Classes from JCatalyst.
    "Centrality_0-1/", "Centrality_1-2/", "Centrality_2-5/", "Centrality_5-10/",
    "Centrality_10-20/", "Centrality_20-30/", "Centrality_30-40/", "Centrality_40-50/",
    "Centrality_50-60/"
  };

  ClassDefNV(JAC2hHistManager, 1);  
};
} // namespace o2::analysis::PWGCF

#endif // JAC2HHISTMANAGER_H