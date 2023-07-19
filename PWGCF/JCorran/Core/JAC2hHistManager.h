
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
  void SetEtaGap(int myGap)
  {
    mEtaGap = myGap;
    LOGF(info, "Eta gap set to %.1f.", mEtaGap);
  }
  int GetEtaGap() const {return mEtaGap;}

  void CreateHistos();

  // The template functions are defined in the header to prevent compilation issues.
  template <int cBin, typename T>
  void FillEventQA(const T& coll, int multi, int sample)
  {
    if (!mHistoRegistry) {
      LOGF(error, "No histogram manager provided. Quit.");
      return;
    }

    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("QA/histCent"), coll.cent());
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("QA/histZvtx"), coll.posZ());
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("QA/histMulti"), multi);
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("QA/histSamples"), sample);
  }

  template <int cBin, typename T>
  void FillTrackQA(const T& track)
  {
    if (!mHistoRegistry) {
      LOGF(error, "No histogram manager provided. Quit.");
      return;
    }

    // NOTE: Crosscheck again that the corrected histograms are filled correctly!!!
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("QA/histPtUncorrected"),
                        track.pt());
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("QA/histPtCorrected"),
                        track.pt(), 1./(track.weightEff()));
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("QA/histEta"),
                        track.eta());
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("QA/histPhiUncorrected"),
                        track.phi());
    mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("QA/histPhiCorrected"),
                        track.phi(), 1./(track.weightNUA()));
  }

private:
  HistogramRegistry *mHistoRegistry = nullptr;  ///< QA + analysis output.
  int mNsamples = 10;   ///< Number of samples for the bootstrap.
  float mEtaGap = 1.0;  ///< Value of the pseudorapidity gap.

  static constexpr std::string_view mCentClasses[] = {  // Classes from JCatalyst.
    "Centrality_00-01/", "Centrality_01-02/", "Centrality_02-05/", "Centrality_05-10/",
    "Centrality_10-20/", "Centrality_20-30/", "Centrality_30-40/", "Centrality_40-50/",
    "Centrality_50-60/"
  };

  ClassDefNV(JAC2hHistManager, 1);  
};
} // namespace o2::analysis::PWGCF

#endif // JAC2HHISTMANAGER_H