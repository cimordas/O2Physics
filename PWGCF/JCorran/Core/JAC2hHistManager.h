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
#include "Framework/RunningWorkflowInfo.h"

#include "PWGCF/JCorran/DataModel/JCatalyst.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis::PWGCF
{
class JAC2hHistManager
{
public:
  JAC2hHistManager() = default;   // Constructor.

  // Setters and getters.
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
  void SetNcombis2h(int myCombis)
  {
    mNcombis2h = myCombis;
    LOGF(info, "Number of pairs of harmos set to %d.", mNcombis2h);
  }
  int GetNcombis2h() const {return mNcombis2h;}

  // Class-specific methods.
  void CreateHistos();

  /// The template functions are defined in the header to prevent compilation issues.
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

    // LOKI: added just to make the compilator happy...
    const int nCentBins = sizeof(jflucCentBins)/sizeof(jflucCentBins[0]);
    printf("Number of centrality classes in JCatalyst: %d\n", nCentBins);
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

  template <int cBin, int sBin, typename T>
  void Fill2pProf(const T& correl2p, double weight2p)
  {
    if (!mHistoRegistry) {
      LOGF(error, "No histogram manager provided. Quit.");
      return;
    }

    for (int iB = 0; iB < 6; iB++) {
      mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("Full/prof2pCorrel"),
                          iB+0.5, correl2p[iB], weight2p);
      mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST(mSamples[sBin])+HIST("prof2pCorrel"),
                          iB+0.5, correl2p[iB], weight2p);
    }
  }

  template <int cBin, int sBin, int pBin, typename T>
  void Fill2hProf(const T& correl2h, const T& weight2h)
  {
    if (!mHistoRegistry) {
      LOGF(error, "No histogram manager provided. Quit.");
      return;
    }

    static constexpr std::string_view combis[] = {"Combi1", "Combi2", "Combi3"};

    for (int iB = 0; iB < 14; iB++) {
      mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("Full/prof2hCorrel")+HIST(combis[pBin]),
                          iB+0.5, correl2h[iB], weight2h[iB]);
      mHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST(mSamples[sBin])+HIST("prof2hCorrel")+HIST(combis[pBin]),
                          iB+0.5, correl2h[iB], weight2h[iB]);
    }
  }

private:
  HistogramRegistry *mHistoRegistry = nullptr;  ///< QA + analysis output.
  int mNsamples = 20;   ///< Number of samples for the bootstrap.
  float mEtaGap = 1.0;  ///< Value of the pseudorapidity gap.
  int mNcombis2h = 3;   ///< Number of pairs of harmonics.

  static constexpr std::string_view mCentClasses[] = {  // Classes from JCatalyst.
    "Centrality_00-01/", "Centrality_01-02/", "Centrality_02-05/", "Centrality_05-10/",
    "Centrality_10-20/", "Centrality_20-30/", "Centrality_30-40/", "Centrality_40-50/",
    "Centrality_50-60/"
  };
  static constexpr std::string_view mSamples[] = {
    "Sample00/", "Sample01/", "Sample02/", "Sample03/", "Sample04/", "Sample05/",
    "Sample06/", "Sample07/", "Sample08/", "Sample09/", "Sample10/", "Sample11/",
    "Sample12/", "Sample13/", "Sample14/", "Sample15/", "Sample16/", "Sample17/",
    "Sample18/", "Sample19/"
  };
  static constexpr std::string_view mPowers[] = {
    "{1,1}", "{2,2}", "{2,0}", "{2,1}", "{3,0}", "{3,1}", "{4,0}", "{4,1}",
                      "{0,2}", "{1,2}", "{0,3}", "{1,3}", "{0,4}", "{1,4}"
  };

  ClassDefNV(JAC2hHistManager, 1);  
};
} // namespace o2::analysis::PWGCF

#endif // JAC2HHISTMANAGER_H