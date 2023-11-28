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

// \brief   Task for the calculation of the 2-harmonic AC with filtered data.
// \author  Cindy Mordasini (cindy.mordasini@cern.ch)

/* Header files. */
//#include <vector>
#include <TRandom3.h>

// O2 headers. //
// The first two are mandatory.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

// O2 Physics headers. //
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGCF/Flow/Core/FlowAC2hHistManager.h"
#include "PWGCF/Flow/Core/FlowAC2hAnalysis.h"

/* Namespaces and definitions. */
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::PWGCF;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                               aod::FT0sCorrected, aod::CentFT0Ms,
                               aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As,
                               aod::CentFDDMs, aod::CentNTPVs>;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;

// ----------------------------------------------------------------------------
// Analysis task starts here.
// ----------------------------------------------------------------------------
struct flowAC2hAnalysisTask
{
  /* Declare the instances of other classes, and histogram registries. */
  // Objects in the registries are saved in alphabetical order (true) and in
  // TDirectoryFile (true). Multiple registries are needed to avoid going over
  // the limit of objects per registry (max: 512).
  FlowAC2hAnalysis acAnalysis;
  FlowAC2hHistManager histManager;
  HistogramRegistry qaHistRegistry{"qaHistRegistry", {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true, true};
  HistogramRegistry anHistRegistry{"anHistRegistry", {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true, true};

  /* Enable debug and additional saving options. */
  Configurable<bool> cfgPrintDebug{"cfgPrintDebug", true, "Print additional log for debug."};
  Configurable<bool> cfgSaveQA{"cfgSaveQA", true, "Save all QA objects for debug."};
  
  /* Set the event quality cuts. */
  Configurable<int> cfgCentEsti{"cfgCentEsti", 2,
                                "Centrality estimator (Run3), 0: FT0M, 1: FT0A, 2: FT0C, 3: FV0A, 4: FDDM, 5: NTPV."};  
  Configurable<int> cfgMultMin{"cfgMultMin", 10, "Minimum multiplicity for a collision to be analyzed."};
  Configurable<float> cfgZvtxMax{"cfgCutZvtx", 10.0f, "Maximum value for Z-vtx."};

  /* Set the track quality cuts. */
  Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum pT for tracks"};
  Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "Maximum pT for tracks"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8f, "Maximum eta range for tracks"};

  /* Set the analysis details. */
  Configurable<int> cfgNcombis2h{"cfgNcombis2h", 3, "Number of harmonics pairs."};
  Configurable<int> cfgNsamples{"cfgNsamples", 20, "Number of bootstrap samples."};
  Configurable<float> cfgEtaGap{"cfgEtaGap", 1.0, "Minimum eta gap value."};
  Configurable<std::vector<int>> cfg2hHarmos{"cfg2hHarmos", {23, 24, 34}, "List of 2h pairs"};

  /* Define the filters to apply to the data. */
  // The analysis assumes the data has been subjected to a QA of its selection,
  // and thus only the final distributions of the data for analysis are saved.
  Filter collFilter = (nabs(aod::collision::posZ) < cfgZvtxMax);
  Filter trackFilter =  (aod::track::pt > cfgPtMin) && (aod::track::pt < cfgPtMax)
                     && (nabs(aod::track::eta) < cfgEtaMax);
                     //&& ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  /// \brief O2 function executed at the beginning of the task workflow.
  void init(InitContext const&)
  {
    /* Reset gRandom to ensure a new sequence is used each time. */
    delete gRandom;
    gRandom = new TRandom3(0);

    /* Initialize the histogram manager with the values it requires to work. */
    histManager.SetHistRegistryQA(&qaHistRegistry);
    histManager.SetHistRegistryAN(&anHistRegistry);
    histManager.SetNcombis2h(cfgNcombis2h);
    histManager.SetNsamples(cfgNsamples);
    histManager.SetEtaGap(cfgEtaGap);
    if (cfgSaveQA) {histManager.CreateHistQA();}
    histManager.CreateHistAN();

    /* Initialize the analysis class instance with the values it requires. */
    acAnalysis.SetAC2hHistManager(histManager);
    acAnalysis.SetDebugPrint(cfgPrintDebug);
    acAnalysis.Set2hPairs(cfg2hHarmos);
  }

  /// \brief O2 function executed for each collision.
  /// \note The collisions and tracks are filtered here by default.
  //        Therefore, no "Before" QA can be done.
  void process(soa::Filtered<MyCollisions>::iterator const& coll, soa::Filtered<MyTracks> const& tracks)
  {
    /* Reject invalid events: no tracks or invalid centrality. */
    int nTracks = tracks.size();
    if (nTracks == 0) {return;}

    float centAllEstim[6] = {
      coll.centFT0M(), coll.centFT0A(), coll.centFT0C(),
      coll.centFV0A(), coll.centFDDM(), coll.centNTPV()
    };
    float cent = centAllEstim[cfgCentEsti];
    if (cent < 0. || cent > 70.) {return;}  // We don't go higher than 70% centrality.
    int cBin = histManager.GetCentBin(cent);

    if (cfgPrintDebug) {
      LOGF(info, "Number of tracks in this collision: %d", nTracks);
      LOGF(info, "Centrality value: %.2f Centrality bin: %d", cent, cBin);
    }

    /* Apply additional event and track selection not applied in the filter. */
    // The minimum cut on the multiplicity ensures that the correlators are
    // properly defined (must be the last cut applied).
    if (!coll.sel8()) {return;}
    if (nTracks < cfgMultMin) {return;}

    /* Give a sample ID to this collision based on a uniform distribution. */
    int sampleID = int(gRandom->Uniform(0, cfgNsamples));
    if (cfgPrintDebug) {
      LOGF(info, "Sample ID of this collision: %d", sampleID);
    }

    /* Save the azimuthal angles and weights for the correlators. */
    // The QA histograms are filled if enabled, the EventQA is filled only for
    // the first track of the collision.
    bool isFirstTrack = true;
    std::vector<float> trackPhi;
    std::vector<float> trackWeight;

    for (auto& track : tracks) {
      trackPhi.push_back(track.phi());
      trackWeight.push_back(1.);  // TODO: Add gestions of non-unit NUE and NUA weights.

      if (cfgSaveQA) {
        switch (cBin) {
        case 0:
          if (isFirstTrack) {histManager.FillEventQA<0>(coll, cent, nTracks);}
          histManager.FillTrackQA<0>(track);
          break;
        case 1:
          if (isFirstTrack) {histManager.FillEventQA<1>(coll, cent, nTracks);}
          histManager.FillTrackQA<1>(track);
          break;
        case 2:
          if (isFirstTrack) {histManager.FillEventQA<2>(coll, cent, nTracks);}
          histManager.FillTrackQA<2>(track);
          break;
        case 3:
          if (isFirstTrack) {histManager.FillEventQA<3>(coll, cent, nTracks);}
          histManager.FillTrackQA<3>(track);
          break;
        case 4:
          if (isFirstTrack) {histManager.FillEventQA<4>(coll, cent, nTracks);}
          histManager.FillTrackQA<4>(track);
          break;
        case 5:
          if (isFirstTrack) {histManager.FillEventQA<5>(coll, cent, nTracks);}
          histManager.FillTrackQA<5>(track);
          break;
        case 6:
          if (isFirstTrack) {histManager.FillEventQA<6>(coll, cent, nTracks);}
          histManager.FillTrackQA<6>(track);
          break;
        case 7:
          if (isFirstTrack) {histManager.FillEventQA<7>(coll, cent, nTracks);}
          histManager.FillTrackQA<7>(track);
          break;
        case 8:
          if (isFirstTrack) {histManager.FillEventQA<8>(coll, cent, nTracks);}
          histManager.FillTrackQA<8>(track);
          break;
        case 9:
          if (isFirstTrack) {histManager.FillEventQA<9>(coll, cent, nTracks);}
          histManager.FillTrackQA<9>(track);
          break;
        default:
          if (cfgPrintDebug) {
            LOGF(info, "Centrality percentile not included in analysis. Next...");
          }
          break;
        }
      }

      // Indicate when the first track of the collision has been analysed.
      if (isFirstTrack) {isFirstTrack = false;}
    }

    /* Compute the Q-vectors, then the multiparticle correlators. */
    acAnalysis.CalculateQvectors(trackPhi, trackWeight);
    acAnalysis.ComputeAllCorrel(cBin, sampleID);

    /* Reset the variables for the next collision. */
    // This ensures no mixing between collision can happen accidentally.
    trackPhi.clear();
    trackWeight.clear();

    LOGF(info, "Collision analysed. Next...");
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowAC2hAnalysisTask>(cfgc)
  };
}
