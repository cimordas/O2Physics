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

// Standard headers.
#include <chrono>
#include <string>
#include <vector>
#include <TRandom3.h>

// O2 headers. //
// The first two are mandatory.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "CCDB/BasicCCDBManager.h"
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

// Namespaces and definitions.
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
  // Declare the histogram registries and instance of their manager.
  // Objects in the registries are saved in alphabetical order (true) and in
  // TDirectoryFile (true). Multiple registries are needed to avoid going over
  // the limit of objects per registry (max: 512).
  FlowAC2hHistManager histManager;
  HistogramRegistry qaHistRegistry{"qaHistRegistry", {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true, true};
  HistogramRegistry anHistRegistry{"anHistRegistry", {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true, true};

  // Enable the general features of the analysis task.
  FlowAC2hAnalysis acAnalysis;
  Configurable<bool> cfgDebugLog{"cfgDebugLog", true, "Enable log for debug."};
  Configurable<bool> cfgSaveQA{"cfgSaveQA", true, "Enable the QA histograms."};
  Configurable<bool> cfgUseEtaGap{"cfgUseEtaGap", true,
                                  "Calculate the 2-particle terms with eta gap."};
  Configurable<int> cfgNcombis2h{"cfgNcombis2h", 3, "Number of harmonics pairs."};
  Configurable<int> cfgNsamples{"cfgNsamples", 20, "Number of bootstrap samples."};
  Configurable<float> cfgEtaGap{"cfgEtaGap", 1.0, "Minimum eta gap value."};
  Configurable<std::vector<int>> cfg2hHarmos{"cfg2hHarmos", {23, 24, 34}, "List of 2h pairs"}; 

  // Set the event quality cuts.
  // The centrality estimators are the ones available for Run 3.
  enum centEstimators {FT0M, FT0A, FT0C, FDDM, NTPV};
  Configurable<int> cfgCentEsti{"cfgCentEsti", 2, "Centrality estimator."};
  Configurable<int> cfgMultMin{"cfgMultMin", 10,
                               "Strict minimum number of tracks per collision."};
  Configurable<float> cfgZvtxMax{"cfgCutZvtx", 10.0f, "Maximum range for Zvtx."};

  // Set the track quality cuts. 
  //Configurable<bool> cfgUseNUA{"cfgUseNUA", true, "Enable the use of NUA weights."};
  //Configurable<bool> cfgUseNUE{"cfgUseNUE", true, "Enable the use of NUE weights."};
  Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum pT for tracks"};
  Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "Maximum pT for tracks"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8f, "Maximum eta range."};

  // Set the access to the CCDB for the NUA/NUE weights.
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch",
                                     "Address of the CCDB to get the NUA/NUE."};
    Configurable<long> cfgTime{"ccdb-no-later-than",
                              std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),
                              "Latest acceptable timestamp of creation for the object."};
  } cfgCCDB;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // NOTE: Add here variables related to the ccdb.

  /* Set the analysis details. */

  //Configurable<float> cfgEtaGap{"cfgEtaGap", 1.0, "Minimum eta gap value."};
  //Configurable<std::vector<int>> cfg2hHarmos{"cfg2hHarmos", {23, 24, 34}, "List of 2h pairs"};

  // Define here the filters to apply to the received data.
  // The analysis assumes the data has been subjected to a QA of its selection,
  // and thus only the final distributions of the data for analysis are saved.
  Filter collFilter = (nabs(aod::collision::posZ) < cfgZvtxMax)
                   && (aod::evsel::sel8);
  Filter trackFilter =  (aod::track::pt > cfgPtMin) && (aod::track::pt < cfgPtMax)
                     && (nabs(aod::track::eta) < cfgEtaMax);

  // Define the partitions used to distinguish the 'good' collisions and tracks
  // from the 'bad' ones.
  //SliceCache cache;

  /// \brief O2 function executed at the beginning of the task workflow.
  void init(InitContext const&)
  {
    // Reset gRandom to ensure a new sequence is used for each new execution.
    delete gRandom;
    gRandom = new TRandom3(0);

    // Initialise the histogram manager for the QA and AN objects.
    histManager.SetHistRegistryQA(&qaHistRegistry);
    histManager.SetDebugLog(cfgDebugLog);
    histManager.SetSaveAllQA(false);  // Full QA is always disabled by default.
    histManager.SetSaveQABefore(false);   // No QA before selection saved.
    if (cfgSaveQA) {histManager.CreateHistQA();}

    histManager.SetHistRegistryAN(&anHistRegistry);
    histManager.SetNcombis2h(cfgNcombis2h);
    histManager.SetNsamples(cfgNsamples);
    histManager.SetEtaGap(cfgUseEtaGap);
    histManager.CreateHistAN();

    // Initialize the analysis class instance.
    acAnalysis.SetAC2hHistManager(histManager);
    acAnalysis.SetDebugPrint(cfgDebugLog);
    acAnalysis.Set2hPairs(cfg2hHarmos);
    acAnalysis.SetEtaGap(cfgUseEtaGap, cfgEtaGap);
  
    // Setup the access to the CCDB objects.
    ccdb->setURL(cfgCCDB.cfgURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(cfgCCDB.cfgTime.value);
    // NOTE: Add here the access to specific objects.
  }

  /// \brief O2 function executed for each collision.
  /// \note The collisions and tracks are filtered here by default.
  //        Therefore, no "Before" QA can be done in this analysis task.
  void process(soa::Filtered<MyCollisions>::iterator const& coll, soa::Filtered<MyTracks> const& tracks)
  {
    // We don't process further the collisions without enough tracks or with
    // a centrality percentile outside of the range 0-70%.
    if (tracks.size() <= cfgMultMin) {return;}
    int nTracks = tracks.size();

    float cent = -1.;
    switch (cfgCentEsti) {
      case FT0M : cent = coll.centFT0M(); break;
      case FT0A : cent = coll.centFT0A(); break;
      case FT0C : cent = coll.centFT0C(); break;
      case FDDM : cent = coll.centFDDM(); break;
      case NTPV : cent = coll.centNTPV(); break;
    }
    if (cent < 0. || cent > 70.) {return;}
    int cBin = histManager.GetCentBin(cent);

    if (cfgDebugLog) {
      LOGF(info, "Number of tracks in this collision: %d", nTracks);
      LOGF(info, "Centrality value: %.2f Centrality bin: %d", cent, cBin);
    }

    // Give a sample ID to this collision based on a uniform distribution.
    int sampleID = int(gRandom->Uniform(0, cfgNsamples));
    if (cfgDebugLog) {
      LOGF(info, "Sample ID of this collision: %d", sampleID);
    }

    // Save the azimuthal angles and weights for all correlators, and the
    // pseudorapidities for the 2-particle terms with eta gap. The QA histograms
    // are filled if enabled, the EventQA is filled only for the first track
    // of the collision.
    bool isFirstTrack = true;
    std::vector<float> trackPhi;
    std::vector<float> trackWeight;
    std::vector<float> trackEta;

    for (auto& track : tracks) {
      trackPhi.push_back(track.phi());
      trackWeight.push_back(1.);  // TODO: Add non-unit NUE and NUA weights.
      // Make it so that the vector of weights contain NUE*NUA to avoid having two vectors.
      if (cfgUseEtaGap) {trackEta.push_back(track.eta());}

      if (cfgSaveQA) {
        switch (cBin) {
        case 0:
          if (isFirstTrack) {histManager.FillEventQA<0>(coll, cent, nTracks);}
          histManager.FillTrackQA<0,1>(track);
          break;
        case 1:
          if (isFirstTrack) {histManager.FillEventQA<1>(coll, cent, nTracks);}
          histManager.FillTrackQA<1,1>(track);
          break;
        case 2:
          if (isFirstTrack) {histManager.FillEventQA<2>(coll, cent, nTracks);}
          histManager.FillTrackQA<2,1>(track);
          break;
        case 3:
          if (isFirstTrack) {histManager.FillEventQA<3>(coll, cent, nTracks);}
          histManager.FillTrackQA<3,1>(track);
          break;
        case 4:
          if (isFirstTrack) {histManager.FillEventQA<4>(coll, cent, nTracks);}
          histManager.FillTrackQA<4,1>(track);
          break;
        case 5:
          if (isFirstTrack) {histManager.FillEventQA<5>(coll, cent, nTracks);}
          histManager.FillTrackQA<5,1>(track);
          break;
        case 6:
          if (isFirstTrack) {histManager.FillEventQA<6>(coll, cent, nTracks);}
          histManager.FillTrackQA<6,1>(track);
          break;
        case 7:
          if (isFirstTrack) {histManager.FillEventQA<7>(coll, cent, nTracks);}
          histManager.FillTrackQA<7,1>(track);
          break;
        case 8:
          if (isFirstTrack) {histManager.FillEventQA<8>(coll, cent, nTracks);}
          histManager.FillTrackQA<8,1>(track);
          break;
        case 9:
          if (isFirstTrack) {histManager.FillEventQA<9>(coll, cent, nTracks);}
          histManager.FillTrackQA<9,1>(track);
          break;
        default:
          if (cfgDebugLog) {
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
    if (cfgUseEtaGap) {acAnalysis.ComputeCorrelEtaGap(trackPhi, trackWeight, trackEta);}    
    acAnalysis.ComputeAllCorrel(cBin, sampleID);

    /* Reset the variables for the next collision. */
    // This ensures no mixing between collision can happen accidentally.
    trackPhi.clear();
    trackWeight.clear();
    trackEta.clear();

    LOGF(info, "Collision analysed. Next...");
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowAC2hAnalysisTask>(cfgc)
  };
}
