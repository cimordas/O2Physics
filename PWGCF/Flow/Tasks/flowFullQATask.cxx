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

// \brief   .
// \author  Cindy Mordasini (cindy.mordasini@cern.ch)

// Standard headers.
#include <chrono>
#include <string>
#include <vector>

// O2 headers.
// The first two are mandatory.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"

// O2 Physics headers.
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGCF/Flow/Core/FlowAC2hHistManager.h"

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
// QA task starts here.
// ----------------------------------------------------------------------------
struct flowFullQATask
{
  // Declare the histogram registry and instance of its manager.
  // Objects in the registries are saved in alphabetical order (true) and in
  // TDirectoryFile (true).
  FlowAC2hHistManager histManager;
  HistogramRegistry qaHistRegistry{"qaHistRegistry", {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true, true};

  // Enable the general features of the analysis task.
  Configurable<bool> cfgDebugLog{"cfgDebugLog", true, "Enable log for debug."};
  Configurable<bool> cfgObtainNUA{"cfgObtainNUA", false, "Enable the 3D histogram for NUA."};
  Configurable<bool> cfgSaveQABefore{"cfgSaveQABefore", true,
                                     "Enable the QA before any selection."};

  // Set the event quality cuts.
  // The centrality estimators are the ones available for Run 3.
  enum centEstimators {FT0M, FT0A, FT0C, FDDM, NTPV};
  Configurable<int> cfgCentEsti{"cfgCentEsti", 2, "Centrality estimator."};
  Configurable<int> cfgMultMin{"cfgMultMin", 10,
                               "Strict minimum number of tracks per collision."};
  Configurable<float> cfgZvtxMax{"cfgCutZvtx", 10.0f, "Maximum range for Zvtx."};

  // Set the track quality cuts. 
  Configurable<bool> cfgUseNUA{"cfgUseNUA", true, "Enable the use of NUA weights."};
  Configurable<bool> cfgUseNUE{"cfgUseNUE", true, "Enable the use of NUE weights."};
  Configurable<bool> cfgUseVariablePtBins{"cfgUseVariablePtBins", false,
                                          "Enable the use of variable pT bins."};
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

  /// \brief O2 function executed at the beginning of the task workflow.
  void init(InitContext const&)
  {
    // Initialise the histogram manager for the QA objects.
    histManager.SetHistRegistryQA(&qaHistRegistry);
    histManager.SetDebugLog(cfgDebugLog);
    histManager.SetObtainNUA(cfgObtainNUA);
    histManager.SetSaveAllQA(true);   // Full QA is always saved by default.
    histManager.SetSaveQABefore(cfgSaveQABefore);
    histManager.SetUseVariablePtBins(cfgUseVariablePtBins);
    histManager.CreateHistQA();

    // Setup the access to the CCDB objects.
    ccdb->setURL(cfgCCDB.cfgURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(cfgCCDB.cfgTime.value);
  }

  /// \brief O2 function executed for each collision.
  void process(MyCollisions::iterator const& coll, MyTracks const& tracks)
  {
    // We don't process further the collisions without tracks or with
    // a centrality percentile outside of the range 0-70%.
    if (tracks.size() == 0) {return;}
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

    // We now perform the event selection to keep only the 'good' collisions
    // for further processing. The EventQA is filled before and after the
    // selection.
    if (cfgSaveQABefore) {histManager.FillEventQA<0>(coll, cBin, cent, nTracks);}
    if (abs(coll.posZ()) > cfgZvtxMax) {return;}
    if (!coll.sel8()) {return;} // TODO: Make it configurable.
    histManager.FillEventQA<1>(coll, cBin, cent, nTracks);

    for (auto& track : tracks) {
      // We look at the track QA only for the 'good' collisions only. First, we
      // fill the TrackQA before selection.
      if (cfgSaveQABefore) {histManager.FillTrackQA<0>(track, cBin);}

      // We now apply the full track selection on the current track. It is marked as
      // 'good' only if is passes all of it. In this case, the TrackQA is also
      // filled for After/ once the weights have been obtained if needed..
      bool isTrackGood = (track.pt() > cfgPtMin) && (track.pt() < cfgPtMax)
                      && (abs(track.eta()) < cfgEtaMax);

      float wNUE = 1.;  // NUE weight straight from weight distribution.
                        // The inverse will be applied to pT.
      if (cfgUseNUE) {
        wNUE = 2.;  // TODO: Implement fetching of the weights.
      }

      float wNUA = 1.;  // NUA weight straight from weight distribution.
                        // The inverse will be applied to phi.
      if (cfgUseNUA) {
        wNUA = 2.;  // TODO: Implement fetching of the weights.
      }

      if (isTrackGood) {histManager.FillTrackQA<1>(track, cBin, wNUE, wNUA, coll.posZ());}
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowFullQATask>(cfgc)
  };
}
