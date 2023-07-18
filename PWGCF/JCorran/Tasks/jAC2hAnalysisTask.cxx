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

/// \brief  Analysis task for the calculation of the 2-harmonic AC.
/// \author Cindy Mordasini (cindy.mordasini@cern.ch)
/// \since  July 2023

// Includes.
#include "TRandom3.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include <CCDB/BasicCCDBManager.h>
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/V0.h"

#include "PWGCF/JCorran/DataModel/JCatalyst.h"
#include "PWGCF/JCorran/Core/JAC2hHistManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::analysis::PWGCF;

using MyCollision = soa::Join<aod::Collisions, aod::CollisionData>::iterator;

struct jAC2hAnalysisTask {
  // Configurables.
  Configurable<int> cfgNsamples{"cfgNsamples", 10, "Number of bootstrap samples"};
  Configurable<int> cfgMultiMin{"cfgMultiMin", 10, "Low multiplicity cut"};

  // Filters.

  // Histogram registry for the output objects. They are not saved in
  // alphabetical order (false), and the names are not saved as TDirectoryFile
  // (false).
  HistogramRegistry myHistoRegistry{"myHistoRegistry", {},
                                  OutputObjHandlingPolicy::AnalysisObject,
                                  false, false};

  // Variables.
  JAC2hHistManager myHistoManager;

  void init(InitContext const&)
  {
    // Initialize the histograms within the registry itself.
    myHistoManager.SetHistManager(&myHistoRegistry);
    myHistoManager.SetNsamples(cfgNsamples);

    myHistoManager.CreateHistosQA();
    myHistoManager.CreateHistosAN();

    // LOKI: added just to make the compilator happy...
    const int nCentBins = sizeof(jflucCentBins)/sizeof(jflucCentBins[0]);
    printf("Number of centrality classes in JCatalyst: %d\n", nCentBins);
  }

  void process(MyCollision const& collision, aod::ParticleTrack const& tracks)
  {
    // Reject the invalid events: no catalyst tracks or invalid centrality.
    int nTracks = tracks.size();

    if (nTracks == 0) {return;}
    if ((collision.cbin() < 0) || (collision.cbin() > 8)) {return;}
    printf("Number of tracks in this event: %d\n", nTracks);
    printf("Centrality bin: %d\n", collision.cbin());

    // Fill the centrality distribution received from the catalyst.
    // No personal cut has been applied at this point.
    myHistoRegistry.fill(HIST("histCentCatalyst"), collision.cent());

    // Apply any additional selection criteria: low multiplicity cut applied
    // to ensure properly defined correlators.
    if (nTracks < cfgMultiMin) {return;}

    // Attribute a bootstrap sample to this collision from a uniform distribution.
    int mySample = static_cast<int>(gRandom->Uniform(0, cfgNsamples));
    printf("Sample: %d\n", mySample);

    // Fill the remaining QA histograms and 
    // The EventQA is filled only for the first track of the loop.
    bool isFirstTrack = true;
    for (auto& track : tracks) {

      switch (collision.cbin()) {
      case 0:
        if (isFirstTrack) {myHistoManager.FillEventQA<0>(collision, nTracks);}
        myHistoManager.FillTrackQA<0>(track);
        break;
      case 1:
        if (isFirstTrack) {myHistoManager.FillEventQA<1>(collision, nTracks);}
        myHistoManager.FillTrackQA<1>(track);
        break;
      case 2:
        if (isFirstTrack) {myHistoManager.FillEventQA<2>(collision, nTracks);}
        myHistoManager.FillTrackQA<2>(track);
        break;
      case 3:
        if (isFirstTrack) {myHistoManager.FillEventQA<3>(collision, nTracks);}
        myHistoManager.FillTrackQA<3>(track);
        break;
      case 4:
        if (isFirstTrack) {myHistoManager.FillEventQA<4>(collision, nTracks);}
        myHistoManager.FillTrackQA<4>(track);
        break;
      case 5:
        if (isFirstTrack) {myHistoManager.FillEventQA<5>(collision, nTracks);}
        myHistoManager.FillTrackQA<5>(track);
        break;
      case 6:
        if (isFirstTrack) {myHistoManager.FillEventQA<6>(collision, nTracks);}
        myHistoManager.FillTrackQA<6>(track);
        break;
      case 7:
        if (isFirstTrack) {myHistoManager.FillEventQA<7>(collision, nTracks);}
        myHistoManager.FillTrackQA<7>(track);
        break;
      case 8:
        if (isFirstTrack) {myHistoManager.FillEventQA<8>(collision, nTracks);}
        myHistoManager.FillTrackQA<8>(track);
        break;
      }

      // Indicate when the first track of the collision has been analysed.
      if (isFirstTrack) {isFirstTrack = false;}
    } // Go to the next track.


  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jAC2hAnalysisTask>(cfgc, TaskName{"jac2h-analysis"})
  };
}