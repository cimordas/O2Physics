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

// Header files.
#include <TRandom3.h>

/* O2 headers. */
/// Those two headers files are mandatory.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include <CCDB/BasicCCDBManager.h>
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"

/* O2Physics headers. */

/* JCorran headers. */
#include "PWGCF/JCorran/DataModel/JCatalyst.h"

#include "PWGCF/JCorran/Core/JAC2hAnalysis.h"
#include "PWGCF/JCorran/Core/JAC2hHistManager.h"

// Namespaces.
using namespace o2;
using namespace o2::framework;
using namespace o2::analysis::PWGCF;

using JCollision = soa::Join<aod::Collisions, aod::CollisionData>::iterator;

/// #include "Framework/ASoAHelpers.h"
/// #include "Framework/RunningWorkflowInfo.h"

// #include "Common/DataModel/Centrality.h"
// #include "Common/DataModel/EventSelection.h"
// #include "Common/Core/TrackSelection.h"
// #include "Common/DataModel/TrackSelectionTables.h"
// #include "ReconstructionDataFormats/V0.h"

// ----------------------------------------------------------------------------
// The analysis task for the 2h AC starts here.
// ----------------------------------------------------------------------------
struct jAC2hAnalysisTask
{
    // Variables related to the histogram manager.
    Configurable<int> cfgMultiMin{"cfgMultiMin", 10, "Low multiplicity cut"};
    Configurable<int> cfgNsamples{"cfgNsamples", 20, "Number of bootstrap samples"};
    Configurable<float> cfgEtaGap{"cfgEtaGap", 1.0, "Minimum eta gap value"};

    /* Registry for the output objects. They are savev in alphabetical order
    (true) and not in TDirectoryFile (false). */
    HistogramRegistry myHistoRegistry{"myHistoRegistry", {},
        OutputObjHandlingPolicy::AnalysisObject, true, false};

    JAC2hHistManager myHistoManager;

    // Variables related to the analysis calculations.
    Configurable<int> cfgDebugLvl{"cfgDebugLvl", 2, "Verbosity (0: none, 1: low, 2: high)"};

    JAC2hAnalysis myACanalysis;
  
/*
  static constexpr int list2h[3][2] = {{2, 3}, {2, 4}, {3, 4}};
  Configurable<int> cfgNcombis2h{"cfgNcombis2h", 3, "Number of pairs of harmos"};
  Configurable<Array2D<int>> cfgList2h{"cfgList2h", {&list2h[0][0], 3, 2}, "List of 2h pairs"};
*/

    /// \brief Framework method where the two managers are initialized.
  void init(InitContext const&)
  {
    // Initialize the histograms within the registry itself.
    myHistoManager.SetHistManager(&myHistoRegistry);
    myHistoManager.SetNsamples(cfgNsamples);
    myHistoManager.SetEtaGap(cfgEtaGap);
    //myHistoManager.SetNcombis2h(cfgNcombis2h);
    myHistoManager.CreateHistos();

    // Initialize the analysis task.
    myACanalysis.SetAC2hHistManager(myHistoManager);
    myACanalysis.SetDebugLevel(cfgDebugLvl);
    myACanalysis.InitArrays();

    if (cfgDebugLvl > 2) {
      const int nCentBins = sizeof(jflucCentBins)/sizeof(jflucCentBins[0]); //LOKI
      printf("Number of centrality classes in JCatalyst: %d\n", nCentBins);
    }
  }

  void process(JCollision const& collision, aod::ParticleTrack const& tracks)
  {
    // Reject the invalid events: no catalyst tracks or invalid centrality.
    int nTracks = tracks.size();

    if (nTracks == 0) {return;}
    if ((collision.cbin() < 0) || (collision.cbin() > 8)) {return;}
/*    if (cfgDebugLvl == 1) {
      printf("Number of tracks in this event: %d\n", nTracks);
      printf("Centrality bin: %d\n", collision.cbin());
    }

    // Fill the centrality distribution received from the catalyst.
    // No personal cut has been applied at this point.
    myHistoRegistry.fill(HIST("histCentCatalyst"), collision.cent());

    // Apply any additional selection criteria: low multiplicity cut applied
    // to ensure properly defined correlators.
    if (nTracks < cfgMultiMin) {return;}

    // Attribute a bootstrap sample to this collision from a uniform distribution.
    int mySample = static_cast<int>(gRandom->Uniform(0, cfgNsamples));
    if (cfgDebugLvl == 1) {printf("Sample: %d\n", mySample);}

    // Save the azimuthal angles and particle weights in dynamic arrays, then
    // fill all the profiles and QA histograms. 
    // The EventQA is filled only for the first track of the loop.
    double *myPhi = new double[nTracks]();
    double *myPartWeights = new double[nTracks]();
    int iTrack = 0;
    bool isFirstTrack = true;

    for (auto& track : tracks) {
      myPhi[iTrack] = track.phi();
      myPartWeights[iTrack] = 1./(track.weightEff() * track.weightNUA());
      if (cfgDebugLvl == 2) {
        printf("iTrack: %d Phi: %.2f Weight: %.2f\n",
                iTrack, myPhi[iTrack], myPartWeights[iTrack]);
      }

      switch (collision.cbin()) {
      case 0:
        if (isFirstTrack) {myHistoManager.FillEventQA<0>(collision, nTracks, mySample);}
        myHistoManager.FillTrackQA<0>(track);
        break;
      case 1:
        if (isFirstTrack) {myHistoManager.FillEventQA<1>(collision, nTracks, mySample);}
        myHistoManager.FillTrackQA<1>(track);
        break;
      case 2:
        if (isFirstTrack) {myHistoManager.FillEventQA<2>(collision, nTracks, mySample);}
        myHistoManager.FillTrackQA<2>(track);
        break;
      case 3:
        if (isFirstTrack) {myHistoManager.FillEventQA<3>(collision, nTracks, mySample);}
        myHistoManager.FillTrackQA<3>(track);
        break;
      case 4:
        if (isFirstTrack) {myHistoManager.FillEventQA<4>(collision, nTracks, mySample);}
        myHistoManager.FillTrackQA<4>(track);
        break;
      case 5:
        if (isFirstTrack) {myHistoManager.FillEventQA<5>(collision, nTracks, mySample);}
        myHistoManager.FillTrackQA<5>(track);
        break;
      case 6:
        if (isFirstTrack) {myHistoManager.FillEventQA<6>(collision, nTracks, mySample);}
        myHistoManager.FillTrackQA<6>(track);
        break;
      case 7:
        if (isFirstTrack) {myHistoManager.FillEventQA<7>(collision, nTracks, mySample);}
        myHistoManager.FillTrackQA<7>(track);
        break;
      case 8:
        if (isFirstTrack) {myHistoManager.FillEventQA<8>(collision, nTracks, mySample);}
        myHistoManager.FillTrackQA<8>(track);
        break;
      }

      // Indicate when the first track of the collision has been analysed.
      if (isFirstTrack) {isFirstTrack = false;}
      iTrack++;
    } // Go to the next track.

    // Compute the multiparticle correlators, starting with the Q-vectors.
    myACanalysis.CalculateQvectors(nTracks, myPhi, myPartWeights);
    myACanalysis.ComputeAllTerms(collision.cbin(), mySample);

    // Reset all event-related quantities for the next collision.
    delete[] myPhi;
    delete[] myPartWeights;
    LOGF(info, "Collision analysed. Next...");
*/
  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jAC2hAnalysisTask>(cfgc, TaskName{"jac2h-analysis"})
  };
}