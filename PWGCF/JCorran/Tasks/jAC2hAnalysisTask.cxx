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

/* Header files. */
#include <vector>
#include <TRandom3.h>

// O2 headers.
/// These two headers files are mandatory.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include <CCDB/BasicCCDBManager.h>
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"

// O2Physics headers.

// JCorran headers.
#include "PWGCF/JCorran/DataModel/JCatalyst.h"

#include "PWGCF/JCorran/Core/JAC2hAnalysis.h"
#include "PWGCF/JCorran/Core/JAC2hHistManager.h"

/* Namespaces. */
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
    /* Variables needed first in the histogram manager. */
    Configurable<int> cfgNcombis2h{"cfgNcombis2h", 3, "Number of pairs of harmos"};
    Configurable<int> cfgNsamples{"cfgNsamples", 20, "Number of bootstrap samples"};
    Configurable<float> cfgEtaGap{"cfgEtaGap", 1.0, "Minimum eta gap value"};

    // Registries for the output objects. They are saved in alphabetical order
    // (true) and in TDirectoryFile (true). Multiple registries are needed
    // to avoid going over the limit of histograms per registry (512).
    HistogramRegistry qaHistoRegistry{"qaHistoRegistry", {},
        OutputObjHandlingPolicy::AnalysisObject, true, true};
    HistogramRegistry acHistoRegistry{"acHistoRegistry", {},
        OutputObjHandlingPolicy::AnalysisObject, true, true};
    JAC2hHistManager histoManager;

    /* Variables needed for the analysis calculations. */
    Configurable<int> cfgDebugLvl{"cfgDebugLvl", 2, "Verbosity (none/low/high)"};
    Configurable<int> cfgMultiMin{"cfgMultiMin", 10, "Low multiplicity cut"};
    Configurable<std::vector<int>> cfg2hHarmos{"cfg2hHarmos", {23, 24, 34}, "List of 2h pairs"};
    JAC2hAnalysis acAnalysis;

    /// \brief O2 function executed at the beginning of the task workflow.
    void init(InitContext const&)
    {
        /* Reinitize gRandom to ensure a new sequence is used. */
        delete gRandom;
        gRandom = new TRandom3(0);

        /* Initialize the histogram manager with the values it needs to work. */
        histoManager.SetQAHistRegistry(&qaHistoRegistry);
        histoManager.SetACHistRegistry(&acHistoRegistry);
        histoManager.SetNcombis2h(cfgNcombis2h);
        histoManager.SetNsamples(cfgNsamples);
        histoManager.SetEtaGap(cfgEtaGap);
        histoManager.CreateQAHistos();
        histoManager.CreateACHistos();

        /* Initialize the analysis task with the values it needs to work. */
        acAnalysis.SetAC2hHistManager(histoManager);
        acAnalysis.SetDebugLevel(cfgDebugLvl);
        acAnalysis.Set2hPairs(cfg2hHarmos);

        if (cfgDebugLvl == 5) {     // High value to prevent any printing of this debug.
            const int nCentBins = sizeof(jflucCentBins)/sizeof(jflucCentBins[0]);
            printf("Number of centrality classes in JCatalyst: %d\n", nCentBins);
        }
    }

    /// \brief O2 function executed for each collision.
    void process(JCollision const& collision, aod::ParticleTrack const& tracks)
    {
        /* Reject the invalid events: no catalyst tracks or invalid centrality. */
        int nTracks = tracks.size();

        if (nTracks == 0) {return;}
        if ((collision.cbin() < 0) || (collision.cbin() > 8)) {return;}
        if (cfgDebugLvl >= 1) {
            printf("Number of tracks in this event: %d\n", nTracks);
            printf("Centrality bin: %d\n", collision.cbin());
        }

        /* Fill the centrality distribution straight out of the catalyst. */
        // No personal cut has been applied at this point, this is to
        // cross-check the catalyst output.
        qaHistoRegistry.fill(HIST("histCentCatalyst"), collision.cent());

        /* Apply any additional selection criteria specific to AC analyses. */
        // Low multiplicity cut applied to ensure properly defined correlators.
        if (nTracks < cfgMultiMin) {return;}

        /* Give a sample ID to this collision based on a uniform distribution. */
        int sampleID = static_cast<int>(gRandom->Uniform(0, cfgNsamples));
        if (cfgDebugLvl >= 1) {printf("Sample: %d\n", sampleID);}

        /* Save the angles and weights for the correlators, and fill the QA. */
        // The EventQA is filled for the first track only.
        bool isFirstTrack = true;
        std::vector<double> trackPhi;
        std::vector<double> trackWeight;

        for (auto& track : tracks) {
            trackPhi.push_back(track.phi());
            trackWeight.push_back(1./(track.weightEff() * track.weightNUA()));

            switch (collision.cbin()) {
            case 0:
                if (isFirstTrack) {histoManager.FillEventQA<0>(collision, nTracks, sampleID);}
                histoManager.FillTrackQA<0>(track);
                break;
            case 1:
                if (isFirstTrack) {histoManager.FillEventQA<1>(collision, nTracks, sampleID);}
                histoManager.FillTrackQA<1>(track);
                break;
            case 2:
                if (isFirstTrack) {histoManager.FillEventQA<2>(collision, nTracks, sampleID);}
                histoManager.FillTrackQA<2>(track);
                break;
            case 3:
                if (isFirstTrack) {histoManager.FillEventQA<3>(collision, nTracks, sampleID);}
                histoManager.FillTrackQA<3>(track);
                break;
            case 4:
                if (isFirstTrack) {histoManager.FillEventQA<4>(collision, nTracks, sampleID);}
                histoManager.FillTrackQA<4>(track);
                break;
            case 5:
                if (isFirstTrack) {histoManager.FillEventQA<5>(collision, nTracks, sampleID);}
                histoManager.FillTrackQA<5>(track);
                break;
            case 6:
                if (isFirstTrack) {histoManager.FillEventQA<6>(collision, nTracks, sampleID);}
                histoManager.FillTrackQA<6>(track);
                break;
            case 7:
                if (isFirstTrack) {histoManager.FillEventQA<7>(collision, nTracks, sampleID);}
                histoManager.FillTrackQA<7>(track);
                break;
            case 8:
                if (isFirstTrack) {histoManager.FillEventQA<8>(collision, nTracks, sampleID);}
                histoManager.FillTrackQA<8>(track);
                break;
            }

        // Indicate when the first track of the collision has been analysed.
        if (isFirstTrack) {isFirstTrack = false;}
        } // Go to the next track.

        /* Compute the multiparticle correlators, starting with the Q-vectors. */
        acAnalysis.CalculateQvectors(trackPhi, trackWeight);
        acAnalysis.ComputeAllTerms(collision.cbin(), sampleID);

        LOGF(info, "Collision analysed. Next...");
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jAC2hAnalysisTask>(cfgc, TaskName{"jac2h-analysis"})
  };
}