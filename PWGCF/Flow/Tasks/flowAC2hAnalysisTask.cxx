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

// \brief   Task for the calculation of the 2-harmonic asymmetric cumulants.
// \author  Cindy Mordasini (cindy.mordasini@cern.ch)

/* Header files. */
//#include <vector>
#include <TRandom3.h>

// O2 headers. //
// The first two are mandatory.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"

// O2 Physics headers. //
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGCF/Flow/Core/FlowAC2hHistManager.h"
//#include "PWGCF/Flow/Core/FlowAC2hAnalysis.h"

/* Namespaces and definitions. */
using namespace o2;
using namespace o2::framework;
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
  /* Configurables and variables for the histogram manager. */
  // Objects in the registries are saved in alphabetical order (true) and in
  // TDirectoryFile (true). Multiple registries are needed to avoid going over
  // the limit of objects per registry (max: 512).
  Configurable<int> cfgNcombis2h{"cfgNcombis2h", 3, "Number of harmonics pairs."};
  Configurable<int> cfgNsamples{"cfgNsamples", 20, "Number of bootstrap samples"};

  HistogramRegistry qaHistoRegistry{"qaHistoRegistry", {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true, true};
  //HistogramRegistry anHistoRegistry{"anHistoRegistry", {},
  //                                  OutputObjHandlingPolicy::AnalysisObject,
  //                                  true, true};
  FlowAC2hHistManager histManager;

  /* Configurables and variables for the data analysis. */
  Configurable<bool> cfgDebugLvl{"cfgDebugLvl", true, "Debug verbosity (true/false)"};
  Configurable<int> cfgCentEsti{"cfgCentEsti", 2,
                                "Centrality estimator (Run3), 0: FT0M, 1: FT0A, 2: FT0C, 3: FV0A, 4: FDDM, 5: NTPV"};
  Configurable<float> cfgZvtxMax{"cfgCutZvtx", 10.0f, "Maximum value for Z-vtx"};


/*
   Variables needed for the analysis calculations. 

  Configurable<int> cfgMultiMin{"cfgMultiMin", 10, "Low multiplicity cut"};
  Configurable<std::vector<int>> cfg2hHarmos{"cfg2hHarmos", {23, 24, 34}, "List of 2h pairs"};
  FlowAC2hAnalysis acAnalysis;
*/
  /// \brief O2 function executed at the beginning of the task workflow.
  void init(InitContext const&)
  {
    /* Reset gRandom to ensure a new sequence is used each time. */
    delete gRandom;
    gRandom = new TRandom3(0);

    /* Initialize the histogram manager with the values it requires to work. */
    histManager.SetHistRegistryQA(&qaHistoRegistry);
    //histoManager.SetANHistRegistry(&anHistoRegistry);
    histManager.SetNcombis2h(cfgNcombis2h);
    histManager.SetNsamples(cfgNsamples);
    //histoManager.SetEtaGap(cfgEtaGap);
    histManager.CreateHistsQA();
    //histoManager.CreateANHistos();
/*
    Initialize the analysis task with the variables it requires to work.
    acAnalysis.SetAC2hHistManager(histoManager);
    acAnalysis.SetDebugPrint(cfgDebugLvl);
    acAnalysis.Set2hPairs(cfg2hHarmos);
*/
  }

  /// \brief O2 function executed for each collision.
  void process(MyCollisions::iterator const& coll, MyTracks const& tracks)
  {
    /* Reject invalid events: no tracks or invalid centrality. */
    int nTracks = tracks.size();
    if (nTracks == 0) {return;}

    float centAllEstim[6] = {
      coll.centFT0M(), coll.centFT0A(), coll.centFT0C(),
      coll.centFV0A(), coll.centFDDM(), coll.centNTPV()
    };
    float cent = centAllEstim[cfgCentEsti];
    if (cent < 0. || cent > 100.) {return;}

    if (cfgDebugLvl) {
      printf("Number of tracks in this event: %d\n", nTracks);
      printf("Centrality value: %.2f\n", cent);
    }

    /* Apply the event selection to the collision. */
    

    /* Fill the event QA histograms before the event cuts.*/
    int cBin = histManager.GetCentBin(cent);
    switch (cBin) {
    case 0:
      histManager.FillEventQA<0,0>(coll, cent, nTracks);
      //if (isFirstTrack) {histoManager.FillEventQA<0>(collision, nTracks, sampleID);}
      //histoManager.FillTrackQA<0>(track);
      break;
    case 1:
      histManager.FillEventQA<1,0>(coll, cent, nTracks);
      //if (isFirstTrack) {histoManager.FillEventQA<1>(collision, nTracks, sampleID);}
      //histoManager.FillTrackQA<1>(track);
      break;
    case 2:
      histManager.FillEventQA<2,0>(coll, cent, nTracks);
      //if (isFirstTrack) {histoManager.FillEventQA<2>(collision, nTracks, sampleID);}
      //histoManager.FillTrackQA<2>(track);
      break;
    case 3:
      histManager.FillEventQA<3,0>(coll, cent, nTracks);
      //if (isFirstTrack) {histoManager.FillEventQA<3>(collision, nTracks, sampleID);}
      //histoManager.FillTrackQA<3>(track);
      break;
    case 4:
      histManager.FillEventQA<4,0>(coll, cent, nTracks);
      //if (isFirstTrack) {histoManager.FillEventQA<4>(collision, nTracks, sampleID);}
      //histoManager.FillTrackQA<4>(track);
      break;
    case 5:
      histManager.FillEventQA<5,0>(coll, cent, nTracks);
      //if (isFirstTrack) {histoManager.FillEventQA<5>(collision, nTracks, sampleID);}
      //histoManager.FillTrackQA<5>(track);
      break;
    case 6:
      histManager.FillEventQA<6,0>(coll, cent, nTracks);
      //if (isFirstTrack) {histoManager.FillEventQA<6>(collision, nTracks, sampleID);}
      //histoManager.FillTrackQA<6>(track);
      break;
    case 7:
      histManager.FillEventQA<7,0>(coll, cent, nTracks);
      //if (isFirstTrack) {histoManager.FillEventQA<7>(collision, nTracks, sampleID);}
      //histoManager.FillTrackQA<7>(track);
      break;
    case 8:
      histManager.FillEventQA<8,0>(coll, cent, nTracks);
      //if (isFirstTrack) {histoManager.FillEventQA<8>(collision, nTracks, sampleID);}
      //histoManager.FillTrackQA<8>(track);
      break;
    case 9:
      histManager.FillEventQA<9,0>(coll, cent, nTracks);
      //if (isFirstTrack) {histoManager.FillEventQA<8>(collision, nTracks, sampleID);}
      //histoManager.FillTrackQA<8>(track);
      break;
    default:
      LOGF(info, "Centrality percentile not included in analysis. Next...");
      break;
    }


  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowAC2hAnalysisTask>(cfgc)
  };
}
