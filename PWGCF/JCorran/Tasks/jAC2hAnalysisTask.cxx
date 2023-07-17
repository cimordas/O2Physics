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
    myHistoManager.CreateHistosQA();

    // LOKI: added just to make the compilator happy...
    const int nCentBins = sizeof(jflucCentBins)/sizeof(jflucCentBins[0]);
    printf("Number of centrality classes in JCatalyst: %d\n", nCentBins);
  }

  void process(MyCollision const& collision, aod::ParticleTrack const& tracks)
  {
    // Reject the invalid events: no catalyst tracks or invalid centrality.
    if (tracks.size() == 0) {return;}
    if ((collision.cbin() < 0) || (collision.cbin() > 8)) {return;}

    // Fill the centrality distribution received from the catalyst.
    myHistoRegistry.fill(HIST("histCentBefore"), collision.cent());
  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jAC2hAnalysisTask>(cfgc, TaskName{"jac2h-analysis"})
  };
}