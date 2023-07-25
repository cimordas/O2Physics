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

///
/// \file   evtPlanesResolution.cxx
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  ...
///

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/EvtPlanes.h"
#include "Common/Core/EventPlaneHelper.h"

// o2 includes.

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TProfile.h>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;

namespace ep
{
  static constexpr std::string_view centClasses[] = {
    "Centrality_0-5/", "Centrality_5-10/", "Centrality_10-20/", "Centrality_20-30/",
    "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/", "Centrality_60-80/"
  };

  static constexpr std::string_view detNames[] = {
    "FT0A", "FT0C", "FV0A", "BPos", "BNeg"
  };
} // namespace ep

struct evtPlanesResolution {
  // Configurables.

  // Histogram registry for the output QA figures and list of centrality classes for it.
  // Objects are NOT saved in alphabetical orders, and registry names are NOT saved
  // as TDirectoryFile.
  HistogramRegistry histosQA{"histosQA", {},
                            OutputObjHandlingPolicy::AnalysisObject,
                            false, false};

  // Helper variables.
  EventPlaneHelper helperEP;


  void init(InitContext const&)
  {
    // Fill the registry with the needed objects.
    const AxisSpec axisEP{3, 0., 3.}; // Bins: A-B, A-C, B-C
    histosQA.add("Centrality_0-5/profCosFT0A",
      ("Average cosines for " + (std::string)ep::detNames[0]).c_str(),
      HistType::kTProfile, {axisEP});
    histosQA.add("Centrality_0-5/profCosFT0C",
      ("Average cosines for " + (std::string)ep::detNames[1]).c_str(),
      HistType::kTProfile, {axisEP});
    histosQA.add("Centrality_0-5/profCosFV0A",
      ("Average cosines for " + (std::string)ep::detNames[2]).c_str(),
      HistType::kTProfile, {axisEP});

    for (int iBin = 1; iBin < 8; iBin++) {
      histosQA.addClone("Centrality_0-5/", ep::centClasses[iBin].data());
    }
  }   // End void init(InitContext const&)

  // Definition of all the needed template functions.
  template <int bin, int fDet, typename T>
  void fillProfCos(const T& psiA, const T& psiB, const T& psiC)
  {
    // Fill the different cosines for the 3-subevents method.
    /// NOTE: FV0 (and FT0C?) are not fully implemented yet
    /// --> Values are just dummy placeholders.
    histosQA.fill(HIST(ep::centClasses[bin])+HIST("profCos")+HIST(ep::detNames[fDet]),
      0.5, TMath::Cos(2.*(psiA-psiB)));
    histosQA.fill(HIST(ep::centClasses[bin])+HIST("profCos")+HIST(ep::detNames[fDet]),
      1.5, TMath::Cos(2.*(psiA-psiC)));
    histosQA.fill(HIST(ep::centClasses[bin])+HIST("profCos")+HIST(ep::detNames[fDet]),
      2.5, TMath::Cos(2.*(psiB-psiC)));

    LOGF(info, "The profile has been filled with a new event.");
  }

  void process(aod::EvtPlanes const& evPls)
  {
    // Iterate over the table and fill the profiles with the combinations of evPl.
    for (auto &evPl : evPls) {
      // Get the centrality bin.
      int centBin = helperEP.GetCentBin(evPl.cent());
      LOGF(info, "Centrality percentile = %.1f Centrality bin: %d", evPl.cent(), centBin);

      // For each FIT subdetector, fill the profiles for the resolution.
      float epFIT[3] = {evPl.evtPlFT0A(), evPl.evtPlFT0C(), evPl.evtPlFV0A()};
      static_for<0, 2>([&](auto iFIT) {
        constexpr int detFIT = iFIT.value;
        switch(centBin) {
        case 0:
          fillProfCos<0, detFIT>(epFIT[detFIT], evPl.evtPlBPos(), evPl.evtPlBNeg());
          break;
        case 1:
          fillProfCos<1, detFIT>(epFIT[detFIT], evPl.evtPlBPos(), evPl.evtPlBNeg());
          break;
        case 2:
          fillProfCos<2, detFIT>(epFIT[detFIT], evPl.evtPlBPos(), evPl.evtPlBNeg());
          break;
        case 3:
          fillProfCos<3, detFIT>(epFIT[detFIT], evPl.evtPlBPos(), evPl.evtPlBNeg());
          break;
        case 4:
          fillProfCos<4, detFIT>(epFIT[detFIT], evPl.evtPlBPos(), evPl.evtPlBNeg());
          break;
        case 5:
          fillProfCos<5, detFIT>(epFIT[detFIT], evPl.evtPlBPos(), evPl.evtPlBNeg());
          break;
        case 6:
          fillProfCos<6, detFIT>(epFIT[detFIT], evPl.evtPlBPos(), evPl.evtPlBNeg());
          break;
        case 7:
          fillProfCos<7, detFIT>(epFIT[detFIT], evPl.evtPlBPos(), evPl.evtPlBNeg());
          break;
        } // End switch(centBin)
      });

    } // Go to the next evPl.

  }   // End void process(...)
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<evtPlanesResolution>(cfgc)
  };
}