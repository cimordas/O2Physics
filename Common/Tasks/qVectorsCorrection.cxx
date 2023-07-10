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
/// \file   qVectorsCorrection.cxx
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

#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"

// o2 includes.

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;

namespace qV
{
  static constexpr std::string_view centClasses[] = {
    "Centrality_0-5/", "Centrality_5-10/", "Centrality_10-20/", "Centrality_20-30/",
    "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/", "Centrality_60-80/"
  };

  static constexpr std::string_view detNames[] = {"FT0-A", "FT0-C", "FV0"};
  static constexpr std::string_view histId[] = {"Before", "After"};
} // namespace qV

struct qVectorsCorrection {
  // Configurables.
  Configurable<std::string> cfgCentEsti{"cfgCentEsti",  // List from qVectorsTable.cxx
      "FT0-C", "Centrality estimator (Run3): 0 = FT0-M, 1 = FT0-A, 2 = FT0-C, 3 = FV0-A"};
  Configurable<std::string> cfgCorrStep{"cfgCorrStep",  // Used in the plotting.
      "Recenter", "Correction step: Recenter, Twist, Rescale"}; // No correction = recenter.

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
    const AxisSpec axisCent{100, 0., 100., fmt::format("Centrality percentile ({})",
                                                      (std::string)cfgCentEsti)};
    const AxisSpec axisQvec{400, -0.2, 0.2};
    const AxisSpec axisConst{8, 0., 8.};

    histosQA.add("histCentFull", "Centrality distribution for valid events",
        HistType::kTH1F, {axisCent});
    histosQA.add("Centrality_0-5/histCent", "Centrality distribution",
        HistType::kTH1F, {axisCent});

    histosQA.add("Centrality_0-5/histQvecFT0ABefore",
        ("(Qx,Qy) for " + (std::string)qV::detNames[0] + " before " + (std::string)cfgCorrStep).c_str(),
        HistType::kTH2D, {axisQvec, axisQvec});
    histosQA.add("Centrality_0-5/histQvecFT0CBefore",
        ("(Qx,Qy) for " + (std::string)qV::detNames[1] + " before " + (std::string)cfgCorrStep).c_str(),
          HistType::kTH2D, {axisQvec, axisQvec});
    histosQA.add("Centrality_0-5/histQvecFV0Before",
        ("(Qx,Qy) for " + (std::string)qV::detNames[2] + " before " + (std::string)cfgCorrStep).c_str(),
        HistType::kTH2D, {axisQvec, axisQvec});

    histosQA.add("Centrality_0-5/histQvecFT0AConst", "Correction constants for FT0A",
      HistType::kTH1F, {axisConst});
    histosQA.add("Centrality_0-5/histQvecFT0CConst", "Correction constants for FT0C",
      HistType::kTH1F, {axisConst});
    histosQA.add("Centrality_0-5/histQvecFV0Const", "Correction constants for FV0",
      HistType::kTH1F, {axisConst});

    histosQA.add("Centrality_0-5/histQvecFT0AAfter",
        ("(Qx,Qy) for " + (std::string)qV::detNames[0] + " after " + (std::string)cfgCorrStep).c_str(),
        HistType::kTH2D, {axisQvec, axisQvec});
    histosQA.add("Centrality_0-5/histQvecFT0CAfter",
        ("(Qx,Qy) for " + (std::string)qV::detNames[1] + " after " + (std::string)cfgCorrStep).c_str(),
          HistType::kTH2D, {axisQvec, axisQvec});
    histosQA.add("Centrality_0-5/histQvecFV0After",
        ("(Qx,Qy) for " + (std::string)qV::detNames[2] + " after " + (std::string)cfgCorrStep).c_str(),
        HistType::kTH2D, {axisQvec, axisQvec});

    for (int iBin = 1; iBin < 8; iBin++) {
      histosQA.addClone("Centrality_0-5/", qV::centClasses[iBin].data());
    }

  }   // End void init(InitContext const&)

  // Definition of all the needed functions.
  template <int bin, int id, typename T>
  void fillHistosQA(const T& vec)
  {
    // Fill the centrality distribution per class for the given bin.
    histosQA.fill(HIST(qV::centClasses[bin])+HIST("histCent"), vec.cent());

    // Fill the (Qx,Qy) distributions for each detector, after removing dummy values.
    /// NOTE: FV0 (and FT0C?) are not fully implemented yet
    /// --> Values are just dummy placeholders.
    if (TMath::Abs(vec.qvecFT0ARe()) < 900 && TMath::Abs(vec.qvecFT0AIm()) < 900) {
      histosQA.fill(HIST(qV::centClasses[bin])+HIST("histQvecFT0A")+HIST(qV::histId[id]),
          vec.qvecFT0ARe(), vec.qvecFT0AIm());
    }
    if (TMath::Abs(vec.qvecFT0CRe()) < 900 && TMath::Abs(vec.qvecFT0CIm()) < 900) {
      histosQA.fill(HIST(qV::centClasses[bin])+HIST("histQvecFT0C")+HIST(qV::histId[id]),
          vec.qvecFT0CRe(), vec.qvecFT0CIm());
    }
    if (TMath::Abs(vec.qvecFV0Re()) < 900 && TMath::Abs(vec.qvecFV0Im()) < 900) {
      histosQA.fill(HIST(qV::centClasses[bin])+HIST("histQvecFV0")+HIST(qV::histId[id]),
          vec.qvecFV0Re(), vec.qvecFV0Im());
    }
    LOGF(info, "QA has been filled");
  }

  //void process(aod::Qvector const& qVec)  // LOKI: Changed to all Qvectors as we need the full histo for the corrections.
  void process(aod::Qvectors const& qVecs)
  {
    // Iterate over the Qvectors table and fill the QA histograms with the uncorrected qVec.
    for (auto &qVec : qVecs) {
      // Get the right centrality bin, and fill the centrality QA histograms.
      int centBin = helperEP.GetCentBin(qVec.cent());
      LOGF(info, "Centrality percentile = %.0f Centrality bin: %d", qVec.cent(), centBin);
      histosQA.fill(HIST("histCentFull"), qVec.cent());
      
      if (centBin < 0 || centBin > 8) {continue;}
      switch(centBin) {
      case 0:
        fillHistosQA<0,0>(qVec);
        break;
      case 1:
        fillHistosQA<1,0>(qVec);
        break;
      case 2:
        fillHistosQA<2,0>(qVec);
        break;
      case 3:
        fillHistosQA<3,0>(qVec);
        break;
      case 4:
        fillHistosQA<4,0>(qVec);
        break;
      case 5:
        fillHistosQA<5,0>(qVec);
        break;
      case 6:
        fillHistosQA<6,0>(qVec);
        break;
      case 7:
        fillHistosQA<7,0>(qVec);
        break;
      } // End switch(centBin)
    
    } // Go to the next uncorrected qVec.

  // The uncorrected QA histograms can now be used to obtain the correction
  // constants for each centrality class and detector, before filling the
  // corresponding TH1 histograms.
  if ((std::string)cfgCorrStep == "Recenter") {  // Get the constants for the recentering.
    float meanX = 0.;
    float meanY = 0.;
    //const std::shared_ptr<TH2> hist = histosQA.get<TH2D>(HIST(qV::centClasses[0])+HIST("histQvecFT0ABefore"));
    const std::shared_ptr<TH2> hist = histosQA.get<TH2>(HIST("Centrality_0-5/histQvecFT0ABefore"));
    helperEP.GetCorrRecentering(hist, meanX, meanY);
    histosQA.get<TH1>(HIST("Centrality_0-5/histQvecFT0AConst"))->SetBinContent(1, meanX);
    histosQA.get<TH1>(HIST("Centrality_0-5/histQvecFT0AConst"))->GetXaxis()->SetBinLabel(1,"meanX");
    histosQA.get<TH1>(HIST("Centrality_0-5/histQvecFT0AConst"))->SetBinContent(2, meanY);
    histosQA.get<TH1>(HIST("Centrality_0-5/histQvecFT0AConst"))->GetXaxis()->SetBinLabel(2,"meanY");
    printf("MeanX: %.6f, MeanY: %.6f\n", meanX, meanY);
  }

  }   // End void process(...)

};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsCorrection>(cfgc)
  };
}