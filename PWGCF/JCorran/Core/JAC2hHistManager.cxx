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

/// \brief  Histogram manager for the AC analysis.
/// \author Cindy Mordasini (cindy.mordasini@cern.ch)
/// \since  July 2023

// Includes.
#include <vector>

#include "PWGCF/JCorran/Core/JAC2hHistManager.h"

using namespace o2;
using namespace o2::framework;
namespace o2::analysis::PWGCF
{
/*
static constexpr std::string_view centClasses[] = {  // Classes from JCatalyst.
  "Centrality_0-1/", "Centrality_1-2/", "Centrality_2-5/", "Centrality_5-10/",
  "Centrality_10-20/", "Centrality_20-30/", "Centrality_30-40/", "Centrality_40-50/",
  "Centrality_50-60/"
};*/

/// \brief .
/// .
void JAC2hHistManager::CreateHistosQA()
{
  if (!mHistoRegistry) {
    LOGF(error, "No histogram manager provided. Quit.");
    return;
  }

  const int nCentBins = sizeof(jflucCentBins)/sizeof(jflucCentBins[0]); //LOKI
  printf("Number of centrality classes in JCatalyst: %d\n", nCentBins);

  // All histograms are defined for the first centrality class in details, with
  // callSumw2 set to true for all. They will be cloned later for all the other
  // classes.
  /// Centrality distributions.
  const AxisSpec axisCent{100, 0., 100., "Centrality percentile"};
  mHistoRegistry->add("histCentCatalyst", "Centrality from the JCatalyst",
                      HistType::kTH1F, {axisCent}, true);
  mHistoRegistry->add("Centrality_0-1/histCent", "Centrality after own cuts",
                      HistType::kTH1F, {axisCent}, true);

  /// Event distributions.
  const AxisSpec axisZvtx{75, -15., 15., "Z_{vtx} [cm]"};
  mHistoRegistry->add("Centrality_0-1/histZvtx", "Z_{vtx} after own cuts",
                      HistType::kTH1F, {axisZvtx}, true);

  const AxisSpec axisMulti{1000, 0., 5000., "N_{tracks}"};
  mHistoRegistry->add("Centrality_0-1/histMulti", "Multiplicity after own cuts",
                      HistType::kTH1I, {axisMulti}, true);

  /// Track distributions.
  std::vector<double> ptBinning = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                                  0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                  0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4,
                                  1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6,
                                  2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.5, 5., 6.};
  const AxisSpec axisPt = {ptBinning, "#it{p}_{T} [GeV/#it{c}]"};
  mHistoRegistry->add("Centrality_0-1/histPtUncorrected", "#it{p}_{T} (not NUE corrected)",
                      HistType::kTH1F, {axisPt}, true);
  mHistoRegistry->add("Centrality_0-1/histPtCorrected", "#it{p}_{T} (NUE corrected)",
                      HistType::kTH1F, {axisPt}, true);

  const AxisSpec axisEta = {20, -1., 1., "#eta"};
  mHistoRegistry->add("Centrality_0-1/histEta", "Pseudorapidity",
                      HistType::kTH1F, {axisEta}, true);

  const AxisSpec axisPhi = {100, 0., 2.*M_PI, "#varphi"};
  mHistoRegistry->add("Centrality_0-1/histPhiUncorrected", "Azimuthal angles (not NUA corrected)",
                      HistType::kTH1F, {axisPhi}, true);
  mHistoRegistry->add("Centrality_0-1/histPhiCorrected", "Azimuthal angles (NUA corrected)",
                      HistType::kTH1F, {axisPhi}, true);

/*  Commented till the charge is implemented in the catalyst.
  const AxisSpec axisCharge = {2, -1., 1., "Charge"};
  mHistoRegistry->add("Centrality_0-1/histCharge", "Electric charge",
                      HistType::kTH1I, {axisCharge}, true);
*/

  // Clone the content of Centrality_0-1 into the other centrality classes.
  for (int iBin = 1; iBin < 9; iBin++) {
    mHistoRegistry->addClone("Centrality_0-1/", mCentClasses[iBin].data());
  }
}

/// \brief .
/// .
void JAC2hHistManager::CreateHistosAN()
{
  if (!mHistoRegistry) {
    LOGF(error, "No histogram manager provided. Quit.");
    return;
  }

  // All objects are defined for the first centrality class in details, with
  // callSumw2 set to true for all. They will be cloned later for all the other
  // classes.
  /// Number of events for each bootstrap sample.
  const AxisSpec axisSamples{mNsamples, 0., (double)mNsamples, "Sample #"};
  mHistoRegistry->add("Centrality_0-1/histSamples", "N_{events} in each bootstrap sample",
                      HistType::kTH1I, {axisSamples}, true);

  /// Profiles for the 2-particle and 2-harmonic terms.
  const AxisSpec axis2pCorrel{8, 0., 8., "n"};
  mHistoRegistry->add("Centrality_0-1/hist2pCorrel_Full", "<v_{n}^{2}>",
                      HistType::kTProfile, {axis2pCorrel}, true);

  const AxisSpec axis2hCorrel{14, 0., 14., "{a,b}"};
  mHistoRegistry->add("Centrality_0-1/hist2hCorrel_Full", "<v_{m}^{2a}v_{n}^{2b}>",
                      HistType::kTProfile, {axis2hCorrel}, true);

  // Clone the full profiles for all the samples.
  for (int iS = 0; iS < mNsamples; iS++) {
    std::string strSample = Form("_Sample%d", iS);
    mHistoRegistry->addClone("Full", strSample.data());
  } // Go to the next sample.

  // Clone the content of Centrality_0-1 into the other centrality classes.
  for (int iBin = 1; iBin < 9; iBin++) {
    mHistoRegistry->addClone("Centrality_0-1/", mCentClasses[iBin].data());
  }

}

} // namespace o2::analysis::PWGCF