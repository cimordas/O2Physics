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
#include <string>
#include <vector>

#include "PWGCF/JCorran/Core/JAC2hHistManager.h"

using namespace o2;
using namespace o2::framework;
namespace o2::analysis::PWGCF
{
/// \brief .
/// .
void JAC2hHistManager::CreateHistos()
{
  if (!mHistoRegistry) {
    LOGF(error, "No histogram manager provided. Quit.");
    return;
  }

  // All histograms are defined for the first centrality class in details, with
  // callSumw2 set to true for all. They will be cloned later for all the other
  // classes.
  /// Centrality distributions.
  const AxisSpec axisCent{100, 0., 100., "Centrality percentile"};
  mHistoRegistry->add("histCentCatalyst", "Centrality from the JCatalyst",
                      HistType::kTH1F, {axisCent}, true);
  mHistoRegistry->add("Centrality_00-01/QA/histCent", "Centrality after own cuts",
                      HistType::kTH1F, {axisCent}, true);

  /// Event distributions.
  const AxisSpec axisZvtx{75, -15., 15., "Z_{vtx} [cm]"};
  mHistoRegistry->add("Centrality_00-01/QA/histZvtx", "Z_{vtx} after own cuts",
                      HistType::kTH1F, {axisZvtx}, true);

  const AxisSpec axisMulti{1000, 0., 5000., "N_{tracks}"};
  mHistoRegistry->add("Centrality_00-01/QA/histMulti", "Multiplicity after own cuts",
                      HistType::kTH1I, {axisMulti}, true);

  /// Number of events for each bootstrap sample.
  const AxisSpec axisSamples{mNsamples, 0., (double)mNsamples, "Sample"};
  mHistoRegistry->add("Centrality_00-01/QA/histSamples", "N_{events} in each bootstrap sample",
                      HistType::kTH1I, {axisSamples}, true);
  for (int i = 1; i <= mNsamples; i++) {
    mHistoRegistry->get<TH1>(HIST("Centrality_00-01/QA/histSamples"))->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
  }

  /// Track distributions.
  std::vector<double> ptBinning = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                                  0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                  0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4,
                                  1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6,
                                  2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.5, 5., 6.};
  const AxisSpec axisPt = {ptBinning, "#it{p}_{T} [GeV/#it{c}]"};
  mHistoRegistry->add("Centrality_00-01/QA/histPtUncorrected", "#it{p}_{T} (not NUE corrected)",
                      HistType::kTH1F, {axisPt}, true);
  mHistoRegistry->add("Centrality_00-01/QA/histPtCorrected", "#it{p}_{T} (NUE corrected)",
                      HistType::kTH1F, {axisPt}, true);

  const AxisSpec axisEta = {20, -1., 1., "#eta"};
  mHistoRegistry->add("Centrality_00-01/QA/histEta", "Pseudorapidity",
                      HistType::kTH1F, {axisEta}, true);

  const AxisSpec axisPhi = {100, 0., 2.*M_PI, "#varphi"};
  mHistoRegistry->add("Centrality_00-01/QA/histPhiUncorrected", "Azimuthal angles (not NUA corrected)",
                      HistType::kTH1F, {axisPhi}, true);
  mHistoRegistry->add("Centrality_00-01/QA/histPhiCorrected", "Azimuthal angles (NUA corrected)",
                      HistType::kTH1F, {axisPhi}, true);

/*  Commented till the charge is implemented in the catalyst.
  const AxisSpec axisCharge = {2, -1., 1., "Charge"};
  mHistoRegistry->add("Centrality_00-01/QA/histCharge", "Electric charge",
                      HistType::kTH1I, {axisCharge}, true);
*/

  /// Profiles for the 2-particle and 2-harmonic terms.
  const AxisSpec axis2pCorrel{6, 0., 6., "n"};
  mHistoRegistry->add("Centrality_00-01/Full/prof2pCorrel", "<v_{n}^{2}>",
                      HistType::kTProfile, {axis2pCorrel}, true);
  mHistoRegistry->add("Centrality_00-01/Full/prof2pCorrelEta",
                      ("<v_{n}^{2}> |#it{#Delta #eta}| > " + std::to_string(mEtaGap)).c_str(),
                      HistType::kTProfile, {axis2pCorrel}, true);
  for (int i = 1; i <= 6; i++) {
    mHistoRegistry->get<TProfile>(HIST("Centrality_00-01/Full/prof2pCorrel"))->GetXaxis()->SetBinLabel(i, Form("v_{%d}", i));
    mHistoRegistry->get<TProfile>(HIST("Centrality_00-01/Full/prof2pCorrelEta"))->GetXaxis()->SetBinLabel(i, Form("v_{%d}", i));
  }

  const AxisSpec axis2hCorrel{14, 0., 14., "{a,b}"};
  mHistoRegistry->add("Centrality_00-01/Full/prof2hCorrel", "<v_{m}^{2a}v_{n}^{2b}>",
                      HistType::kTProfile, {axis2hCorrel}, true);

  // Clone the full profiles for all the samples.
  for (int iS = 0; iS < mNsamples; iS++) {
    std::string strSample = Form("Centrality_00-01/Sample%d/", iS);
    if (iS < 10) {strSample = Form("Centrality_00-01/Sample0%d/", iS);}
    mHistoRegistry->addClone("Centrality_00-01/Full/", strSample.data());
  } // Go to the next sample.

  // Clone the content of Centrality_00-01 into the other centrality classes.
  for (int iBin = 1; iBin < 9; iBin++) {
    mHistoRegistry->addClone("Centrality_00-01/", mCentClasses[iBin].data());
  }

}

} // namespace o2::analysis::PWGCF