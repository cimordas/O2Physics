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

// Header files.

// O2 headers.

// O2 Physics headers.
#include "PWGCF/Flow/Core/FlowJAC2hHistManager.h"

// Namespaces.
using namespace o2;
using namespace o2::framework;

namespace o2::analysis::PWGCF
{
/// \brief Create the histograms in the analysis registry.
void FlowJAC2hHistManager::CreateHistAN()
{
  /* Security checks. */
  if (!mHistRegistryAN) {
    LOGF(error, "AN histogram registry missing. Quitting...");
    return;
  }

  /* Define the histogram bookkeeping the individual harmonics from each pair. */
  const AxisSpec axisCombis{2*mNcombis2h, 0., 2.*(double)mNcombis2h, "Pairs (m,n)"};
  mHistRegistryAN->add("histCombis", "Configured pairs of harmonics",
                        HistType::kTProfile, {axisCombis}, true);
  for (int i = 0; i < mNcombis2h; i++) {
    mHistRegistryAN->get<TProfile>(HIST("histCombis"))->GetXaxis()
      ->SetBinLabel(2*i+1, Form("Combi%d-m", i));
    mHistRegistryAN->get<TProfile>(HIST("histCombis"))->GetXaxis()
      ->SetBinLabel(2*(i+1), Form("Combi%d-n", i));
  }

  /* Define the histogram bookkeeping the number of events per samples. */
  const AxisSpec axisSamples{mNsamples+1, 0., (double)mNsamples+1., "Sample ID"};
  mHistRegistryAN->add("Centrality_00-01/histSamples", "N_{events} in each bootstrap sample",
                       HistType::kTH1I, {axisSamples}, true);
  mHistRegistryAN->get<TH1>(HIST("Centrality_00-01/histSamples"))
    ->GetXaxis()->SetBinLabel(1, "Full");
  for (int i = 2; i <= mNsamples+1; i++) {
    mHistRegistryAN->get<TH1>(HIST("Centrality_00-01/histSamples"))
      ->GetXaxis()->SetBinLabel(i, Form("%d", i-2));
  }

  /* Define the TProfile2D used to save the multiparticle correlators. */
  // The x-axis keeps the correlators, the y-axis the full data and samples.
  // All profiles are defined for the first centrality class in details,
  // then cloned for the other combinations and centrality classes.
  const AxisSpec axis2pCorrelX{12, 0., 12., "n"};
  const AxisSpec axis2hCorrelX{14, 0., 14., "Powers"};
  mHistRegistryAN->add("Centrality_00-01/prof2pCorrel", "<v_{n}^{2}>",
                       HistType::kTProfile2D, {axis2pCorrelX, axisSamples}, true);
  mHistRegistryAN->add("Centrality_00-01/prof2hCorrelCombi0",
                       "<v_{m}^{2a}v_{n}^{2b}>",
                       HistType::kTProfile2D, {axis2hCorrelX, axisSamples}, true);

  // The x-axis of 'prof2pCorrel' follows <n, n with gap> for n = 1...6.
  for (int i = 1; i <= 6; i++) {
    mHistRegistryAN->get<TProfile2D>(HIST("Centrality_00-01/prof2pCorrel"))
      ->GetXaxis()->SetBinLabel(2*i-1, Form("v_{%d}", i));
    mHistRegistryAN->get<TProfile2D>(HIST("Centrality_00-01/prof2pCorrel"))
      ->GetXaxis()->SetBinLabel(2*i, Form("v_{%d, #Delta#eta}", i));
  }

  for (int i = 1; i <= 14; i++) {
    mHistRegistryAN->get<TProfile2D>(HIST("Centrality_00-01/prof2hCorrelCombi0"))
      ->GetXaxis()->SetBinLabel(i, Form("%s", mPowers[i-1].data()));
  }

  // The first y-bin is always the full dataset.
  mHistRegistryAN->get<TProfile2D>(HIST("Centrality_00-01/prof2pCorrel"))
    ->GetYaxis()->SetBinLabel(1, "Full");
  mHistRegistryAN->get<TProfile2D>(HIST("Centrality_00-01/prof2hCorrelCombi0"))
    ->GetYaxis()->SetBinLabel(1, "Full");
  for (int i = 2; i <= mNsamples+1; i++) {
    mHistRegistryAN->get<TProfile2D>(HIST("Centrality_00-01/prof2pCorrel"))
      ->GetYaxis()->SetBinLabel(i, Form("%d", i-2));
    mHistRegistryAN->get<TProfile2D>(HIST("Centrality_00-01/prof2hCorrelCombi0"))
      ->GetYaxis()->SetBinLabel(i, Form("%d", i-2));
  }

  for (int iH = 1; iH < mNcombis2h; iH++) {
    std::string strCombi = Form("Centrality_00-01/prof2hCorrelCombi%d", iH);
    mHistRegistryAN->addClone("Centrality_00-01/prof2hCorrelCombi0",
                              strCombi.data());
  }

 for (int iBin = 1; iBin < mNcentBins; iBin++) {
    mHistRegistryAN->addClone("Centrality_00-01/", mCentClasses[iBin].data());
  }

  LOGF(info, "AN histograms created.");
}

/// \brief Fill the profile with the pairs of harmonics for the provided array.
/// \param hBin Bin to fill in the histogram.
/// \param pairHarmo Decomposed pair of harmonics.
void FlowJAC2hHistManager::FillPairProf(const int hBin, const std::array<int, 2>& pairHarmo)
{
    /* Security checks. */
    if (!mHistRegistryAN) {
        LOGF(error, "AN histogram registry missing. Quitting...");
        return;
    }

    // Fill the profile itself.
    mHistRegistryAN->fill(HIST("histCombis"), 2.*(double)hBin+0.5, pairHarmo[0]);
    mHistRegistryAN->fill(HIST("histCombis"), 2.*(double)hBin+1.5, pairHarmo[1]);
}

} // namespace o2::analysis::PWGCF
