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

/* Header files. */

// O2 headers.

// O2 Physics headers.
#include "PWGCF/Flow/Core/FlowAC2hHistManager.h"

/* Namespaces. */
using namespace o2;
using namespace o2::framework;

namespace o2::analysis::PWGCF
{
/// \brief Create the histograms in the QA registry.
/// \note QA and AN are created separately to allow getting the QA only.
void FlowAC2hHistManager::CreateHistQA()
{
  /* Security checks. */
  if (!mHistRegistryQA) {
    LOGF(error, "QA histogram registry missing. Quitting...");
    return;
  }

  /* Definition of the QA histograms. */
  // All the histograms are defined in details for the first centrality
  // class before additional cuts, then cloned for the other classes.

  const AxisSpec axisCent{100, 0., 100., "Centrality percentile"};
  mHistRegistryQA->add("Centrality_00-01/histCent", "Centrality",
                        HistType::kTH1F, {axisCent}, true);

  const AxisSpec axisMulti{1000, 0., 25000., "N_{tracks}"};
  mHistRegistryQA->add("Centrality_00-01/histMulti", "Multiplicity",
                        HistType::kTH1I, {axisMulti}, true);

  const AxisSpec axisZvtx{24, -12., 12., "Z_{vtx} [cm]"};
  mHistRegistryQA->add("Centrality_00-01/histZvtx", "Z_{vtx}",
                        HistType::kTH1F, {axisZvtx}, true);

/*If a variable binning is needed.
  std::vector<double> ptBinning = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                                  0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                  0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4,
                                  1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6,
                                  2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.5, 5., 6.};
  const AxisSpec axisPt = {ptBinning, "#it{p}_{T} [GeV/#it{c}]"};
*/
  const AxisSpec axisPt = {60, 0., 6., "#it{p}_{T} [GeV/#it{c}]"};
  mHistRegistryQA->add("Centrality_00-01/histPt", "#it{p}_{T} (not NUE corrected)",
                        HistType::kTH1F, {axisPt}, true);
  mHistRegistryQA->add("Centrality_00-01/histPtCorrected", "#it{p}_{T} (NUE corrected)",
                        HistType::kTH1F, {axisPt}, true);

  const AxisSpec axisEta = {20, -1., 1., "#eta"};
  mHistRegistryQA->add("Centrality_00-01/histEta", "Pseudorapidity",
                        HistType::kTH1F, {axisEta}, true);

  const AxisSpec axisPhi = {100, 0., 2.*M_PI, "#varphi"};
  mHistRegistryQA->add("Centrality_00-01/histPhi", "Azimuthal angles (not NUA corrected)",
                        HistType::kTH1F, {axisPhi}, true);
  mHistRegistryQA->add("Centrality_00-01/histPhiCorrected", "Azimuthal angles (NUA corrected)",
                        HistType::kTH1F, {axisPhi}, true);

  const AxisSpec axisCharge = {2, -2., 2., "Charge"};
  mHistRegistryQA->add("Centrality_00-01/histCharge", "Electric charge",
                        HistType::kTH1I, {axisCharge}, true);

  // Clone the first centrality class.
  for (int iBin = 1; iBin < mNcentBins; iBin++) {
    mHistRegistryQA->addClone("Centrality_00-01/", mCentClasses[iBin].data());
  }

  LOGF(info, "QA histograms created.");
}

/// \brief Create the histograms in the analysis registry.
void FlowAC2hHistManager::CreateHistAN()
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

/// \brief Get the centrality bin value corresponding to the percentile.
/// \param Centrality percentile of the collision.
/// \return Bin for the histograms,...
int FlowAC2hHistManager::GetCentBin(float cValue)
{
  const float centClasses[] = {0., 1., 2., 5., 10., 20., 30., 40., 50., 60., 70.};

  for (int i = 0; i < mNcentBins+1; i++) {
    if (cValue >= centClasses[i]) {continue;}
    else {return i-1;}
  }

  // We went through all centrality edges without returning at all.
  // --> The measured percentile is larger than the final class we consider.
  return -1;
}

/// \brief Fill the profile with the pairs of harmonics for the provided array.
/// \param hBin Bin to fill in the histogram.
/// \param pairHarmo Decomposed pair of harmonics.
void FlowAC2hHistManager::FillPairProf(const int hBin, const std::array<int, 2>& pairHarmo)
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
