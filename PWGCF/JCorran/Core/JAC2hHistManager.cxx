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

/* Header files. */
#include <string>
#include <vector>

// O2 headers.

// O2Physics headers.

// JCorran headers.
#include "PWGCF/JCorran/Core/JAC2hHistManager.h"

/* Namespaces. */
using namespace o2;
using namespace o2::framework;
namespace o2::analysis::PWGCF
{
/// \brief Create the histograms in the QA registry.
/// \note QA and AC are created separately to allow to get only the QA if needed.
void JAC2hHistManager::CreateQAHistos()
{
    /* Security checks. */
    if (!mQAHistoRegistry) {
        LOGF(error, "QA histogram registry missing. Quitting...");
        return;
    }

    /* Define the histograms for the QA and for the AC analysis. */
    // All histograms are defined for the first centrality class in details,
    // then cloned for the other centrality classes.
    const AxisSpec axisCent{100, 0., 100., "Centrality percentile"};
    mQAHistoRegistry->add("histCentCatalyst", "Centrality from the JCatalyst",
        HistType::kTH1F, {axisCent}, true);
    mQAHistoRegistry->add("Centrality_00-01/histCent", "Centrality after AC cuts",
        HistType::kTH1F, {axisCent}, true);

    const AxisSpec axisZvtx{75, -15., 15., "Z_{vtx} [cm]"};
    mQAHistoRegistry->add("Centrality_00-01/histZvtx", "Z_{vtx} after AC cuts",
        HistType::kTH1F, {axisZvtx}, true);

    const AxisSpec axisMulti{1000, 0., 5000., "N_{tracks}"};
    mQAHistoRegistry->add("Centrality_00-01/histMulti", "Multiplicity after AC cuts",
        HistType::kTH1I, {axisMulti}, true);

    const AxisSpec axisSamples{mNsamples, 0., (double)mNsamples, "Sample ID"};
    mQAHistoRegistry->add("Centrality_00-01/histSamples", "N_{events} in each bootstrap sample",
        HistType::kTH1I, {axisSamples}, true);
    for (int i = 1; i <= mNsamples; i++) {
        mQAHistoRegistry->get<TH1>(HIST("Centrality_00-01/histSamples"))
            ->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
    }

/*  If a variable binning is needed.
    std::vector<double> ptBinning = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                                    0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                    0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, 1.4,
                                    1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6,
                                    2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.5, 5., 6.};
    const AxisSpec axisPt = {ptBinning, "#it{p}_{T} [GeV/#it{c}]"};
*/
    const AxisSpec axisPt = {60, 0., 6., "#it{p}_{T} [GeV/#it{c}]"};
    mQAHistoRegistry->add("Centrality_00-01/histPtUncorrected", "#it{p}_{T} (not NUE corrected)",
        HistType::kTH1F, {axisPt}, true);
    mQAHistoRegistry->add("Centrality_00-01/histPtCorrected", "#it{p}_{T} (NUE corrected)",
        HistType::kTH1F, {axisPt}, true);

    const AxisSpec axisEta = {20, -1., 1., "#eta"};
    mQAHistoRegistry->add("Centrality_00-01/histEta", "Pseudorapidity",
        HistType::kTH1F, {axisEta}, true);

    const AxisSpec axisPhi = {100, 0., 2.*M_PI, "#varphi"};
    mQAHistoRegistry->add("Centrality_00-01/histPhiUncorrected", "Azimuthal angles (not NUA corrected)",
        HistType::kTH1F, {axisPhi}, true);
    mQAHistoRegistry->add("Centrality_00-01/histPhiCorrected", "Azimuthal angles (NUA corrected)",
        HistType::kTH1F, {axisPhi}, true);

    const AxisSpec axisCharge = {2, -1., 1., "Charge"};
    mQAHistoRegistry->add("Centrality_00-01/histCharge", "Electric charge",
        HistType::kTH1I, {axisCharge}, true);

    for (int iBin = 1; iBin < 9; iBin++) {
        mQAHistoRegistry->addClone("Centrality_00-01/", mCentClasses[iBin].data());
    }

    // LOKI: Debug trick to stop the compilation to complain.
    const int nCentBins = sizeof(jflucCentBins)/sizeof(jflucCentBins[0]);
    printf("Number of centrality classes in JCatalyst: %d\n", nCentBins);
}

/// \brief Create the histograms/profiles in the AC registry.
void JAC2hHistManager::CreateACHistos()
{
    /* Security checks. */
    if (!mACHistoRegistry) {
        LOGF(error, "AC histogram registry missing. Quitting...");
        return;
    }

    /* Define the histogram bookkeeping the pairs of harmonics as integers. */
    const AxisSpec axisCombis{2*mNcombis2h, 0., 2.*(double)mNcombis2h, "Pairs (m,n)"};
    mACHistoRegistry->add("histCombis", "Configured pairs of harmonics",
        HistType::kTProfile, {axisCombis}, true);
    for (int i = 0; i < mNcombis2h; i++) {
        mACHistoRegistry->get<TProfile>(HIST("histCombis"))
            ->GetXaxis()->SetBinLabel(2*i+1, Form("Combi%d-m", i));
        mACHistoRegistry->get<TProfile>(HIST("histCombis"))
            ->GetXaxis()->SetBinLabel(2*(i+1), Form("Combi%d-n", i));
    }

    /* Define the TProfile2D used to save the multiparticle correlators. */
    // The x-axis keeps the correlators, the y-axis the full data and samples.
    // All profiles are defined for the first centrality class in details,
    // then cloned for the other centrality classes.
    const AxisSpec axis2pCorrelX{12, 0., 12., "n"};
    const AxisSpec axis2hCorrelX{14, 0., 14., "Powers"};
    const AxisSpec axis2pCorrelY{mNsamples+1, 0., (double)mNsamples+1., "Sample ID"};
    mACHistoRegistry->add("Centrality_00-01/prof2pCorrel", "<v_{n}^{2}>",
        HistType::kTProfile2D, {axis2pCorrelX, axis2pCorrelY}, true);
    mACHistoRegistry->add("Centrality_00-01/prof2hCorrelCombi0",
        "<v_{m}^{2a}v_{n}^{2b}>",
        HistType::kTProfile2D, {axis2hCorrelX, axis2pCorrelY}, true);

    // The x-axis of 'prof2pCorrel' follows <n, n with gap> for n = 1...6.
    for (int i = 1; i <= 6; i++) {
        mACHistoRegistry->get<TProfile2D>(HIST("Centrality_00-01/prof2pCorrel"))
            ->GetXaxis()->SetBinLabel(2*i-1, Form("v_{%d}", i));
        mACHistoRegistry->get<TProfile2D>(HIST("Centrality_00-01/prof2pCorrel"))
            ->GetXaxis()->SetBinLabel(2*i, Form("v_{%d, #Delta#eta}", i));
    }
    for (int i = 1; i <= 14; i++) {
        mACHistoRegistry->get<TProfile2D>(HIST("Centrality_00-01/prof2hCorrelCombi0"))
            ->GetXaxis()->SetBinLabel(i, Form("%s", mPowers[i-1].data()));
    }

    // The first y-bin is always the full dataset.
    mACHistoRegistry->get<TProfile2D>(HIST("Centrality_00-01/prof2pCorrel"))
        ->GetYaxis()->SetBinLabel(1, "Full");
    mACHistoRegistry->get<TProfile2D>(HIST("Centrality_00-01/prof2hCorrelCombi0"))
        ->GetYaxis()->SetBinLabel(1, "Full");
    for (int i = 2; i <= mNsamples+1; i++) {
        mACHistoRegistry->get<TProfile2D>(HIST("Centrality_00-01/prof2pCorrel"))
            ->GetYaxis()->SetBinLabel(i, Form("%d", i-2));
        mACHistoRegistry->get<TProfile2D>(HIST("Centrality_00-01/prof2hCorrelCombi0"))
            ->GetYaxis()->SetBinLabel(i, Form("%d", i-2));
    }

    for (int iH = 1; iH < mNcombis2h; iH++) {
        std::string strCombi = Form("Centrality_00-01/prof2hCorrelCombi%d", iH);
        mACHistoRegistry->addClone("Centrality_00-01/prof2hCorrelCombi0",
            strCombi.data());
    }

    for (int iBin = 1; iBin < 9; iBin++) {
        mACHistoRegistry->addClone("Centrality_00-01/", mCentClasses[iBin].data());
    }
}

/// \brief Fill the profile with the pairs of harmonics for the provided array.
/// \param .
/// \param .
void JAC2hHistManager::FillPairProf(const int hBin, const std::array<int, 2>& pairHarmo)
{
    /* Security checks. */
    if (!mACHistoRegistry) {
        LOGF(error, "AC histogram registry missing. Quitting...");
        return;
    }

    // Fill the profile itself.
    mACHistoRegistry->fill(HIST("histCombis"), 2.*(double)hBin+0.5, pairHarmo[0]);
    mACHistoRegistry->fill(HIST("histCombis"), 2.*(double)hBin+1.5, pairHarmo[1]);
}

} // namespace o2::analysis::PWGCF