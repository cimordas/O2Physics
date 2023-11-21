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
void FlowAC2hHistManager::CreateHistsQA()
{
  /* Security checks. */
  if (!mQAHistRegistry) {
    LOGF(error, "QA histogram registry missing. Quitting...");
    return;
  }

  /* Definition of the QA histograms. */
  // All the histograms are defined in details for the first centrality
  // class before additional cuts, then cloned for the other classes.

  const AxisSpec axisCent{100, 0., 100., "Centrality percentile"};
  mQAHistRegistry->add("Centrality_00-01/Before/histCent", "Centrality",
                        HistType::kTH1F, {axisCent}, true);

  const AxisSpec axisMulti{1000, 0., 25000., "N_{tracks}"};
  mQAHistRegistry->add("Centrality_00-01/Before/histMulti", "Multiplicity",
                        HistType::kTH1I, {axisMulti}, true);

  const AxisSpec axisZvtx{24, -12., 12., "Z_{vtx} [cm]"};
  mQAHistRegistry->add("Centrality_00-01/Before/histZvtx", "Z_{vtx}",
                        HistType::kTH1F, {axisZvtx}, true);

  mQAHistRegistry->addClone("Centrality_00-01/Before", "Centrality_00-01/After");

  // Clone the first centrality class.
  for (int iBin = 1; iBin < mNcentBins; iBin++) {
    mQAHistRegistry->addClone("Centrality_00-01/", mCentClasses[iBin].data());
  }

  LOGF(info, "QA histograms created.");
}

/// \brief Get the centrality bin value corresponding to the percentile.
/// \param Centrality percentile of the collision.
/// \return Bin for the histograms,...
int FlowAC2hHistManager::GetCentBin(float cValue)
{
  const float centClasses[] = {0., 1., 2., 5., 10., 20., 30., 40., 50., 60., 70.};

  for (int i = 0; i < mNcentBins+1; i++) {
    if (cValue >= centClasses[i]) {continue;}
    else {return i - 1;}
  }

  // We went through all centrality edges without returning at all.
  // --> The measured percentile is larger than the final class we consider.
  return -1;
}

} // namespace o2::analysis::PWGCF
