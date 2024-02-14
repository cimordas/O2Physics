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

// \brief   Histogram manager for the AC-related analyses.
// \author  Cindy Mordasini (cindy.mordasini@cern.ch)

#ifndef PWGCF_FLOW_CORE_FLOWJAC2HHISTMANAGER_H
#define PWGCF_FLOW_CORE_FLOWJAC2HHISTMANAGER_H

/* Header files. */
#include <array>
#include <string>
#include <string_view>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
//#include <THnSparse.h> // TODO: See if more advantageous for more events.
#include <TProfile.h>
#include <TProfile2D.h>

// O2 headers. //
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"

// O2 Physics headers.
#include "PWGCF/Flow/Core/FlowJHistManager.h"

/* Namespaces. */
using namespace o2;
using namespace o2::framework;

// ----------------------------------------------------------------------------
// Histogram manager to fill both QA and AN registries at the same time.
// ----------------------------------------------------------------------------
namespace o2::analysis::PWGCF
{
class FlowJAC2hHistManager
{
public:
  FlowJAC2hHistManager() = default;

  // Setters and getters, in the same order as the data members.
  void SetHistRegistryAN(HistogramRegistry *myRegistry)
  {
    mHistRegistryAN = myRegistry;
    LOGF(info, "AN histogram registry successfully set.");
  }
  HistogramRegistry *GetHistRegistryAN() const {return mHistRegistryAN;}

  void SetEtaGap(bool useGap)
  {
    mUseEtaGap = useGap;
    LOGF(info, "Use eta gap: %d.", mUseEtaGap);
  }
  bool GetUseEtaGap() const {return mUseEtaGap;}

  void SetNcombis2h(int myCombis)
  {
    mNcombis2h = myCombis;
    LOGF(info, "Number of pairs of harmonics: %d.", mNcombis2h);
  }
  int GetNcombis2h() const {return mNcombis2h;}

  void SetNsamples(int mySamples)
  {
    mNsamples = mySamples;
    LOGF(info, "Number of samples: %d.", mNsamples);
  }
  int GetNsamples() const {return mNsamples;}

  /* Methods specific to this class. */
  // The template functions are defined down here to prevent compilation errors.
  void CreateHistAN();
  void FillPairProf(const int hBin, const std::array<int, 2>& pairHarmo);

  /// \brief Fill the number of events per bootstrap sample and the full dataset.
  /// \tparam cBin Centrality class bin of the collision.
  /// \param sBin SampleID attributed to the collision.
  template <int cBin>
  void FillHistSamples(int sBin)
  {
    if (!mHistRegistryAN) {
      LOGF(fatal, "AN histogram registry missing. Quitting...");
      return;
    }

    /* Add the collision to the full and sample bins. */
    // 1st bin: full dataset.
    // sBin+1: bin associated with sBin (so sBin = 0 -> bin = 1).
    mHistRegistryAN->fill(HIST(mCentClasses[cBin])+HIST("histSamples"), 0);
    mHistRegistryAN->fill(HIST(mCentClasses[cBin])+HIST("histSamples"), sBin+1);
    LOGF(info, "Sample histogram successfully filled.");
  }

  /// \brief Fill the profile2D for the 2-particle correlators.
  /// \tparam cBin Centrality class bin of the collision.
  /// \tparam T Template standing for std::array<double>.
  /// \tparam U Template standing for std::array<std::array<double>>.
  /// \param correl2p srd::array<double> for v_1 to v_6 no eta gap.
  /// \param weight2p Associated event weight.
  /// \param sBin SampleID attributed to the collision.
  /// \param correl2pEtaGap std::array<std::array<double>> for correlator and weight with eta gap.
  template <int cBin, typename T, typename U>
  void FillProf2p(const T& correl2p, const double weight2p, const int sBin,
                  const U& correl2pEtaGap)
  {
    if (!mHistRegistryAN) {
      LOGF(fatal, "AN histogram registry missing. Quitting...");
      return;
    }

    /* Fill the values for v_n without eta gap. */
    // x-bin: values without gap on odd bins, with gap on even bins.
    // y-bin: full data on bin 0.5, sample on bin sBin+1.5.
    for (int iB = 0; iB < 6; iB++) {
      mHistRegistryAN->fill(HIST(mCentClasses[cBin])+HIST("prof2pCorrel"),
        2.*iB+0.5, 0.5, correl2p[iB], weight2p);
      mHistRegistryAN->fill(HIST(mCentClasses[cBin])+HIST("prof2pCorrel"),
        2.*iB+0.5, sBin+1.5, correl2p[iB], weight2p);

      if (mUseEtaGap) {
        mHistRegistryAN->fill(HIST(mCentClasses[cBin])+HIST("prof2pCorrel"),
          2.*iB+1.5, 0.5, correl2pEtaGap[iB][0], correl2pEtaGap[iB][1]);
        mHistRegistryAN->fill(HIST(mCentClasses[cBin])+HIST("prof2pCorrel"),
          2.*iB+1.5, sBin+1.5, correl2pEtaGap[iB][0], correl2pEtaGap[iB][1]);
      }
    }
     LOGF(info, "2p profile successfully filled.");
  }

  /// \brief Fill the profile2D for the 2-harmonic correlators.
  /// \tparam T Template standing for std::array<double>.
  /// \tparam cBin Centrality class bin of the collision.
  /// \tparam hBin Index to identify which combination it is.
  /// \param correl2h srd::array<double> for the different powers.
  /// \param weight2h Associated event weights.
  /// \param sBin SampleID attributed to the collision.
  template <int cBin, int hBin, typename T>
  void FillProf2h(const T& correl2h, const T& weight2h, const int sBin)
  {
    if (!mHistRegistryAN) {
      LOGF(fatal, "AN histogram registry missing. Quitting...");
      return;
    }

    /* Fill the values for the full and samples for all powers and each pair. */
    // The number of combinations needs to be hard-coded as we need to loop over it.
    // TODO: Try to make it compatible with the configurable...
    for (int iB = 0; iB < 14; iB++) {
      mHistRegistryAN->fill(HIST(mCentClasses[cBin])+HIST("prof2hCorrel")+HIST(mCombi[hBin]),
        iB+0.5, 0.5, correl2h.at(hBin).at(iB), weight2h.at(hBin).at(iB));
      mHistRegistryAN->fill(HIST(mCentClasses[cBin])+HIST("prof2hCorrel")+HIST(mCombi[hBin]),
        iB+0.5, sBin+1.5, correl2h.at(hBin).at(iB), weight2h.at(hBin).at(iB));
    }
     LOGF(info, "2h profile successfully filled.");
  }

private:
  HistogramRegistry *mHistRegistryQA = nullptr;   ///< For the QA output.
  HistogramRegistry *mHistRegistryAN = nullptr;   ///< For the analysis output.

  bool mUseEtaGap = true;   ///< Enable the eta gap of the 2-particle correlators.
  int mNcombis2h = 3;   ///< Number of pairs of harmonics.
  int mNsamples = 20;   ///< Number of samples for the bootstrap.

  static const int mNcentBins = 10;   ///< Number of centrality classes.
  static constexpr std::string_view mCentClasses[] = {  ///< Centrality classes.
    "Centrality_00-01/", "Centrality_01-02/", "Centrality_02-05/",
    "Centrality_05-10/", "Centrality_10-20/", "Centrality_20-30/",
    "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/",
    "Centrality_60-70/"
  };

  static constexpr std::string_view mPowers[] = {   ///< Labels for 2h profile.
    "{1,1}", "{2,2}", "{2,0}", "{2,1}", "{3,0}", "{3,1}", "{4,0}", "{4,1}",
                      "{0,2}", "{1,2}", "{0,3}", "{1,3}", "{0,4}", "{1,4}"
  };

  static constexpr std::string_view mCombi[] = {     ///< 
      "Combi0", "Combi1", "Combi2"
  };

  ClassDefNV(FlowJAC2hHistManager, 1);
};
} // namespace o2::analysis::PWGCF

#endif  // PWGCF_FLOW_CORE_FLOWJAC2HHISTMANAGER_H
