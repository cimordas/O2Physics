
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

#ifndef JAC2HHISTMANAGER_H
#define JAC2HHISTMANAGER_H

// Includes.
#include <string>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
//#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "PWGCF/JCorran/DataModel/JCatalyst.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis::PWGCF
{
class JAC2hHistManager
{
public:
  JAC2hHistManager() = default;
  //~JAC2hHistManager();

  void SetHistManager(HistogramRegistry *myRegistry) {
    mHistoRegistry = myRegistry;
    LOGF(info, "Histogram registry has been set.\n");
  }
  HistogramRegistry *GetHistManager() const {return mHistoRegistry;}

  void CreateHistosQA();
  //void CreateHistosAN();

  template <int cBin, typename T>
  void fillEventQA(const T& coll);

  //template <int cBin, typename T>
  //void fillTrackQA(const T& track);

private:
  HistogramRegistry *mHistoRegistry = nullptr;  ///< QA + analysis output.

  ClassDefNV(JAC2hHistManager, 1);  
};
} // namespace o2::analysis::PWGCF

#endif // JAC2HHISTMANAGER_H