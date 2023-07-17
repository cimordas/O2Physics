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

#include <vector>
#include <string>
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

struct : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<std::vector<float>> cfgFT0ACentBin0{"cfgFT0ACentBin0", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0A, cent bin 0"};
  o2::framework::Configurable<std::vector<float>> cfgFT0ACentBin1{"cfgFT0ACentBin1", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0A, cent bin 1"};
  o2::framework::Configurable<std::vector<float>> cfgFT0ACentBin2{"cfgFT0ACentBin2", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0A, cent bin 2"};
  o2::framework::Configurable<std::vector<float>> cfgFT0ACentBin3{"cfgFT0ACentBin3", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0A, cent bin 3"};
  o2::framework::Configurable<std::vector<float>> cfgFT0ACentBin4{"cfgFT0ACentBin4", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0A, cent bin 4"};
  o2::framework::Configurable<std::vector<float>> cfgFT0ACentBin5{"cfgFT0ACentBin5", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0A, cent bin 5"};
  o2::framework::Configurable<std::vector<float>> cfgFT0ACentBin6{"cfgFT0ACentBin6", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0A, cent bin 6"};
  o2::framework::Configurable<std::vector<float>> cfgFT0ACentBin7{"cfgFT0ACentBin7", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0A, cent bin 7"};
} cfgCorrConstFT0A;

struct : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<std::vector<float>> cfgFT0CCentBin0{"cfgFT0CCentBin0", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0C, cent bin 0"};
  o2::framework::Configurable<std::vector<float>> cfgFT0CCentBin1{"cfgFT0CCentBin1", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0C, cent bin 1"};
  o2::framework::Configurable<std::vector<float>> cfgFT0CCentBin2{"cfgFT0CCentBin2", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0C, cent bin 2"};
  o2::framework::Configurable<std::vector<float>> cfgFT0CCentBin3{"cfgFT0CCentBin3", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0C, cent bin 3"};
  o2::framework::Configurable<std::vector<float>> cfgFT0CCentBin4{"cfgFT0CCentBin4", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0C, cent bin 4"};
  o2::framework::Configurable<std::vector<float>> cfgFT0CCentBin5{"cfgFT0CCentBin5", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0C, cent bin 5"};
  o2::framework::Configurable<std::vector<float>> cfgFT0CCentBin6{"cfgFT0CCentBin6", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0C, cent bin 6"};
  o2::framework::Configurable<std::vector<float>> cfgFT0CCentBin7{"cfgFT0CCentBin7", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0C, cent bin 7"};
} cfgCorrConstFT0C;

struct : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<std::vector<float>> cfgFV0ACentBin0{"cfgFV0ACentBin0", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FV0A, cent bin 0"};
  o2::framework::Configurable<std::vector<float>> cfgFV0ACentBin1{"cfgFV0ACentBin1", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FV0A, cent bin 1"};
  o2::framework::Configurable<std::vector<float>> cfgFV0ACentBin2{"cfgFV0ACentBin2", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FV0A, cent bin 2"};
  o2::framework::Configurable<std::vector<float>> cfgFV0ACentBin3{"cfgFV0ACentBin3", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FV0A, cent bin 3"};
  o2::framework::Configurable<std::vector<float>> cfgFV0ACentBin4{"cfgFV0ACentBin4", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FV0A, cent bin 4"};
  o2::framework::Configurable<std::vector<float>> cfgFV0ACentBin5{"cfgFV0ACentBin5", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FV0A, cent bin 5"};
  o2::framework::Configurable<std::vector<float>> cfgFV0ACentBin6{"cfgFV0ACentBin6", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FV0A, cent bin 6"};
  o2::framework::Configurable<std::vector<float>> cfgFV0ACentBin7{"cfgFV0ACentBin7", {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FV0A, cent bin 7"};
} cfgCorrConstFV0A;