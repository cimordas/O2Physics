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
/// \file   Qvectors.h
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Declaration of the table for the (un)corrected Q-vectors for the event plane
/// determination.
///

#ifndef O2_ANALYSIS_QVECTORS_H_
#define O2_ANALYSIS_QVECTORS_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace qvec
{
DECLARE_SOA_COLUMN(Cent, cent, float);                //! Centrality percentile.
DECLARE_SOA_COLUMN(QvecFT0ARe, qvecFT0ARe, float);    //! Real part of Qvec in FT0-A.
DECLARE_SOA_COLUMN(QvecFT0AIm, qvecFT0AIm, float);    //! Imaginary part for FT0-A.
DECLARE_SOA_COLUMN(QvecFT0CRe, qvecFT0CRe, float);    //! Real part of Qvec in FT0-c.
DECLARE_SOA_COLUMN(QvecFT0CIm, qvecFT0CIm, float);    //! Imaginary part for FT0-c.
DECLARE_SOA_COLUMN(QvecFV0Re, qvecFV0Re, float);      //! Real part of Qvec in FV0.
DECLARE_SOA_COLUMN(QvecFV0Im, qvecFV0Im, float);      //! Imaginary part for FV0.
/// NOTE: Add here Qx,Qy for other systems.
} // namespace qvec

DECLARE_SOA_TABLE(Qvectors, "AOD", "QVECTORS",          //! Table with all Qvectors.
                  qvec::Cent,
                  qvec::QvecFT0ARe, qvec::QvecFT0AIm,
                  qvec::QvecFT0CRe, qvec::QvecFT0CIm,
                  qvec::QvecFV0Re, qvec::QvecFV0Im);

using Qvector = Qvectors::iterator;
} // namespace o2::aod

#endif // O2_ANALYSIS_QVECTORS_H_
