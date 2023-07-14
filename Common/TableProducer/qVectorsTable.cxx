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
/// \file   qVectorsTable.cxx
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Task calculating the Q-vectors for each collision in a bunch crossing
///         (with or without corrections) and save the results in a dedicated table.
///

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

#include "Common/DataModel/Qvectors.h"

// o2 includes.
#include <CCDB/BasicCCDBManager.h>
#include "DetectorsCommonDataFormats/AlignParam.h"

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;
//using namespace o2::expressions;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected,
  aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;

struct qVectorsTable {
  // Configurables.
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL",
      "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<long> nolaterthan{"ccdb-no-later-than",
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),
      "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  Configurable<int> cfgCentEsti{"cfgCentEsti",
      2, "Centrality estimator (Run3): 0 = FT0M, 1 = FT0A, 2 = FT0C, 3 = FV0A"};
    // LOKI: We have here all centrality estimators for Run 3 (except FDDM and NTPV),
    // but the Q-vectors are calculated only for some of them.

  struct : ConfigurableGroup {
    Configurable<std::vector<float>> cfgCorrConstFT0A{"cfgCorrConstFT0A",
        {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0A"};
    Configurable<std::vector<float>> cfgCorrConstFT0C{"cfgCorrConstFT0C",
        {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0C"};
    Configurable<std::vector<float>> cfgCorrConstFV0{"cfgCorrConstFV0",
        {0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FV0"};
  } cfgCorrConstAll;

  // Table.
  Produces<aod::Qvectors> qVector;

  // Enable access to the CCDB for the offset and correction constants and save them
  // in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::vector<o2::detectors::AlignParam> *offsetFT0;
  std::vector<o2::detectors::AlignParam> *offsetFV0;

  // Variables for other classes.
  EventPlaneHelper helperEP;


  void init(InitContext const&)
  {
    // Setup the access to the CCDB objects of interest.
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(cfgCcdbParam.nolaterthan.value);

    LOGF(info, "Getting alignment offsets from the CCDB...");
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align",
      cfgCcdbParam.nolaterthan.value);
    offsetFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align",
      cfgCcdbParam.nolaterthan.value);

    // Get the offset values for the different parts of FIT.
    if (offsetFT0 != nullptr) {
      // FT0 has vector size 2: one element for A side, one for C side.
      helperEP.SetOffsetFT0A((*offsetFT0)[0].getX(), (*offsetFT0)[0].getY());
      helperEP.SetOffsetFT0C((*offsetFT0)[1].getX(), (*offsetFT0)[1].getY());
    }
    else {LOGF(fatal, "Could not get the alignment parameters for FT0.");}

    if (offsetFV0 != nullptr) {
      // FV0 has vector size 2: one element for left side, one for right side.
      helperEP.SetOffsetFV0left((*offsetFV0)[0].getX(), (*offsetFV0)[0].getY());
      helperEP.SetOffsetFV0right((*offsetFV0)[1].getX(), (*offsetFV0)[1].getY());
    }
    else {LOGF(fatal, "Could not get the alignment parameters for FV0.");}
    /*  // Debug printing.
      printf("Offset for FT0A: x = %.3f y = %.3f\n", (*offsetFT0)[0].getX(), (*offsetFT0)[0].getY());
      printf("Offset for FT0C: x = %.3f y = %.3f\n", (*offsetFT0)[1].getX(), (*offsetFT0)[1].getY());
      printf("Offset for FV0-left: x = %.3f y = %.3f\n", (*offsetFV0)[0].getX(), (*offsetFV0)[0].getY());
      printf("Offset for FV0-right: x = %.3f y = %.3f\n", (*offsetFV0)[1].getX(), (*offsetFV0)[1].getY());
    */

   // LOKI: If we need to access the corrections from the CCDB, insert that here.
  }

  void process(MyCollisions::iterator const& coll, aod::FT0s const& ft0s, aod::FV0As const& fv0s)//, aod::FV0Cs const&)
  {
    // Get the centrality value for all subscribed estimators and takes the one
    // corresponding to cfgCentEsti. Reject also the events with invalid values.
    // NOTE: centFDDM and centNTPV not implemented as it makes the compilation crashes...
    float centAllEstim[4] = {
      coll.centFT0M(), coll.centFT0A(), coll.centFT0C(),
      coll.centFV0A()
    };
    float cent = centAllEstim[cfgCentEsti];
    LOG(info) << "COLLISION INDEX: " << coll.globalIndex()
      << " Centrality percentile: " << cent;
    if (cent < 0. || cent > 100.) {
      LOGF(info, "Invalid centrality value. Skipping this event.");
      return;
    }

    // Calculate the Q-vectors values for this event.
    // TODO: Add here qVect for other detectors,...
    float qVectFT0A[2] = {0.};    // Real and imaginary parts of the Q-vector in FT0A.
    float qVectFT0C[2] = {0.};    // Real and imaginary parts of the Q-vector in FT0C.
    float qVectFV0[2] = {0.};     // Real and imaginary parts of the Q-vector in FV0A.

    TComplex QvecDet(0);      // Complex value of the Q-vector for any detector.
    double sumAmplDet = 0.;   // Sum of the amplitudes of all non-dead channels in any detector.

    /// First check if the collision has a found FT0. If yes, calculate the
    /// Q-vectors for FT0A and FT0C (both real and imaginary parts). If no,
    /// attribute dummy values to the corresponding qVect.
    if (coll.has_foundFT0()) {
      LOGF(info, "A FT0 has been found. Getting Q-vectors for it...");
      auto ft0 = coll.foundFT0();

      // Iterate over the non-dead channels for FT0-A to get the total Q-vector
      // and sum of amplitudes.
      for (std::size_t iChA = 0; iChA < ft0.channelA().size(); iChA++) {
        // Get first the corresponding amplitude.
        float ampl = ft0.amplitudeA()[iChA];

        // Update the Q-vector and sum of amplitudes using the helper function.
        // LOKI: Note this assumes nHarmo = 2!! Likely generalise in the future.
        helperEP.SumQvectors(0, iChA, ampl, QvecDet, sumAmplDet);
      } // Go to the next channel iChA.

      // Set the Qvectors for FT0A with the normalised Q-vector values if the sum of
      // amplitudes is non-zero. Otherwise, set it to a dummy 999.
      if (sumAmplDet != 0) {
        QvecDet /= sumAmplDet;
        qVectFT0A[0] = QvecDet.Re();
        qVectFT0A[1] = QvecDet.Im();
        //printf("qVectFT0A[0] = %.2f ; qVectFT0A[1] = %.2f \n", qVectFT0A[0], qVectFT0A[1]); // Debug printing.
      }
      else {
        qVectFT0A[0] = 999.; qVectFT0A[1] = 999.;
      }

      // Repeat the procedure with FT0-C for the found FT0.
      // Start by resetting to zero the intermediate quantities.
      QvecDet = TComplex(0.,0.); sumAmplDet = 0;
      for (std::size_t iChC = 0; iChC < ft0.channelC().size(); iChC++) {
        // iChC ranging from 0 to max 112. We need to add 96 (= max channels in FT0-A)
        // to ensure a proper channel number in FT0 as a whole.
        float ampl = ft0.amplitudeC()[iChC];
        helperEP.SumQvectors(0, iChC+96, ampl, QvecDet, sumAmplDet);
      }

      if (sumAmplDet != 0) {
        QvecDet /= sumAmplDet;
        qVectFT0C[0] = QvecDet.Re();
        qVectFT0C[1] = QvecDet.Im();
        //printf("qVectFT0C[0] = %.2f ; qVectFT0C[1] = %.2f \n", qVectFT0C[0], qVectFT0C[1]); // Debug printing.
      }
      else {
        qVectFT0C[0] = 999.; qVectFT0C[1] = 999.;
      }
    }
    else {
      LOGF(info, "No FT0 has been found. Setting Q-vectors to -999.");
      qVectFT0A[0] = -999.; qVectFT0A[1] = -999.;
      qVectFT0C[0] = -999.; qVectFT0C[1] = -999.;
    }

    /// Repeat the procedure for FV0 if one has been found for this collision.
    /// Again reset the intermediate quantities to zero.
    QvecDet = TComplex(0.,0.); sumAmplDet = 0;
    if (coll.has_foundFV0()) {
      LOGF(info, "A FV0 has been found. Calculating Q-vectors for FV0A...");
      auto fv0 = coll.foundFV0();

      for (std::size_t iCh = 0; iCh < fv0.channel().size(); iCh++) {
        float ampl = fv0.amplitude()[iCh];
        helperEP.SumQvectors(1, iCh, ampl, QvecDet, sumAmplDet);
      }

      if (sumAmplDet != 0) {
        QvecDet /= sumAmplDet;
        qVectFV0[0] = QvecDet.Re();
        qVectFV0[1] = QvecDet.Im();
        //printf("qVectFV0[0] = %.2f ; qVectFV0[1] = %.2f \n", qVectFV0[0], qVectFV0[1]); // Debug printing.
      }
      else {
        qVectFV0[0] = 999.; qVectFV0[1] = 999.;
      }
    }
    else {
      LOGF(info, "No FV0 has been found. Setting Q-vectors to -999.");
      qVectFV0[0] = -999.; qVectFV0[1] = -999.;
    }

    /// TODO: Repeat here the procedure for any other Qvector columns.
    /// Do not forget to add the configurable for the correction constants.

    // Apply the correction constants (configurable) to the obtained Q-vectors.
    // The function needs to be called for each detector/set separately.
    // A correction constant set to zero means this correction is not applied.
    // LOKI: Each detector must have their own vector of correction constants.
    helperEP.DoCorrections(qVectFT0A[0], qVectFT0A[1], cfgCorrConstAll.cfgCorrConstFT0A);
    helperEP.DoCorrections(qVectFT0C[0], qVectFT0C[1], cfgCorrConstAll.cfgCorrConstFT0C);
    helperEP.DoCorrections(qVectFV0[0], qVectFV0[1], cfgCorrConstAll.cfgCorrConstFV0);

    // Fill the columns of the Qvectors table.
    ///qVector(coll.globalIndex(), qVectFT0A[0], qVectFT0A[1],
    qVector(cent,
            qVectFT0A[0], qVectFT0A[1],
            qVectFT0C[0], qVectFT0C[1],
            qVectFV0[0], qVectFV0[1]);
  } // End process.

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsTable>(cfgc)
  };
}