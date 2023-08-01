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

/* Header files. */
#include <array>
#include <string>
#include <vector>
#include <TH1.h>
#include <TProfile2D.h>

// O2 headers.
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"

// O2Physics headers.

// JCorran headers.
#include "PWGCF/JCorran/DataModel/JCatalyst.h"

/*
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
*/

/* Namespaces. */
using namespace o2;
using namespace o2::framework;

// ----------------------------------------------------------------------------
// Histogram manager to fill both QA and analysis histogram registries at the
// same time.
// ----------------------------------------------------------------------------
namespace o2::analysis::PWGCF
{
class JAC2hHistManager
{
public:
    JAC2hHistManager() = default;   // Constructor.

    /* Setters and getters. */
    void SetQAHistRegistry(HistogramRegistry *myRegistry)
    {
        mQAHistoRegistry = myRegistry;
        LOGF(info, "QA histogram registry successfully set.");
    }
    HistogramRegistry *GetQAHistRegistry() const {return mQAHistoRegistry;}

    void SetACHistRegistry(HistogramRegistry *myRegistry)
    {
        mACHistoRegistry = myRegistry;
        LOGF(info, "AC histogram registry successfully set.");
    }
    HistogramRegistry *GetACHistRegistry() const {return mACHistoRegistry;}

    void SetNcombis2h(int myCombis)
    {
        mNcombis2h = myCombis;
        LOGF(info, "Number of pairs of harmos: %d.", mNcombis2h);
    }
    int GetNcombis2h() const {return mNcombis2h;}

    void SetNsamples(int mySamples)
    {
        mNsamples = mySamples;
        LOGF(info, "Number of samples: %d.", mNsamples);

        // LOKI: Debug trick to stop the compilation to complain.
        const int nCentBins = sizeof(jflucCentBins)/sizeof(jflucCentBins[0]);
        printf("Number of centrality classes in JCatalyst: %d\n", nCentBins);
    }
    int GetNsamples() const {return mNsamples;}

    void SetEtaGap(int myGap)
    {
        mEtaGap = myGap;
        LOGF(info, "Eta gap: %.1f.", mEtaGap);
    }
    int GetEtaGap() const {return mEtaGap;}

    /* Methods related to this class. */
    // The template functions are defined here to prevent compilation errors.
    void CreateQAHistos();
    void CreateACHistos();

    /// \brief Fill the event QA histograms in the centrality class 'cBin'.
    /// \param coll Collision per defined in the JCatalyst data model.
    /// \param multi Number of tracks in the collision.
    /// \param sample ID of the bootstrap sample.
    template <int cBin, typename T>
    void FillEventQA(const T& coll, int multi, int sample)
    {
        if (!mQAHistoRegistry) {
            LOGF(error, "QA histogram registry missing. Quitting...");
            return;
        }

        mQAHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histCent"), coll.cent());
        mQAHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histZvtx"), coll.posZ());
        mQAHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histMulti"), multi);
        mQAHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histSamples"), sample);
    }

    /// \brief Fill the track QA histograms in the centrality class 'cBin'.
    /// \param track Track per defined in the JCatalyst data model.
    template <int cBin, typename T>
    void FillTrackQA(const T& track)
    {
        if (!mQAHistoRegistry) {
            LOGF(error, "QA histogram registry missing. Quitting...");
            return;
        }

        // NOTE: Crosscheck again that the corrected histograms are filled correctly!!!
        mQAHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histPtUncorrected"),
            track.pt());
        mQAHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histPtCorrected"),
            track.pt(), 1./(track.weightEff()));
        mQAHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histEta"),
            track.eta());
        mQAHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histPhiUncorrected"),
            track.phi());
        mQAHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histPhiCorrected"),
            track.phi(), 1./(track.weightNUA()));
        // mQAHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("histCharge"),
        //    track.charge());
    }

    /// \brief Fill 2-particle profile2D for 'cBin' class.
    /// \param
    template <int cBin, typename T>
    void Fill2pProf(const T& correl2p, const int sBin, const double weight2p)
    {
        if (!mACHistoRegistry) {
            LOGF(error, "AC histogram registry missing. Quitting...");
            return;
        }

        /* Fill the values for the full and sample for vn no gap. */
        // x-bin: no gap values on odd bins.
        // y-bin for full: 0, y-bin for sample: sampleID+1.5 (eg. 0 -> bin 1.5).
        for (int iB = 0; iB < 6; iB++) {
            mACHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("prof2pCorrel"),
                2.*iB+0.5, 0.5, correl2p[iB], weight2p);
            mACHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("prof2pCorrel"),
                2.*iB+0.5, sBin+1.5, correl2p[iB], weight2p);
        }
    }

    /// \brief Fill the 2-harmonic profile2D for 'cBin' class.
    /// \param 
    template <int cBin, int pBin, typename T>
    void Fill2hProf(const T& correl2h, const T& weight2h, const int sBin)
    {
        if (!mACHistoRegistry) {
            LOGF(error, "AC histogram registry missing. Quitting...");
            return;
        }

        /* Fill the values for the full and samples for all powers and each pair. */
        // The number of combinations needs to be hard-coded as we need to loop over it.
        // TODO: Try to make it compatible with the configurable...
        for (int iB = 0; iB < 14; iB++) {
            mACHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("prof2hCorrel")+HIST(mCombi[pBin]),
                iB+0.5, 0.5, correl2h.at(pBin).at(iB), weight2h.at(pBin).at(iB));
            mACHistoRegistry->fill(HIST(mCentClasses[cBin])+HIST("prof2hCorrel")+HIST(mCombi[pBin]),
                iB+0.5, sBin+1.5, correl2h.at(pBin).at(iB), weight2h.at(pBin).at(iB));
        }
    }

private:
    HistogramRegistry *mQAHistoRegistry = nullptr;  ///< For the QA output objects.
    HistogramRegistry *mACHistoRegistry = nullptr;  ///< For the analysis output objects.

    int mNcombis2h = 3;   ///< Number of pairs of harmonics.
    int mNsamples = 20;   ///< Number of samples for the bootstrap.
    float mEtaGap = 1.0;  ///< Value of the pseudorapidity gap.

    static constexpr std::string_view mCentClasses[] = {    ///< JCatalyst classes.
        "Centrality_00-01/", "Centrality_01-02/", "Centrality_02-05/",
        "Centrality_05-10/", "Centrality_10-20/", "Centrality_20-30/",
        "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/"
    };
    static constexpr std::string_view mPowers[] = {     ///< Labels for 2h profile.
        "{1,1}", "{2,2}", "{2,0}", "{2,1}", "{3,0}", "{3,1}", "{4,0}", "{4,1}",
                        "{0,2}", "{1,2}", "{0,3}", "{1,3}", "{0,4}", "{1,4}"
    };
    static constexpr std::string_view mCombi[] = {     ///< 
        "Combi0", "Combi1", "Combi2"
    };

  ClassDefNV(JAC2hHistManager, 1);  
};
} // namespace o2::analysis::PWGCF

#endif // JAC2HHISTMANAGER_H