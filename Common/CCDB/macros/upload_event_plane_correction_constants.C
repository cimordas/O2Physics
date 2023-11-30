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

#include <map>
#include <string>
#include <vector>

#include "CCDB/CcdbApi.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"

using std::map;
using std::string;

void GetRecenter(TH2* h, std::vector<float>& corr)
{
    corr.push_back(h->GetMean(1));
    corr.push_back(h->GetMean(2));
}

double CalcB( double rho, double sigmax, double sigmay ){
    return rho * sigmax * sigmay * TMath::Sqrt(2.0 * (sigmax * sigmax + sigmay * sigmay - 2.0 * sigmax * sigmay * TMath::Sqrt(1.0 - rho * rho)) / ((sigmax * sigmax - sigmay * sigmay) * (sigmax * sigmax - sigmay * sigmay) + 4.0 * (sigmax * sigmay * rho) * (sigmax * sigmay * rho)));
}

void GetTwist(TH2* h, std::vector<float>& corr)
{
    double b = CalcB( h->GetCorrelationFactor(), h->GetStdDev(1), h->GetStdDev(2) );

    float aPlus = (float)TMath::Sqrt(2. * TMath::Power(h->GetStdDev(1), 2.) - TMath::Power(b, 2.));
    float aMinus = (float)TMath::Sqrt(2. * TMath::Power(h->GetStdDev(2), 2.) - TMath::Power(b, 2.));

    corr.push_back(b / aPlus);
    corr.push_back(b / aMinus);
}

void GetRescale(TH2* h, std::vector<float>& corr)
{
    double b = CalcB( h->GetCorrelationFactor(), h->GetStdDev(1), h->GetStdDev(2) );

    float aPlus = (float)TMath::Sqrt(2. * TMath::Power(h->GetStdDev(1), 2.) - TMath::Power(b, 2.));
    float aMinus = (float)TMath::Sqrt(2. * TMath::Power(h->GetStdDev(2), 2.) - TMath::Power(b, 2.));

    corr.push_back(aPlus);
    corr.push_back(aMinus);
}

// \param fInput Input AnalysisResult.root where the histograms Qx-Qy can be found.
// NOTE: period and pass are some of the metadata for the ccdb.
// NOTE: run and harmonic disabled for now. They will be renabled if a run/harmonic dependence is seen.
void upload_event_plane_correction_constants(const string fInput = "/home/cindy/cernbox/MyProjects/EP_O2-Tests/Train135707_AnalysisResults.root",
                                             const string ccdbPath = "http://ccdb-test.cern.ch:8080",
                                            //const string ccdbPath = "http://alice-ccdb.cern.ch",
                                             const string period = "LHC22s",
                                             const string pass = "pass5")
                                            //const string run = "",
                                            //const string detector = "",
                                            //const string harmonic = "2")
{
    // Init the access to the CCDB.
    o2::ccdb::CcdbApi ccdb;
    map<string, string> metadata;//, metadataRCT, header;   // NOTE: Re-enable the other two if timing information is needed.
    ccdb.init(Form("%s", ccdbPath.data()));

    // Format the vector of correction constants for each centrality class and detector.
    const int nDet = 3; //6;    // FIXME: Update when the final structure with all detectors is in.
    const char detNames[nDet][1000] = {
        "FT0C",
        //"FT0A",
        //"FT0M",
        //"FV0A",
        "BPos",
        "BNeg"
    };

    const int nBins = 8;
    const char dirNames[nBins][1000] = {
        "Centrality_0-5",
        "Centrality_5-10",
        "Centrality_10-20",
        "Centrality_20-30",
        "Centrality_30-40",
        "Centrality_40-50",
        "Centrality_50-60",
        "Centrality_60-80"
    };

    TFile* fileInput = new TFile(Form("%s", fInput.data()), "read");
    TDirectory* dirInput;
    TH2D* hQv;

    std::vector<float> corrConst;

    for (int j = 0; j < nDet; j++) {
        for (int i = 0; i < nBins; i++) {
            dirInput = (TDirectoryFile*)fileInput->GetDirectory(Form("q-vectors-correction/%s", dirNames[i]));
            if (j == 0) {hQv = (TH2D*)dirInput->Get("histQvecUncor");}    // FIXME: Update when the final structure with all detectors is in.
            else {hQv = (TH2D*)dirInput->Get(Form("histQvec%sUncor", detNames[j]));}
            if (!hQv) {
                printf("Histogram not found\n");
                return;
            }

            // Get all the constants for this centrality class, in the order required by 'qVectorTable.cxx'.
            GetRecenter(hQv, corrConst);
            GetTwist(hQv, corrConst);
            GetRescale(hQv, corrConst);
        }
        /* 
        printf("Vector: ");
        for (int i = 0; i < nBins*6; i++) {
            printf("%e ,", corrConst[i]);
        }
        printf("\n\n");
        */

        // Prepare the metadata information for this set of correction constants.
        metadata["period"] = period;
        metadata["pass"] = pass; 
        metadata["detector"] = detNames[j];
        //metadata["harmonics"] = harmonic;
        //metadata["runNumber"] = run;

        //header = ccdb.retrieveHeaders(Form("RCT/Info/RunInformation/%i", run), metadataRCT, -1);
        //ULong64_t sor = atol(header["SOR"].c_str());
        //ULong64_t eor = atol(header["EOR"].c_str());
        ULong64_t sor = 1672531200000;
        ULong64_t eor = 1893456000000;

        ccdb.storeAsTFileAny(&corrConst, Form("EventPlaneCalib/%s", detNames[j]), metadata, sor, eor);

        // Reset the information for the next detector.
        corrConst.clear();
        metadata.clear();
    }   // Go to the next detector.

    return;
}
