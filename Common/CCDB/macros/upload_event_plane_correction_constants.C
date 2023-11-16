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

// NOTE: run and harmonic disabled for now. They will be renabled if a run/harmonic dependence is seen.
void upload_event_plane_correction_constants(const string ccdbPath = "http://ccdb-test.cern.ch:8080",
                                            const string fInput = "/home/cindy/cernbox/MyProjects/EP_O2-Tests/AnalysisResults.root",
                                            const string period = "LHC22s",
                                            const string pass = "pass5")//,
                                            //const string run = "",
                                            //const string detector = "",
                                            //const string harmonic = "2")
{
    // Open the connexion to the CCDB.
    o2::ccdb::CcdbApi ccdb;
    ccdb.init(Form("%s", ccdbPath.data()));

    // Format the vector of correction constants for each centrality class.
    const int nbins = 8;
    const char dirNames[nbins][1000] = {
        "Centrality_0-5",
        "Centrality_5-10",
        "Centrality_10-20",
        "Centrality_20-30",
        "Centrality_30-40",
        "Centrality_40-50",
        "Centrality_50-60",
        "Centrality_60-80" };    

    TFile* fileInput = new TFile(Form("%s", fInput.data()), "read");
    TDirectory* dirInput;
    TH2D* hQv;

    std::vector<float> corrConst;

    for (int i = 0; i < nbins; i++) {
        dirInput = (TDirectoryFile*)fileInput->GetDirectory(Form("q-vectors-correction/%s", dirNames[i]));
        hQv = (TH2D*)dirInput->Get("histQvecUncor");

        // Get all the constants for this centrality class, in the order required by 'qVectorTable.cxx'.
        GetRecenter(hQv, corrConst);
        GetTwist(hQv, corrConst);
        GetRescale(hQv, corrConst);
    }
    printf("Vector: ");
    for (int i = 0; i < nbins*6; i++) {
        printf("%e ,", corrConst[i]);
    }
    printf("\n");

    // Prepare the metadata information for this set of correction constants.
    map<string, string> metadata, metadataRCT, header;
    metadata["period"] = period;
    metadata["pass"] = pass;
    //metadata["runNumber"] = run;
    //metadata["detector"] = detector;
    //metadata["harmonics"] = harmonic;

    //header = ccdb.retrieveHeaders(Form("RCT/Info/RunInformation/%i", run), metadataRCT, -1);
    //ULong64_t sor = atol(header["SOR"].c_str());
    //ULong64_t eor = atol(header["EOR"].c_str());
    ULong64_t sor = 1672531200000;
    ULong64_t eor = 1893456000000;

    ccdb.storeAsTFileAny(&corrConst, "EventPlaneCalib", metadata, sor, eor);

    return;
}
