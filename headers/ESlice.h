#include "TGraphErrors.h"
#include "TVector3.h"
//#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "util.h"
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
//#include "./Unfold.h"

//#include "BetheBloch.h"

//Basic config. -----------------------------------------------------//
std::string fOutputFileName;
TFile *outputFile;
void SetOutputFileName(std::string name){fOutputFileName = name;};
void BookHistograms();
void SaveHistograms();

//XS histograms ------------------------//
//reco inc
TH1D *h_recosliceid_allevts_cuts;
//reco st inc
TH1D *h_reco_st_sliceid_allevts_cuts;
//int
TH1D *h_recosliceid_recoinelastic_cuts;
//--------------------------------------//

void BookHistograms() { //BookHistograms
	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); //open output file

	//XS histograms -------------------------------------------------------------------------------------------------------------------------------------------------------//
	h_recosliceid_allevts_cuts = new TH1D("h_recosliceid_allevts_cuts","h_recosliceid_allevts_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_reco_st_sliceid_allevts_cuts = new TH1D("h_reco_st_sliceid_allevts_cuts","h_reco_st_sliceid_allevts_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_recoinelastic_cuts = new TH1D("h_recosliceid_recoinelastic_cuts","h_recosliceid_recoinelastic_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);

	h_recosliceid_allevts_cuts->Sumw2();
	h_reco_st_sliceid_allevts_cuts->Sumw2();
	h_recosliceid_recoinelastic_cuts->Sumw2();

} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();

} //SaveHistograms

