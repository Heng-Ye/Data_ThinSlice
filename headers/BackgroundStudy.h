#include "TGraphErrors.h"
#include "TVector3.h"
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

//Basic config. -----------------------------------------------------//
std::string fOutputFileName;
TFile *outputFile;
void SetOutputFileName(std::string name){fOutputFileName = name;};
void BookHistograms();
void SaveHistograms();

//reco sliceID histograms
//misID:p-rich sample
TH1D *h_recosliceid_cosLE09;
TH1D *h_recosliceid_cosLE08;
TH1D *h_recosliceid_cosLE07;
TH1D *h_recosliceid_cosGT09;

TH1D *h_cosTheta_Pos[nthinslices+2];
TH1D *h_cosTheta_Pos_all;
TH1D *h_cosTheta_Pos_all_nosliceidOne;

TH1D *h_chi2_BQ[nthinslices+2];

TH2D *trklen_cosineTheta_Pos;

//dxy vs cosine
TH2D *h2d_dxy_cosine_BQ;
TH1D *h1d_dxy_BQ;
TH2D *h2d_dxy_cosine_Pos;
TH1D *h1d_dxy_Pos;

//pos+recoinel
TH2D *h2d_trklen_cos_PosRecoInel;

void BookHistograms() { //BookHistograms
	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); //open output file

	h_recosliceid_cosLE09 = new TH1D("h_recosliceid_cosLE09","h_recosliceid_cosLE09;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE09->Sumw2();
	h_recosliceid_cosLE08 = new TH1D("h_recosliceid_cosLE08","h_recosliceid_cosLE08;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE08->Sumw2();
	h_recosliceid_cosLE07 = new TH1D("h_recosliceid_cosLE07","h_recosliceid_cosLE07;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE07->Sumw2();

	h_recosliceid_cosGT09 = new TH1D("h_recosliceid_cosGT09","h_recosliceid_cosGT09;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosGT09->Sumw2();

        //cosTheta dists.
        int n_cosine=100;
        double cosine_min=0;
        double cosine_max=1.0;

        //int nn_cosine=2; //
	//double cosbins[3];
	//cosbins[0]=0;
	//cosbins[1]=0.9;
	//cosbins[2]=1;

        int nn_cosine=3; //
	double cosbins[4];
	cosbins[0] = 0;
	cosbins[1] = 0.9;
	cosbins[2] = 0.96;
	cosbins[3] = 1.;

        for (int i=0; i<nthinslices+2; ++i) {
                h_cosTheta_Pos[i]=new TH1D(Form("h_cosTheta_Pos_recosliceid_%d", i-1), Form("Pos; Reco SliceID %d",i-1), nn_cosine, cosbins);
                h_cosTheta_Pos[i]->Sumw2();
        }
        h_cosTheta_Pos_all=new TH1D("h_cosTheta_Pos_all","h_cosTheta_Pos_all", n_cosine, cosine_min, cosine_max);
        h_cosTheta_Pos_all->Sumw2();
        h_cosTheta_Pos_all_nosliceidOne=new TH1D("h_cosTheta_Pos_all_nosliceidOne","h_cosTheta_Pos_all_nosliceidOne",n_cosine, cosine_min, cosine_max);
        h_cosTheta_Pos_all_nosliceidOne->Sumw2();
     
        int nn_chi2=2; //
	double chi2bins[3];
	chi2bins[0] = 0;
	chi2bins[1] = 10;
	chi2bins[2] = 150;
	for (int i=0; i<nthinslices+2; ++i) {
		h_chi2_BQ[i]=new TH1D(Form("h_chi2_BQ_recosliceid_%d", i-1), Form("BQ; Reco SliceID %d",i-1), nn_chi2, chi2bins);
		h_chi2_BQ[i]->Sumw2();
	}

	int n_b=300;
	double b_min=0;
	double b_max=150;

	trklen_cosineTheta_Pos=new TH2D(Form("trklen_cosineTheta_Pos"), Form(""), n_b, b_min, b_max, n_cosine, cosine_min, cosine_max);


	int n_dxy=500;
	double dxy_min=0;
	double dxy_max=50;
	int n_cos1=100;
	double cos1_min=0.;
	double cos1_max=1;
	h2d_dxy_cosine_BQ=new TH2D("h2d_dxy_cosine_BQ","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);
	h2d_dxy_cosine_Pos=new TH2D("h2d_dxy_cosine_Pos","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);

	h1d_dxy_BQ=new TH1D("h1d_dxy_BQ","",n_dxy,dxy_min,dxy_max);
	h1d_dxy_Pos=new TH1D("h1d_dxy_Pos","",n_dxy,dxy_min,dxy_max);

	h2d_dxy_cosine_BQ->Sumw2();
	h2d_dxy_cosine_Pos->Sumw2();
	h1d_dxy_BQ->Sumw2();
	h1d_dxy_Pos->Sumw2();

	h2d_trklen_cos_PosRecoInel=new TH2D(Form("h2d_trklen_cos_PosRecoInel"), Form("h2d_trklen_cos_PosRecoInel"), n_b, b_min, b_max, n_cosine, cosine_min, cosine_max);
        h2d_trklen_cos_PosRecoInel->Sumw2();

} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();

} //SaveHistograms


