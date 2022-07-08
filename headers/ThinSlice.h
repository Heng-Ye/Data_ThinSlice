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

//Basic config. -----------------------------------------------------//
std::string fOutputFileName;
TFile *outputFile;
void SetOutputFileName(std::string name){fOutputFileName = name;};
void BookHistograms();
void SaveHistograms();

//XS histograms ----------------------//

//Histograms for basic parameters ---------//
//reco x, y, z [after SCE corr]
TH1D *reco_startX_sce; 
TH1D *reco_startY_sce;
TH1D *reco_startZ_sce;

TH1D *hdeltaX;
TH1D *hdeltaY;
TH1D *hdeltaZ;
TH1D *hdeltaXY; //similiar to XY cut

TH1D *reco_cosineTheta;
TH1D *reco_cosineTheta_Pos;

//dE/dx vs rr related histograms ------------//
TH2D *rr_dedx_recostop;

//KE calc using reco stopping protons ----------------------------------------------------//
TH1D *KE_calo_recostop;
TH1D *KE_rrange_recostop;
TH1D *KE_rrange2_recostop;
TH1D *KE_range_recostop;
TH2D *KE_range_calo_recostop;

//Track length histograms -----------------------------------------------------------------//
//reco range
TH1D *trklen_reco_NoCut; //no cut
TH1D *trklen_reco_PanS; //pandora slice cu
TH1D *trklen_reco_CaloSz; //calosz cut
TH1D *trklen_reco_BQ; //beam quality cut
TH1D *trklen_reco_RecoInel; //RecoInel cut

//ntrklen histograms ----------//
TH1D *ntrklen_BQ; //beam quality cut

//xs slice histograms
TH1D *h_recosliceid_allevts_cuts; //inc
TH1D *h_recosliceid_inelastic_cuts; //int

void BookHistograms() { //BookHistograms
	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); //open output file

	//XS histograms -------------------------------------------------------------------------------------------------------------------------------------------------------//

	//Histograms for basic parameters ----------------------------------------------------------------------------------------------------//
	reco_startX_sce = new TH1D(Form("reco_startX_sce"), Form("reco_startX_sce"), 100, -80, 20);  reco_startX_sce->Sumw2();
	reco_startY_sce = new TH1D(Form("reco_startY_sce"), Form("reco_startY_sce"), 100, 350, 500); reco_startY_sce->Sumw2();
	reco_startZ_sce = new TH1D(Form("reco_startZ_sce"), Form("reco_startZ_sce"), 100, -5, 10);   reco_startZ_sce->Sumw2();

	hdeltaX = new TH1D(Form("hdeltaX"), Form("#Deltax/#sigma_{x}"), 100, -10, 10);  hdeltaX->Sumw2();
	hdeltaY = new TH1D(Form("hdeltaY"), Form("#Deltay/#sigma_{y}"), 100, -10, 10);  hdeltaY->Sumw2();
	hdeltaZ = new TH1D(Form("hdeltaZ"), Form("#Deltaz/#sigma_{z}"), 100, -10, 10);  hdeltaZ->Sumw2();
	hdeltaXY = new TH1D(Form("hdeltaXY"), Form("Sqrt((#Deltax/#sigma_{x})^2+(#Deltay/#sigma_{y})^2)"), 100, -10, 10);  hdeltaXY->Sumw2();

        int n_cosine=100;
        //double cosine_min=0.9;
        double cosine_min=0;
        double cosine_max=1.0;
	reco_cosineTheta = new TH1D("reco_cosineTheta","", n_cosine, cosine_min, cosine_max);	reco_cosineTheta->Sumw2();
	reco_cosineTheta_Pos = new TH1D("reco_cosineTheta_Pos","", n_cosine, cosine_min, cosine_max);	reco_cosineTheta_Pos->Sumw2();

	//dE/dx vs rr related histograms ---------------------------------------------------------//
	rr_dedx_recostop=new TH2D("rr_dedx_recostop","", 240,0,120, 300,0, 30);

	//KE calc using reco stopping protons ----------------------------------------------------//
	int n_ke=140;
	float ke_min=0;
	float ke_max=700;
	KE_calo_recostop=new TH1D("KE_calo_recostop","",n_ke, ke_min, ke_max); KE_calo_recostop->Sumw2();
	KE_rrange_recostop=new TH1D("KE_rrange_recostop", "", n_ke, ke_min, ke_max); KE_rrange_recostop->Sumw2();
	KE_rrange2_recostop=new TH1D("KE_rrange2_recostop", "", n_ke, ke_min, ke_max); KE_rrange2_recostop->Sumw2();
	KE_range_recostop=new TH1D("KE_range_recostop", "", n_ke, ke_min, ke_max); KE_range_recostop->Sumw2();
	KE_range_calo_recostop=new TH2D("KE_range_calo_recostop","",n_ke, ke_min, ke_max, n_ke, ke_min, ke_max); 
	KE_range_calo_recostop->GetXaxis()->SetTitle("KE_{range} [MeV]"); KE_range_calo_recostop->GetYaxis()->SetTitle("KE_{calo} [MeV]"); 
	KE_range_calo_recostop->Sumw2();

	//trklen ------------------------------------------------------------------------------------------------------------------//
        int n_trklen=34;
        double trklen_min=-4;
        double trklen_max=132;

	//reco range
	//no cut
	trklen_reco_NoCut = new TH1D("trklen_reco_NoCut","", n_trklen, trklen_min, trklen_max); trklen_reco_NoCut->Sumw2();

	//pandora cut
	trklen_reco_PanS = new TH1D("trklen_reco_PanS","", n_trklen, trklen_min, trklen_max); trklen_reco_PanS->Sumw2();

	//CaloSz
	trklen_reco_CaloSz = new TH1D("trklen_reco_CaloSz","", n_trklen, trklen_min, trklen_max); trklen_reco_CaloSz->Sumw2();

	//beam quality
	trklen_reco_BQ = new TH1D("trklen_reco_BQ","",n_trklen, trklen_min, trklen_max); trklen_reco_BQ->Sumw2();

	//reco inel cut
	trklen_reco_RecoInel = new TH1D("trklen_reco_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_reco_RecoInel->Sumw2();

	//ntrklen ------------------------------------------------------------------------------------//
	int n_ntrklen=61;
	double st_ntrklen=-0.2;
	double ed_ntrklen=1.2;
	//bq cut
	ntrklen_BQ=new TH1D("ntrklen_BQ", "", n_ntrklen, st_ntrklen, ed_ntrklen);

	//xs sliceID
	h_recosliceid_allevts_cuts = new TH1D("h_recosliceid_allevts_cuts","h_recosliceid_allevts_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_inelastic_cuts = new TH1D("h_recosliceid_inelastic_cuts","h_recosliceid_inelastic_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts->Sumw2();
	h_recosliceid_inelastic_cuts->Sumw2();

} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();
} //SaveHistograms

//void CalcXS(const Unfold & uf) { //CalcXS

//} //CalcXS

