#define ProtonMomentumReweight_run5387_cxx
#include "ProtonMomentumReweight_run5387.h"

#include <stdexcept>      // std::out_of_range
#include <vector>
#include <iostream>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TProfile2D.h>
//#include <iostream>
#include <fstream>
#include <string>
#include "TCanvas.h" 
#include "TGraph.h" 
#include "TGraphSmooth.h" 
#include "TVectorD.h"
#include "TParameter.h"
#include "TGraphErrors.h"
#include <TProfile3D.h>
#include <TProfile2D.h>
#include "TVector3.h"

#include <TMath.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>

//#include "./cali/dedx_function_r5387.h"
#include "./headers/BasicParameters.h"
#include "./headers/BasicAnaFunc.h"
#include "./headers/util.h"
#include "./headers/BetheBloch.h"

#include <cassert>

using namespace std;
using namespace ROOT::Math;

void ProtonMomentumReweight_run5387::Loop() {
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	//time info ---------------------------------------------------------//
	Long64_t nnbytes = 0, nnb = 0;
	int n_stopping_protons=0;
	Double_t t0=999999999999999999999999999.;
	Double_t tmax=-1.;
	for (Long64_t kentry=0; kentry<nentries;kentry++) { //evt loop to get t0
		Long64_t kkentry = LoadTree(kentry);
		if (kkentry < 0) break;
		nnb = fChain->GetEntry(kentry);   nnbytes += nnb;
		if (evttime<t0) {
			t0=evttime;
		}
		if (evttime>tmax) {
			tmax=evttime;
		}
	} //evt loop to get t0

	//summary of total time of data taking --------------//
	//Double_t tot_sec=(tmax-t0);
	//Double_t tot_min=tot_sec/60.;
	//Double_t tot_hr=tot_min/60.;
	//std::cout<<"t0:"<<t0-t0<<std::endl;
	//std::cout<<"tmax:"<<tmax-t0<<std::endl;
	//std::cout<<"tot:"<<tot_hr<<" [hrs]"<<std::endl;
	//-------------------------------------------------------------------//

	//Basic configure ------//
	BetheBloch BB;
	BB.SetPdgCode(2212);
	//----------------------//

	//book histograms -------------------------------------------------------------------------//
	//trklen
	int n_b=30;
	//int n_b=150;
	double b_min=0;
	double b_max=150;
	TH1D *h1d_dist_stop=new TH1D(Form("h1d_dist_stop"),Form("reco stop"),n_b,b_min,b_max);
	TH1D *h1d_trklen_stop=new TH1D(Form("h1d_trklen_stop"),Form("reco stop"),n_b,b_min,b_max);
	TH1D *h1d_trklen_stop_XY=new TH1D(Form("h1d_trklen_stop_XY"),Form("reco stop with xy cut"),n_b,b_min,b_max);
	TH1D *h1d_trklen_XY=new TH1D(Form("h1d_trklen_XY"),Form("reco"),n_b,b_min,b_max);

	TH1D *h1d_trklen_CaloSz=new TH1D(Form("h1d_trklen_CaloSz"),Form("calosz"),n_b,b_min,b_max);
	TH1D *h1d_trklen_Pos=new TH1D(Form("h1d_trklen_Pos"),Form("Pos"),n_b,b_min,b_max);
	TH1D *h1d_trklen_BQ=new TH1D(Form("h1d_trklen_BQ"),Form("bq"),n_b,b_min,b_max);
	TH1D *h1d_trklen_RecoInel=new TH1D(Form("h1d_trklen_RecoInel"), Form("reco inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich=new TH1D(Form("h1d_trklen_MidpRich"), Form("misidp_rich"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich=new TH1D(Form("h1d_trklen_ElRich"), Form("el rich"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich=new TH1D(Form("h1d_trklen_CaloSz_ElRich"), Form("el rich"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich=new TH1D(Form("h1d_trklen_Pos_ElRich"), Form("el rich"), n_b, b_min, b_max);

	//KEs
	int ny_edept=450;
	double ymin_edept=-100;
	double ymax_edept=800;
	TH1D *h1d_kebeam=new TH1D("h1d_kebeam","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_kebeam_stop=new TH1D("h1d_kebeam_stop","",ny_edept,ymin_edept,ymax_edept);

	TH1D *h1d_kehy=new TH1D("h1d_kehy","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_kehy_stop=new TH1D("h1d_kehy_stop","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_kehy_inel=new TH1D("h1d_kehy_inel","",ny_edept,ymin_edept,ymax_edept);

	TH1D *h1d_kerange_stop=new TH1D("h1d_kerange_stop","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kecalo_stop=new TH1D("h1d_kecalo_stop","", ny_edept, ymin_edept, ymax_edept);

	TH1D *h1d_kend_calo_stop=new TH1D("h1d_kend_calo_stop","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_calo_el=new TH1D("h1d_kend_calo_el","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_calo_inel=new TH1D("h1d_kend_calo_inel","", ny_edept, ymin_edept, ymax_edept);

	TH1D *h1d_kend_bb_stop=new TH1D("h1d_kend_bb_stop","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_bb_el=new TH1D("h1d_kend_bb_el","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_bb_inel=new TH1D("h1d_kend_bb_inel","", ny_edept, ymin_edept, ymax_edept);



	//mom
	int nx=250;	
	double xmin=0.; //pmin [MeV/c]
	double xmax=2000.; //pmax [MeV/c]
	TH1D *h1d_pbeam=new TH1D("h1d_pbeam","",nx,xmin,xmax);
	TH1D *h1d_pbeam_stop=new TH1D("h1d_pbeam_stop","",nx,xmin,xmax);

	TH1D *h1d_phy=new TH1D("h1d_phy","",nx,xmin,xmax);
	TH1D *h1d_phy_stop=new TH1D("h1d_phy_stop","",nx,xmin,xmax);

	TH1D *h1d_prange_stop=new TH1D("h1d_prange_stop","",nx,xmin,xmax);
	TH1D *h1d_pcalo_stop=new TH1D("h1d_pcalo_stop","",nx,xmin,xmax);

	//dedx_rr
	TH2D *h2d_rr_dedx_recoSTOP=new TH2D("h2d_rr_dedx_recoSTOP","",240,0,120,90,0,30);

	//chi2_pid
	int n_ntrklen=61;
	double st_ntrklen=-0.02;
	double ed_ntrklen=1.2;

	int n_chi2=750;
	double st_chi2=0;
	double ed_chi2=150;
	TH2D *ntrklen_chi2pid_BQ=new TH2D("ntrklen_chi2pid_BQ","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_BQ=new TH1D("h1d_chi2pid_BQ","", n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_ntrklen_BQ=new TH1D("h1d_ntrklen_BQ","", n_ntrklen, st_ntrklen, ed_ntrklen);

	TH1D *chi2pid_recostop=new TH1D("chi2pid_recostop","",500,0,100);
	TH1D *chi2pid_recoinel=new TH1D("chi2pid_recoinel","",500,0,100);

	//xy-dist
	TH2D *h2d_xy_noSCE=new TH2D("h2d_xy_noSCE","", 70,-60,10,60,390,450); //nosce
	TH2D *h2d_xy_SCE=new TH2D("h2d_xy_SCE","", 70,-60,10,60,390,450); //after sce

	//z-st
	TH1D *h1d_zst_noSCE_stop=new TH1D("h1d_zst_noSCE_stop","",220,-10,100);
	TH1D *h1d_zst_stop=new TH1D("h1d_zst_stop","",220,-10,100);

	//zend
	int dz=500;
	float z_st=0;
	float z_end=250;
	TH1D *h1d_zend_stop = new TH1D("h1d_zend_stop", "", dz, z_st, z_end);
	TH1D *h1d_zend_noSCE_stop=new TH1D("h1d_zend_noSCE_stop","",110,-10,100);
	//TH1D *h1d_zend_XY = new TH1D("h1d_zend_XY", "reco+BQ+XY", dz, z_st, z_end);	

	//median dedx
	TH1D *h1d_mediandedx_BQ=new TH1D("h1d_mediandedx_BQ","",100,0,10);

	//Edept/L
	TH1D *h1d_dEdL_BQ=new TH1D("h1d_dEdL_BQ","",100,0,10);	

	//Time dist. for
	int n_t=24*60;
	double t_min=0;
	double t_max=24; 
	TH2D *h2d_time_pcalo_stop=new TH2D("h2d_time_pcalo_stop","",n_t,t_min,t_max,nx,xmin,xmax); //x in min
	TH2D *h2d_time_prange_stop=new TH2D("h2d_time_prange_stop","",n_t,t_min,t_max,nx,xmin,xmax); //x in min
	TH2D *h2d_time_pcaloOverprange_stop=new TH2D("h2d_time_pcaloOverprange_stop","",n_t,t_min,t_max, 25000,-5,20); //x in min
	TH2D *h2d_time_zst_stop=new TH2D("h2d_time_zst_stop","",n_t,t_min,t_max,220,-10,100); //x in min
	TH2D *h2d_time_zst_noSCE_stop=new TH2D("h2d_time_zst_noSCE_stop","",n_t,t_min,t_max,1100,-10,100); //x in min

	//Fitted Data Beam Momentum
	double m1=1013.71; //Data prod4 reco2
	double s1=68.6327; //Data prod4 reco2

	//momentum cut range	
	double mu_min=m1-3.*s1;
	double mu_max=m1+3.*s1;


	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) { //evt loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

		//size_t j=0;
		//if (Cut(jentry) < 0) continue;

		//define good beam trigger -----------------------------------------------------//
		bool IsBeamInstValid=false;
		bool IsBeamInst_nMomenta=false;
		bool IsBeamInst_nTracks=false;
		bool IsBeamOK=false;

		if (beam_inst_valid) IsBeamInstValid=true;
		if (beam_inst_nMomenta==1) IsBeamInst_nMomenta=true;
		if (beam_inst_nTracks==1) IsBeamInst_nTracks=true;	
		if (IsBeamInstValid&&IsBeamInst_nMomenta&&IsBeamInst_nTracks) IsBeamOK=true;

		if (!IsBeamOK) continue; //skip evt if beam not okay

		//Incoming candidate event ------------------------------------------------------------------------//	
		bool IsProton=false;
		if (!beam_inst_PDG_candidates->empty()&&IsBeamInstValid) {
			//n_can++;
			for (int h=0; h<(int)beam_inst_PDG_candidates->size(); ++h) { 
				//cout<<"beam_inst_PDG_candidates["<<h<<"]:"<<beam_inst_PDG_candidates->at(h)<<endl;
				if (beam_inst_PDG_candidates->at(h)==2212) {
					//n_proton++;
					IsProton=true;
				}
			}
		}
		if (!IsProton) continue; //only protons, not other candidate particles
		//(the ntuple actually only save proton evt for the beam track)

		//Pandora slice cut ------------------------------------------------------------------//
		bool IsPandoraSlice=false; //pandora slice
		bool IsCaloSize=false; //if calo size not empty

		if (isprimarytrack==1&&isprimaryshower==0) { //pandora slice cut
			IsPandoraSlice=true;
		}

		if (!primtrk_hitz->empty()) { 
			IsCaloSize=true;
		}

		//intersection cut
		bool IsIntersection=false; //if any track intersect with our reco track
		if (timeintersection->size()) IsIntersection=true;


		//reco pos info & cut --------------------------------------------------------//
		double reco_stx=-99, reco_sty=-99, reco_stz=-99;
		double reco_endx=-99, reco_endy=-99, reco_endz=-99;
		bool IsPos=false;
		if (IsCaloSize) { //calosz
			if (primtrk_hitz->at(0)>primtrk_hitz->at(primtrk_hitz->size()-1)) {
				reco_endx=primtrk_hitx->at(0);
				reco_endy=primtrk_hity->at(0);
				reco_endz=primtrk_hitz->at(0);

				reco_stx=primtrk_hitx->at(primtrk_dedx->size()-1);
				reco_sty=primtrk_hity->at(primtrk_dedx->size()-1);
				reco_stz=primtrk_hitz->at(primtrk_dedx->size()-1);
			}
			else {
				reco_stx=primtrk_hitx->at(0);
				reco_sty=primtrk_hity->at(0);
				reco_stz=primtrk_hitz->at(0);

				reco_endx=primtrk_hitx->at(primtrk_dedx->size()-1);
				reco_endy=primtrk_hity->at(primtrk_dedx->size()-1);
				reco_endz=primtrk_hitz->at(primtrk_dedx->size()-1);
			}

			double beam_dx=(reco_stx-mean_StartX)/sigma_StartX;
			double beam_dy=(reco_sty-mean_StartY)/sigma_StartY;
			double beam_dz=(reco_stz-mean_StartZ)/sigma_StartZ;
			double beam_dxy=sqrt(pow(beam_dx,2)+pow(beam_dy,2));

			if (beam_dx>=dx_min&&beam_dx<=dx_max) { //dx
				if (beam_dy>=dy_min&&beam_dy<=dy_max) { //dy
					if (beam_dz>=dz_min&&beam_dz<=dz_max) { //dz
						if (beam_dxy>=dxy_min&&beam_dxy<=dxy_max) { //dxy
							IsPos=true;
						} //dxy
					} //dz
				} //dy
			} //dx
		} //calosz

		//cosine_theta/cut ---------------------------------------------------------------------------------//
		bool IsCosine=false;
		double cosine_beam_spec_primtrk=-999;
		//cosine_beam_spec_primtrk=cosine_beam_primtrk; //cosine between beam_spec and primary trk direction
		TVector3 dir;
		if (IsCaloSize) {
			//trk direction after SCE corr.
			TVector3 pt0(primtrk_hitx->at(0), primtrk_hity->at(0), primtrk_hitz->at(0));
			TVector3 pt1(primtrk_hitx->at(-1+primtrk_hitx->size()), primtrk_hity->at(-1+primtrk_hity->size()), primtrk_hitz->at(-1+primtrk_hitz->size()));
			dir = pt1 - pt0;
			dir = dir.Unit();

			//beam direction
			//TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180), cos(beam_angleY_mc*TMath::Pi()/180), cos(beam_angleZ_mc*TMath::Pi()/180));
			TVector3 beamdir(beamDirx->at(0), beamDiry->at(0), beamDirz->at(0));
			beamdir = beamdir.Unit();
			cosine_beam_spec_primtrk=dir.Dot(beamdir);
		}
		if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
		//if (cosine_beam_spec_primtrk>min_cosine) { IsCosine=true; }
		if (cosine_beam_spec_primtrk>costh_min&&cosine_beam_spec_primtrk<costh_max) { IsCosine=true; }		


		//beam quality cut -----------------------------//
		bool IsBQ=false;
		if (IsCosine&&IsPos) IsBQ=true;
		//if (IsBQ&&IsCaloSize&&IsPandoraSlice) n_bq++;

		//beam info //
		double bx=beamPosx->at(0);
		double by=beamPosy->at(0);
		double bz=beamPosz->at(0);
		double p_beam=beamMomentum->at(0);
		double ke_beam_MeV=1000.*p2ke(p_beam); //unit:MeV

		//Beam XY cut
		bool IsBeamXY=false;
		if ((pow(((bx-meanX_data)/(1.5*rmsX_data)),2)+pow(((by-meanY_data)/(1.5*rmsY_data)),2))<=1.) IsBeamXY=true;

		//beam-mom cut (within 3-sigma)
		bool IsBeamMom=false;
		if ((1000.*p_beam)>=mu_min&&(1000.*p_beam)<=mu_max) IsBeamMom=true;

		//Get calo info -----------------------------------------------------------------------------------------------//
		double range_reco=-99;
		vector<double> reco_trklen_accum;
		reco_trklen_accum.reserve(primtrk_hitz->size());
		vector<double> EDept;
		double pid=-99;

		if (IsCaloSize) { //if calo size not empty
			vector<double> trkdedx;
			vector<double> trkres;
			for (size_t h=0; h<primtrk_dedx->size(); ++h) { //loop over reco hits in a track
				double hitx_reco=primtrk_hitx->at(h);
				double hity_reco=primtrk_hity->at(h);
				double hitz_reco=primtrk_hitz->at(h);
				double resrange_reco=primtrk_resrange->at(h);

				double dqdx=primtrk_dqdx->at(h);
				double dedx=primtrk_dedx->at(h); //prod4-reco 2 has dedx already
				double pitch=primtrk_pitch->at(h);

				int wid_reco=primtrk_wid->at(-1+primtrk_wid->size()-h);
				double pt_reco=primtrk_pt->at(-1+primtrk_wid->size()-h);

				//double cali_dedx=0.;
				//cali_dedx=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);
				EDept.push_back(dedx*pitch);

				if (h==1) range_reco=0;
				if (h>=1) {
					range_reco += sqrt( pow(primtrk_hitx->at(h)-primtrk_hitx->at(h-1), 2)+
							pow(primtrk_hity->at(h)-primtrk_hity->at(h-1), 2)+
							pow(primtrk_hitz->at(h)-primtrk_hitz->at(h-1), 2) );
					reco_trklen_accum[h] = range_reco;
				}

				trkdedx.push_back(dedx);
				trkres.push_back(resrange_reco);
			} //loop over reco hits in a track
			//range_reco=primtrk_range->at(0); 

			pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis

		} //if calo size not empty	


		double csda_val=csda_range_vs_mom_sm->Eval(p_beam); //expected csda value if proton stops; unit: cm
		double ke_range=ke_vs_csda_range_sm->Eval(range_reco); //unit:GeV
		double ke_range_MeV=1000.*ke_range; //unit:MeV
		double p_range_MeV=1000.*ke2p(ke_range_MeV/1000.); //MeV/c
		double norm_trklen=range_reco/csda_val; //trklen/csda_val
		//cout<<"csda_val:"<<csda_val<<" range_reco:"<<range_reco<<" norm_trklen:"<<norm_trklen<<" min_norm_trklen_csda:"<<min_norm_trklen_csda<<endl;

		bool IsRecoStop=false;
		bool IsRecoInel=false;
		bool IsRecoEl=false;
		//if (norm_trklen>min_norm_trklen_csda) IsRecoStop=true; //stopping proton cut
		//if (norm_trklen<=min_norm_trklen_csda) IsRecoInel=true; //inelastic proton cut
		//if (IsRecoInel&&IsBQ&&IsCaloSize&&IsPandoraSlice) n_recoinel++;

		if ((range_reco/csda_val)>=min_norm_trklen_csda&&(range_reco/csda_val)<max_norm_trklen_csda) IsRecoStop=true;

		if ((range_reco/csda_val)<min_norm_trklen_csda) { //inel region
			if (pid>pid_1) IsRecoInel=true;
			if (pid<=pid_1) IsRecoEl=true;
		} //inel region
		if ((range_reco/csda_val)>=min_norm_trklen_csda&&(range_reco/csda_val)<max_norm_trklen_csda) { //stopping p region
			if (pid>pid_2) IsRecoInel=true;
			if (pid<=pid_2) IsRecoEl=true;
		} //stopping p region

		//XY Cut ---------------------------------------------------------------------------------------------------------------//
		bool IsXY=false;
		//start(x,y,z) without SCE corr. 
		double xst_nosce=0;
		double yst_nosce=0;
		double zst_nosce=0;
		double xend_nosce=0;
		double yend_nosce=0;
		double zend_nosce=0;
		if (IsCaloSize&&IsPandoraSlice) { //if calo size not empty
			if ((primtrk_startz->at(-1+primtrk_startz->size()))>(primtrk_startz->at(0))) { //check if Pandora flip the sign
				xst_nosce=primtrk_startx->at(0);
				yst_nosce=primtrk_starty->at(0);
				zst_nosce=primtrk_startz->at(0);

				xend_nosce=primtrk_endx->at(-1+primtrk_startx->size());
				yend_nosce=primtrk_endy->at(-1+primtrk_starty->size());
				zend_nosce=primtrk_endz->at(-1+primtrk_startz->size());
			} //check if Pandora flip the sign
			else {
				xst_nosce=primtrk_startx->at(-1+primtrk_startx->size());
				yst_nosce=primtrk_starty->at(-1+primtrk_starty->size());
				zst_nosce=primtrk_startz->at(-1+primtrk_startz->size());

				xend_nosce=primtrk_endx->at(0);
				yend_nosce=primtrk_endy->at(0);
				zend_nosce=primtrk_endz->at(0);
			}
		} //if calo size not empty
		if ((pow(((xst_nosce-mean_x)/dev_x),2)+pow(((yst_nosce-mean_y)/dev_y),2))<=1.) IsXY=true;

		//calo info
		double ke_calo_MeV=0;
		double xst_sce=0;
		double yst_sce=0;
		double zst_sce=0;
		double xend_sce=0;
		double yend_sce=0;
		double zend_sce=0;
		double d_sce=0;
		vector<double> trkdedx; //dedx 
		vector<double> trkres; //rr
		//if (IsBQ&&IsCaloSize&&IsPandoraSlice) { //if calo size not empty
		//if (IsCaloSize) { //if calo size not empty
		if (IsBeamMom&&IsBeamXY&&IsBQ&&IsCaloSize&&IsPandoraSlice) {

			h2d_xy_noSCE->Fill(xst_nosce, yst_nosce);

			//start(x,y,z) after SCE corr. ----------------------------------------------------------------------------//
			if ((primtrk_hitz->at(-1+primtrk_hitz->size()))>(primtrk_hitz->at(0))) { //check if Pandora flip the sign
				xst_sce=primtrk_hitx->at(0);
				yst_sce=primtrk_hity->at(0);
				zst_sce=primtrk_hitz->at(0);

				xend_sce=primtrk_hitx->at(-1+primtrk_hitx->size());
				yend_sce=primtrk_hity->at(-1+primtrk_hity->size());
				zend_sce=primtrk_hitz->at(-1+primtrk_hitz->size());
			} //check if Pandora flip the sign
			else {
				xst_sce=primtrk_hitx->at(-1+primtrk_hitx->size());
				yst_sce=primtrk_hity->at(-1+primtrk_hity->size());
				zst_sce=primtrk_hitz->at(-1+primtrk_hitz->size());

				xend_sce=primtrk_hitx->at(0);
				yend_sce=primtrk_hity->at(0);
				zend_sce=primtrk_hitz->at(0);
			}
			h2d_xy_SCE->Fill(xst_sce,yst_sce);
			d_sce=sqrt(pow(xst_sce-xend_sce,2)+pow(yst_sce-yend_sce,2)+pow(zst_sce-zend_sce,2));

			//calo info --------------------------------------------------------------------------//
			for (size_t h=0; h<primtrk_dqdx->size(); ++h) { //loop over hits
				double dqdx=primtrk_dqdx->at(h);
				double dedx=primtrk_dedx->at(h);		
				double resrange_reco=primtrk_resrange->at(h);

				double hitx_reco=primtrk_hitx->at(h);
				double hity_reco=primtrk_hity->at(h);
				double hitz_reco=primtrk_hitz->at(h);
				double pitch=primtrk_pitch->at(h);
				int ch=primtrk_ch->at(h);

				ke_calo_MeV+=dedx*pitch;

				trkdedx.push_back(dedx); //vector for pid	
				trkres.push_back(resrange_reco); //vector for pid

				if (IsRecoStop) { //reco_stop
					h2d_rr_dedx_recoSTOP->Fill(resrange_reco, dedx);
				} //reco_stop	
			} //loop over hits

			double pid=-99; pid=chi2pid(trkdedx,trkres);
			double median_dedx=-99; median_dedx=TMath::Median(trkdedx.size(), &trkdedx.at(0));

			ntrklen_chi2pid_BQ->Fill(norm_trklen, pid);
			Fill1DHist(h1d_chi2pid_BQ, pid);
			Fill1DHist(h1d_ntrklen_BQ, norm_trklen);

			Fill1DHist(h1d_mediandedx_BQ, median_dedx);
			Fill1DHist(h1d_dEdL_BQ, ke_calo_MeV/range_reco);			

			if (IsRecoStop) { //reco_stop
				chi2pid_recostop->Fill(pid);
			} //reco_stop	
			if (IsRecoInel) { //reco_inel
				chi2pid_recoinel->Fill(pid);
			} //reco_inel
		} //if calo size not empty

		//hypothetical length -------------------------------------------------------------------------------------//
		double fitted_length=-1;
                std::reverse(trkdedx.begin(),trkdedx.end());  
		std::reverse(trkres.begin(),trkres.end()); 
		double tmp_fitted_length=BB.Fit_dEdx_Residual_Length(trkdedx, trkres, 2212, false);
		if (tmp_fitted_length>0) fitted_length=tmp_fitted_length;
		double fitted_KE=-50; 
		if (fitted_length>0) { 
			fitted_KE=BB.KEFromRangeSpline(fitted_length);
			//cout<<"event:"<<event<<" evttime:"<<evttime<<" fitted_KE:"<<endl;
		}

		//ke at end point ---------------------------------------------------------------------//
		double kebb=-50; if (fitted_KE>0) kebb=BB.KEAtLength(fitted_KE, range_reco);
		double kecalo=-50; kecalo=fitted_KE-ke_calo_MeV;


		if (IsCaloSize&&IsPandoraSlice) { //CaloSz
			h1d_trklen_CaloSz->Fill(range_reco);
		} //CaloSz

		if (IsPos&&IsCaloSize&&IsPandoraSlice) { //Pos
			h1d_trklen_Pos->Fill(range_reco);
		} //Pos


		//if (IsBQ&&IsCaloSize&&IsPandoraSlice) {
		//if (IsBeamXY&&IsBQ&&IsCaloSize&&IsPandoraSlice) {
		if (IsBeamMom&&IsBeamXY&&IsBQ&&IsCaloSize&&IsPandoraSlice) {
			//if (IsIntersection==false&&IsBeamMom&&IsBeamXY&&IsBQ&&IsCaloSize&&IsPandoraSlice) {
			h1d_kebeam->Fill(ke_beam_MeV);
			h1d_pbeam->Fill(1000.*p_beam);

			h1d_kehy->Fill(fitted_KE);

			h1d_phy->Fill(1000.*ke2p(fitted_KE/1000.));


			h1d_trklen_BQ->Fill(range_reco);
			if (IsXY) { 
				h1d_trklen_XY->Fill(range_reco);
				//h1d_zend_XY->Fill(zend_sce);
			}

			if (IsRecoStop) {
				h1d_zst_noSCE_stop->Fill(zst_nosce);
				h1d_zst_stop->Fill(zst_sce);

				double DT=(evttime-t0)/(60.*60.); //in hour
				h2d_time_zst_stop->Fill(DT,zst_sce);
				h2d_time_zst_noSCE_stop->Fill(DT,zst_nosce);

				h1d_zend_noSCE_stop->Fill(zend_nosce);
				h1d_zend_stop->Fill(zend_sce);

				h1d_kebeam_stop->Fill(ke_beam_MeV);
				h1d_pbeam_stop->Fill(1000.*p_beam);

				h1d_phy_stop->Fill(1000.*ke2p(fitted_KE/1000.));

				h1d_kehy_stop->Fill(fitted_KE);

				h1d_dist_stop->Fill(d_sce);
				h1d_trklen_stop->Fill(range_reco);
				h1d_kerange_stop->Fill(ke_range_MeV);
				h1d_kecalo_stop->Fill(ke_calo_MeV);


				h1d_kend_calo_stop->Fill(kecalo);
				h1d_kend_bb_stop->Fill(kebb);


				double tmp_prange=1000.*ke2p(ke_range);
				double tmp_pcalo=1000.*ke2p(ke_calo_MeV/1000.);
				double tmp_ratio=tmp_pcalo/tmp_prange;

				h1d_prange_stop->Fill(tmp_prange);
				h1d_pcalo_stop->Fill(tmp_pcalo);

				h2d_time_pcalo_stop->Fill(DT,tmp_pcalo);
				h2d_time_prange_stop->Fill(DT,tmp_prange);
				h2d_time_pcaloOverprange_stop->Fill(DT,tmp_ratio);

				if (IsXY) h1d_trklen_stop_XY->Fill(range_reco);

				Fill1DHist(h1d_trklen_ElRich, range_reco);
			}
			if (IsRecoInel) { //reco_inel
				//h1d_trklen_RecoInel->Fill(range_reco);
				Fill1DHist(h1d_trklen_RecoInel,range_reco);

				h1d_kehy_inel->Fill(fitted_KE);

				h1d_kend_calo_inel->Fill(kecalo);
				h1d_kend_bb_inel->Fill(kebb);

			} //reco inel
			if (IsRecoEl) { //reco_el
				h1d_kend_calo_el->Fill(kecalo);
				h1d_kend_bb_el->Fill(kebb);
			} //reco_el
		}

		//misid:p-rich
		if (cosine_beam_spec_primtrk<=0.9&&IsPos&&IsRecoInel&&IsPandoraSlice&&IsCaloSize) { //cos<=0.9+pos+RecoInel
			Fill1DHist(h1d_trklen_MidpRich, range_reco); 
		} //cos<=0.9+pos+RecoInel

		if (IsCaloSize&&IsPandoraSlice) { //calosz
			if (IsRecoStop) {
				Fill1DHist(h1d_trklen_CaloSz_ElRich, range_reco);
				if (IsPos) {
					Fill1DHist(h1d_trklen_Pos_ElRich, range_reco);
				}
			}
		} //calosz
		} //evt loop

		TF1* kebeam_fit; kebeam_fit=VFit(h1d_kebeam, 1); kebeam_fit->SetName("kebeam_fit");
		TF1* pbeam_fit;  pbeam_fit=VFit(h1d_pbeam, 1); 	 pbeam_fit->SetName("pbeam_fit");
		TF1* pbeam_stop_fit;  pbeam_stop_fit=VFit(h1d_pbeam_stop, 1); 	 pbeam_stop_fit->SetName("pbeam_stop_fit");

		//time info ---------------------------------------------------//
		TParameter<Double_t>* T0=new TParameter<Double_t>("T0",0.);
		T0->SetVal(t0);
		TParameter<Double_t>* Tmax=new TParameter<Double_t>("Tmax",0.);
		Tmax->SetVal(tmax);
		TParameter<Double_t>* Tmed=new TParameter<Double_t>("Tmed",0.);
		Tmed->SetVal(t0+0.5*(tmax-t0));


		//save results...
		//TFile *fout = new TFile(Form("data_proton_bmrw.root"),"RECREATE");
		//TFile *fout = new TFile(Form("data_proton_tightxy_bmrw.root"),"RECREATE");
		//TFile *fout = new TFile(Form("data_proton_bmrw_HD.root"),"RECREATE");
		//TFile *fout = new TFile(Form("data_proton_bmrw.root"),"RECREATE");
		//TFile *fout = new TFile(Form("data_proton_beamxy_bmrw.root"),"RECREATE");
		//TFile *fout = new TFile(Form("data_proton_beamxy_beammom_bmrw.root"),"RECREATE");
		//TFile *fout = new TFile(Form("data_proton_beamxy_beammom_rmintersection_bmrw.root"),"RECREATE");
		//TFile *fout = new TFile(Form("data_proton_beamxy_beammom_rmintersection_bmrw2.root"),"RECREATE");
		//TFile *fout = new TFile(Form("/dune/data2/users/hyliao/protonana/v09_39_01/KEstudy/proton_beamxy_beammom_run%d.root",run),"RECREATE");
		TFile *fout = new TFile(Form("/dune/data2/users/hyliao/protonana/v09_39_01/KEHY/proton_beamxy_beammom_run%d_new_stepsz0.105.root",run),"RECREATE");
		//TFile *fout = new TFile(Form("/dune/data2/users/hyliao/protonana/v09_39_01/KEstudy/proton_beamxy_beammom_All.root"),"RECREATE");
		//TFile *fout = new TFile(Form("data_proton_bmrw2_usedefault_range_calc.root"),"RECREATE");
		T0->Write();
		Tmed->Write();
		Tmax->Write();

		h2d_rr_dedx_recoSTOP->Write();
		gr_predict_dedx_resrange->Write();

		h1d_pbeam->Write();
		h1d_pbeam_stop->Write();

		h1d_phy->Write();
		h1d_phy_stop->Write();

		h1d_kebeam->Write();
		h1d_kebeam_stop->Write();

		h1d_kehy->Write();
		h1d_kehy_stop->Write();
		h1d_kehy_inel->Write();


		h1d_prange_stop->Write();
		h1d_pcalo_stop->Write();

		kebeam_fit->Write();
		pbeam_fit->Write();
		pbeam_stop_fit->Write();


		h1d_kend_calo_stop->Write();
		h1d_kend_calo_el->Write();
		h1d_kend_calo_inel->Write();

		h1d_kend_bb_stop->Write();
		h1d_kend_bb_el->Write();
		h1d_kend_bb_inel->Write();


		h1d_kerange_stop->Write();
		h1d_kecalo_stop->Write();

		h1d_trklen_CaloSz->Write();
		h1d_trklen_Pos->Write();
		h1d_trklen_BQ->Write();
		h1d_trklen_RecoInel->Write();
		h1d_trklen_MidpRich->Write();
		h1d_trklen_ElRich->Write();
		h1d_trklen_CaloSz_ElRich->Write();
		h1d_trklen_Pos_ElRich->Write();

		h1d_dist_stop->Write();
		h1d_trklen_XY->Write();
		h1d_trklen_stop->Write();
		h1d_trklen_stop_XY->Write();


		h2d_xy_noSCE->Write();
		h2d_xy_SCE->Write();

		h1d_zst_noSCE_stop->Write();
		h1d_zst_stop->Write();
		h2d_time_zst_stop->Write();
		h2d_time_zst_noSCE_stop->Write();

		h1d_zend_noSCE_stop->Write();
		h1d_zend_stop->Write();
		//h1d_zend_XY->Write();

		chi2pid_recostop->Write();
		chi2pid_recoinel->Write();

		ntrklen_chi2pid_BQ->Write();
		h1d_chi2pid_BQ->Write();
		h1d_ntrklen_BQ->Write();

		h1d_mediandedx_BQ->Write();
		h1d_dEdL_BQ->Write();

		h2d_time_pcalo_stop->Write();
		h2d_time_prange_stop->Write();
		h2d_time_pcaloOverprange_stop->Write();

		fout->Close();

	}
