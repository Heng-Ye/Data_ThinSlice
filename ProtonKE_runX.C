#define ProtonKE_run5387_cxx
#include "ProtonKE_run5387.h"

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

#include "./cali/dedx_function_r5387.h"
#include "./headers/BasicParameters.h"
#include "./headers/BasicAnaFunc.h"
#include "./headers/util.h"
#include "./headers/BetheBloch.h"

#include <cassert>

using namespace std;
using namespace ROOT::Math;

void ProtonKE_run5387::Loop() {
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	//book histograms -------------------------------------------------------------------------//
	int nke=100;
	double kemin=0;
	double kemax=600;

	//beam
	TH1D *h1d_kebeam=new TH1D("h1d_kebeam","", nke, kemin, kemax); //ke from beamline inst.
	TH1D *h1d_kebeam_stop=new TH1D("h1d_kebeam_stop","", nke, kemin, kemax); //ke from beamline inst. (stopping protons)

	//ke of stopping protons
	TH1D *h1d_kerange_stop=new TH1D("h1d_kerange_stop","", nke, kemin, kemax); //range-based calc. (stopping protons)
	TH1D *h1d_kecalo_stop=new TH1D("h1d_kecalo_stop","", nke, kemin, kemax); //calorimetric-based calc. (stopping protons)
	TH1D *h1d_kecalo_recoinel=new TH1D("h1d_kecalo_recoinel","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)
	TH1D *h1d_kecalo_recoel=new TH1D("h1d_kecalo_recoel","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)

	//upstream energy loss
	//int ndke=300;
	//double dkemin=0;
	//double dkemax=600;

	//int ndke=800;
	int ndke=320;
	double dkemin=-800;
	double dkemax=800;

	//reco_label
	TH1D *h1d_dKE=new TH1D("h1d_dKE","", ndke, dkemin, dkemax); //ke_beam-E_dept
	TH1D *h1d_dKE_stop=new TH1D("h1d_dKE_stop","", ndke, dkemin, dkemax); //ke_beam-E_dept
	TH1D *h1d_dKE_recoinel=new TH1D("h1d_dKE_recoinel","", ndke, dkemin, dkemax); //ke_beam-E_dept
	TH1D *h1d_dKE_recoel=new TH1D("h1d_dKE_recoel","", ndke, dkemin, dkemax); //ke_beam-E_dept

	TH1D *h1d_KEbb_reco=new TH1D("h1d_KEbb_reco","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_RecoInel=new TH1D("h1d_KEbb_reco_RecoInel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_RecoEl=new TH1D("h1d_KEbb_reco_RecoEl","", ndke, dkemin, dkemax);

	TH1D *h1d_KEcalo=new TH1D("h1d_KEcalo","", ndke, dkemin, dkemax); //calorimetric-based calc. (all reco. protons)
	TH1D *h1d_KEcalo_RecoInel=new TH1D("h1d_KEcalo_RecoInel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEcalo_RecoEl=new TH1D("h1d_KEcalo_RecoEl","", ndke, dkemin, dkemax);

	//KE vs trklen
	TH2D *h2d_trklen_KEcalo=new TH2D("h2d_trklen_KEcalo","", 140, 0,140, ndke, dkemin, dkemax);
	TH2D *h2d_trklen_KEcalo_RecoInel=new TH2D("h2d_trklen_KEcalo_RecoInel","", 140, 0,140, ndke, dkemin, dkemax);
	TH2D *h2d_trklen_KEcalo_RecoEl=new TH2D("h2d_trklen_KEcalo_RecoEl","", 140, 0,140, ndke, dkemin, dkemax);
	//------------------------------------------------------------------------------------------------------------------------//

	//Energy loss -------------------------------//
	double mean_kecalo_stop=3.78066e+02;
	double err_mean_kecalo_stop=9.81156e-01;
	double sigma_kecalo_stop=5.60675e+01;
	double err_sigma_kecalo_stop=8.12221e-01;

	double mean_kerange_stop=4.00446e+02;
	double err_mean_kerange_stop=1.07727e+00;
	double sigma_kerange_stop=5.12091e+01;
	double err_sigma_kerange_stop=1.06785e+00;

	double mean_kebeam=4.41392e+02;
	double err_mean_kebeam=5.95236e-01;
	double sigma_kebeam=5.15066e+01;
	double err_sigma_kebeam=4.28244e-01;

	double mean_Elossrange_stop=mean_kebeam-mean_kerange_stop;
	double mean_Elosscalo_stop=mean_kebeam-mean_kecalo_stop;

	//-------------------------------------------//

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

		//Basic configure ------//
		BetheBloch BB;
		BB.SetPdgCode(2212);
		//----------------------//

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
		bool IsRecoInEL=false;
		bool IsRecoEL=false;
		if (norm_trklen>min_norm_trklen_csda) IsRecoStop=true; //stopping proton cut
		//if (norm_trklen<=min_norm_trklen_csda) IsRecoInel=true; //inelastic proton cut
		//if (IsRecoInel&&IsBQ&&IsCaloSize&&IsPandoraSlice) n_recoinel++;

                if ((range_reco/csda_val)<min_norm_trklen_csda) { //inel region
                        if (pid>pid_1) IsRecoInEL=true;
                        if (pid<=pid_1) IsRecoEL=true;
                } //inel region
                if ((range_reco/csda_val)>=min_norm_trklen_csda&&(range_reco/csda_val)<max_norm_trklen_csda) { //stopping p region
                        if (pid>pid_2) IsRecoInEL=true;
                        if (pid<=pid_2) IsRecoEL=true;
                } //stopping p region

		//XY Cut ---------------------------------------------------------------------------------------------------------------//
		bool IsXY=false;
		//start(x,y,z) without SCE corr. 
		double xst_nosce=0;
		double yst_nosce=0;
		double zst_nosce=0;
		if (IsCaloSize&&IsPandoraSlice) { //if calo size not empty
			if ((primtrk_startz->at(-1+primtrk_startz->size()))>(primtrk_startz->at(0))) { //check if Pandora flip the sign
				xst_nosce=primtrk_startx->at(0);
				yst_nosce=primtrk_starty->at(0);
				zst_nosce=primtrk_startz->at(0);
			} //check if Pandora flip the sign
			else {
				xst_nosce=primtrk_startx->at(-1+primtrk_startx->size());
				yst_nosce=primtrk_starty->at(-1+primtrk_starty->size());
				zst_nosce=primtrk_startz->at(-1+primtrk_startz->size());
			}
		} //if calo size not empty
		if ((pow(((xst_nosce-mean_x)/dev_x),2)+pow(((yst_nosce-mean_y)/dev_y),2))<=1.) IsXY=true;


		//-----------------------------------------------------------------------//

		//calo info
		double ke_calo_MeV=0;
		double xst_sce=0;
		double yst_sce=0;
		double zst_sce=0;
		double xend_sce=0;
		double yend_sce=0;
		double zend_sce=0;
		if (IsBQ&&IsCaloSize&&IsPandoraSlice) { //if calo size not empty
			//h2d_xy_noSCE->Fill(xst_nosce, yst_nosce);
			//h1d_zst_noSCE->Fill(zst_nosce);

			//Fill1DHist(h1d_KEbb_reco, KEbb_reco2);
			//Fill1DHist(h1d_KEcalo, keff_reco4-ke_calo_MeV);

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
			//h2d_xy_SCE->Fill(xst_sce,yst_sce);
			//h1d_zst_SCE->Fill(zst_sce);

			//calo info --------------------------------------------------------------------------//
			vector<double> trkdedx; //dedx 
			vector<double> trkres; //rr
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

				//if (IsRecoStop) { //reco_stop
					//h2d_rr_dedx_recoSTOP->Fill(resrange_reco, dedx);
				//} //reco_stop	
			} //loop over hits

			double pid=-99; pid=chi2pid(trkdedx,trkres);
			double median_dedx=-99; median_dedx=TMath::Median(trkdedx.size(), &trkdedx.at(0));

			//ntrklen_chi2pid_BQ->Fill(norm_trklen, pid);
                        //Fill1DHist(h1d_chi2pid_BQ, pid);
                        //Fill1DHist(h1d_ntrklen_BQ, norm_trklen);

			//Fill1DHist(h1d_mediandedx_BQ, median_dedx);
			//Fill1DHist(h1d_dEdL_BQ, ke_calo_MeV/range_reco);			

			//if (IsRecoStop) { //reco_stop
				//chi2pid_recostop->Fill(pid);
			//} //reco_stop	
			//if (IsRecoInel) { //reco_inel
				//chi2pid_recoinel->Fill(pid);
			//} //reco_inel
		} //if calo size not empty

		//KEff, KEend -----------------------------------------------------------------------------------------//
		double KEcsda=1000.*ke_vs_csda_range_sm->Eval(range_reco);
		double KEcsda_inel=ke_beam_MeV-mean_Elossrange_stop;

		double KEbb_reco_constErange=-1; KEbb_reco_constErange=BB.KEAtLength(KEcsda_inel, range_reco);
		double KEbb_reco_Edept_range=-1; KEbb_reco_Edept_range=BB.KEAtLength(KEcsda, range_reco);

		double KEcalo_reco_constEcalo=-1; KEcalo_reco_constEcalo=ke_beam_MeV-mean_Elosscalo_stop-ke_calo_MeV;

		//if (IsCaloSize&&IsPandoraSlice) { //CaloSz
			//h1d_trklen_CaloSz->Fill(range_reco);
		//} //CaloSz

		//if (IsPos&&IsCaloSize&&IsPandoraSlice) { //Pos
			//h1d_trklen_Pos->Fill(range_reco);
		//} //Pos


		if (IsBQ&&IsCaloSize&&IsPandoraSlice) {
			//h1d_pbeam->Fill(1000.*p_beam);		
			h1d_kebeam->Fill(ke_beam_MeV);
			h1d_dKE->Fill(ke_beam_MeV-ke_calo_MeV);
			Fill1DHist(h1d_KEcalo, KEcalo_reco_constEcalo);	
			h2d_trklen_KEcalo->Fill(range_reco, KEcalo_reco_constEcalo);

			//h1d_trklen_BQ->Fill(range_reco);
			//h1d_zend->Fill(zend_sce);

			//if (IsXY) { 
				//h1d_trklen_XY->Fill(range_reco);
				//h1d_zend_XY->Fill(zend_sce);
			//}

			if (IsRecoStop) {
				h1d_kebeam_stop->Fill(ke_beam_MeV);
				//h1d_trklen_stop->Fill(range_reco);
				h1d_kerange_stop->Fill(ke_range_MeV);
				h1d_kecalo_stop->Fill(ke_calo_MeV);
				h1d_dKE_stop->Fill(ke_beam_MeV-ke_calo_MeV);

				//h1d_pbeam_stop->Fill(1000.*ke2p(ke_range));
				//if (IsXY) h1d_trklen_stop_XY->Fill(range_reco);
				//Fill1DHist(h1d_trklen_ElRich, range_reco);
			}

			if (IsRecoInEL) { //reco_inel
				h1d_kecalo_recoinel->Fill(ke_calo_MeV);
				h1d_dKE_recoinel->Fill(ke_beam_MeV-ke_calo_MeV);

				Fill1DHist(h1d_KEbb_reco, KEbb_reco_constErange);
				Fill1DHist(h1d_KEbb_reco_RecoInel, KEbb_reco_constErange);

				Fill1DHist(h1d_KEcalo_RecoInel, KEcalo_reco_constEcalo);	
				h2d_trklen_KEcalo_RecoInel->Fill(range_reco, KEcalo_reco_constEcalo);

				//Fill1DHist(h1d_KEcalo_RecoInel, keff_reco4-ke_calo_MeV);
				//h1d_trklen_RecoInel->Fill(range_reco);
				//Fill1DHist(h1d_trklen_RecoInel,range_reco);
			} //reco inel
			if (IsRecoEL) { //reco_el
				h1d_kecalo_recoel->Fill(ke_calo_MeV);
				h1d_dKE_recoel->Fill(ke_beam_MeV-ke_calo_MeV);

				Fill1DHist(h1d_KEbb_reco, KEbb_reco_Edept_range);
				Fill1DHist(h1d_KEbb_reco_RecoEl, KEbb_reco_Edept_range);

				Fill1DHist(h1d_KEcalo_RecoEl, KEcalo_reco_constEcalo);
				h2d_trklen_KEcalo_RecoEl->Fill(range_reco, KEcalo_reco_constEcalo);

				//Fill1DHist(h1d_KEcalo_RecoEl, keff_reco4-ke_calo_MeV);
			} //reco_el
		}

		//misid:p-rich
		//if (cosine_beam_spec_primtrk<=0.9&&IsPos&&IsRecoInel&&IsPandoraSlice&&IsCaloSize) { //cos<=0.9+pos+RecoInel
				//Fill1DHist(h1d_trklen_MidpRich, range_reco); 
		//} //cos<=0.9+pos+RecoInel

		//if (IsCaloSize&&IsPandoraSlice) { //calosz
			//if (IsRecoStop) {
				//Fill1DHist(h1d_trklen_CaloSz_ElRich, range_reco);
				//if (IsPos) {
					//Fill1DHist(h1d_trklen_Pos_ElRich, range_reco);
				//}
			//}
		//} //calosz
	} //evt loop

	TF1* kebeam_fit; kebeam_fit=VFit(h1d_kebeam, 1); kebeam_fit->SetName("kebeam_fit");
	//TF1* pbeam_fit;  pbeam_fit=VFit(h1d_pbeam, 1); 	 pbeam_fit->SetName("pbeam_fit");
	//TF1* pbeam_stop_fit;  pbeam_stop_fit=VFit(h1d_pbeam_stop, 1); 	 pbeam_stop_fit->SetName("pbeam_stop_fit");

	TF1* kebeam_stop_fit; kebeam_stop_fit=VFit(h1d_kebeam_stop, 2);
	kebeam_stop_fit->SetName("kebeam_stop_fit");

	TF1* kecalo_stop_fit; kecalo_stop_fit=VFit(h1d_kecalo_stop, 2);
	kecalo_stop_fit->SetName("kecalo_stop_fit");

	TF1* kerange_stop_fit; kerange_stop_fit=VFit(h1d_kerange_stop, 2);
	kerange_stop_fit->SetName("kerange_stop_fit");

	//save results...
	TFile *fout = new TFile(Form("data_ke_AdvancedKEff.root"),"RECREATE");
		h1d_kebeam->Write();
		h1d_kebeam_stop->Write();

		h1d_kerange_stop->Write();
		h1d_kecalo_stop->Write();

		h1d_kecalo_recoinel->Write();
		h1d_kecalo_recoel->Write();


		kebeam_fit->Write();
		kebeam_stop_fit->Write();
		kecalo_stop_fit->Write();
		kerange_stop_fit->Write();

		h1d_dKE->Write();
		h1d_dKE_stop->Write();
		h1d_dKE_recoinel->Write();
		h1d_dKE_recoel->Write();

		h1d_KEbb_reco->Write();
		h1d_KEbb_reco_RecoInel->Write();
		h1d_KEbb_reco_RecoEl->Write();

		h1d_KEcalo->Write();
		h1d_KEcalo_RecoInel->Write();
		h1d_KEcalo_RecoEl->Write();

		h2d_trklen_KEcalo->Write();
		h2d_trklen_KEcalo_RecoInel->Write();
		h2d_trklen_KEcalo_RecoEl->Write();

		//h1d_KEcalo->Write();
		//h1d_KEcalo_RecoInel->Write();
		//h1d_KEcalo_RecoEl->Write();



/*
	h2d_rr_dedx_recoSTOP->Write();
	gr_predict_dedx_resrange->Write();

	h1d_pbeam->Write();
	h1d_kebeam->Write();
	h1d_pbeam_stop->Write();

	kebeam_fit->Write();
	pbeam_fit->Write();
	pbeam_stop_fit->Write();

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

	h1d_trklen_XY->Write();
	h1d_trklen_stop->Write();
	h1d_trklen_stop_XY->Write();


	h2d_xy_noSCE->Write();
	h2d_xy_SCE->Write();

	h1d_zst_noSCE->Write();
	h1d_zst_SCE->Write();

	h1d_zend->Write();
	h1d_zend_XY->Write();

	chi2pid_recostop->Write();
	chi2pid_recoinel->Write();

	ntrklen_chi2pid_BQ->Write();
        h1d_chi2pid_BQ->Write();
        h1d_ntrklen_BQ->Write();

	h1d_mediandedx_BQ->Write();
	h1d_dEdL_BQ->Write();
*/

		
	fout->Close();

}
