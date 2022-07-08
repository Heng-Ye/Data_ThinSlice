#define ProtonEvtDisplay_run5387_cxx
#include "ProtonEvtDisplay_run5387.h"

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
#include "./headers/SliceParams.h"
#include "./headers/util.h"
#include "./headers/BackgroundStudy.h"
//#include "./headers/ThinSlice.h"

#include <cassert>

using namespace std;
using namespace ROOT::Math;

void ProtonEvtDisplay_run5387::Loop() {
	if (fChain == 0) return;

	//ThinSlice config. ---------------------------------------------------------------------------------------------------//
        //SetOutputFileName(Form("data_run5387_prod4reco2_bkgstudy_dx%dcm_%dslcs.root", name_thinslicewidth, nthinslices)); //output file name
        //SetOutputFileName(Form("data_run5387_prod4reco2_bkgstudy_dx%dcm_%dslcs_largerbin.root", name_thinslicewidth, nthinslices)); //output file name

        //book histograms --//
        BookHistograms();

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	//histograms

	//countings ------------------//
	Long64_t nbytes = 0, nb = 0;
	int n_tot=0, n_pan=0, n_calosz=0, n_bq=0, n_recoinel=0;
	int n_pan_shower=0, n_pan_other=0;
	int n_beam=0, n_beam_inst_tof=0, n_beam_inst_valid=0, n_beam_mom=0;
	int n_can=0, n_proton=0;

	for (Long64_t jentry=0; jentry<nentries;jentry++) { //evt loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

		size_t j=0;
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
		n_beam++;

		//Incoming candidate event ------------------------------------------------------------------------//	
		bool IsProton=false;
		if (!beam_inst_PDG_candidates->empty()&&IsBeamInstValid) {
			n_can++;
			for (int h=0; h<(int)beam_inst_PDG_candidates->size(); ++h) { 
				if (beam_inst_PDG_candidates->at(h)==2212) {
					n_proton++;
					IsProton=true;
				}
			}
		}
		if (!IsProton) continue; //only protons, not other candidate particles
		//(the ntuple actually only save proton evt for the beam track)

		//Pandora slice cut ------------------------------------------------------------------//
                bool IsPandoraSlice=false; //pandora slice
		bool IsCaloSize=false; //if calo size not empty

		//cout<<"isprimarytrack:"<<isprimarytrack<<" isprimaryshower:"<<isprimaryshower<<endl;
		if (isprimarytrack==1&&isprimaryshower==0) { //pandora slice cut
			IsPandoraSlice=true;
			n_pan++;
		}

		if (!primtrk_hitz->empty()) { 
			IsCaloSize=true;
			if (IsPandoraSlice) n_calosz++;
		}

		//other beam info -----------------------------------------------------------------------------------//
/*
		if (!beam_inst_TOF->empty()) {
			for (int h=0; h<(int)beam_inst_TOF->size(); ++h) { 
				cout<<"beam_inst_TOF["<<h<<"]:"<<beam_inst_TOF->at(h)<<" | "
				    <<"beam_inst_TOF_Chan["<<h<<"]:"<<beam_inst_TOF_Chan->at(h)<<endl;
			}
		}

		if (!tofs->empty()) {
			for (int h=0; h<(int)tofs->size(); ++h) { 
				cout<<"tofs["<<h<<"]:"<<tofs->at(h)<<" | ch_tofs["<<h<<"]:"<<ch_tofs->at(h)<<endl;
			}
		}
*/

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


                //if (reco_stz>min1_z&&reco_stz<min2_z) {
                	//if (reco_sty>min1_y&&reco_sty<min2_y) {
                        	//if (reco_stx>min1_x&&reco_stx<min2_x) {
                                	//IsPos=true;
                                //}
                        //}
                //}
                
                //cosine_theta/cut ---------------------------------------------------------------------------------//
                bool IsCosine=false;
                double cosine_beam_spec_primtrk=-999;
                //cosine_beam_spec_primtrk=cosine_beam_primtrk; //cosine between beam_spec and primary trk direction(trk before sce corr.)
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
		if (IsBQ&&IsCaloSize&&IsPandoraSlice) n_bq++;

		//if (IsPos&&IsCaloSize&&IsPandoraSlice) Fill1DHist(reco_cosineTheta_Pos, cosine_beam_spec_primtrk);
		


		//Get calo info ---------------------------------------------------------------------------------//
                int index_reco_endz=0;
                double wid_reco_max=-9999;
                double range_reco=-99;
                vector<double> reco_trklen_accum;
                reco_trklen_accum.reserve(primtrk_hitz->size());
                double kereco_calo=0;
                double kereco_range=0;
                double kereco_range2=0;
                vector<double> EDept;
		double pid=-99;
		//Fill1DHist(trklen_reco_NoCut, range_reco);
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

                        	//if (wid_reco==-9999) continue; //outside TPC
                        	if (wid_reco>wid_reco_max) {
                                	wid_reco_max=wid_reco;
                                	index_reco_endz=(int)-1+primtrk_wid->size()-h;
                        	}

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

                        	kereco_calo+=dedx*pitch;
                        	kereco_range+=pitch*dedx_predict(resrange_reco);
                        	kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);

                                trkdedx.push_back(dedx);
                                trkres.push_back(resrange_reco);

                  	} //loop over reco hits in a track

                        pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis

		} //if calo size not empty	

		//if (IsCaloSize) range_reco=primtrk_range->at(0); //unit:cm

		//beam info ------------------------------------------------------------------------------------------------------------------------------------//
		double bx=beamPosx->at(0);
		double by=beamPosy->at(0);
		double bz=beamPosz->at(0);
		double p_beam=beamMomentum->at(0);
		double ke_beam_MeV=1000.*p2ke(p_beam); //unit:MeV

		double csda_val=csda_range_vs_mom_sm->Eval(p_beam); //expected csda value if proton stops; unit: cm
		double ke_range_MeV=1000.*ke_vs_csda_range_sm->Eval(range_reco); //unit:MeV
		double p_range_MeV=1000.*ke2p(ke_range_MeV/1000.); //MeV/c
		double norm_trklen=range_reco/csda_val; //trklen/csda_val
		//cout<<"csda_val:"<<csda_val<<" range_reco:"<<range_reco<<" norm_trklen:"<<norm_trklen<<" min_norm_trklen_csda:"<<min_norm_trklen_csda<<endl;

		//reco Inel cut ------------------------------------------------------------------------------------//
		bool IsRecoStop=false;
		bool IsRecoInel=false;
		//if (norm_trklen>min_norm_trklen_csda) IsRecoStop=true; //old stopping proton cut
		//if (norm_trklen<=min_norm_trklen_csda) IsRecoInel=true; //old inelastic proton cut

                if ((range_reco/csda_val)<min_norm_trklen_csda) { //inel region
                        if (pid>pid_1) IsRecoInel=true;
                        if (pid<=pid_1) IsRecoStop=true;
                } //inel region
                if ((range_reco/csda_val)>=min_norm_trklen_csda&&(range_reco/csda_val)<max_norm_trklen_csda) { //stopping p region
                        if (pid>pid_2) IsRecoInel=true;
                        if (pid<=pid_2) IsRecoStop=true;
                } //stopping p region
		

                //Slice ID definition ----------------------------------------------//
                //reco slice ID
                int reco_sliceID=-1;
                if (IsPandoraSlice&&IsCaloSize) {
                        reco_sliceID = int(range_reco/thinslicewidth);
                        if (reco_sliceID < 0) reco_sliceID = -1;
                        if (reco_endz<0) reco_sliceID = -1;
                        if (reco_sliceID >= nthinslices) reco_sliceID = nthinslices;
                }
                //------------------------------------------------------------------//

		//[1]misID:p rich sample ---------------------------------------------//
                //if (IsPos&&IsPandoraSlice&&IsCaloSize) { //if Pos
                //if (IsBQ&&IsPandoraSlice&&IsCaloSize) { //if Pos
                if (IsRecoStop&&IsBQ&&IsPandoraSlice&&IsCaloSize) { //if Pos
			//if (range_reco<=22&&range_reco>=9) { //only care about slice id<=5
			if (range_reco>=100) { //only care about slice id<=5
				cout<<"size of beampos:"<<beamPosz->size()<<endl;
				cout<<"range_reco:"<<range_reco<<endl;
				for (size_t jb=0; jb<beamPosz->size(); ++jb ) {
					cout<<"beamZ["<<jb<<"]:"<<beamPosz->at(jb)<<endl;
				}

				TString tit; tit=Form("run%d subrun:%d evt:%d|beamSz:%zu", run, subrun, event, beamPosz->size() );
				TString fig_out;
				cout<<tit<<endl;
				//fig_out=Form("/dune/data2/users/hyliao/evt_display/proton_cosTheta/sliceIDlt5/data/run%d_event%d_evt%d_trklen%.2fcm_costheta%.2f.png",run,subrun,event,range_reco,cosine_beam_spec_primtrk);
				//fig_out=Form("/dune/data2/users/hyliao/evt_display/trklen922/data/run%d_event%d_evt%d_trklen%.2fcm_costheta%.2f.png",run,subrun,event,range_reco,cosine_beam_spec_primtrk);
				fig_out=Form("/dune/data2/users/hyliao/evt_display/trklen100/data/run%d_event%d_evt%d_trklen%.2fcm_costheta%.2f.png",run,subrun,event,range_reco,cosine_beam_spec_primtrk);
				
				gStyle->SetOptStat(0); //no statistics box
				gStyle->SetFrameBorderSize(3);

				TCanvas *c1 = new TCanvas(Form("c1"),tit.Data(),1600,800);
				//cnt_canvas++;
				//c1->SetName(Form("c_%d",cnt_canvas));
				c1->SetTitle(tit.Data());
				c1->Divide(2,1);

				//draw track in y-z view
				//after entering TPC
				TGraph *gr_trk_yz=new TGraph(primtrk_hitz->size(), &primtrk_hitz->at(0), &primtrk_hity->at(0));
				TGraph *gr_trk_xz=new TGraph(primtrk_hitz->size(), &primtrk_hitz->at(0), &primtrk_hitx->at(0));
				gr_trk_yz->SetMarkerStyle(20);
				gr_trk_xz->SetMarkerStyle(20);
				gr_trk_yz->SetMarkerColor(2);
				gr_trk_xz->SetMarkerColor(2);

				//after entering TPC, trkend
				TGraph *gr_trkend_yz=new TGraph(1, &primtrk_hitz->at(-1+primtrk_hitz->size()), &primtrk_hity->at(-1+primtrk_hitz->size()));
				TGraph *gr_trkend_xz=new TGraph(1, &primtrk_hitz->at(-1+primtrk_hitz->size()), &primtrk_hitx->at(-1+primtrk_hitz->size()));
				gr_trkend_yz->SetMarkerStyle(20);
				gr_trkend_xz->SetMarkerStyle(20);
				gr_trkend_yz->SetMarkerColor(3);
				gr_trkend_xz->SetMarkerColor(3);


				//before tpc
				TGraph *gr_btrk_yz=new TGraph(beamPosz->size(), &beamPosz->at(0), &beamPosy->at(0));
				TGraph *gr_btrk_xz=new TGraph(beamPosz->size(), &beamPosz->at(0), &beamPosx->at(0));
				gr_btrk_yz->SetMarkerStyle(24);
				gr_btrk_xz->SetMarkerStyle(24);
				gr_btrk_yz->SetMarkerColor(4);
				gr_btrk_xz->SetMarkerColor(4);

				TH2D *h2d_yz=new TH2D("h2d_yz","",170,-50,120,120,350,450); //bkg
				c1->cd(1);
				h2d_yz->GetXaxis()->SetTitle("Z [cm]");
				h2d_yz->GetYaxis()->SetTitle("Y [cm]");
				h2d_yz->SetTitle(tit.Data());
				h2d_yz->Draw();
				gr_trk_yz->Draw("p same");
				gr_btrk_yz->Draw("p same");
				gr_trkend_yz->Draw("p same");

				TH2D *h2d_xz=new TH2D("h2d_xz","",170,-50,120,120,-120,-10);
				c1->cd(2);
				h2d_xz->GetXaxis()->SetTitle("Z [cm]");
				h2d_xz->GetYaxis()->SetTitle("X [cm]");
				TString str_right;
					str_right=Form("reco_len:%.1f cm| SID:%d| cos#Theta:%.2f| #chi^{2}:%.2f| ntrklen:%.2f", range_reco, reco_sliceID, cosine_beam_spec_primtrk,pid,range_reco/csda_val);
				h2d_xz->SetTitle(str_right.Data());
				h2d_xz->Draw();
				gr_trk_xz->Draw("p same");
				gr_btrk_xz->Draw("p same");
				gr_trkend_xz->Draw("p same");


				c1->Modified(); c1->Update();
				c1->Print(fig_out.Data());


			} //only care about slice id<=5             

                } //if Pos
		//--------------------------------------------------------------------//

		//[2]El rich sample
		if (IsBQ&&IsCaloSize&&IsPandoraSlice) { //if calo size not empty
			for (int j=0; j<nthinslices+2; ++j) { //slice loop
				if (j==(reco_sliceID+1)) {
					//h_chi2_BQ[j]->Fill(pid);
				}
			} //slice loop
		}



	} //evt loop

	//save results ---------//
	//SaveHistograms();







}
