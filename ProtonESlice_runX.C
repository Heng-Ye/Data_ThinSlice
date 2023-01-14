#define ProtonESlice_run5387_cxx
#include "ProtonESlice_run5387.h"

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
#include "./headers/ESliceParams.h"
#include "./headers/ESlice.h"
#include "./headers/BetheBloch.h"

#include <cassert>

using namespace std;
using namespace ROOT::Math;

void ProtonESlice_run5387::Loop() {
	if (fChain == 0) return;

	//ThinSlice config. --------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        //SetOutputFileName(Form("data_run5387_prod4a_thinslice_dx%dcm_%dslcs.root", name_thinslicewidth, nthinslices)); //output file name
        //SetOutputFileName(Form("data_run5387_prod4a_eslice_dx%dcm_%dslcs_stid+0.5.root", name_thinslicewidth, nthinslices)); //output file name
        //SetOutputFileName(Form("/dune/data2/users/hyliao/protonana/v09_39_01/XS/prod4a_Eslice_dE%dMeV_%dslcs_beamxy_run%d_v09_39_01.root", name_thinslicewidth, nthinslices,5387)); //output file name
        SetOutputFileName(Form("/dune/data2/users/hyliao/protonana/v09_39_01/XS_newslcid/prod4a_Eslice_dE%dMeV_%dslcs_beamxy_run%d_v09_39_01_newslcid.root", name_thinslicewidth, nthinslices,5387)); //output file name

	//Basic configure ------//
	BetheBloch BB;
	BB.SetPdgCode(2212);
	//----------------------//
	
	//Fitted Data Beam Momentum ---------------//
	double m1=1013.71; //Data prod4 reco2
	double s1=68.6327; //Data prod4 reco2

	//momentum cut range	
	double mu_min=m1-3.*s1;
	double mu_max=m1+3.*s1;


        //const. E-loss assumption -----------------------------------------------//     
        double const_eloss_data=45.6084/0.99943; //const E-loss from fit (calo)
	//p[0]:45.6084; err_p[0]:0.296889; p[1]:-0.99943 err_p[1]:0.00617258
        
	//hy -------------------------------------------//
	double Eloss_data_hy_stop=25.1911/1.00063;
	double R_fit_hy=0.9994965057702763;
	double er_R_fit_hy=0.04516809864770481;

	//p[0]:25.1911
	//err_p[0]:0.114666
	//p[1]:-1.00063
	//err_p[1]:0.00439028

        //book histograms --//
        BookHistograms();

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	//countings ------------------//
	Long64_t nbytes = 0, nb = 0;
	int reco_sliceID = -1;
	int reco_st_sliceID = -1;
	for (Long64_t jentry=0; jentry<nentries;jentry++) { //evt loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

		reco_sliceID = -1;
		reco_st_sliceID = -1;

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

		//cout<<"beam_inst_trigger:"<<beam_inst_trigger<<endl;
		//cout<<"\n\nrun/subrun/event:"<<run<<"/"<<subrun<<"/"<<event<<endl;
		//cout<<"IsBeamInstValid:"<<IsBeamInstValid<<endl;
		//cout<<"beam_inst_nMomenta:"<<beam_inst_nMomenta<<endl;
		//cout<<"beam_inst_nTracks:"<<beam_inst_nTracks<<endl;
		//cout<<"IsBeamOK:"<<IsBeamOK<<endl;
		//cout<<"hp:"<<high_pressure_status<<" lp:"<<low_pressure_status<<endl;
		if (!IsBeamOK) continue; //skip evt if beam not okay
		//n_beam++;

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

		//cout<<"isprimarytrack:"<<isprimarytrack<<" isprimaryshower:"<<isprimaryshower<<endl;
		if (isprimarytrack==1&&isprimaryshower==0) { //pandora slice cut
			IsPandoraSlice=true;
			//n_pan++;
		}

		if (!primtrk_hitz->empty()) { 
			IsCaloSize=true;
			//if (IsPandoraSlice) n_calosz++;
		}
		//cout<<"IsProton:"<<IsProton<<" IsPandoraSlice:"<<IsPandoraSlice<<endl;

		//other beam info -----------------------------------------------------------------------------------//
		//if (!beam_inst_TOF->empty()) {
			//for (int h=0; h<(int)beam_inst_TOF->size(); ++h) { 
				//cout<<"beam_inst_TOF["<<h<<"]:"<<beam_inst_TOF->at(h)<<" | "
				   // <<"beam_inst_TOF_Chan["<<h<<"]:"<<beam_inst_TOF_Chan->at(h)<<endl;
			//}
		//}

		//if (!tofs->empty()) {
			//for (int h=0; h<(int)tofs->size(); ++h) { 
				//cout<<"tofs["<<h<<"]:"<<tofs->at(h)<<" | ch_tofs["<<h<<"]:"<<ch_tofs->at(h)<<endl;
			//}
		//}

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

			//Fill1DHist(reco_startX_sce, reco_stx);
			//Fill1DHist(reco_startY_sce, reco_sty);
			//Fill1DHist(reco_startZ_sce, reco_stz);

                        double beam_dx=(reco_stx-mean_StartX)/sigma_StartX;
                        double beam_dy=(reco_sty-mean_StartY)/sigma_StartY;
                        double beam_dz=(reco_stz-mean_StartZ)/sigma_StartZ;
                        double beam_dxy=sqrt(pow(beam_dx,2)+pow(beam_dy,2));

			//Fill1DHist(hdeltaX, beam_dx);
			//Fill1DHist(hdeltaY, beam_dy);
			//Fill1DHist(hdeltaZ, beam_dz);
			//Fill1DHist(hdeltaXY, beam_dxy);

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
		//if (IsCaloSize) Fill1DHist(reco_cosineTheta, cosine_beam_spec_primtrk);

                //beam quality cut -----------------------------//
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


                bool IsBQ=false;
                if (IsBeamXY&&IsBeamMom&&IsCosine&&IsPos) IsBQ=true;
		//if (IsBQ&&IsCaloSize&&IsPandoraSlice) n_bq++;
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
                vector<double> trkdedx;
                vector<double> trkres;

		//Fill1DHist(trklen_reco_NoCut, range_reco);
		double reco_calo_MeV=0;
                if (IsCaloSize) { //if calo size not empty
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

                        	//kereco_calo+=dedx*pitch;
                        	//kereco_range+=pitch*dedx_predict(resrange_reco);
                        	//kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);
				reco_calo_MeV+=dedx*pitch;

                                trkdedx.push_back(dedx);
                                trkres.push_back(resrange_reco);

                  	} //loop over reco hits in a track

                        pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis

		} //if calo size not empty
			

		//if (IsCaloSize) range_reco=primtrk_range->at(0); //unit:cm

		//beam info ------------------------------------------------------------------------------------------------------------------------------------//
		//double bx=beamPosx->at(0);
		//double by=beamPosy->at(0);
		//double bz=beamPosz->at(0);
		//double p_beam=beamMomentum->at(0);
		//double ke_beam_MeV=1000.*p2ke(p_beam); //unit:MeV

		double csda_val=csda_range_vs_mom_sm->Eval(p_beam); //expected csda value if proton stops; unit: cm
		double ke_range_MeV=1000.*ke_vs_csda_range_sm->Eval(range_reco); //unit:MeV
		double p_range_MeV=1000.*ke2p(ke_range_MeV/1000.); //MeV/c
		double norm_trklen=range_reco/csda_val; //trklen/csda_val
		//cout<<"csda_val:"<<csda_val<<" range_reco:"<<range_reco<<" norm_trklen:"<<norm_trklen<<" min_norm_trklen_csda:"<<min_norm_trklen_csda<<endl;

		//reco Inel cut ------------------------------------------------------------------------------------//
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		bool IsRecoEL=false;
		if (norm_trklen>min_norm_trklen_csda) IsRecoStop=true; //old stopping proton cut
		//if (norm_trklen<=min_norm_trklen_csda) IsRecoInEL=true; //old inelastic proton cut

                if ((range_reco/csda_val)<min_norm_trklen_csda) { //inel region
                        if (pid>pid_1) IsRecoInEL=true;
                        if (pid<=pid_1) IsRecoEL=true;
                } //inel region
                if ((range_reco/csda_val)>=min_norm_trklen_csda&&(range_reco/csda_val)<max_norm_trklen_csda) { //stopping p region
                        if (pid>pid_2) IsRecoInEL=true;
                        if (pid<=pid_2) IsRecoEL=true;
                } //stopping p region
		

		//Energy loss -----------------------------------------------------------------//
		double mean_Elosscalo_stop=(4.85990e+01)/(1.00263e+00); //from fit result

		//fit result ==========================================================//
   		//p0           4.85990e+01   8.01697e-01  -2.71531e-08   8.79542e-12
   		//p1          -1.00263e+00   1.67201e-02   1.67201e-02   1.66581e-08
		//fit result ==========================================================//

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

		//const E-loss asump ----------------------------------------
		//double keffbeam=ke_beam_MeV-const_eloss_data;
		double keffbeam=(ke_beam_MeV-Eloss_data_hy_stop)*R_fit_hy;
		//double keffbeam=fitted_KE;

		//ke at end point ---------------------------------------------------------------------//
		//double kebb=-50; if (fitted_KE>0) kebb=BB.KEAtLength(fitted_KE, range_reco);
		double kebb=-9999.; kebb=BB.KEAtLength(keffbeam, range_reco);
		double kecalo=-9999.; kecalo=keffbeam-reco_calo_MeV;

		//KE definition using const E-loss assumption -----------------------------------------------------------//
		//double KEcalo_reco_constEcalo=-1; KEcalo_reco_constEcalo=ke_beam_MeV-mean_Elosscalo_stop-reco_calo_MeV;
		//double KEff_reco=ke_beam_MeV-mean_Elosscalo_stop;
		//double KEend_reco=KEff_reco-reco_calo_MeV;
		double KEff_reco=keffbeam;
		double KEend_reco=kebb;

		//evt selection cuts
		bool PassCuts_INT=false; //all bq cut+reco inel cut
		bool PassCuts_INC=false; //all bq cut
		if (IsPandoraSlice&&IsCaloSize) {

			//SliceID definitions
			//end
			//reco_sliceID = int((Emax-KEend_reco)/thinslicewidth);
			reco_sliceID = int(floor((Emax-KEend_reco)/thinslicewidth));			
			if (reco_sliceID < 0) reco_sliceID = -1;
			if (reco_endz<0) reco_sliceID = -1;
			if (reco_sliceID >= nthinslices) reco_sliceID = nthinslices;

			//st
			//reco_st_sliceID=int((Emax-KEff_reco)/thinslicewidth+0.5);
			//reco_st_sliceID=int(ceil((Emax-KEff_reco)/thinslicewidth));
			reco_st_sliceID=int(ceil((Emax-KEff_reco)/thinslicewidth));

      			if (reco_st_sliceID<0) reco_st_sliceID=-1; //KE higher than Emax
			if (reco_endz < 0) reco_st_sliceID = -1; 
			if (reco_st_sliceID >= nthinslices) reco_st_sliceID = nthinslices;

			if (IsBQ&&IsRecoInEL) {
				PassCuts_INT=true; //for INT 
				h_recosliceid_recoinelastic_cuts->Fill(reco_sliceID);
			}
			if (IsBQ) {
				PassCuts_INC=true; //for INC

				h_recosliceid_allevts_cuts->Fill(reco_sliceID);
				h_reco_st_sliceid_allevts_cuts->Fill(reco_st_sliceID);
			}
		}		

	} //evt loop

	//save results ---------//
	SaveHistograms();









}
