#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <cmath>
#//include <iostream>


#include <vector>
#include <iostream>
#include "THStack.h"

#include "../headers/BasicParameters.h"
#include "../headers/BasicAnaFunc.h"
#include "../headers/TemplateFitter.h"

void plot_BeamXY_All_and_Others() {

	//load data & MC
	TString fmc=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_kebkg_bmrw.root");
	TString fdata=Form("../data_kebkg.root");

	//output folder
	TString fout="./plot_BeamXY/beamxy_all_others_data_mc.eps";

	//plot style --------------------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(0);

	//read MC -----------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());
	TH2D *bxy_mc=(TH2D*)f_mc->Get("bx_by_RecoAll");
	TH2D *bxy_mc_midp=(TH2D*)f_mc->Get("bx_by_RecoMidP");
	TH2D *bxy_mc_inel=(TH2D*)f_mc->Get("bx_by_RecoInEl");

	//read data -----------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fdata.Data());
	TH2D *bxy_data=(TH2D*)f_data->Get("bx_by_RecoAll");
	TH2D *bxy_data_midp=(TH2D*)f_data->Get("bx_by_RecoMidP");
	TH2D *bxy_data_inel=(TH2D*)f_data->Get("bx_by_RecoInEl");

        float bx_min=-50;
        float bx_max=-10;
        int n_bx=40;
        float by_min=405;
        float by_max=440;
        int n_by=35;

	float kemin=0;
	float ke_max=600;
	TCanvas *c1=new TCanvas("c1","", 4600, 2000);
	c1->Divide(3,2);

	float meanX_data=bxy_data->GetMean(1);
	float rmsX_data=bxy_data->GetRMS(1);
	float meanY_data=bxy_data->GetMean(2);
	float rmsY_data=bxy_data->GetRMS(2);
	cout<<"double meanX_data="<<meanX_data<<endl;
	cout<<"double rmsX_data="<<rmsX_data<<endl;
	cout<<"double meanY_data="<<meanY_data<<endl;
	cout<<"double rmsY_data="<<rmsY_data<<endl;		

   	TEllipse *all_data = new TEllipse(meanX_data,meanY_data,1.5*rmsX_data,1.5*rmsY_data);
   	all_data->SetLineWidth(1);
	all_data->SetLineColor(2);
	all_data->SetLineStyle(2);
	all_data->SetFillStyle(0);

	float meanX_mc=bxy_mc->GetMean(1);
	float rmsX_mc=bxy_mc->GetRMS(1);
	float meanY_mc=bxy_mc->GetMean(2);
	float rmsY_mc=bxy_mc->GetRMS(2);
   	TEllipse *all_mc = new TEllipse(meanX_mc,meanY_mc,1.5*rmsX_mc,1.5*rmsY_mc);
   	all_mc->SetLineWidth(1);
	all_mc->SetLineColor(2);
	all_mc->SetLineStyle(2);
	all_mc->SetFillStyle(0);

	cout<<"double meanX_mc="<<meanX_mc<<endl;
	cout<<"double rmsX_mc="<<rmsX_mc<<endl;
	cout<<"double meanY_mc="<<meanY_mc<<endl;
	cout<<"double rmsY_mc="<<rmsY_mc<<endl;

	//el
	double meanX_el_data=-31.3139;
	double rmsX_el_data=3.79366;
	double meanY_el_data=422.116;
	double rmsY_el_data=3.48005;
	double meanX_el_mc=-29.1637;
	double rmsX_el_mc=4.50311;
	double meanY_el_mc=421.76;
	double rmsY_el_mc=3.83908;

   	TEllipse *el_data = new TEllipse(meanX_el_data, meanY_el_data,1.5*rmsX_el_data,1.5*rmsY_el_data);
   	el_data->SetLineWidth(1);
	el_data->SetLineColor(3);
	el_data->SetLineStyle(2);
	el_data->SetFillStyle(0);

   	TEllipse *el_mc = new TEllipse(meanX_el_mc, meanY_el_mc,1.5*rmsX_el_mc,1.5*rmsY_el_mc);
   	el_mc->SetLineWidth(1);
	el_mc->SetLineColor(3);
	el_mc->SetLineStyle(2);
	el_mc->SetFillStyle(0);

	c1->cd(1);
	TH2D *f2d_data=new TH2D("f2d_data","",n_bx,bx_min,bx_max,n_by,by_min,by_max);
	f2d_data->SetTitle("All Protons [Data]; Beam X-position [cm]; Beam Y-position [cm]");
	f2d_data->Draw();
	bxy_data->Draw("colz same");
	all_data->Draw("same");
	el_data->Draw("same");

	c1->cd(2);
	TH2D *f2d_data2=new TH2D("f2d_data2","",n_bx,bx_min,bx_max,n_by,by_min,by_max);
	f2d_data2->SetTitle("MisID:P-rich Sample [Data]; Beam X-position [cm]; Beam Y-position [cm]");
	f2d_data2->Draw();
	bxy_data_midp->Draw("colz same");
	all_data->Draw("same");
	el_data->Draw("same");


	c1->cd(3);
	TH2D *f2d_data3=new TH2D("f2d_data3","",n_bx,bx_min,bx_max,n_by,by_min,by_max);
	f2d_data3->SetTitle("Reco Inelastic [Data]; Beam X-position [cm]; Beam Y-position [cm]");
	f2d_data3->Draw();
	bxy_data_inel->Draw("colz same");
	all_data->Draw("same");
	el_data->Draw("same");



	//
	c1->cd(4);
	TH2D *f2d_mc=new TH2D("f2d_mc","",n_bx,bx_min,bx_max,n_by,by_min,by_max);
	f2d_mc->SetTitle("All Protons [MC]; Beam X-position [cm]; Beam Y-position [cm]");
	f2d_mc->Draw();
	bxy_mc->Draw("colz same");
	all_mc->Draw("same");
	el_mc->Draw("same");


	c1->cd(5);
	TH2D *f2d_mc2=new TH2D("f2d_mc2","",n_bx,bx_min,bx_max,n_by,by_min,by_max);
	f2d_mc2->SetTitle("MisID:P-rich Sample [MC]; Beam X-position [cm]; Beam Y-position [cm]");
	f2d_mc2->Draw();
	bxy_mc_midp->Draw("colz same");
	all_mc->Draw("same");
	el_mc->Draw("same");


	c1->cd(6);
	TH2D *f2d_mc3=new TH2D("f2d_mc3","",n_bx,bx_min,bx_max,n_by,by_min,by_max);
	f2d_mc3->SetTitle("Reco Inelastic [MC]; Beam X-position [cm]; Beam Y-position [cm]");
	f2d_mc3->Draw();
	bxy_mc_inel->Draw("colz same");
	all_mc->Draw("same");
	el_mc->Draw("same");







	c1->Print(Form("%s",fout.Data()));






}
