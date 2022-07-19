#include "../headers/BasicParameters.h"
#include "../headers/BasicAnaFunc.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "TTimeStamp.h" 

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <memory>
#include <ctime>

#include <time.h>

#include "TTimeStamp.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFrame.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TH1D.h"
#include "TH1.h"
#include "TLatex.h"
#include "TText.h"

void plot_beamXY_distributions() {
	//TString rep="trklen";
	//TString x_axis_label="Proton Track Length [cm]";

	//plot style -------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(1);
	//------------------------------------------------------------//

	//data runs ---------------//
	const int n=26;
	vector<int> RUN(n);
	RUN[0]=5219;
	RUN[1]=5225;
	RUN[2]=5235;
	RUN[3]=5240;
	RUN[4]=5244;
	RUN[5]=5308;
	RUN[6]=5311;
	RUN[7]=5315;
	RUN[8]=5338;
	RUN[9]=5387;
	RUN[10]=5423;
	RUN[11]=5424;
	RUN[12]=5426;
	RUN[13]=5455;
	RUN[14]=5456;
	RUN[15]=5457;
	RUN[16]=5458;
	RUN[17]=5460;
	RUN[18]=5809;
	RUN[19]=5810;
	RUN[20]=5814;
	RUN[21]=5816;
	RUN[22]=5817;
	RUN[23]=5842;
	RUN[24]=5843;
	RUN[25]=5844;
	//-------------------------//

	vector<double> t0;
	vector<double> tmed;
	vector<double> tmax;
	vector<double> zero(n,0.);
	vector<double> ex;
	char *run_label[n];

	vector<std::string> str_evt_time;

	TFile **file=new TFile*[n];

	//TH1D *reco_startZ_sce_all;
	//TH1D *reco_startY_sce_all;
	//TH1D *reco_startX_sce_all;

	TH2D **bx_by=new TH2D*[n];
	TH2D *bx_by_all;



	
/*
	TF1 **fit_zst=new TF1*[n];
	TF1 **fit_yst=new TF1*[n];
	TF1 **fit_xst=new TF1*[n];

	vector<double> mu_zst;
	vector<double> err_mu_zst;
	vector<double> sigma_zst;
	vector<double> err_sigma_zst;

	vector<double> mu_yst;
	vector<double> err_mu_yst;
	vector<double> sigma_yst;
	vector<double> err_sigma_yst;

	vector<double> mu_xst;
	vector<double> err_mu_xst;
	vector<double> sigma_xst;
	vector<double> err_sigma_xst;
*/

	//c1->cd(1)->SetLogy();

/*
	TH2D *f2d_time_ratio=new TH2D("f2d_time_ratio","",17,0,17,10,0.8,1.3);
	TH2D *f2d_time_prange=new TH2D("f2d_time_prange","",17,0,17,100,0,1800);
	TH2D *f2d_time_pcalo=new TH2D("f2d_time_pcalo","",17,0,17,100,0,1800);

	f2d_time_ratio->SetTitle("Stopping Protons; Time since Data Taking[hr]; P_{calo}/P_{range}");
	f2d_time_prange->SetTitle("Stopping Protons; Time since Data Taking[hr]; P_{range} [MeV/c]");
	f2d_time_pcalo->SetTitle("Stopping Protons; Time since Data Taking[hr]; P_{calo} [MeV/c]");
	c1->cd(1);
	f2d_time_ratio->Draw();
	c1->cd(2);
	f2d_time_prange->Draw();
	c1->cd(3);
	f2d_time_pcalo->Draw();
*/



  	//ofstream myfile_output;
  	//myfile_output.open (Form("./poscut/pos_cut.csv"));
	for (int i=0; i<n; ++i) {
		//int i=9;
		//Read data -----------------------------------------------------------------------------------------------------//
		TString str_name;
		str_name=Form("/dune/data2/users/hyliao/protonana/v09_39_01/Basic_Cuts/proton_ApplybasicCuts_run%d.root", RUN.at(i));
		std::cout<<"Processing data:  "<<str_name.Data()<<std::endl;

		file[i]=new TFile(str_name.Data());

		bx_by[i]=(TH2D *)file[i]->Get(Form("bx_by"));

		TCanvas *c_xy=new TCanvas("c_xy","", 600, 600);
		c_xy->Divide(1, 1);

		TH2D *f2d_bx_by=new TH2D("f2d_bx_by","", 60, -50, -10, 35, 405, 440);
		f2d_bx_by->GetXaxis()->SetTitle("X [cm]");
		f2d_bx_by->GetYaxis()->SetTitle("Y [cm]");

		TString str_xy;
		if (i==0) str_xy=Form("./beamxy_all.pdf(");
		if (i>0&&i<n-1) str_xy=Form("./beamxy_all.pdf");
		if(i==n-1) str_xy=Form("./beamxy_all.pdf)");

		if (i==0) {
			bx_by_all=(TH2D*)bx_by[i]->Clone();
		}
		if (i>0) { 
			bx_by_all->Add(bx_by[i]);
		}


		//Time Info
		TParameter<double> *tmp_t0=(TParameter<double>*)file[i]->Get("T0");
		TParameter<double> *tmp_tmed=(TParameter<double>*)file[i]->Get("Tmed");
		TParameter<double> *tmp_tmax=(TParameter<double>*)file[i]->Get("Tmax");

		auto ttmp_t0=tmp_t0->GetVal();
		auto ttmp_tmed=tmp_tmed->GetVal();
		auto ttmp_tmax=tmp_tmax->GetVal();
		auto dt=(ttmp_tmax-ttmp_t0)/(3600.);

		t0.push_back(ttmp_t0);
		tmed.push_back(ttmp_tmed);
		tmax.push_back(ttmp_tmax);
		ex.push_back(0.5*(ttmp_tmax-ttmp_t0));

		//evt time conversion
		//epochtime to string
		struct tm tm;
		char buf_time[255];
		memset(&tm, 0, sizeof(struct tm));
		strptime(Form("%f",ttmp_t0), "%s", &tm);
		strftime(buf_time, sizeof(buf_time), "%b %d %H:%M %Y", &tm);
		str_evt_time.push_back(buf_time);
		//puts(buf_time); //frotmat:: Nov 12 00:46 2018
		std::cout<<"  EvtTime: "<<ttmp_tmed<<" ("<<buf_time<<")"<<endl;
		//txt_zst_stop[j]=new TText(5, 200, Form("(%s)",substr_evt_time[j].c_str()));

		//float plot_ymax_z=reco_startZ_sce[i]->GetBinContent(reco_startZ_sce[i]->GetMaximumBin());
		//float plot_ymax_y=reco_startY_sce[i]->GetBinContent(reco_startY_sce[i]->GetMaximumBin());
		//float plot_ymax_x=reco_startX_sce[i]->GetBinContent(reco_startX_sce[i]->GetMaximumBin());

		//TH2D *f2d_zst=new TH2D("f2d_zst","",15,-5,10,plot_ymax_z,0,plot_ymax_z*1.1);
		//TH2D *f2d_yst=new TH2D("f2d_yst","",150,350,500,plot_ymax_y,0,plot_ymax_y*1.1);
		//TH2D *f2d_xst=new TH2D("f2d_xst","",100,-80,20,plot_ymax_x,0,plot_ymax_x*1.1);

		bx_by_all->SetTitle(Form("Run %d: t_{0}:%s (%.2f hrs); Z [cm]; Counts",RUN.at(i),buf_time, dt));
		//f2d_yst->SetTitle(Form("; Y [cm]; Counts"));
		//f2d_xst->SetTitle(Form("; X [cm]; Counts"));


		bx_by[i]->GetXaxis()->SetTitle("X [cm]");
		bx_by[i]->GetYaxis()->SetTitle("Y [cm]");
		bx_by[i]->Draw("colz");	

		c_xy->Print(str_xy.Data());

		delete c_xy;
		delete f2d_bx_by;



/*
		TCanvas *c_zyx_st=new TCanvas("c_zyx_st","", 1200, 1800);
		c_zyx_st->Divide(1, 3);

		TString str_zyx_st_out;
		if (i==0) str_zyx_st_out=Form("./zyx_st_all.pdf(");
		if (i>0&&i<n-1) str_zyx_st_out=Form("./zyx_st_all.pdf");
		if(i==n-1) str_zyx_st_out=Form("./zyx_st_all.pdf)");

		c_zyx_st->cd(1);
		f2d_zst->Draw();
		reco_startZ_sce[i]->Draw("same");
		if (i==10||i==8) fit_zst[i]=VNFit(reco_startZ_sce[i], 4.5, 6);
		else fit_zst[i]=VNFit(reco_startZ_sce[i], 3.4, 6);
		fit_zst[i]->SetLineColor(2);
		fit_zst[i]->SetLineStyle(2);
		fit_zst[i]->Draw("same");
		double tmp_mu_zst=fit_zst[i]->GetParameter(0);
		double tmp_sigma_zst=fit_zst[i]->GetParameter(1);
		mu_zst.push_back(tmp_mu_zst);
		//err_mu_zst.push_back(fit_zst[i]->GetParErrors(0));
		sigma_zst.push_back(tmp_sigma_zst);
		//err_sigma_zst.push_back(fit_zst[i]->GetParErrors(1));

		c_zyx_st->cd(2);
		f2d_yst->Draw();
		reco_startY_sce[i]->Draw("same");
		fit_yst[i]=VNFit(reco_startY_sce[i], 425, 3);
		fit_yst[i]->SetLineColor(2);
		fit_yst[i]->SetLineStyle(2);
		fit_yst[i]->Draw("same");
		double tmp_mu_yst=fit_yst[i]->GetParameter(0);
		double tmp_sigma_yst=fit_yst[i]->GetParameter(1);
		mu_yst.push_back(tmp_mu_yst);
		//err_mu_yst.push_back(fit_yst[i]->GetParErrors(0));
		sigma_yst.push_back(tmp_sigma_yst);
		//err_sigma_yst.push_back(fit_yst[i]->GetParErrors(1));

		c_zyx_st->cd(3);
		f2d_xst->Draw();
		reco_startX_sce[i]->Draw("same");
		fit_xst[i]=VNFit(reco_startX_sce[i], -30, 3);
		fit_xst[i]->SetLineColor(2);
		fit_xst[i]->SetLineStyle(2);
		fit_xst[i]->Draw("same");
		double tmp_mu_xst=fit_xst[i]->GetParameter(0);
		double tmp_sigma_xst=fit_xst[i]->GetParameter(1);
		mu_xst.push_back(tmp_mu_xst);
		//err_mu_xst.push_back(fit_xst[i]->GetParErrors(0));
		sigma_xst.push_back(tmp_sigma_xst);
		//err_sigma_xst.push_back(fit_xst[i]->GetParErrors(1));

		//output file ----------------------------------------------//
		myfile_output << RUN.at(i) << ","<< tmp_mu_zst<<","<< tmp_sigma_zst<<","<<tmp_mu_yst<<","<<tmp_sigma_yst<<","<<tmp_mu_xst<<","<<tmp_sigma_xst<<"\n";
		c_zyx_st->Print(str_zyx_st_out.Data());
		delete f2d_zst;
		delete f2d_yst;
		delete f2d_xst;
		delete c_zyx_st;
*/


	}
  	//myfile_output.close();


/*
	TLegend *leg = new TLegend(0.14,0.65,.6,0.85);
	leg->SetFillStyle(0);
	leg->AddEntry(h1d_dz_all, "Fit then Sum", "ep");
	leg->AddEntry(h1d_dz_global_all, "Sum then Fit", "ep");

	TCanvas *c1_xyz=new TCanvas("c1_xyz","", 1200, 1800);
	c1_xyz->Divide(1,4);

	c1_xyz->cd(1);
	h1d_dz_all->SetMarkerColor(2);
	h1d_dz_all->SetLineColor(2);
	h1d_dz_all->GetXaxis()->SetTitle("#DeltaZ/#sigma_{z}");
	h1d_dz_all->GetXaxis()->SetTitleOffset(1.1);
	h1d_dz_all->GetXaxis()->CenterTitle(true);
	h1d_dz_all->Draw();
	h1d_dz_global_all->SetLineColor(1);
	h1d_dz_global_all->SetMarkerColor(1);
	h1d_dz_global_all->Draw("same");
	leg->Draw();
	
	c1_xyz->cd(2);
	h1d_dy_all->SetMarkerColor(2);
	h1d_dy_all->SetLineColor(2);
	h1d_dy_all->GetXaxis()->SetTitle("#DeltaY/#sigma_{Y}");
	h1d_dy_all->GetXaxis()->SetTitleOffset(1.1);
	h1d_dy_all->GetXaxis()->CenterTitle(true);
	h1d_dy_all->Draw();
	h1d_dy_global_all->SetLineColor(1);
	h1d_dy_global_all->SetMarkerColor(1);
	h1d_dy_global_all->Draw("same");
	leg->Draw();

	c1_xyz->cd(3);
	h1d_dx_all->SetMarkerColor(2);
	h1d_dx_all->SetLineColor(2);
	h1d_dx_all->GetXaxis()->SetTitle("#DeltaX/#sigma_{X}");
	h1d_dx_all->GetXaxis()->SetTitleOffset(1.1);
	h1d_dx_all->GetXaxis()->CenterTitle(true);
	h1d_dx_all->Draw();
	h1d_dx_global_all->SetLineColor(1);
	h1d_dx_global_all->SetMarkerColor(1);
	h1d_dx_global_all->Draw("same");
	leg->Draw();


	c1_xyz->cd(4);
	h1d_dxy_all->SetMarkerColor(2);
	h1d_dxy_all->SetLineColor(2);
	h1d_dxy_all->GetXaxis()->SetTitle("#sqrt{(#DeltaX/#sigma_{X})^{2}+(#DeltaY/#sigma_{Y})^{2}}");
	h1d_dxy_all->GetXaxis()->SetTitleOffset(1.1);
	h1d_dxy_all->GetXaxis()->CenterTitle(true);
	h1d_dxy_all->Draw();
	h1d_dxy_global_all->SetLineColor(1);
	h1d_dxy_global_all->SetMarkerColor(1);
	h1d_dxy_global_all->Draw("same");
	leg->Draw();

	c1_xyz->Print("poscut.eps");


	//All distributions -------------------------------------------------//
  	ofstream myfile_all_output;
  	myfile_all_output.open(Form("./poscut/pos_cut_all.csv"));

	TCanvas *c1_zyx_st_all=new TCanvas("c1_zyx_st_all","", 1200, 1800);
	c1_zyx_st_all->Divide(1,3);
	c1_zyx_st_all->cd(1);
	reco_startZ_sce_all->Draw();
	TF1 *fit_zst_all=VNFit(reco_startZ_sce_all, 3.4, 5);
	fit_zst_all->SetLineColor(2);
	fit_zst_all->SetLineStyle(2);
	fit_zst_all->Draw("same");
	double mu_zst_all=fit_zst_all->GetParameter(0);
	double sigma_zst_all=fit_zst_all->GetParameter(1);

	c1_zyx_st_all->cd(2);
	reco_startY_sce_all->Draw();
	TF1 *fit_yst_all=VNFit(reco_startY_sce_all, 425, 3);
	fit_yst_all->SetLineColor(2);
	fit_yst_all->SetLineStyle(2);
	fit_yst_all->Draw("same");
	double mu_yst_all=fit_yst_all->GetParameter(0);
	double sigma_yst_all=fit_yst_all->GetParameter(1);

	c1_zyx_st_all->cd(3);
	reco_startX_sce_all->Draw();
	TF1 *fit_xst_all=VNFit(reco_startX_sce_all, -24, 3);
	fit_xst_all->SetLineColor(2);
	fit_xst_all->SetLineStyle(2);
	fit_xst_all->Draw("same");
	c1_zyx_st_all->Print("zyx_st_all.eps");
	double mu_xst_all=fit_xst_all->GetParameter(0);
	double sigma_xst_all=fit_xst_all->GetParameter(1);
	myfile_all_output << mu_zst_all<<","<< sigma_zst_all<<","<<mu_yst_all<<","<<sigma_yst_all<<","<<mu_xst_all<<","<<sigma_xst_all<<"\n";
  	myfile_all_output.close();


	TCanvas *c1_t_xyz=new TCanvas("c1_t_xyz","", 1200, 1800);
	c1_t_xyz->Divide(1,3);
	c1_t_xyz->cd(1);
	TGraphErrors *t_zst=new TGraphErrors(tmed.size(),&tmed.at(0),&mu_zst.at(0),&ex.at(0),&sigma_zst.at(0));
	t_zst->GetXaxis()->SetTimeDisplay(1);
	t_zst->GetXaxis()->SetTimeFormat("%b/%d");
	t_zst->GetXaxis()->SetTimeOffset(0,"gmt");
	t_zst->SetTitle(Form("; Time; Z Position [cm]"));
	t_zst->Draw("ap| ");
	TLine* l_zall=new TLine(tmed.at(0), mu_zst_all, tmed.at(-1+tmed.size()), mu_zst_all);
	TLine* lup_zall=new TLine(tmed.at(0), mu_zst_all+sigma_zst_all, tmed.at(-1+tmed.size()), mu_zst_all+sigma_zst_all);
	TLine* ldn_zall=new TLine(tmed.at(0), mu_zst_all-sigma_zst_all, tmed.at(-1+tmed.size()), mu_zst_all-sigma_zst_all);
	l_zall->SetLineColor(2);
	lup_zall->SetLineColor(2);
	ldn_zall->SetLineColor(2);
	l_zall->SetLineWidth(2);
	lup_zall->SetLineWidth(2);
	ldn_zall->SetLineWidth(2);
	lup_zall->SetLineStyle(2);
	ldn_zall->SetLineStyle(2);
	l_zall->Draw();
	lup_zall->Draw();
	ldn_zall->Draw();

	c1_t_xyz->cd(2);
	TGraphErrors *t_yst=new TGraphErrors(tmed.size(),&tmed.at(0),&mu_yst.at(0),&ex.at(0),&sigma_yst.at(0));
	t_yst->GetXaxis()->SetTimeDisplay(1);
	t_yst->GetXaxis()->SetTimeFormat("%b/%d");
	t_yst->GetXaxis()->SetTimeOffset(0,"gmt");
	t_yst->SetTitle(Form("; Time; Y Position [cm]"));
	t_yst->Draw("ap| ");
	TLine* l_yall=new TLine(tmed.at(0), mu_yst_all, tmed.at(-1+tmed.size()), mu_yst_all);
	TLine* lup_yall=new TLine(tmed.at(0), mu_yst_all+sigma_yst_all, tmed.at(-1+tmed.size()), mu_yst_all+sigma_yst_all);
	TLine* ldn_yall=new TLine(tmed.at(0), mu_yst_all-sigma_yst_all, tmed.at(-1+tmed.size()), mu_yst_all-sigma_yst_all);
	l_yall->SetLineColor(2);
	lup_yall->SetLineColor(2);
	ldn_yall->SetLineColor(2);
	l_yall->SetLineWidth(2);
	lup_yall->SetLineWidth(2);
	ldn_yall->SetLineWidth(2);
	lup_yall->SetLineStyle(2);
	ldn_yall->SetLineStyle(2);
	l_yall->Draw();
	lup_yall->Draw();
	ldn_yall->Draw();

	c1_t_xyz->cd(3);
	TGraphErrors *t_xst=new TGraphErrors(tmed.size(),&tmed.at(0),&mu_xst.at(0),&ex.at(0),&sigma_xst.at(0));
	t_xst->GetXaxis()->SetTimeDisplay(1);
	t_xst->GetXaxis()->SetTimeFormat("%b/%d");
	t_xst->GetXaxis()->SetTimeOffset(0,"gmt");
	t_xst->SetTitle(Form("; Time; X Position [cm]"));
	t_xst->Draw("ap| ");
	TLine* l_xall=new TLine(tmed.at(0), mu_xst_all, tmed.at(-1+tmed.size()), mu_xst_all);
	TLine* lup_xall=new TLine(tmed.at(0), mu_xst_all+sigma_xst_all, tmed.at(-1+tmed.size()), mu_xst_all+sigma_xst_all);
	TLine* ldn_xall=new TLine(tmed.at(0), mu_xst_all-sigma_xst_all, tmed.at(-1+tmed.size()), mu_xst_all-sigma_xst_all);
	l_xall->SetLineColor(2);
	lup_xall->SetLineColor(2);
	ldn_xall->SetLineColor(2);
	l_xall->SetLineWidth(2);
	lup_xall->SetLineWidth(2);
	ldn_xall->SetLineWidth(2);
	lup_xall->SetLineStyle(2);
	ldn_xall->SetLineStyle(2);
	l_xall->Draw();
	lup_xall->Draw();
	ldn_xall->Draw();

	c1_t_xyz->Print(Form("time_vs_ZYX.eps"));
*/



	// Vertical alignment.
	//auto 
	//auto *tv1 = new TText(0.66,0.165,"Bottom adjusted");
	//tv1->SetTextAlign(11); tv1->SetTextSize(0.12);
	//tv1->Draw();

	// Draw labels on the y axis
	//TText *t = new TText();
	//t->SetTextAlign(32);
	//t->SetTextSize(0.035);
	//t->SetTextFont(72);
	//char *labels[6] = {"Jan98","Feb98","Mar98","Apr98","May98","Jun98"};
	//TFile *fout = new TFile(Form("/dune/data2/users/hyliao/protonana/v09_39_01/KEstudy/proton_beamxy_beammom_run%d.root",run),"RECREATE");












}
