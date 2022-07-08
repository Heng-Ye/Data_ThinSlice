#include "../headers/BasicParameters.h"
#include "../headers/BasicAnaFunc.h"
#include "../headers/Eff.h"
#include "TEfficiency.h"
#include "TGraphSmooth.h" 

#include <iostream>
#include <fstream>
#include <string>

void plotRecoEff_Data_MC() {

	//plot style -----------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(0);
	//---------------------------------------------------------------//

	//data/mc config. -------------------//
	//data
	//vector<int> data{0, 5, 8, 9, 10, 12, 15, 20, 25, 30};
	//vector<int> color_data{11, 6, 4, 38, 46, 7, 2, 3, 8, 1};
	vector<int> data{0, 5, 8, 9, 10, 12, 15, 20, 30};
	vector<int> color_data{11, 6, 4, 38, 46, 7, 2, 3, 1};
	int nn=(int)data.size();
	const int n_data=(const int)nn;
	TH1D *h_data[n_data];
	TFile *f_data[n_data];
	std::cout<<"n_data:"<<n_data<<std::endl;

	//mc
	//vector<int> mc{0, 5, 8, 9, 10, 12, 15, 20, 25, 30}; //50 is weird, 17 eff> 1
	//vector<int> color_mc{11, 6, 4, 38, 46, 7, 2, 3, 8, 1};
	vector<int> mc{0, 5, 8, 9, 10, 12, 15, 20, 30}; //50 is weird, 17 eff> 1
	vector<int> color_mc{11, 6, 4, 38, 46, 7, 2, 3, 1};
	//vector<int> mc{0, 5, 10, 20, 30};
	//vector<int> color_mc{11, 4, 2, 3, 1};
	int nn_mc=(int)mc.size();
	const int n_mc=(const int)nn_mc;
	TH1D *h_mc[n_mc];
	TFile *f_mc[n_mc];
	std::cout<<"n_mc:"<<n_mc<<std::endl;

	double plot_min=0;
	double plot_max=40;

	//Truncated histograms -------------------------------------//
        TCanvas *ct=new TCanvas("ct","");
        ct->Divide(1,1);
	//TH2D *f2dt=new TH2D("f2dt","",80,0,80,100,0.1,5000);
	TH2D *f2dt=new TH2D("f2dt","",70,plot_min,plot_max,100,0.1, 1000000);
	f2dt->SetTitle("; Reco Track Length after SCE Correction [cm]; Counts");
        //ct->cd(1)->SetLogy();
        //ct->cd(1)->SetLogx();
      
	TPad *grid0 = new TPad("grid0","",0,0,1,1);
   	grid0->Draw();
   	grid0->cd();
   	grid0->SetGrid();
        grid0->SetLogy();

  
	f2dt->Draw();
        //TLegend *legt = new TLegend(0.15616,0.716387,0.820917,0.87395);
        TLegend *legt = new TLegend(0.1,0.716387,0.9,0.87395);
	legt->SetLineColor(1);
	legt->SetFillColor(0);	
   	legt->SetBorderSize(3);
	legt->SetNColumns(8);

	//----------------------------------------------------------//

	//data
	double n0_data=-99;
	vector<double> x_data;
	vector<double> ex_data;
	vector<double> y_data;
	vector<double> ey_up_data;
	vector<double> ey_dn_data;
	for (int j=0; j<n_data; ++j) {
		//read data ---------------------------------------------------------------------------------//
		TString str_data=Form("/dune/data2/users/hyliao/reco_eff_study/data_%dcm.root",data.at(j));
		f_data[j] = new TFile(str_data.Data());
		h_data[j]=(TH1D *)f_data[j]->Get("h1d_trklen");
		h_data[j]->SetName(Form("h1d_%d",j));
		h_data[j]->SetLineColor(color_data.at(j));
		h_data[j]->SetMarkerColor(color_data.at(j));
		
		double xi_data=h_data[j]->GetMean();
		double ni_data=h_data[j]->GetEntries();
		double exi_data=h_data[j]->GetStdDev();
		//double exi_data=0.5;
		if (j==0) n0_data=ni_data;

		std::cout<<"data::"<<data.at(j)<<" cm"<<std::endl;
		std::cout<<"x:"<<xi_data<<" +-"<<exi_data<<std::endl;		
		std::cout<<ni_data<<std::endl;		

		if (j>0) {
			x_data.push_back(xi_data);
			ex_data.push_back(exi_data);

			double ri_data=ni_data/n0_data;
			double eyi_up_data = TEfficiency::ClopperPearson(n0_data, ni_data, 0.682689492137, true) - ri_data;
			double eyi_dn_data = -TEfficiency::ClopperPearson(n0_data, ni_data, 0.682689492137, false) + ri_data;
			y_data.push_back(ri_data);
			ey_up_data.push_back(eyi_up_data);
			ey_dn_data.push_back(eyi_dn_data);

			std::cout<<"eff_data:("<<ni_data<<"/"<<n0_data<<")="<<ri_data<<"\t+"<<eyi_up_data<<"\t-"<<eyi_dn_data<<"\n\n"<<std::endl;

			//add histogram to the plot ----------------------------------------//	
			legt->AddEntry(h_data[j], Form("Data: %d cm", data.at(j)), "ep");
        		//ct->cd(1);
			h_data[j]->Draw("ep same");
		}
		std::cout<<"\n\n"<<std::endl;

		std::cout<<"test::"<<h_data[j]->GetEntries()<<std::endl;		
	}

	//mc
	double n0_mc=-99;
	vector<double> x_mc;
	vector<double> ex_mc;
	vector<double> y_mc;
	vector<double> ey_up_mc;
	vector<double> ey_dn_mc;
	for (int j=0; j<n_mc; ++j) {
		//read mc ---------------------------------------------------------------------------------//
		TString str_mc=Form("/dune/data2/users/hyliao/reco_eff_study/mc_%dcm.root",mc.at(j));
		f_mc[j] = new TFile(str_mc.Data());
		h_mc[j]=(TH1D *)f_mc[j]->Get("h1d_trklen");
		h_mc[j]->SetName(Form("h1d_%d",j));
		h_mc[j]->SetLineColor(color_mc.at(j));
		h_mc[j]->SetMarkerColor(color_mc.at(j));
		
		//remove outliers
		h_mc[j]->GetXaxis()->SetRangeUser(plot_min, plot_max); //only focus on the range from 0 to 45, ignore outliers in long track length
		double xi_mc=h_mc[j]->GetMean();
		double ni_mc=h_mc[j]->GetEntries();
		double exi_mc=h_mc[j]->GetStdDev();
		//double exi_mc=0.5;
		if (j==0) n0_mc=ni_mc;

		std::cout<<"mc::"<<mc.at(j)<<" cm"<<std::endl;
		std::cout<<"x:"<<xi_mc<<" +-"<<exi_mc<<std::endl;		
		std::cout<<ni_mc<<std::endl;		

		if (j>0) {
			x_mc.push_back(xi_mc);
			ex_mc.push_back(exi_mc);

			double ri_mc=ni_mc/n0_mc;
   			//calculates the boundaries for the frequentist Clopper-Pearson interval
			//TEfficiency::ClopperPearson(Int_t total,Int_t passed,Double_t level,Bool_t bUpper)

   			//This interval is recommended by the PDG.
   			//Input: - total : number of total events
   			//       - passed: 0 <= number of passed events <= total
   			//       - level : confidence level
   			//       - bUpper: true  - upper boundary is returned
   			//                 false - lower boundary is returned
   			//Calculation reference: https://web.pa.msu.edu/people/brock/file_sharing/ATLAS/root/hist/hist/src/TEfficiency.cxx

			double eyi_up_mc = TEfficiency::ClopperPearson(n0_mc, ni_mc, 0.682689492137, true) - ri_mc;
			double eyi_dn_mc = -TEfficiency::ClopperPearson(n0_mc, ni_mc, 0.682689492137, false) + ri_mc;
			y_mc.push_back(ri_mc);
			ey_up_mc.push_back(eyi_up_mc);
			ey_dn_mc.push_back(eyi_dn_mc);

			std::cout<<"eff_mc:("<<ni_mc<<"/"<<n0_mc<<")="<<ri_mc<<"\t+"<<eyi_up_mc<<"\t-"<<eyi_dn_mc<<"\n\n"<<std::endl;

			//add histogram to the plot ----------------------------------------//	
			legt->AddEntry(h_mc[j], Form("MC: %d cm", mc.at(j)), "l");
        		//ct->cd(1);
			h_mc[j]->Scale(n0_data/n0_mc);
			h_mc[j]->Draw("hist same");
		}
		std::cout<<"\n\n"<<std::endl;

	}



       	//ct->cd(1);
	legt->Draw();
	ct->Print("./plots_reco_eff/dist.eps");	



	TGraphAsymmErrors *gr_data = new TGraphAsymmErrors(y_data.size(), &x_data.at(0), &y_data.at(0), &ex_data.at(0), &ex_data.at(0), &ey_dn_data.at(0), &ey_up_data.at(0));
	TGraphAsymmErrors *gr_mc = new TGraphAsymmErrors(y_mc.size(), &x_mc.at(0), &y_mc.at(0), &ex_mc.at(0), &ex_mc.at(0), &ey_dn_mc.at(0), &ey_up_mc.at(0));
	gr_data->SetName("gr_data");
	gr_mc->SetName("gr_mc");

	gr_mc->SetMarkerStyle(21);
	gr_mc->SetMarkerColor(2);
	gr_mc->SetLineColor(2);

        TCanvas *c0=new TCanvas("c0","");
        c0->Divide(1,1);
        c0->cd(1);
	//c0->SetGrid();

	TPad *grid = new TPad("grid","",0,0,1,1);
   	grid->Draw();
   	grid->cd();
   	grid->SetGrid();

	//TH2D *f2d=new TH2D("f2d","",80,0,80,10,0,1);
	TH2D *f2d=new TH2D("f2d","",70,plot_min,plot_max,10,0,1.2);
	f2d->SetTitle("; Reco Track Length after SCE Correction [cm]; Efficiency");
	f2d->Draw();
	gr_data->Draw("lp same");
	gr_mc->Draw("lp same");

        TLegend *leg = new TLegend(0.147681,0.705797,0.296796,0.87515);
	//leg->SetFillStyle(0);
	leg->SetLineColor(1);
	leg->SetFillColor(0);
   	leg->SetBorderSize(3);
	
	leg->AddEntry(gr_data, "Data", "ep");
	leg->AddEntry(gr_mc, "MC", "ep");
	leg->Draw();

	//try super-smoother
/*
  	TGraph *gr_mc_sm=new TGraph(y_mc.size(), &x_mc.at(0), &y_mc.at(0));
  	TGraphSmooth *gs_mc = new TGraphSmooth("supsmu");
  	gr_mc_sm= gs_mc->Approx(gr_mc_sm,"linear",9000);
  	gr_mc_sm->SetLineColor(4);
  	gr_mc_sm->Draw("same");

	//best-fit
	//data
	double fit_min=0;
	double fit_max=40;




	TF1 *eff_fit_data=new TF1("eff_fit_data",tanh_fit,fit_min,fit_max,2);
  	eff_fit_data->SetParameter(0, 12);
	eff_fit_data->SetParameter(1, 0.8/40);
	//func->SetParLimits(0,-1,1);
	eff_fit_data->SetLineColor(4);
	gr_data->Fit("eff_fit_data","remn");
	gr_data->Draw("same");


	//mc
	TF1 *eff_fit_mc=new TF1("eff_fit_mc",eff_fit,fit_min,fit_max,3);
  	eff_fit_mc->SetParameter(0, 0.96);
  	eff_fit_mc->SetParLimits(0, 0, 1.1);
	eff_fit_mc->SetParameter(1, 12);
  	eff_fit_mc->SetParLimits(1, 0, 100);
	eff_fit_mc->SetParameter(2, 100);
  	eff_fit_mc->SetParLimits(2, 0, 1000000);
	//func->SetParLimits(0,-1,1);
	eff_fit_mc->SetLineColor(4);
	gr_mc->Fit("eff_fit_mc","rem+");
	gr_mc->Draw("same");
*/

	c0->Print("./plots_reco_eff/eff.eps");

/*
	//Save output eff curve
	TString str_out=Form("/dune/data2/users/hyliao/reco_eff_study/eff.root");
	TFile *f_out = new TFile(str_out.Data(), "recreate");
		gr_data->Write();
		gr_mc->Write();
	f_out->Close();
*/


/*
	//read data ---------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fdata.Data());
	TH1D *data=(TH1D*)f_data->Get(Form("h1d_%s",rep.Data()));
	int n_data=data->Integral(); 

	data->SetLineColor(1); data->SetMarkerColor(1);

	
	//read MC -----------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *mc=(TH1D*)f_mc->Get(Form("h1d_%s",rep.Data()));
	TH1D *mc_inel=(TH1D*)f_mc->Get(Form("h1d_%s_inel",rep.Data()));
	TH1D *mc_el=(TH1D*)f_mc->Get(Form("h1d_%s_el",rep.Data()));
	TH1D *mc_midcosmic=(TH1D*)f_mc->Get(Form("h1d_%s_midcosmic",rep.Data()));
	TH1D *mc_midpi=(TH1D*)f_mc->Get(Form("h1d_%s_midpi",rep.Data()));
	TH1D *mc_midp=(TH1D*)f_mc->Get(Form("h1d_%s_midp",rep.Data()));
	TH1D *mc_midmu=(TH1D*)f_mc->Get(Form("h1d_%s_midmu",rep.Data()));
	TH1D *mc_mideg=(TH1D*)f_mc->Get(Form("h1d_%s_mideg",rep.Data()));
	TH1D *mc_midother=(TH1D*)f_mc->Get(Form("h1d_%s_midother",rep.Data()));

	mc->SetLineColor(15);
	mc->SetLineStyle(2);
	mc_inel->SetFillColor(2); mc_inel->SetLineColor(2);
	mc_el->SetFillColor(4); mc_el->SetLineColor(4);
	mc_midp->SetFillColor(3); mc_midp->SetLineColor(3);
	mc_midcosmic->SetFillColor(5); mc_midcosmic->SetLineColor(5);
	mc_midpi->SetFillColor(6); mc_midpi->SetLineColor(6);
	mc_midmu->SetFillColor(28); mc_midmu->SetLineColor(28);
	mc_mideg->SetFillColor(30); mc_mideg->SetLineColor(30);
	mc_midother->SetFillColor(15); mc_midother->SetLineColor(15);

        int n_mc_inel=mc_inel->Integral();
        int n_mc_el=mc_el->Integral();
        int n_midcosmic=mc_midcosmic->Integral();
        int n_midpi=mc_midpi->Integral();
        int n_midp=mc_midp->Integral();
        int n_midmu=mc_midmu->Integral();
        int n_mideg=mc_mideg->Integral();
        int n_midother=mc_midother->Integral();
        int n_mc=n_mc_inel+n_mc_el+n_midcosmic+n_midpi+n_midp+n_midmu+n_mideg+n_midother;
	double scale_mc=(double)n_data/(double)n_mc;
	//double scale_mc=(double)16298/(double)124117;


	TH1D *mc_rw=(TH1D*)f_mc->Get(Form("h1d_%s",rep1.Data()));
	TH1D *mc_rw_inel=(TH1D*)f_mc->Get(Form("h1d_%s_inel",rep1.Data()));
	TH1D *mc_rw_el=(TH1D*)f_mc->Get(Form("h1d_%s_el",rep1.Data()));
	TH1D *mc_rw_midcosmic=(TH1D*)f_mc->Get(Form("h1d_%s_midcosmic",rep1.Data()));
	TH1D *mc_rw_midpi=(TH1D*)f_mc->Get(Form("h1d_%s_midpi",rep1.Data()));
	TH1D *mc_rw_midp=(TH1D*)f_mc->Get(Form("h1d_%s_midp",rep1.Data()));
	TH1D *mc_rw_midmu=(TH1D*)f_mc->Get(Form("h1d_%s_midmu",rep1.Data()));
	TH1D *mc_rw_mideg=(TH1D*)f_mc->Get(Form("h1d_%s_mideg",rep1.Data()));
	TH1D *mc_rw_midother=(TH1D*)f_mc->Get(Form("h1d_%s_midother",rep1.Data()));

	mc_rw_inel->SetFillColor(2); mc_rw_inel->SetLineColor(2);
	mc_rw_el->SetFillColor(4); mc_rw_el->SetLineColor(4);
	mc_rw_midp->SetFillColor(3); mc_rw_midp->SetLineColor(3);
	mc_rw_midcosmic->SetFillColor(5); mc_rw_midcosmic->SetLineColor(5);
	mc_rw_midpi->SetFillColor(6); mc_rw_midpi->SetLineColor(6);
	mc_rw_midmu->SetFillColor(28); mc_rw_midmu->SetLineColor(28);
	mc_rw_mideg->SetFillColor(30); mc_rw_mideg->SetLineColor(30);
	mc_rw_midother->SetFillColor(15); mc_rw_midother->SetLineColor(15);

        int n_mc_rw_inel=mc_rw_inel->Integral();
        int n_mc_rw_el=mc_rw_el->Integral();
        int n_mc_rw_midcosmic=mc_rw_midcosmic->Integral();
        int n_mc_rw_midpi=mc_rw_midpi->Integral();
        int n_mc_rw_midp=mc_rw_midp->Integral();
        int n_mc_rw_midmu=mc_rw_midmu->Integral();
        int n_mc_rw_mideg=mc_rw_mideg->Integral();
        int n_mc_rw_midother=mc_rw_midother->Integral();
        int n_mc_rw=n_mc_rw_inel+n_mc_rw_el+n_mc_rw_midcosmic+n_mc_rw_midpi+n_mc_rw_midp+n_mc_rw_midmu+n_mc_rw_mideg+n_mc_rw_midother;
	double scale_mc_rw=(double)n_data/(double)n_mc_rw;

        //chi2 calc. -----------------------------------------------------------------//
        vector<double> D; //data
        vector<double> er_D; //error of data
        //h1d->Sumw2();
        for (int k=1; k<=data->GetNbinsX(); k++){
                D.push_back(data->GetBinContent(k));
                er_D.push_back(sqrt(data->GetBinContent(k)));
        }
        vector<double> MC; //MC
        vector<double> er_MC; //error of MC
        vector<double> MC_RW;
        vector<double> er_MC_RW;
        //double norm=(double)n_data/(double)(mc->Integral());
        //double norm_rw=(double)n_data/(double)(mc_rw->Integral());
        for (int k=1; k<=mc->GetNbinsX(); k++) {
                        MC.push_back(scale_mc*mc->GetBinContent(k));
                        er_MC.push_back(scale_mc*sqrt(mc->GetBinContent(k)));

                        MC_RW.push_back(scale_mc_rw*mc_rw->GetBinContent(k));
                        er_MC_RW.push_back(scale_mc_rw*sqrt(mc_rw->GetBinContent(k)));
	}

        double chi2=ml_data_mc(D, er_D, MC, er_MC);
        double chi2_rw=ml_data_mc(D, er_D, MC_RW, er_MC_RW);

	mc->Scale(scale_mc);
	mc_inel->Scale(scale_mc);
	mc_el->Scale(scale_mc);
	mc_midp->Scale(scale_mc);
	mc_midcosmic->Scale(scale_mc);
	mc_midpi->Scale(scale_mc);
	mc_midmu->Scale(scale_mc);
	mc_mideg->Scale(scale_mc);
	mc_midother->Scale(scale_mc);


	mc_rw->Scale(scale_mc_rw);
	mc_rw_inel->Scale(scale_mc_rw);
	mc_rw_el->Scale(scale_mc_rw);
	mc_rw_midp->Scale(scale_mc_rw);
	mc_rw_midcosmic->Scale(scale_mc_rw);
	mc_rw_midpi->Scale(scale_mc_rw);
	mc_rw_midmu->Scale(scale_mc_rw);
	mc_rw_mideg->Scale(scale_mc_rw);
	mc_rw_midother->Scale(scale_mc_rw);


	//mc->SetLineColor(4); mc->SetMarkerColor(4);
	//mc_rw->SetLineColor(2); mc_rw->SetMarkerColor(2);

	THStack* hs_rw=new THStack("hs_rw","");
	hs_rw->Add(mc_rw_inel);
	hs_rw->Add(mc_rw_el);
	hs_rw->Add(mc_rw_midp);
	hs_rw->Add(mc_rw_midcosmic);
	hs_rw->Add(mc_rw_midpi);
	hs_rw->Add(mc_rw_midmu);
	hs_rw->Add(mc_rw_mideg);
	hs_rw->Add(mc_rw_midother);

        //data/MC -----------------------------------------//
        TH1D *R=(TH1D*)data->Clone();
        R->Divide(data, mc);
	R->SetLineColor(15);
	R->SetMarkerColor(15);
	//R->SetLineStyle(2);

        TH1D *R_rw=(TH1D*)data->Clone();
        R_rw->Divide(data, mc_rw);
	R_rw->SetLineColor(1);
	//R_rw->SetLineStyle(2);
        //-------------------------------------------------//

	//Proton Momentum --------------------------------------------------------------//
	//TCanvas *c0=new TCanvas("c0","");
	TCanvas *c0=new TCanvas("c0","",1200,1500);
	c0->Divide(1,2);
	c0->cd(1);

	float xmin=0;
	//float xmax=40;
	float xmax=140;
	//float xmax=100;

	//float xmax=40;

	//float xmax=100;

	//float ymax=1200;
	//float ymax=1400;
	//float ymax=300;
	//float ymax=25;
	//float ymax=80;
	//float ymax=200;
	//float ymax=40;
	float ymax=900;
	

	//TH2D* frame2d=new TH2D("frame2d","", 140, xmin, xmax, 1000, 0, 1000); //trklen
	TH2D* frame2d=new TH2D("frame2d","", 140, xmin, xmax, 1000, 0, ymax); //trklen
	//TH2D* frame2d=new TH2D("frame2d","", 140, 0, 140, 200, 0, 200); //endz
	frame2d->SetTitle(Form(";%s;",x_axis_label.Data()));
	frame2d->GetXaxis()->CenterTitle();
	frame2d->Draw();
	hs_rw->Draw("hist same");
	mc_rw->Draw("hist same");
	mc->Draw("hist same");
	data->Draw("ep same");


	TLegend *leg0 = new TLegend(0.14,0.7,.75,0.9);
        TLegend *leg = new TLegend(0.14,0.54,.75,0.69);

	//TLegend *leg0 = new TLegend(0.34,0.7,.85,0.9);
        //TLegend *leg = new TLegend(0.34,0.54,.85,0.69);


	//TLegend *leg0 = new TLegend(0.34,0.7,.85,0.9);
        //TLegend *leg = new TLegend(0.34,0.54,.85,0.69);


	leg0->SetFillStyle(0);
	leg0->AddEntry(data, "Data", "ep");
	leg0->AddEntry(mc,Form("MC (before reweighting) - #chi^{2}/dof:%.2f/%d",chi2,data->GetNbinsX()), "l");
	leg0->AddEntry(mc_rw,Form("MC (after reweighting) - #chi^{2}/dof:%.2f/%d",chi2_rw,data->GetNbinsX()-2), "l");
	leg0->Draw();

        //TLegend *leg = new TLegend(0.14,0.54,.75,0.69);
        leg->AddEntry(mc_rw_inel, "Inel","f");
        leg->AddEntry(mc_rw_el, "El","f");

        leg->AddEntry(mc_rw_midcosmic,"misID:cosmic","f");
        leg->AddEntry(mc_rw_midp, "misID:p","f");
        leg->AddEntry(mc_rw_midpi, "misID:#pi","f");

        leg->AddEntry(mc_rw_midmu,"misID:#mu","f");
        leg->AddEntry(mc_rw_mideg, "misID:e/#gamma","f");
        leg->AddEntry(mc_rw_midother, "misID:other","f");

        leg->SetNColumns(3);
        leg->Draw();


        //pDUNE Logo
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(0.002, ymax+8, Form("#bf{DUNE:ProtoDUNE-SP}"));
        //txt_pdune1[0]=new TLatex(0.002, ymax+2, Form("#bf{DUNE:ProtoDUNE-SP}"));
        //txt_pdune1[0]=new TLatex(0.002, ymax+.5, Form("#bf{DUNE:ProtoDUNE-SP}"));

        txt_pdune1[0]->SetTextColor(1);
        txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        txt_p1[0]=new TLatex(100,ymax+8, Form("Protons (1 GeV/c)"));
        //txt_p1[0]=new TLatex(100,ymax+2, Form("Protons (1 GeV/c)"));
        //txt_p1[0]=new TLatex(70,ymax+.5, Form("Protons (1 GeV/c)"));

        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        txt_p1[0]->Draw();

        c0->cd(2);

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0.02);
        pad2->SetBottomMargin(0.25);
        pad2->SetGridx();
        pad2->SetGridy();
        pad2->Draw();
        pad2->cd();


        TH2D *f2d2=new TH2D("f2d2",Form("%s","Pos"),100,xmin,xmax,64,0,6.4); //liny
        //TH2D *f2d2=new TH2D("f2d2",Form("%s",""),150,xmin,xmax,100,0,5); //liny
        f2d2->SetTitle(Form(";%s;Data/MC",x_axis_label.Data()));
        //f2d2->GetXaxis()->SetLabelSize(0.1);
        //f2d2->GetYaxis()->SetLabelSize(0.1);

        //f2d2->GetXaxis()->SetTitleSize(0.1);
        //f2d2->GetYaxis()->SetTitleSize(0.1);
        //f2d2->GetYaxis()->SetTitleOffset(.5);
        f2d2->Draw();
        TLine *line1=new TLine(xmin,1,xmax,1);
        line1->SetLineColor(7);
        line1->SetLineWidth(2);
        line1->Draw("same");
        R->Draw("ep same");
	R_rw->Draw("ep same");

        TLegend *leg2 = new TLegend(0.14,0.74,.6,0.86);
        //TLegend *leg2 = new TLegend(0.54,0.74,.9,0.86);

        leg2->AddEntry(R, "Before reweighting","l");
        leg2->AddEntry(R_rw, "After reweighting","l");
        leg2->Draw();


	c0->Print(Form("%s",outpath.Data()));




	pbeam_data->Draw("ep same");
	pbeam_stop_data->Draw("ep same");
	pbeam_stop_fit_data->Draw("same");
	//pbeam_fit_mc->Draw("same");

	pbeam_mc->Draw("hist same");
        pbeam_stop_mc->Draw("hist same");

        TF1* fit_pbeam_mc=VFit(pbeam_mc, 2);
	fit_pbeam_mc->Draw("same");
        TF1* fit_pbeam_stop_mc=VFit(pbeam_stop_mc, 3);
	fit_pbeam_stop_mc->Draw("same");

	//pbeam_stop_fit_mc->Draw("same");
	//read rw histograms -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

	TLegend *leg0 = new TLegend(0.14,0.65,.6,0.85);
	leg0->SetFillStyle(0);
	leg0->AddEntry(pbeam_data, "Data (Protons before entering TPC)", "ep");
	leg0->AddEntry(pbeam_stop_data, "Data (Stopping Protons)", "ep");
	leg0->AddEntry(pbeam_mc, "MC (Protons before entering TPC)", "l");
	leg0->AddEntry(pbeam_stop_mc, "MC (Stopping Protons)", "l");
	leg0->Draw();

        //pDUNE Logo
        TLatex **txt_pdune=new TLatex*[1];
        txt_pdune[0]=new TLatex(600.002, 706, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune[0]->SetTextColor(1);
        txt_pdune[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p=new TLatex*[1];
        txt_p[0]=new TLatex(1145,706, Form("Protons (1 GeV/c)"));
        txt_p[0]->SetTextColor(1);
        txt_pdune[0]->SetTextSize(0.07);
        txt_p[0]->Draw();
        //
	c0->Print(Form("%s/pbeam_data_mc.eps",outpath.Data()));


	//Beam (mu, sigma) -------------------------------------------------------//
	vector<float> mu_data;
	vector<float> err_mu_data;
	vector<float> sigma_data;
	vector<float> err_sigma_data;
	mu_data.push_back(pbeam_fit_data->GetParameter(0));		err_mu_data.push_back(pbeam_fit_data->GetParError(0));
	sigma_data.push_back(pbeam_fit_data->GetParameter(1));		err_sigma_data.push_back(pbeam_fit_data->GetParError(1));
	mu_data.push_back(pbeam_stop_fit_data->GetParameter(0));	err_mu_data.push_back(pbeam_stop_fit_data->GetParError(0));
	sigma_data.push_back(pbeam_stop_fit_data->GetParameter(1));	err_sigma_data.push_back(pbeam_stop_fit_data->GetParError(1));

	vector<float> mu_mc;
	vector<float> err_mu_mc;
	vector<float> sigma_mc;
	vector<float> err_sigma_mc;
	mu_mc.push_back(pbeam_fit_mc->GetParameter(0));		err_mu_mc.push_back(pbeam_fit_mc->GetParError(0));
	sigma_mc.push_back(pbeam_fit_mc->GetParameter(1));	err_sigma_mc.push_back(pbeam_fit_mc->GetParError(1));

	mu_mc.push_back(pff_stop_fit_mc->GetParameter(0));	err_mu_mc.push_back(pff_stop_fit_mc->GetParError(0));
	sigma_mc.push_back(pff_stop_fit_mc->GetParameter(1));	err_sigma_mc.push_back(pff_stop_fit_mc->GetParError(1));

	mu_mc.push_back(pbeam_stop_fit_mc->GetParameter(0));	err_mu_mc.push_back(pbeam_stop_fit_mc->GetParError(0));
	sigma_mc.push_back(pbeam_stop_fit_mc->GetParameter(1));	err_sigma_mc.push_back(pbeam_stop_fit_mc->GetParError(1));

	TGraphErrors *ms_data=new TGraphErrors(mu_data.size(), &mu_data.at(0), &sigma_data.at(0), &err_mu_data.at(0), &err_sigma_data.at(0));
	TGraphErrors *ms_mc=new TGraphErrors(mu_mc.size(), &mu_mc.at(0), &sigma_mc.at(0), &err_mu_mc.at(0), &err_sigma_mc.at(0));
	
	ms_mc->SetMarkerColor(2);
	ms_mc->SetLineColor(2);

	TCanvas *c1=new TCanvas("c1","");
	c1->Divide(1,1);
	c1->cd(1);
	TH2D* frame2d2=new TH2D("frame2d2","", 200, 930, 1040, 42, 48, 90);
	frame2d2->SetTitle(";#mu [MeV/c];#sigma [MeV/c]");
	frame2d2->GetXaxis()->CenterTitle();
	frame2d2->Draw();
	ms_data->Draw("p same");
	ms_mc->Draw("p same");

	float emin=930.;
	float emax=1040.;
	TF1 *fit_data = new TF1("fit_data", "[0]+[1]*x", emin, emax);
	fit_data->SetLineColor(1);
	fit_data->SetLineStyle(2);
	ms_data->Fit(fit_data,"remn");
	fit_data->Draw("same");
	
        TF1 *fit_mc = new TF1("fit_mc", "[0]+[1]*x", emin, emax);
        fit_mc->SetLineColor(2);
        fit_mc->SetLineStyle(2);
        ms_mc->Fit(fit_mc,"remn");
        fit_mc->Draw("same");

	TLegend *leg1 = new TLegend(0.14,0.65,.6,0.85);
	leg1->SetFillStyle(0);
	leg1->AddEntry(ms_data, Form("Data: #sigma=%.2f+%.2f*#mu",fit_data->GetParameter(0),fit_data->GetParameter(1)), "ep");
	leg1->AddEntry(ms_mc, Form("MC: #sigma=%.2f+%.2f*#mu",fit_mc->GetParameter(0),fit_mc->GetParameter(1)), "ep");
	leg1->Draw();

        //pDUNE Logo
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(930.002, 90.6, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune1[0]->SetTextColor(1);
        txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        txt_p1[0]=new TLatex(1002,90.6, Form("Protons (1 GeV/c)"));
        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        txt_p1[0]->Draw();


	c1->Print(Form("%s/pbeam_mus_igma_data_mc.eps",outpath.Data()));

*/



}
