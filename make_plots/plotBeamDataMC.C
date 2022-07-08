#include "../headers/BasicParameters.h"
#include "../headers/BasicAnaFunc.h"

void plotBeamDataMC(TString fdata, TString fmc, TString outpath) {
	//TString rep="trklen";
	//TString x_axis_label="Proton Track Length [cm]";

	//plot style --------------------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(0);

	//read data ---------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fdata.Data());
	TH1D *pbeam_data=(TH1D*)f_data->Get("h1d_pbeam");
	TF1 *pbeam_fit_data=(TF1*)f_data->Get("pbeam_fit");
	//TH1D* h1d=(TH1D*)f_data->Get(Form("h1d_%s_stop",rep.Data()));

	TH1D *pbeam_stop_data=(TH1D*)f_data->Get("h1d_pbeam_stop");
	TF1 *pbeam_stop_fit_data=(TF1*)f_data->Get("pbeam_stop_fit");
	int n_data=pbeam_data->Integral(); 
	int n_stop_data=pbeam_stop_data->Integral(); 
	cout<<"n_data:"<<n_data<<endl;
	cout<<"n_stop_data:"<<n_stop_data<<endl;

	pbeam_data->SetLineColor(1); 	      pbeam_data->SetMarkerColor(1);
	pbeam_fit_data->SetLineColor(1);      pbeam_fit_data->SetMarkerColor(1);
	pbeam_stop_data->SetLineColor(4);     pbeam_stop_data->SetMarkerColor(4);
	pbeam_stop_fit_data->SetLineColor(4); pbeam_stop_fit_data->SetMarkerColor(4);

	//read MC -----------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *pbeam_mc=(TH1D*)f_mc->Get("h1d_pbeam");
	TF1 *pbeam_fit_mc=(TF1*)f_mc->Get("pbeam_fit");

	TH1D *pbeam_stop_mc=(TH1D*)f_mc->Get("h1d_pbeam_stop");
	TF1 *pbeam_stop_fit_mc=(TF1*)f_mc->Get("pbeam_stop_fit");

	TH1D *pff_stop_mc=(TH1D*)f_mc->Get("h1d_pff_stop");
	TF1 *pff_stop_fit_mc=(TF1*)f_mc->Get("pff_stop_fit");

	int n_mc=pbeam_mc->Integral(); 
	int n_stop_mc=pbeam_stop_mc->Integral(); 
	int n_ff_mc=pff_stop_mc->Integral(); 
	cout<<"n_mc:"<<n_mc<<endl;
	cout<<"n_stop_mc:"<<n_stop_mc<<endl;

	pbeam_mc->Scale((double)n_data/(double)n_mc);
	pbeam_stop_mc->Scale((double)n_data/(double)n_mc);
	//pbeam_stop_mc->Scale((double)(n_stop_data)/(double)n_stop_mc);
	pff_stop_mc->Scale((double)n_data/(double)n_mc);
	//pbeam_fit_mc->Scale((double)n_data/(double)n_mc);
	//cout<<"pbeam_mc:"<<pbeam_mc->Integral()<<endl;

	pbeam_mc->SetLineColor(2); 	     pbeam_mc->SetMarkerColor(2);
	pbeam_fit_mc->SetLineColor(2);       pbeam_fit_mc->SetMarkerColor(2);
	pbeam_stop_mc->SetLineColor(3);      pbeam_stop_mc->SetMarkerColor(3);
	pbeam_stop_fit_mc->SetLineColor(2);  pbeam_stop_fit_mc->SetMarkerColor(2);
	pff_stop_mc->SetLineColor(6);        pff_stop_mc->SetMarkerColor(6);
	pff_stop_fit_mc->SetLineColor(6);    pff_stop_fit_mc->SetMarkerColor(6);

	//Proton Momentum --------------------------------------------------------------//
	TCanvas *c0=new TCanvas("c0","");
	c0->Divide(1,1);
	c0->cd(1);
	TH2D* frame2d=new TH2D("frame2d","", 600, 600, 1400, 700, 0, 700); //zend_2d
	frame2d->SetTitle(";Proton Momentum [MeV/c];");
	frame2d->GetXaxis()->CenterTitle();
	frame2d->Draw();
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



	//Proton Momentum [only MC] --------------------------------------------------------------//
	TCanvas *c01=new TCanvas("c01","");
	c01->Divide(1,1);
	c01->cd(1);
	TH2D* frame2d01=new TH2D("frame2d01","", 600, 600, 1400, 700, 0, 700); //zend_2d
	frame2d01->SetTitle(";Proton Momentum [MeV/c];");
	frame2d01->GetXaxis()->CenterTitle();
	frame2d01->Draw();
	pbeam_mc->Draw("hist same");
	fit_pbeam_mc->Draw("hist same");

        pbeam_stop_mc->Draw("hist same");
	fit_pbeam_stop_mc->Draw("same");

	pff_stop_mc->Draw("hist same");
        TF1* fit_pff_stop_mc=VFit(pff_stop_mc, 6);
	fit_pff_stop_mc->Draw("same");


	TLegend *leg01 = new TLegend(0.14,0.65,.6,0.85);
	leg01->SetFillStyle(0);
	leg01->AddEntry(pbeam_mc, "MC (Protons before entering TPC)", "l");
	leg01->AddEntry(pff_stop_mc, "MC (Protons at TPC Front Face)", "l");
	leg01->AddEntry(pbeam_stop_mc, "MC (Stopping Protons)", "l");
	leg01->Draw();

        //pDUNE Logo
        //TLatex **txt_pdune=new TLatex*[1];
        //txt_pdune[0]=new TLatex(600.002, 706, Form("#bf{DUNE:ProtoDUNE-SP}"));
        //txt_pdune[0]->SetTextColor(1);
        txt_pdune[0]->Draw();
        //
        //Beam Logo
        //TLatex **txt_p=new TLatex*[1];
        //txt_p[0]=new TLatex(1145,706, Form("Protons (1 GeV/c)"));
        //txt_p[0]->SetTextColor(1);
        //txt_pdune[0]->SetTextSize(0.07);
        txt_p[0]->Draw();

	c01->Print(Form("%s/pbeam_only_mc.eps",outpath.Data()));




	//Beam (mu, sigma) -------------------------------------------------------------------------------------------------------------------//
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

	//mu_mc.push_back(pff_stop_fit_mc->GetParameter(0));	err_mu_mc.push_back(pff_stop_fit_mc->GetParError(0));
	//sigma_mc.push_back(pff_stop_fit_mc->GetParameter(1));	err_sigma_mc.push_back(pff_stop_fit_mc->GetParError(1));

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
        //txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        txt_p1[0]=new TLatex(1002,90.6, Form("Protons (1 GeV/c)"));
        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        //txt_p1[0]->Draw();


	c1->Print(Form("%s/pbeam_mus_igma_data_mc.eps",outpath.Data()));





}
