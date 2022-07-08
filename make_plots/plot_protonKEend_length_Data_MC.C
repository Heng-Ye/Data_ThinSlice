#include <vector>
#include <iostream>
#include "THStack.h"
#include "../headers/BasicParameters.h"
#include "../headers/BasicAnaFunc.h"
#include "../headers/TemplateFitter.h"

void plot_protonKEend_Data_MC() {

	//MisID:P scaling factor --------//
        //double scal_MisIDP=0.662056;
	//double err_scal_MisIDP=0.186766;	
	//-------------------------------//

	//TString x_axis_label="Proton Track Length [cm]";
	
	TString fmc=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_kecalo_nobmrw.root");
	TString fmc_bmrw=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_kecalo_bmrw.root");
	TString fdata=Form("../data_kecalo.root");

	//plot style --------------------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(0);

	//read data ---------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fdata.Data());
	TH1D *h2d_data=(TH1D*)f_data->Get(Form("h2d_trklen_KEcalo"));

	//read MC ---------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *h2d_mc=(TH1D*)f_mc->Get(Form("h2d_trklen_KEcalo"));





	TH1D *data=(TH1D*)f_data->Get(Form("h1d_%s",rep_data.Data()));
	TH1D *data_inel=(TH1D*)f_data->Get(Form("h1d_%s_RecoInel",rep_data.Data()));
	TH1D *data_el=(TH1D*)f_data->Get(Form("h1d_%s_RecoEl",rep_data.Data()));

	int n_data=data->Integral(); 
	data->SetLineColor(1); data->SetMarkerColor(1);
	data_inel->SetLineColor(1); data_inel->SetMarkerColor(1);
	data_el->SetLineColor(1); data_el->SetMarkerColor(1);

	//read MC -----------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *mc=(TH1D*)f_mc->Get(Form("h1d_%s",rep_mc.Data()));
	TH1D *mc_inel=(TH1D*)f_mc->Get(Form("h1d_%s_inel",rep_mc.Data()));
	TH1D *mc_el=(TH1D*)f_mc->Get(Form("h1d_%s_el",rep_mc.Data()));
	TH1D *mc_midcosmic=(TH1D*)f_mc->Get(Form("h1d_%s_midcosmic",rep_mc.Data()));
	TH1D *mc_midpi=(TH1D*)f_mc->Get(Form("h1d_%s_midpi",rep_mc.Data()));
	TH1D *mc_midp=(TH1D*)f_mc->Get(Form("h1d_%s_midp",rep_mc.Data()));
	TH1D *mc_midmu=(TH1D*)f_mc->Get(Form("h1d_%s_midmu",rep_mc.Data()));
	TH1D *mc_mideg=(TH1D*)f_mc->Get(Form("h1d_%s_mideg",rep_mc.Data()));
	TH1D *mc_midother=(TH1D*)f_mc->Get(Form("h1d_%s_midother",rep_mc.Data()));

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

	//bmrw
	TFile *f_mc_bmrw = TFile::Open(fmc_bmrw.Data());	
	TH1D *mc_rw=(TH1D*)f_mc_bmrw->Get(Form("h1d_%s",rep_mc.Data()));
	TH1D *mc_rw_inel=(TH1D*)f_mc_bmrw->Get(Form("h1d_%s_inel",rep_mc.Data()));
	TH1D *mc_rw_el=(TH1D*)f_mc_bmrw->Get(Form("h1d_%s_el",rep_mc.Data()));
	TH1D *mc_rw_midcosmic=(TH1D*)f_mc_bmrw->Get(Form("h1d_%s_midcosmic",rep_mc.Data()));
	TH1D *mc_rw_midpi=(TH1D*)f_mc_bmrw->Get(Form("h1d_%s_midpi",rep_mc.Data()));
	TH1D *mc_rw_midp=(TH1D*)f_mc_bmrw->Get(Form("h1d_%s_midp",rep_mc.Data()));
	TH1D *mc_rw_midmu=(TH1D*)f_mc_bmrw->Get(Form("h1d_%s_midmu",rep_mc.Data()));
	TH1D *mc_rw_mideg=(TH1D*)f_mc_bmrw->Get(Form("h1d_%s_mideg",rep_mc.Data()));
	TH1D *mc_rw_midother=(TH1D*)f_mc_bmrw->Get(Form("h1d_%s_midother",rep_mc.Data()));

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

	cout<<"data->GetNbinsX():"<<data->GetNbinsX()<<endl;
	cout<<"mc->GetNbinsX():"<<mc->GetNbinsX()<<endl;


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
        //chi2 calc. -----------------------------------------------------------------//



	//normaliation ----------------------------------------//
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

	THStack* hs_rw=new THStack("hs_rw","");
	hs_rw->Add(mc_rw_inel);
	hs_rw->Add(mc_rw_el);
	hs_rw->Add(mc_rw_midp);
	hs_rw->Add(mc_rw_midcosmic);
	hs_rw->Add(mc_rw_midpi);
	hs_rw->Add(mc_rw_midmu);
	hs_rw->Add(mc_rw_mideg);
	hs_rw->Add(mc_rw_midother);
	//------------------------------------------------------//





	float xmin=0;
	//float xmax=40;
	float xmax=500;
	//float xmax=100;

	//float xmax=40;

	//float xmax=100;

	//float ymax=1200;
	//float ymax=1400;
	//float ymax=300;
	//float ymax=25;
	//float ymax=80;
	float ymax=1850;
	//float ymax=40;
	//float ymax=900;


	//TH2D* frame2d=new TH2D("frame2d","", 140, xmin, xmax, 1000, 0, 1000); //trklen
	TH2D* frame2d=new TH2D("frame2d","", 140, xmin, xmax, 1000, 0, ymax); //trklen
	//TH2D* frame2d=new TH2D("frame2d","", 140, 0, 140, 200, 0, 200); //endz
	frame2d->SetTitle(Form(";Proton Kinetic Energy MeV]; Counts"));
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








/*
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
	float ymax=850;
	//float ymax=40;
	//float ymax=900;
	

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

*/

/*
        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0.02);
        pad2->SetBottomMargin(0.25);
        pad2->SetGridx();
        pad2->SetGridy();
        pad2->Draw();
        pad2->cd();
*/

/*
        TH2D *f2d2=new TH2D("f2d2",Form("%s",""),100,xmin,xmax,64,0,6.4); //liny
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

	double xfit_min=30.;
	double xfit_max=110.;
        double n_sigma=1.;

	//BKG-Fit:: MisID:P ------------------------------------------------------//
	//[1]prepare histograms for fitting, ignore short track length area
	double trklen_thr=30.0; //track length threshold
	double trklen_hthr=110.0; //track length threshold

	TH1D *h0_data=new TH1D("h0_data","", data->GetNbinsX(), data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());
	TH1D *h2=new TH1D("h2","", data->GetNbinsX(), data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());
	TH1D *h1_mcother=new TH1D("h1_mcother","", data->GetNbinsX(), data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());

	double av_R_rw=0;
	double err_av_R_rw=0;
	double n_R_rw=0;
	for (int i = 1; i<=data->GetNbinsX(); ++i){
		double x_cen=data->GetBinCenter(i);
		//std::cout<<"x_cen: "<<x_cen<<std::endl;
			
		if (x_cen>xfit_min&&x_cen<xfit_max) { //ignore short track length area
			h0_data->SetBinContent(i, data->GetBinContent(i));
			h0_data->SetBinError(i, data->GetBinError(i));

			double n_mc_el=mc_rw_el->GetBinContent(i);
			double err_n_mc_el=mc_rw_el->GetBinError(i);

			double n_mc_all=mc_rw->GetBinContent(i);
			double err_n_mc_all=mc_rw->GetBinError(i);
			
			double n_mc_other=n_mc_all-n_mc_el;
			double err_n_mc_other=sqrt(err_n_mc_all*err_n_mc_all+err_n_mc_el*err_n_mc_el);

			h1_mcother->SetBinContent(i, n_mc_other);
			h1_mcother->SetBinError(i, err_n_mc_other);

			h2->SetBinContent(i, n_mc_el);
			h2->SetBinError(i, err_n_mc_el);

			//av_R_rw
			n_R_rw++;
			av_R_rw+=R_rw->GetBinContent(i);
			err_av_R_rw+=(R_rw->GetBinError(i))*(R_rw->GetBinError(i));
		} //ignore short track length area	
	}
	av_R_rw/=n_R_rw;
	err_av_R_rw=sqrt(err_av_R_rw)/n_R_rw;

	//TF1 *f_el = new TF1("f_el", "[0]", xfit_min, xfit_max); 
	//f_el->SetLineColor(2);
	//R_rw->Fit("f_el","em0");

	//double bestfit_el=f_el->GetParameter(0);
	//double err_bestfit_el=f_el->GetParError(0);
	double bestfit_el=av_R_rw;
	double err_bestfit_el=err_av_R_rw;
	cout<<"\n\nbestfit_el:"<<bestfit_el<<" +- "<<err_bestfit_el<<endl;

	TLine* l_misidp=new TLine(xfit_min, bestfit_el, xfit_max, bestfit_el);
	TLine* lup_misidp=new TLine(xfit_min, bestfit_el+n_sigma*err_bestfit_el, xfit_max, bestfit_el+n_sigma*err_bestfit_el);
	TLine* ldn_misidp=new TLine(xfit_min, bestfit_el-n_sigma*err_bestfit_el, xfit_max, bestfit_el-n_sigma*err_bestfit_el);
	l_misidp->SetLineColor(2);
	l_misidp->SetLineWidth(3);
	lup_misidp->SetLineColor(2); lup_misidp->SetLineStyle(2);
	ldn_misidp->SetLineColor(2); ldn_misidp->SetLineStyle(2);
	c0->cd(2);
	l_misidp->Draw("same");
	lup_misidp->Draw("same");
	ldn_misidp->Draw("same");

        //TLegend *leg2 = new TLegend(0.14,0.74,.6,0.86);
        //TLegend *leg2 = new TLegend(0.54,0.74,.9,0.86);
        TLegend *leg2 = new TLegend(0.14,0.64,.6,0.86);        
        leg2->AddEntry(R, "Before reweighting","l");
        leg2->AddEntry(R_rw, "After reweighting","l");

	TLegendEntry* text_bestfit = leg2->AddEntry(l_misidp, Form("<Data/MC>: %.2f#pm%.2f",bestfit_el,err_bestfit_el),"l");
	text_bestfit->SetTextColor(2);
        //leg2->AddEntry(l_misidp, Form("Best-fit: %.2f#pm%.2f",scal_midp.at(0),err_scal_midp.at(0))," ");
        leg2->Draw();


	//[2]Min.(h0-h1-s*h2)
        TemplateFitter fitter;
        fitter.SetHistograms(h0_data, h1_mcother, h2);
        //fitter.SetFitRange(4,-1);
        fitter.Fit();

       	vector<double> scal_el;
        vector<double> err_scal_el;
        scal_el.push_back(fitter.GetPar());
        err_scal_el.push_back(fitter.GetParError());

	cout<<"scal_el.size():"<<scal_el.size()<<endl;
        cout<<"scal_el:"<<scal_el.at(0)<<" +- "<<err_scal_el.at(0)<<endl;

	c0->Print(Form("%s",outpath.Data()));

	TH1D *h2_corr=new TH1D("h2_corr","", data->GetNbinsX(), data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());
	for (int i = 1; i<=h2->GetNbinsX(); ++i) {
		//double x_cen=data->GetBinCenter(i);
		//std::cout<<"x_cen: "<<x_cen<<std::endl;
		double n_mc_el_corr=scal_el.at(0)*mc_rw_el->GetBinContent(i);
		double err_n_mc_el_corr=n_mc_el_corr*sqrt(pow((mc_rw_el->GetBinError(i)/mc_rw_el->GetBinContent(i)),2)+pow((err_scal_el.at(0)/scal_el.at(0)),2));
			
		h2_corr->SetBinContent(i, n_mc_el_corr);
		h2_corr->SetBinError(i, err_n_mc_el_corr);
	}

	//[3]Min.(h0-h1-s*h2_corr)
        TemplateFitter fitter_aftermisidp_corr;
        fitter_aftermisidp_corr.SetHistograms(h0_data, h1_mcother, h2_corr);
        fitter_aftermisidp_corr.Fit();

       	vector<double> scal_midp_corr;
        vector<double> err_scal_midp_corr;
        scal_midp_corr.push_back(fitter_aftermisidp_corr.GetPar());
        err_scal_midp_corr.push_back(fitter_aftermisidp_corr.GetParError());

	cout<<"scal_midp_corr.size():"<<scal_midp_corr.size()<<endl;
        cout<<"scal_midp_corr:"<<scal_midp_corr.at(0)<<" +- "<<err_scal_midp_corr.at(0)<<endl;

*/




/*
	TH1D *data=(TH1D*)f_data->Get(Form("h1d_%s",rep.Data()));
	int n_data=data->Integral(); 
	data->SetLineColor(1); data->SetMarkerColor(1);



        TemplateFitter fitter;
  
        //Min.(h0-h1-s*h2)
        //h0:data
        //h1:mc(mc except misID:p)
        //h2:mc(misID:p)

        fitter.SetHistograms(h1d_data, MC_1, MC_2);
        //fitter.SetFitRange(4,-1);
        fitter.Fit();
        scal_midp.push_back(fitter.GetPar());
        err_scal_midp.push_back(fitter.GetParError());

        cout<<"scal_midp:"<<scal_midp.at(0)<<" +- "<<err_scal_midp.at(0)<<endl;
*/














/*

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
