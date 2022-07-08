#include "../headers/BasicParameters.h"
#include "../headers/BasicAnaFunc.h"

void plotBeforeAfterBMRWAfterEffCorr(TString rep, TString rep1, TString x_axis_label, TString fdata, TString fmc, TString outpath) {
	//TString rep="trklen";
	//TString x_axis_label="Proton Track Length [cm]";
	
	//read efficiency curves -------------------------------------------------------//
	TFile *f_eff = TFile::Open("/dune/data2/users/hyliao/reco_eff_study/eff.root");
	TGraphAsymmErrors *gr_data=(TGraphAsymmErrors *)f_eff->Get("gr_data");
	TGraphAsymmErrors *gr_mc=(TGraphAsymmErrors *)f_eff->Get("gr_mc");

	cout<<"\n\n Test_data at 10 cm:"<<gr_data->Eval(10)<<endl;
	cout<<"Test_data at 20 cm:"<<gr_data->Eval(20)<<endl;

	//plot style --------------------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(0);

	//read data ---------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fdata.Data());
	TH1D *data_orig=(TH1D*)f_data->Get(Form("h1d_%s",rep.Data()));
	TH1D *data=new TH1D("data","", data_orig->GetNbinsX(), data_orig->GetXaxis()->GetXmin(), data_orig->GetXaxis()->GetXmax());
        for (int i = 1; i<=data_orig->GetNbinsX(); ++i){
		double n=data_orig->GetBinContent(i);
		double err_n=data_orig->GetBinError(i);
		double cen=data_orig->GetBinCenter(i);
		if (cen<=40) { //make eff corr if trklen < 40 cm
			data->SetBinContent(i, n/(double)gr_data->Eval(cen));
			data->SetBinError(i, err_n/(double)gr_data->Eval(cen));
		} //make eff corr if trklen < 40 cm
	}

	//TH1D *data=(TH1D*)f_data->Get(Form("h1d_%s",rep.Data()));
	int n_data=data->Integral(); 

	data->SetLineColor(1); data->SetMarkerColor(1);

	
	//read MC -----------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());

	TH1D *mc_orig=(TH1D*)f_mc->Get(Form("h1d_%s",rep.Data()));
	TH1D *mc_inel_orig=(TH1D*)f_mc->Get(Form("h1d_%s_inel",rep.Data()));
	TH1D *mc_el_orig=(TH1D*)f_mc->Get(Form("h1d_%s_el",rep.Data()));
	TH1D *mc_midcosmic_orig=(TH1D*)f_mc->Get(Form("h1d_%s_midcosmic",rep.Data()));
	TH1D *mc_midpi_orig=(TH1D*)f_mc->Get(Form("h1d_%s_midpi",rep.Data()));
	TH1D *mc_midp_orig=(TH1D*)f_mc->Get(Form("h1d_%s_midp",rep.Data()));
	TH1D *mc_midmu_orig=(TH1D*)f_mc->Get(Form("h1d_%s_midmu",rep.Data()));
	TH1D *mc_mideg_orig=(TH1D*)f_mc->Get(Form("h1d_%s_mideg",rep.Data()));
	TH1D *mc_midother_orig=(TH1D*)f_mc->Get(Form("h1d_%s_midother",rep.Data()));


	TH1D *mc=new TH1D("mc","", mc_orig->GetNbinsX(), mc_orig->GetXaxis()->GetXmin(), mc_orig->GetXaxis()->GetXmax());
	TH1D *mc_inel=new TH1D("mc_inel","", mc_inel_orig->GetNbinsX(), mc_inel_orig->GetXaxis()->GetXmin(), mc_inel_orig->GetXaxis()->GetXmax());
	TH1D *mc_el=new TH1D("mc_el","", mc_el_orig->GetNbinsX(), mc_el_orig->GetXaxis()->GetXmin(), mc_el_orig->GetXaxis()->GetXmax());
	TH1D *mc_midcosmic=new TH1D("mc_midcosmic","", mc_midcosmic_orig->GetNbinsX(), mc_midcosmic_orig->GetXaxis()->GetXmin(), mc_midcosmic_orig->GetXaxis()->GetXmax());
	TH1D *mc_midpi=new TH1D("mc_midpi","", mc_midpi_orig->GetNbinsX(), mc_midpi_orig->GetXaxis()->GetXmin(), mc_midpi_orig->GetXaxis()->GetXmax());
	TH1D *mc_midp=new TH1D("mc_midp","", mc_midp_orig->GetNbinsX(), mc_midp_orig->GetXaxis()->GetXmin(), mc_midp_orig->GetXaxis()->GetXmax());
	TH1D *mc_midmu=new TH1D("mc_midp","", mc_midmu_orig->GetNbinsX(), mc_midmu_orig->GetXaxis()->GetXmin(), mc_midmu_orig->GetXaxis()->GetXmax());
	TH1D *mc_mideg=new TH1D("mc_mideg","", mc_mideg_orig->GetNbinsX(), mc_mideg_orig->GetXaxis()->GetXmin(), mc_mideg_orig->GetXaxis()->GetXmax());
	TH1D *mc_midother=new TH1D("mc_midother","", mc_midother_orig->GetNbinsX(), mc_midother_orig->GetXaxis()->GetXmin(), mc_midother_orig->GetXaxis()->GetXmax());
        for (int i = 1; i<=mc_orig->GetNbinsX(); ++i){
		double cen=mc_orig->GetBinCenter(i);
		if (cen<=40) { //make eff corr if trklen < 40 cm
			mc->SetBinContent(i, (double)mc_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc->SetBinError(i, (double)mc_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_inel->SetBinContent(i, (double)mc_inel_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_inel->SetBinError(i, (double)mc_inel_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_el->SetBinContent(i, (double)mc_el_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_el->SetBinError(i, (double)mc_el_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_midcosmic->SetBinContent(i, (double)mc_midcosmic_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_midcosmic->SetBinError(i, (double)mc_midcosmic_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_midpi->SetBinContent(i, (double)mc_midpi_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_midpi->SetBinError(i, (double)mc_midpi_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_midp->SetBinContent(i, (double)mc_midp_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_midp->SetBinError(i, (double)mc_midp_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_midmu->SetBinContent(i, (double)mc_midmu_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_midmu->SetBinError(i, (double)mc_midmu_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_mideg->SetBinContent(i, (double)mc_mideg_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_mideg->SetBinError(i, (double)mc_mideg_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_midother->SetBinContent(i, (double)mc_midother_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_midother->SetBinError(i, (double)mc_midother_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

		} //make eff corr if trklen < 40 cm
	}






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


	TH1D *mc_rw_orig=(TH1D*)f_mc->Get(Form("h1d_%s",rep1.Data()));
	TH1D *mc_rw_inel_orig=(TH1D*)f_mc->Get(Form("h1d_%s_inel",rep1.Data()));
	TH1D *mc_rw_el_orig=(TH1D*)f_mc->Get(Form("h1d_%s_el",rep1.Data()));
	TH1D *mc_rw_midcosmic_orig=(TH1D*)f_mc->Get(Form("h1d_%s_midcosmic",rep1.Data()));
	TH1D *mc_rw_midpi_orig=(TH1D*)f_mc->Get(Form("h1d_%s_midpi",rep1.Data()));
	TH1D *mc_rw_midp_orig=(TH1D*)f_mc->Get(Form("h1d_%s_midp",rep1.Data()));
	TH1D *mc_rw_midmu_orig=(TH1D*)f_mc->Get(Form("h1d_%s_midmu",rep1.Data()));
	TH1D *mc_rw_mideg_orig=(TH1D*)f_mc->Get(Form("h1d_%s_mideg",rep1.Data()));
	TH1D *mc_rw_midother_orig=(TH1D*)f_mc->Get(Form("h1d_%s_midother",rep1.Data()));

	TH1D *mc_rw=new TH1D("mc_rw","", mc_rw_orig->GetNbinsX(), mc_rw_orig->GetXaxis()->GetXmin(), mc_rw_orig->GetXaxis()->GetXmax());
	TH1D *mc_rw_inel=new TH1D("mc_rw_inel","", mc_rw_inel_orig->GetNbinsX(), mc_rw_inel_orig->GetXaxis()->GetXmin(), mc_rw_inel_orig->GetXaxis()->GetXmax());
	TH1D *mc_rw_el=new TH1D("mc_rw_el","", mc_rw_el_orig->GetNbinsX(), mc_rw_el_orig->GetXaxis()->GetXmin(), mc_rw_el_orig->GetXaxis()->GetXmax());
	TH1D *mc_rw_midcosmic=new TH1D("mc_rw_midcosmic","", mc_rw_midcosmic_orig->GetNbinsX(), mc_rw_midcosmic_orig->GetXaxis()->GetXmin(), mc_rw_midcosmic_orig->GetXaxis()->GetXmax());
	TH1D *mc_rw_midpi=new TH1D("mc_rw_midpi","", mc_rw_midpi_orig->GetNbinsX(), mc_rw_midpi_orig->GetXaxis()->GetXmin(), mc_rw_midpi_orig->GetXaxis()->GetXmax());
	TH1D *mc_rw_midp=new TH1D("mc_rw_midp","", mc_rw_midp_orig->GetNbinsX(), mc_rw_midp_orig->GetXaxis()->GetXmin(), mc_rw_midp_orig->GetXaxis()->GetXmax());
	TH1D *mc_rw_midmu=new TH1D("mc_rw_midmu","", mc_rw_midmu_orig->GetNbinsX(), mc_rw_midmu_orig->GetXaxis()->GetXmin(), mc_rw_midmu_orig->GetXaxis()->GetXmax());
	TH1D *mc_rw_mideg=new TH1D("mc_rw_mideg","", mc_rw_mideg_orig->GetNbinsX(), mc_rw_mideg_orig->GetXaxis()->GetXmin(), mc_rw_mideg_orig->GetXaxis()->GetXmax());
	TH1D *mc_rw_midother=new TH1D("mc_rw_midother","", mc_rw_midother_orig->GetNbinsX(), mc_rw_midother_orig->GetXaxis()->GetXmin(), mc_rw_midother_orig->GetXaxis()->GetXmax());
        for (int i = 1; i<=mc_rw_orig->GetNbinsX(); ++i){
		double cen=mc_rw_orig->GetBinCenter(i);

		if (cen<=40) { //make eff corr if trklen < 40 cm
			mc_rw->SetBinContent(i, (double)mc_rw_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_rw->SetBinError(i, (double)mc_rw_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_rw_inel->SetBinContent(i, (double)mc_rw_inel_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_rw_inel->SetBinError(i, (double)mc_rw_inel_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_rw_el->SetBinContent(i, (double)mc_rw_el_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_rw_el->SetBinError(i, (double)mc_rw_el_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_rw_midcosmic->SetBinContent(i, (double)mc_rw_midcosmic_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_rw_midcosmic->SetBinError(i, (double)mc_rw_midcosmic_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_rw_midpi->SetBinContent(i, (double)mc_rw_midpi_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_rw_midpi->SetBinError(i, (double)mc_rw_midpi_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_rw_midp->SetBinContent(i, (double)mc_rw_midp_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_rw_midp->SetBinError(i, (double)mc_rw_midp_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_rw_midmu->SetBinContent(i, (double)mc_rw_midmu_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_rw_midmu->SetBinError(i, (double)mc_rw_midmu_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_rw_mideg->SetBinContent(i, (double)mc_rw_mideg_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_rw_mideg->SetBinError(i, (double)mc_rw_mideg_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

			mc_rw_midother->SetBinContent(i, (double)mc_rw_midother_orig->GetBinContent(i)/(double)gr_mc->Eval(cen));
			mc_rw_midother->SetBinError(i, (double)mc_rw_midother_orig->GetBinError(i)/(double)gr_mc->Eval(cen));

		} //make eff corr if trklen < 40 cm
	}









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
	float xmax=50;
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
/*
        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0.02);
        pad2->SetBottomMargin(0.25);
        pad2->SetGridx();
        pad2->SetGridy();
        pad2->Draw();
        pad2->cd();
*/

        TH2D *f2d2=new TH2D("f2d2",Form("%s","Pos"),100,xmin,xmax,64,0,2); //liny
        //TH2D *f2d2=new TH2D("f2d2",Form("%s","Pos"),100,xmin,xmax,64,0,6.4); //liny
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
