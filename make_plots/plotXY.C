#include "../headers/BasicParameters.h"

void plotXY(TString fdata, TString fmc, TString outpath){

	//plot style
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(1);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);

	//read data
	TFile *f_data = TFile::Open(fdata.Data());
        TH2D *h2d_xy_noSCE=(TH2D *)f_data->Get("h2d_xy_noSCE");
        TH2D *h2d_xy_SCE=(TH2D *)f_data->Get("h2d_xy_SCE");
	h2d_xy_noSCE->SetTitle("Data (Prod 4 Reco2): No SCE Correction; X^{1st hit} [cm]; Y^{1st hit} [cm]");
	h2d_xy_SCE->SetTitle("Data (Prod 4 Reco2): With SCE Correction; X^{1st hit} [cm]; Y^{1st hit} [cm]");
	
	//read MC
	TFile *f_mc = TFile::Open(fmc.Data());
        TH2D *hmc2d_xy_noSCE=(TH2D *)f_mc->Get("h2d_xy_noSCE");
        TH2D *hmc2d_xy_SCE=(TH2D *)f_mc->Get("h2d_xy_SCE");
	hmc2d_xy_noSCE->SetTitle("MC (Prod 4a): No SCE Correction; X^{1st hit} [cm]; Y^{1st hit} [cm]");
	hmc2d_xy_SCE->SetTitle("MC (Prod 4a): With SCE Correction; X^{1st hit} [cm]; Y^{1st hit} [cm]");


        //TH2D* frame2d=new TH2D("frame2d","", 400, 200, 600, 400, 0, 400); //zend_2d
        //frame2d->SetTitle("Reco Stopping Protons;Proton KE [MeV];Counts");

	TEllipse *el_data_noSCE=new TEllipse(mean_x, mean_y, dev_x, dev_y);
	el_data_noSCE->SetLineColor(2);
	el_data_noSCE->SetLineWidth(2);
	el_data_noSCE->SetFillStyle(0);
	
        TCanvas *c0=new TCanvas("c0","");
        c0->Divide(1,1);
        c0->cd(1);
        h2d_xy_noSCE->Draw("colz");
	el_data_noSCE->Draw();
        c0->Print(Form("%s/XY_data_noSCE.eps",outpath.Data()));

        TCanvas *c1=new TCanvas("c1","");
        c1->Divide(1,1);
        c1->cd(1);
        h2d_xy_SCE->Draw("colz");
        c1->Print(Form("%s/XY_data_SCE.eps",outpath.Data()));

	TEllipse *el_mc_noSCE=new TEllipse(mean_x_mc, mean_y_mc, dev_x_mc, dev_y_mc);
	el_mc_noSCE->SetLineColor(2);
	el_mc_noSCE->SetLineWidth(2);
	el_mc_noSCE->SetFillStyle(0);

        TCanvas *c2=new TCanvas("c2","");
        c2->Divide(1,1);
        c2->cd(1);
        hmc2d_xy_noSCE->Draw("colz");
	el_mc_noSCE->Draw();
        c2->Print(Form("%s/XY_mc_noSCE.eps",outpath.Data()));

        TCanvas *c3=new TCanvas("c3","");
        c3->Divide(1,1);
        c3->cd(1);
        hmc2d_xy_SCE->Draw("colz");
        c3->Print(Form("%s/XY_mc_SCE.eps",outpath.Data()));




        //TLegend *leg6 = new TLegend(0.15,0.65,0.75,0.85);
        //leg6->SetFillStyle(0);
        //leg6->AddEntry(KE_ff_recostop, "KE_{ff}", "l");
        //leg6->AddEntry(KE_simide_recostop, "KE_{simIDE}", "l");
        //leg6->Draw();


}
