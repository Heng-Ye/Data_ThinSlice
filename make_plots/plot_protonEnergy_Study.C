

void plot_protonEnergy_Study() {

	//TString fmc=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_kecalobkg_bmrw_beamxy_new6.root");
	TString fmc=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_kecalobkg_bmrw_beamxy_new9.root");
	//TString fdata=fmc; //for test
	TString fdata="../data_kebkg_beamxy_new8.root";

	TString fout_path="./plot_KEstudy/";


	//plot style ------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	//gStyle->SetOptStat(1);
	//-----------------------------------------------------------//

	//load data //
	TFile *f_data = TFile::Open(fdata.Data());
	//TH1D *edept_RecoEl_data=(TH1D*)f_data->Get("h1d_Edept_RecoEl");

	TH1D *edept_RecoEl_data=(TH1D*)f_data->Get("edept_RecoEl");

	TH1D *Rf_RecoEl_data=(TH1D*)f_data->Get("Rf_RecoEl");
	TH1D *Re_RecoEl_data=(TH1D*)f_data->Get("Re_RecoEl");
	TH2D *Rf_Re_RecoEl_data=(TH2D*)f_data->Get("Rf_Re_RecoEl");

	TH1D *Rf_RecoInEl_data=(TH1D*)f_data->Get("Rf_RecoInEl");
	TH1D *Re_RecoInEl_data=(TH1D*)f_data->Get("Re_RecoInEl");
	TH2D *Rf_Re_RecoInEl_data=(TH2D*)f_data->Get("Rf_Re_RecoInEl");

	TH2D *trklen_Lkeff_RecoEl_data=(TH2D*)f_data->Get("trklen_Lkeff_RecoEl");
	TH2D *trklen_Ledept_RecoEl_data=(TH2D*)f_data->Get("trklen_Ledept_RecoEl");
	TH2D *trklen_Lkeff_RecoInEl_data=(TH2D*)f_data->Get("trklen_Lkeff_RecoInEl");
	TH2D *trklen_Ledept_RecoInEl_data=(TH2D*)f_data->Get("trklen_Ledept_RecoInEl");

	TH2D *tf_keff_R_RecoEl_data=(TH2D*)f_data->Get("tf_keff_R_RecoEl");
	TH2D *tf_keff_R_RecoInEl_data=(TH2D*)f_data->Get("tf_keff_R_RecoInEl");

	int n_edept_RecoEl_data=edept_RecoEl_data->Integral(); 
	int n_Rf_RecoEl_data=Rf_RecoEl_data->Integral();

	//load MC //
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *edept_RecoEl_mc=(TH1D*)f_mc->Get("edept_RecoEl");

	TH1D *Rf_RecoEl_mc=(TH1D*)f_mc->Get("Rf_RecoEl");
	TH1D *Re_RecoEl_mc=(TH1D*)f_mc->Get("Re_RecoEl");
	TH2D *Rf_Re_RecoEl_mc=(TH2D*)f_mc->Get("Rf_Re_RecoEl");

	TH1D *Rf_RecoInEl_mc=(TH1D*)f_mc->Get("Rf_RecoInEl");
	TH1D *Re_RecoInEl_mc=(TH1D*)f_mc->Get("Re_RecoInEl");
	TH2D *Rf_Re_RecoInEl_mc=(TH2D*)f_mc->Get("Rf_Re_RecoInEl");

	TH2D *trklen_Lkeff_RecoEl_mc=(TH2D*)f_mc->Get("trklen_Lkeff_RecoEl");
	TH2D *trklen_Ledept_RecoEl_mc=(TH2D*)f_mc->Get("trklen_Ledept_RecoEl");
	TH2D *trklen_Lkeff_RecoInEl_mc=(TH2D*)f_mc->Get("trklen_Lkeff_RecoInEl");
	TH2D *trklen_Ledept_RecoInEl_mc=(TH2D*)f_mc->Get("trklen_Ledept_RecoInEl");

	TH2D *tf_keff_R_RecoEl_mc=(TH2D*)f_mc->Get("tf_keff_R_RecoEl");
	TH2D *tf_keff_R_RecoInEl_mc=(TH2D*)f_mc->Get("tf_keff_R_RecoInEl");

	int n_edept_RecoEl_mc=edept_RecoEl_mc->Integral(); 
	int n_Rf_RecoEl_mc=Rf_RecoEl_mc->Integral();


	//edept(El): data/MC ------------------------------------------------------------//
	//data-mc overlap
	TCanvas *c3=new TCanvas("c3","", 1200, 800);
	c3->Divide(1,1);
	c3->cd(1)->SetLogy();

 	float xmin=100;
	float xmax=700;
	TH2D* frame2d=new TH2D("frame2d","", 600, xmin, xmax, 400, 0, 400);
	frame2d->SetTitle("Reco. Elastic; Proton Energy [MeV]; Counts");
	frame2d->Draw();

	edept_RecoEl_data->SetTitle("Reco Elastic; Proton Energy [MeV]");
	edept_RecoEl_data->GetXaxis()->SetTitleOffset(1.2);

	edept_RecoEl_mc->Scale((float)n_edept_RecoEl_data/(float)n_edept_RecoEl_mc);
	edept_RecoEl_data->Draw("ep same");
	edept_RecoEl_mc->SetLineColor(2);
	edept_RecoEl_mc->Draw(" hist same");

        TLegend *leg = new TLegend(0.147681,0.705797,0.296796,0.87515);
	//leg->SetFillStyle(0);
	leg->SetLineColor(1);
	leg->SetFillColor(0);
   	leg->SetBorderSize(3);
	
	leg->AddEntry(edept_RecoEl_data, "Data", "ep");
	leg->AddEntry(edept_RecoEl_mc, "MC", "l");
	leg->Draw();
	c3->Print(Form("%sedept_data_mc_edept_recoEl.eps", fout_path.Data()));

	//ratio (data/MC) of edept
	TCanvas *c4=new TCanvas("c4","", 1200, 800);
	c4->Divide(1,1);
	c4->cd(1);
	TH2D* frame2d_r=new TH2D("frame2d_r","", 600, xmin, xmax, 15, -5, 10);
	frame2d_r->SetTitle("Reco. Elastic; Proton Energy [MeV]; Data/MC");
	frame2d_r->Draw();


	TH1D* R_edept_RecoEl_data=(TH1D*)edept_RecoEl_data->Clone();
	R_edept_RecoEl_data->Divide(edept_RecoEl_mc);
	R_edept_RecoEl_data->SetTitle("Reco Elastic; Proton Energy [MeV]; Data/MC");
	R_edept_RecoEl_data->GetXaxis()->SetTitleOffset(1.2);
	R_edept_RecoEl_data->Draw("ep same");


	TLine *line=new TLine(xmin,1,xmax,1);
	line->SetLineColor(2);
	line->SetLineWidth(2);
	line->SetLineStyle(2);
	line->Draw("same");
	c4->Print(Form("%sRatio_edept_data_mc_edept_recoEl.eps", fout_path.Data()));

	//Elastic -------------------------------------------------------------------------------------------//
	//ff
	TCanvas *c5f_el=new TCanvas("c5f_el","", 1200, 800);
	c5f_el->Divide(1,1);
	c5f_el->cd(1)->SetLogy();

 	float rf_min_el=0;
	float rf_max_el=3;
	TH2D* f2d_rf_el=new TH2D("f2d_rf_el","", 12, rf_min_el, rf_max_el, 50, 0.001, 0.8);
	f2d_rf_el->SetTitle("Reco. Elastic; R_{F}=RangeFromKE(KE_{beam}-#DeltaE)/(track length)");
	f2d_rf_el->GetXaxis()->SetTitleOffset(1.2);
	f2d_rf_el->Draw();

	Rf_RecoEl_data->Scale(1./(float)Rf_RecoEl_data->Integral());
	Rf_RecoEl_data->Draw("ep same");
	Rf_RecoEl_mc->SetLineColor(4);
	Rf_RecoEl_mc->Scale(1./(float)Rf_RecoEl_mc->Integral());
	Rf_RecoEl_mc->Draw("hist same");

        TLegend *leg_rf_el = new TLegend(0.693349,0.663055,0.842407,0.832685);
	//leg->SetFillStyle(0);
	leg_rf_el->SetLineColor(1);
	leg_rf_el->SetFillColor(0);
   	leg_rf_el->SetBorderSize(3);
	leg_rf_el->AddEntry(Rf_RecoEl_data, "Data", "ep");
	leg_rf_el->AddEntry(Rf_RecoEl_mc, "MC", "l");
	leg_rf_el->Draw();
	c5f_el->Print(Form("%sRf_recoEl.eps", fout_path.Data()));

	//end
	TCanvas *c5e_el=new TCanvas("c5e_el","", 1200, 800);
	c5e_el->Divide(1,1);
	c5e_el->cd(1)->SetLogy();

 	float re_min_el=0;
	float re_max_el=3;
	TH2D* f2d_re_el=new TH2D("f2d_re_el","", 12, re_min_el, re_max_el, 50, 0.001, 1);
	f2d_re_el->SetTitle("Reco. Elastic; R_{E}=RangeFromKE(E_{dept.})/(track length)");
	f2d_re_el->GetXaxis()->SetTitleOffset(1.2);
	f2d_re_el->Draw();

	Re_RecoEl_data->Scale(1./(float)Re_RecoEl_data->Integral());
	Re_RecoEl_data->Draw("ep same");
	Re_RecoEl_mc->SetLineColor(4);
	Re_RecoEl_mc->Scale(1./(float)Re_RecoEl_mc->Integral());
	Re_RecoEl_mc->Draw("hist same");

        TLegend *leg_re_el = new TLegend(0.693349,0.663055,0.842407,0.832685);
	//leg->SetFillStyle(0);
	leg_re_el->SetLineColor(1);
	leg_re_el->SetFillColor(0);
   	leg_re_el->SetBorderSize(3);
	leg_re_el->AddEntry(Re_RecoEl_data, "Data", "ep");
	leg_re_el->AddEntry(Re_RecoEl_mc, "MC", "l");
	leg_re_el->Draw();
	c5e_el->Print(Form("%sRe_recoEl.eps", fout_path.Data()));

	//ff(x) vs end(y): el_data
	TCanvas *c5ef_el_data=new TCanvas("c5ef_el_data","", 1200, 800);
	c5ef_el_data->Divide(1,1);

	TH2D* f2d_re_el_data=new TH2D("f2d_re_el_data","", 40,0,4,40,0,4);
	f2d_re_el_data->SetTitle("Reco. Elastic (Data); R_{F}; R_{E}");
	f2d_re_el_data->GetXaxis()->SetTitleOffset(1.2);
	f2d_re_el_data->Draw("");
	Rf_Re_RecoEl_data->Draw("colz same");

	TLine* l_hor_el_data=new TLine(0, 1, 4, 1);
	TLine* l_ver_el_data=new TLine(1, 0, 1, 4);
	l_hor_el_data->SetLineColor(2);
	l_ver_el_data->SetLineColor(2);
	l_hor_el_data->SetLineStyle(2);
	l_ver_el_data->SetLineStyle(2);
	l_hor_el_data->Draw();
	l_ver_el_data->Draw();
	c5ef_el_data->Print(Form("%sRf_Re_recoEl_data.eps", fout_path.Data()));


	//ff(x) vs end(y): el_mc
	TCanvas *c5ef_el_mc=new TCanvas("c5ef_el_mc","", 1200, 800);
	c5ef_el_mc->Divide(1,1);

	TH2D* f2d_re_el_mc=new TH2D("f2d_re_el_mc","", 40,0,4,40,0,4);
	f2d_re_el_mc->SetTitle("Reco. Elastic (MC); R_{F}; R_{E}");
	f2d_re_el_mc->GetXaxis()->SetTitleOffset(1.2);
	f2d_re_el_mc->Draw("");
	Rf_Re_RecoEl_mc->Draw("colz same");

	TLine* l_hor_el_mc=new TLine(0, 1, 4, 1);
	TLine* l_ver_el_mc=new TLine(1, 0, 1, 4);
	l_hor_el_mc->SetLineColor(2);
	l_ver_el_mc->SetLineColor(2);
	l_hor_el_mc->SetLineStyle(2);
	l_ver_el_mc->SetLineStyle(2);
	l_hor_el_mc->Draw();
	l_ver_el_mc->Draw();
	c5ef_el_mc->Print(Form("%sRf_Re_recoEl_mc.eps", fout_path.Data()));



	//mc_el 
	TCanvas *c5_trklen_lkeff_el_mc=new TCanvas("c5_trklen_lkeff_el_mc","", 1200, 800);
	c5_trklen_lkeff_el_mc->Divide(1,1);
	TH2D* f2d_trklen_lkeff_el_mc=new TH2D("f2d_trklen_lkeff_el_mc","", 140, 0, 140, 140, 0, 140);
	f2d_trklen_lkeff_el_mc->SetTitle("Reco. Elastic (MC); Reco Track Length [cm]; RangeFromKE(KE_{beam}-#DeltaE) [cm]");
	f2d_trklen_lkeff_el_mc->GetXaxis()->SetTitleOffset(1.2);
	f2d_trklen_lkeff_el_mc->Draw("");
	trklen_Lkeff_RecoEl_mc->Draw("colz same");

	TLine* l_trklen_dia=new TLine(0,0,140,140);
	l_trklen_dia->SetLineColor(2);
	l_trklen_dia->SetLineStyle(2);
	l_trklen_dia->Draw();
	c5_trklen_lkeff_el_mc->Print(Form("%strklen_lkeff_recoEl_mc.eps", fout_path.Data()));

	TCanvas *c5_trklen_Ledept_RecoEl_mc=new TCanvas("c5_trklen_Ledept_RecoEl_mc","", 1200, 800);
	c5_trklen_Ledept_RecoEl_mc->Divide(1,1);
	TH2D* f2d_trklen_Ledept_RecoEl_mc=new TH2D("f2d_trklen_Ledept_RecoEl_mc","", 140, 0, 140, 140, 0, 140);
	f2d_trklen_Ledept_RecoEl_mc->SetTitle("Reco. Elastic (MC); Reco Track Length [cm]; RangeFromKE(E_{dept.}) [cm]");
	f2d_trklen_Ledept_RecoEl_mc->GetXaxis()->SetTitleOffset(1.2);
	f2d_trklen_Ledept_RecoEl_mc->Draw("");
	trklen_Ledept_RecoEl_mc->Draw("colz same");
	l_trklen_dia->Draw();
	c5_trklen_Ledept_RecoEl_mc->Print(Form("%strklen_ledept_recoEl_mc.eps", fout_path.Data()));


	//data_el 
	TCanvas *c5_trklen_lkeff_el_data=new TCanvas("c5_trklen_lkeff_el_data","", 1200, 800);
	c5_trklen_lkeff_el_data->Divide(1,1);
	TH2D* f2d_trklen_lkeff_el_data=new TH2D("f2d_trklen_lkeff_el_data","", 140, 0, 140, 140, 0, 140);
	f2d_trklen_lkeff_el_data->SetTitle("Reco. Elastic (Data); Reco Track Length [cm]; RangeFromKE(KE_{beam}-#DeltaE) [cm]");
	f2d_trklen_lkeff_el_data->GetXaxis()->SetTitleOffset(1.2);
	f2d_trklen_lkeff_el_data->Draw("");
	trklen_Lkeff_RecoEl_data->Draw("colz same");
	l_trklen_dia->Draw();
	c5_trklen_lkeff_el_data->Print(Form("%strklen_lkeff_recoEl_data.eps", fout_path.Data()));

	TCanvas *c5_trklen_Ledept_RecoEl_data=new TCanvas("c5_trklen_Ledept_RecoEl_data","", 1200, 800);
	c5_trklen_Ledept_RecoEl_data->Divide(1,1);
	TH2D* f2d_trklen_Ledept_RecoEl_data=new TH2D("f2d_trklen_Ledept_RecoEl_data","", 140, 0, 140, 140, 0, 140);
	f2d_trklen_Ledept_RecoEl_data->SetTitle("Reco. Elastic (Data); Reco Track Length [cm]; RangeFromKE(E_{dept.}) [cm]");
	f2d_trklen_Ledept_RecoEl_data->GetXaxis()->SetTitleOffset(1.2);
	f2d_trklen_Ledept_RecoEl_data->Draw("");
	trklen_Ledept_RecoEl_data->Draw("colz same");
	l_trklen_dia->Draw();
	c5_trklen_Ledept_RecoEl_data->Print(Form("%strklen_ledept_recoEl_data.eps", fout_path.Data()));


	//data_el:  
	TCanvas *c5_tf_keff_R_RecoEl_data=new TCanvas("c5_tf_keff_R_RecoEl_data","", 1200, 800);
	c5_tf_keff_R_RecoEl_data->Divide(1,1);
	TH2D* f2d_tf_keff_R_RecoEl_data=new TH2D("f2d_tf_keff_R_RecoEl_data","", 40, 200, 600, 20, 0, 3.5);
	f2d_tf_keff_R_RecoEl_data->SetTitle("Reco. Elastic (Data); KE_{beam}-#DeltaE [MeV]; R_{F}=RangeFromKE(KE_{beam}-#DeltaE)/(track length)");
	f2d_tf_keff_R_RecoEl_data->GetXaxis()->SetTitleOffset(1.2);
	f2d_tf_keff_R_RecoEl_data->Draw("");
	tf_keff_R_RecoEl_data->Draw("colz same");

	TLine* l_hor=new TLine(200,1,600,1);
	l_hor->SetLineColor(2);
	l_hor->SetLineStyle(2);
	l_hor->Draw();

	c5_tf_keff_R_RecoEl_data->Print(Form("%skeff_rf_recoEl_data.eps", fout_path.Data()));


	//mc_el:  
	TCanvas *c5_tf_keff_R_RecoEl_mc=new TCanvas("c5_tf_keff_R_RecoEl_mc","", 1200, 800);
	c5_tf_keff_R_RecoEl_mc->Divide(1,1);
	TH2D* f2d_tf_keff_R_RecoEl_mc=new TH2D("f2d_tf_keff_R_RecoEl_mc","", 40, 200, 600, 20, 0, 3.5);
	f2d_tf_keff_R_RecoEl_mc->SetTitle("Reco. Elastic (MC); KE_{beam}-#DeltaE [MeV]; R_{F}=RangeFromKE(KE_{beam}-#DeltaE)/(track length)");
	f2d_tf_keff_R_RecoEl_mc->GetXaxis()->SetTitleOffset(1.2);
	f2d_tf_keff_R_RecoEl_mc->Draw("");
	tf_keff_R_RecoEl_mc->Draw("colz same");
	l_hor->Draw();
	c5_tf_keff_R_RecoEl_mc->Print(Form("%skeff_rf_recoEl_mc.eps", fout_path.Data()));





	//Inelastic -------------------------------------------------------------------------------------------//
	//ff
	TCanvas *c5f_inel=new TCanvas("c5f_inel","", 1200, 800);
	c5f_inel->Divide(1,1);
	c5f_inel->cd(1)->SetLogy();

 	float rf_min_inel=0;
	float rf_max_inel=10;
	TH2D* f2d_rf_inel=new TH2D("f2d_rf_inel","", 12, rf_min_inel, rf_max_inel, 50, 0.001, 0.5);
	f2d_rf_inel->SetTitle("Reco. Inelastic; R_{F}=RangeFromKE(KE_{beam}-#DeltaE)/(track length)");
	f2d_rf_inel->GetXaxis()->SetTitleOffset(1.2);
	f2d_rf_inel->Draw();

	Rf_RecoInEl_data->Scale(1./(float)Rf_RecoInEl_data->Integral());
	Rf_RecoInEl_data->Draw("ep same");
	Rf_RecoInEl_mc->SetLineColor(2);
	Rf_RecoInEl_mc->Scale(1./(float)Rf_RecoInEl_mc->Integral());
	Rf_RecoInEl_mc->Draw("hist same");

        TLegend *leg_rf_inel = new TLegend(0.693349,0.663055,0.842407,0.832685);
	//leg->SetFillStyle(0);
	leg_rf_inel->SetLineColor(1);
	leg_rf_inel->SetFillColor(0);
   	leg_rf_inel->SetBorderSize(3);
	leg_rf_inel->AddEntry(Rf_RecoInEl_data, "Data", "ep");
	leg_rf_inel->AddEntry(Rf_RecoInEl_mc, "MC", "l");
	leg_rf_inel->Draw();
	c5f_inel->Print(Form("%sRf_recoInEl.eps", fout_path.Data()));


	//end
	TCanvas *c5e_inel=new TCanvas("c5e_inel","", 1200, 800);
	c5e_inel->Divide(1,1);
	c5e_inel->cd(1)->SetLogy();

 	float re_min_inel=0;
	float re_max_inel=10;
	TH2D* f2d_re_inel=new TH2D("f2d_re_inel","", 12, re_min_inel, re_max_inel, 50, 0.001, 1);
	f2d_re_inel->SetTitle("Reco. Inelastic; R_{E}=RangeFromKE(E_{dept.})/(track length)");
	f2d_re_inel->GetXaxis()->SetTitleOffset(1.2);
	f2d_re_inel->Draw();

	Re_RecoInEl_data->Scale(1./(float)Re_RecoInEl_data->Integral());
	Re_RecoInEl_data->Draw("ep same");
	Re_RecoInEl_mc->SetLineColor(2);
	Re_RecoInEl_mc->Scale(1./(float)Re_RecoInEl_mc->Integral());
	Re_RecoInEl_mc->Draw("hist same");

        TLegend *leg_re_inel = new TLegend(0.693349,0.663055,0.842407,0.832685);
	//leg->SetFillStyle(0);
	leg_re_inel->SetLineColor(1);
	leg_re_inel->SetFillColor(0);
   	leg_re_inel->SetBorderSize(3);
	leg_re_inel->AddEntry(Re_RecoEl_data, "Data", "ep");
	leg_re_inel->AddEntry(Re_RecoEl_mc, "MC", "l");
	leg_re_inel->Draw();
	c5e_inel->Print(Form("%sRe_recoInEl.eps", fout_path.Data()));


	//ff(x) vs end(y): inel_data
	TCanvas *c5ef_inel_data=new TCanvas("c5ef_inel_data","", 1200, 800);
	c5ef_inel_data->Divide(1,1);

	TH2D* f2d_re_inel_data=new TH2D("f2d_re_inel_data","", 40,0,4,40,0,4);
	f2d_re_inel_data->SetTitle("Reco. Inelastic (Data); R_{F}; R_{E}");
	f2d_re_inel_data->GetXaxis()->SetTitleOffset(1.2);
	f2d_re_inel_data->Draw("");
	Rf_Re_RecoInEl_data->Draw("colz same");

	TLine* l_hor_inel_data=new TLine(0, 1, 4, 1);
	TLine* l_ver_inel_data=new TLine(1, 0, 1, 4);
	l_hor_inel_data->SetLineColor(2);
	l_ver_inel_data->SetLineColor(2);
	l_hor_inel_data->SetLineStyle(2);
	l_ver_inel_data->SetLineStyle(2);
	l_hor_inel_data->Draw();
	l_ver_inel_data->Draw();
	c5ef_inel_data->Print(Form("%sRf_Re_recoInEl_data.eps", fout_path.Data()));


	//ff(x) vs end(y): inel_mc
	TCanvas *c5ef_inel_mc=new TCanvas("c5ef_inel_mc","", 1200, 800);
	c5ef_inel_mc->Divide(1,1);

	TH2D* f2d_re_inel_mc=new TH2D("f2d_re_inel_mc","", 40,0,4,40,0,4);
	f2d_re_inel_mc->SetTitle("Reco. Inelastic (MC); R_{F}; R_{E}");
	f2d_re_inel_mc->GetXaxis()->SetTitleOffset(1.2);
	f2d_re_inel_mc->Draw("");
	Rf_Re_RecoInEl_mc->Draw("colz same");

	TLine* l_hor_inel_mc=new TLine(0, 1, 4, 1);
	TLine* l_ver_inel_mc=new TLine(1, 0, 1, 4);
	l_hor_inel_mc->SetLineColor(2);
	l_ver_inel_mc->SetLineColor(2);
	l_hor_inel_mc->SetLineStyle(2);
	l_ver_inel_mc->SetLineStyle(2);
	l_hor_inel_mc->Draw();
	l_ver_inel_mc->Draw();
	c5ef_inel_mc->Print(Form("%sRf_Re_recoInEl_mc.eps", fout_path.Data()));




	//mc_inel 
	TCanvas *c5_trklen_lkeff_inel_mc=new TCanvas("c5_trklen_lkeff_inel_mc","", 1200, 800);
	c5_trklen_lkeff_inel_mc->Divide(1,1);
	TH2D* f2d_trklen_lkeff_inel_mc=new TH2D("f2d_trklen_lkeff_inel_mc","", 140, 0, 140, 140, 0, 140);
	f2d_trklen_lkeff_inel_mc->SetTitle("Reco. InElastic (MC); Reco Track Length [cm]; RangeFromKE(KE_{beam}-#DeltaE) [cm]");
	f2d_trklen_lkeff_inel_mc->GetXaxis()->SetTitleOffset(1.2);
	f2d_trklen_lkeff_inel_mc->Draw("");
	trklen_Lkeff_RecoInEl_mc->Draw("colz same");
	l_trklen_dia->Draw();
	c5_trklen_lkeff_inel_mc->Print(Form("%strklen_lkeff_recoInEl_mc.eps", fout_path.Data()));

	TCanvas *c5_trklen_Ledept_RecoInEl_mc=new TCanvas("c5_trklen_Ledept_RecoInEl_mc","", 1200, 800);
	c5_trklen_Ledept_RecoInEl_mc->Divide(1,1);
	TH2D* f2d_trklen_Ledept_RecoInEl_mc=new TH2D("f2d_trklen_Ledept_RecoInEl_mc","", 140, 0, 140, 140, 0, 140);
	f2d_trklen_Ledept_RecoInEl_mc->SetTitle("Reco. Inelastic (MC); Reco Track Length [cm]; RangeFromKE(E_{dept.}) [cm]");
	f2d_trklen_Ledept_RecoInEl_mc->GetXaxis()->SetTitleOffset(1.2);
	f2d_trklen_Ledept_RecoInEl_mc->Draw("");
	trklen_Ledept_RecoInEl_mc->Draw("colz same");
	l_trklen_dia->Draw();
	c5_trklen_Ledept_RecoInEl_mc->Print(Form("%strklen_ledept_recoInEl_mc.eps", fout_path.Data()));


	//data_inel
	TCanvas *c5_trklen_lkeff_inel_data=new TCanvas("c5_trklen_lkeff_inel_data","", 1200, 800);
	c5_trklen_lkeff_inel_data->Divide(1,1);
	TH2D* f2d_trklen_lkeff_inel_data=new TH2D("f2d_trklen_lkeff_inel_data","", 140, 0, 140, 140, 0, 140);
	f2d_trklen_lkeff_inel_data->SetTitle("Reco. InElastic (Data); Reco Track Length [cm]; RangeFromKE(KE_{beam}-#DeltaE) [cm]");
	f2d_trklen_lkeff_inel_data->GetXaxis()->SetTitleOffset(1.2);
	f2d_trklen_lkeff_inel_data->Draw("");
	trklen_Lkeff_RecoInEl_data->Draw("colz same");
	l_trklen_dia->Draw();
	c5_trklen_lkeff_inel_data->Print(Form("%strklen_lkeff_recoInEl_data.eps", fout_path.Data()));

	TCanvas *c5_trklen_Ledept_RecoInEl_data=new TCanvas("c5_trklen_Ledept_RecoInEl_data","", 1200, 800);
	c5_trklen_Ledept_RecoInEl_data->Divide(1,1);
	TH2D* f2d_trklen_Ledept_RecoInEl_data=new TH2D("f2d_trklen_Ledept_RecoInEl_data","", 140, 0, 140, 140, 0, 140);
	f2d_trklen_Ledept_RecoInEl_data->SetTitle("Reco. Inelastic (Data); Reco Track Length [cm]; RangeFromKE(E_{dept.}) [cm]");
	f2d_trklen_Ledept_RecoInEl_data->GetXaxis()->SetTitleOffset(1.2);
	f2d_trklen_Ledept_RecoInEl_data->Draw("");
	trklen_Ledept_RecoInEl_data->Draw("colz same");
	l_trklen_dia->Draw();
	c5_trklen_Ledept_RecoInEl_data->Print(Form("%strklen_ledept_recoInEl_data.eps", fout_path.Data()));


	//data_inel:  
	TCanvas *c5_tf_keff_R_RecoInEl_data=new TCanvas("c5_tf_keff_R_RecoInEl_data","", 1200, 800);
	c5_tf_keff_R_RecoInEl_data->Divide(1,1);
	TH2D* f2d_tf_keff_R_RecoInEl_data=new TH2D("f2d_tf_keff_R_RecoInEl_data","", 40, 200, 600, 20, 0, 10);
	f2d_tf_keff_R_RecoInEl_data->SetTitle("Reco. Inelastic (Data); KE_{beam}-#DeltaE [MeV]; R_{F}=RangeFromKE(KE_{beam}-#DeltaE)/(track length)");
	f2d_tf_keff_R_RecoInEl_data->GetXaxis()->SetTitleOffset(1.2);
	f2d_tf_keff_R_RecoInEl_data->Draw("");
	tf_keff_R_RecoInEl_data->Draw("colz same");
	l_hor->Draw();
	c5_tf_keff_R_RecoInEl_data->Print(Form("%skeff_rf_recoInEl_data.eps", fout_path.Data()));

	//mc_inel:  
	TCanvas *c5_tf_keff_R_RecoInEl_mc=new TCanvas("c5_tf_keff_R_RecoInEl_mc","", 1200, 800);
	c5_tf_keff_R_RecoInEl_mc->Divide(1,1);
	TH2D* f2d_tf_keff_R_RecoInEl_mc=new TH2D("f2d_tf_keff_R_RecoInEl_mc","", 40, 200, 600, 20, 0, 10);
	f2d_tf_keff_R_RecoInEl_mc->SetTitle("Reco. Inelastic (MC); KE_{beam}-#DeltaE [MeV]; R_{F}=RangeFromKE(KE_{beam}-#DeltaE)/(track length)");
	f2d_tf_keff_R_RecoInEl_mc->GetXaxis()->SetTitleOffset(1.2);
	f2d_tf_keff_R_RecoInEl_mc->Draw("");
	tf_keff_R_RecoInEl_mc->Draw("colz same");
	l_hor->Draw();
	c5_tf_keff_R_RecoInEl_mc->Print(Form("%skeff_rf_recoInEl_mc.eps", fout_path.Data()));





}
