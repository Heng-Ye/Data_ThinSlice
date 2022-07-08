
//void RunAna(int run){
void RunAna(TString str){
 //load selector
 //TString load_selector=Form(".L ProtonSelector_run%d.C+",run);
 TString load_selector=Form(".L %s.C+",str.Data());
 gROOT->ProcessLine(load_selector);

 //TString selector=Form("ProtonSelector_run%d ANA",run);
 TString selector=Form("%s ANA",str.Data());

 gROOT->ProcessLine(selector);
 gROOT->ProcessLine("ANA.Loop()"); 
 //gROOT->ProcessLine("ANA.Show(0)"); 





}
