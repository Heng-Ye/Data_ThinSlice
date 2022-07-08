#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>

#include "TString.h"
#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"
#include "TChain.h"
#include "color_text.h"

using namespace std;

int main(int argc, char *argv[], char *envp[]) {
//int main() {
 assert(argc == 3);
 //TString str_in(argv[1]);
 TString fgoodlist(argv[1]);
 //TString runnum(argv[2]);
 TString selector_name(argv[2]);

 //TString fgoodlist="goodfile_list.txt";

 //TString str_out(argv[2]);
 //TApplication theApp("App", &argc, argv); //necessary for tree

 //setup output file names ---------------------//
 //std::cout<<"Input root file path : "<<str_in<<std::endl;
 //std::cout<<"Output root file : "<<str_out<<std::endl;
 
 ifstream f_string(fgoodlist.Data());
 string buffer_string;
 vector<string> filename;
 while (f_string>>buffer_string) { filename.push_back(buffer_string); }
  f_string.close();

 std::cout<<"Looping over the good beam data list : "<<yellow<<fgoodlist.Data()<<std::endl;

 //TChain *chain=new TChain("protonanalysis/beamana");
 TChain *chain=new TChain("protonbeamana/beamana");
 for (size_t ii=0; ii<filename.size();ii++) { //loop over all the selected ana files
   //TString inputfilename(filename[ii]);
   std::cout<<green<<"--> reading beam data:  "<<filename[ii]<<std::endl;

   chain->Add(filename[ii].c_str());
   
 } //loop over all the selected ana files

 std::cout<<reset<<"\n"<<std::endl;


 //---------------------------------------------//
 //TString exe_comm=Form("root -b -q 'HealthCheck_Files.C(\"%s\")'",str_in.Data());
 //std::cout<<exe_comm<<std::endl; 
 //root -b -q 'HealthCheck_Files.C("/pnfs/dune/scratch/users/hyliao/v07_11_00/ana/protodune_1gev_180kv_run5387/")'
 //gSystem->Exec(exe_comm.Data()); 

 //get the related branches from the old tree -------------------------
 //TFile *fin=new TFile(str_in.Data());
 //
 //TChain *chain=new TChain("protonanalysis/beamana");
 //chain->Add("/pnfs/dune/scratch/users/hyliao/v07_11_00/ana/protodune_1gev_180kv_run5387/13784923_1/Beam.root");
 //chain->Add("/pnfs/dune/scratch/users/hyliao/v07_11_00/ana/protodune_1gev_180kv_run5387/13784935_3/Beam.root");
 //chain->GetEntries();
 
 //chain->Print();
 //gDirectory->pwd();
 //HY::Problem while tchain enter directory


 //TDirectory* dir=fin->GetDirectory("protonanalysis");
 //TDirectory* dir=chain->GetDirectory("protonanalysis");
 //dir->ls();

 //get tree
 //TTree *beamana=(TTree*)dir->Get("beamana"); //get basic tree
 //TChain *beamana=(TChain*)dir->Get("beamana"); //get basic tree
 //tr->SetName("Event");
 //tr->MakeSelector("CosmicSNSelector");
 //beamana->Print();
 //chain->MakeClass(Form("ProtonSelector_run%s",runnum.Data()));
 chain->MakeClass(Form("%s",selector_name.Data()));
 //beamana->MakeClass("ProtonSelector");
 //tr->Print();
 //fin->Close();




}                   
